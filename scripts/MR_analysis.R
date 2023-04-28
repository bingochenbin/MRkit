############################################################
#  @author:Bin Chen
#  @created time:2023/4/24
#  @email: chenbin_6901@163.com/a1030539294@gmail.com
#  @comment: Perform the MR analysis to reveal the relationship between gene expression levels and the phenotype
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
############################################################
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table", update = F, ask = F)
if (!requireNamespace("foreach", quietly = TRUE)) install.packages("foreach", update = F, ask = F)
if (!requireNamespace("iterators", quietly = TRUE)) install.packages("iterators", update = F, ask = F)
if (!requireNamespace("doParallel", quietly = TRUE)) install.packages("doParallel", update = F, ask = F)
if (!requireNamespace("optparse", quietly = TRUE)) install.packages("optparse", update = F, ask = F)
library(data.table)
library(foreach)
library(iterators)
library(doParallel)
library(optparse)

option_list <- list(
  make_option(c("-G", "--geno_input"),
              action="store",
              metavar='FILENAME',
              type='character',
              help="The SNP genotyping data file <required>"),
  make_option(c("-E", "--expr_input"),
              action="store",
              metavar='FILENAME',
              type='character',
              help="The gene expression levels data file <required>"),
  make_option(c("-P", "--pheno_input"),
              action="store",
              metavar='FILENAME',
              type='character',
              help="The phenotype data file <required>"),
  make_option(c("-Q", "--eqtl_input"),
              action="store",
              metavar='FILENAME',
              type='character',
              help="The eQTLs data file <required>"),
  make_option(c("-d", "--outdir"),
              action="store",
              metavar='DIRECTORY',
              type='character',
              help="Path to outdir <required>"),
  make_option(c("-p", "--outprefix"),
              action="store",
              metavar='STR',
              type='character',
              help="Prefix to output filename <required>"),
  make_option(c("-n", "--cpu"),
              action="store",
              metavar='INT',
              default=1,
              type="integer",
              help="Number of cpu cores to use [default: %default]")
)

opt <- parse_args(OptionParser(usage = "%prog [-h] -G FILENAME -E FILENAME -P FILENAME -Q FILENAME -d DIRECTORY -p outprefix [-n INT]\n
Perform the MR analysis to reveal the relationship between gene expression levels and the phenotype.
Note: Please see the https://github.com/bingochenbin/MRkit for details.",option_list=option_list))
##################################################
## check the output directory
outdir <- normalizePath(opt$outdir)
if(!dir.exists(outdir)){stop("The directory (", outdir, ") doesn't exist! Please check!\n")}

## define the function to check the outprefix
prefix_checker <- function(prefix, allowed_chars = "[A-Za-z0-9_\\.]") {
  regex <- paste0("^", allowed_chars, "+$")
  
  if(!grepl(regex, prefix)) {
    stop("Invalid outprefix! Only alphanumeric characters and underscores are allowed.")
  }
  prefix
}
outprefix = prefix_checker(opt$outprefix)

# ## check if the files exist
if (file.access(opt$geno_input) == -1){
  stop("The file (", opt$geno_input, ") doesn't exist! Please check!\n")
}else{
  geno_input <- opt$geno_input
}

if (file.access(opt$expr_input) == -1){
  stop("The file (", opt$expr_input, ") doesn't exist! Please check!\n")
}else{
  expr_input <- opt$expr_input
}

if (file.access(opt$pheno_input) == -1){
  stop("The file (", opt$pheno_input, ") doesn't exist! Please check!\n")
}else{
  pheno_input <- opt$pheno_input
}

if (file.access(opt$eqtl_input) == -1){
  stop("The file (", opt$eqtl_input, ") doesn't exist! Please check!\n")
}else{
  eqtl_input <- opt$eqtl_input
}

if (opt$cpu <= 0){
  stop("The number of cpu provided is incorrect! Please check!\n")
}else{
  nchunk = opt$cpu
}

##################################################
##### Define the function to split the eQTLs data into parallel-execution-ready chunks
# ref: https://cran.r-project.org/web/packages/foreach/vignettes/foreach.html
plist_splitor = function(plist, chunks){
  n <- nrow(plist)
  i <- 1
  ichunk = 0  # record the number of chunk
  
  nextElem <- function(){
    if(n <= 0 || chunks <= 0) stop('StopIteration')
    m <- ceiling(n / chunks)
    r <- seq(i, length=m)
    i <<- i + m
    n <<- n - m
    chunks <<- chunks - 1
    ichunk <<- ichunk + 1
    list(ichunk=ichunk, block=plist[r])
  }
  
  structure(list(nextElem=nextElem), class=c('plist_splitor', 'iter'))
}
nextElem.plist_splitor <- function(obj) obj$nextElem()

ptime <- system.time({
##### Read in the data
## Read in the SNP genotyping data
genodat <- fread(file = geno_input, header = T, sep = "\t", drop = c(4,5))
if (nrow(genodat[, .N, by = ID][N > 1]) > 0) {
  stop('Duplicate SNP ID was found in the SNP genotyping data! Please check!\n')
}
# Get the SNPs and their physical positions on the genome
snpspos <- copy(genodat[, .(CHROM, POS, ID)])
# remove the SNPs' physical positions
genodat[, CHROM := NULL]
genodat[, POS := NULL]
# transpose
genodat <- data.table::transpose(genodat, make.names = "ID", keep.names = "ID")
names(genodat)[1] <- "SampleID"
# check the number of samples
nsample = genodat[, .N]
cat('The number of samples detected: ', nsample, '\n')
# Get all SNP IDs
all_snpids <- names(genodat)[-1]
# genodat[1:10, 1:10]

## Read in the eQTLs data
eqtldat <- fread(file = eqtl_input, header = T, sep = "\t")
names(eqtldat) <- c("GeneID", "eQTL", "SNP", "Pvalue")

## Read in the gene expression data
exprdat <- fread(file = expr_input, header = T, sep = "\t")
names(exprdat)[1] <- "SampleID"
all_geneids <- names(exprdat)[-1]

## Read in the phenotype data
phenodat <- fread(file = pheno_input, header = T, sep = "\t")
names(phenodat)[1] <- "SampleID"
phenovec <- phenodat[, 2][[1]]
phenovar <- var(phenovec, na.rm = T)
pheno_NA_index <- which(is.na(phenovec))  # get the index of NAs

## Check if the Sample IDs in all datasets are identical and in same order
if (!identical(genodat[, SampleID], exprdat[, SampleID], attrib.as.set = FALSE) && identical(genodat[, SampleID], phenodat[, SampleID], attrib.as.set = FALSE)){
  stop('The sample IDs in the input datasets are not identical or in the same order!\n')
}

##### 
## Get the top eQTLs most significantly associated with the expression levels of each gene
topeqtldat <- eqtldat[, .SD[which.min(Pvalue)], keyby = GeneID]
# topeqtldat

## check 
if (nchunk > topeqtldat[, .N]){
  cat('>>> Warning: The number of CPU cores provided is greater than the number of eQTLs',
      'most significantly associated with the expression levels of each gene,',
      'so only one CPU core is actually used...\n')
  nchunk = 1
}

##### Run the MR analysis
cl = makeCluster(nchunk, type="FORK", outfile="")
registerDoParallel(cl)
allMRres <- foreach(block=plist_splitor(plist=topeqtldat, chunks=nchunk),
        .combine = 'rbind',
        .packages = c('data.table', 'foreach'),
        .errorhandling = 'stop',
        .verbose = T) %dopar% {
          cat("Processing Block", block[[1]], "\n")
          MRres = foreach(topeqtl=iter(block[[2]], by='row'),
                          .combine = 'rbind',
                          .errorhandling = 'stop',
                          .verbose = F) %do% {
                            # Get the Gene ID and the SNP ID
                            geneid = topeqtl[, GeneID]
                            geneid_topeqtl_snpid = topeqtl[, SNP]
                            
                            ##### MR analysis
                            if(geneid_topeqtl_snpid %in% all_snpids && geneid %in% all_geneids){
                            # Get the SNP's genotypes across all samples
                            geneid_topeqtl_snpid_geno <- genodat[, ..geneid_topeqtl_snpid][[1]]
                            geno_NA_index <- which(is.na(geneid_topeqtl_snpid_geno))  # get the index of NAs
                            
                            # Get the gene expression level of the GeneID
                            geneid_expr <- exprdat[, ..geneid][[1]]
                            expr_NA_index <- which(is.na(geneid_expr))  # get the index of NAs
                            
                            # Get the number of NAs in this calculation
                            NA_count <- length(union(union(pheno_NA_index, geno_NA_index), expr_NA_index))
                            
                            ###
                            ## the least-square estimate of the gene expression level on the SNP
                            lm_expr_snp = lm(geneid_expr ~ geneid_topeqtl_snpid_geno)
                            #summary(lm_expr_snp)
                            beta_expr_snp = coefficients(lm_expr_snp)[2]
                            
                            ## the least-square estimate of the phenotype on the SNP
                            lm_pheno_snp = lm(phenovec ~ geneid_topeqtl_snpid_geno)
                            #summary(lm_pheno_snp)
                            beta_pheno_snp = coefficients(lm_pheno_snp)[2]
                            
                            ## the two-step least-squares (2SLS) estimate of the effect of the gene expression level on the phenotype
                            beta_expr_pheno = beta_pheno_snp/beta_expr_snp
                            
                            ###
                            ## The sampling variance of the 2SLS estimate
                            lm_expr_pheno <- lm(phenovec ~ geneid_expr)
                            r2_expr_pheno <- summary(lm_expr_pheno)$r.square
                            r2_expr_snp <- summary(lm_expr_snp)$r.square
                            
                            upper = phenovar*(1-r2_expr_pheno)
                            down = (nsample-NA_count)*var(geneid_expr, na.rm = T)*r2_expr_snp
                            var_expr_pheno = upper/down
                            
                            ## Use the test statistic TMR to test the significance of the 2SLS estimate of the effect
                            TMR = beta_expr_pheno^2 / var_expr_pheno
                            pvalue = pchisq(TMR, 1, lower.tail = F)
                            
                            ## output
                            res <- data.table(
                              GeneID = geneid,
                              TopSNP = geneid_topeqtl_snpid,
                              Effect = beta_expr_pheno,
                              TMR = TMR,
                              Pvalue = pvalue
                            )
                            
                            } else if(geneid_topeqtl_snpid %in% all_snpids){
                              cat('>>> Warning: No expression level available for the gene:', geneid, '\n')
                              res <- data.table(
                                GeneID = geneid,
                                TopSNP = geneid_topeqtl_snpid,
                                Effect = "NA",
                                TMR = "NA",
                                Pvalue = "NA"
                              )
                              
                            } else if(geneid %in% all_geneids){
                              cat('>>> Warning: No genotypes available for the SNP:', geneid_topeqtl_snpid, '\n')
                              res <- data.table(
                                GeneID = geneid,
                                TopSNP = geneid_topeqtl_snpid,
                                Effect = "NA",
                                TMR = "NA",
                                Pvalue = "NA"
                              )
                              
                            } else {
                              cat('>>> Warning: Either the expression level for the gene', geneid,
                                  'or the genotypes for the SNP', geneid_topeqtl_snpid, ' are not available\n')
                              res <- data.table(
                                GeneID = geneid,
                                TopSNP = geneid_topeqtl_snpid,
                                Effect = "NA",
                                TMR = "NA",
                                Pvalue = "NA"
                              )
                            }
                          }
          cat("Block", block[[1]], "Processed", "\n")
          MRres
        }
stopCluster(cl)

##### write out
# allMRres
MRres_out <- merge.data.table(
  x = allMRres,
  y = snpspos,
  by.x = "TopSNP",
  by.y = "ID",
  all = F,
  sort = F
)
write.table(MRres_out[, .(GeneID, TopSNP, CHROM, POS, Effect, TMR, Pvalue)],
            file = file.path(outdir, paste0(outprefix, ".MR.txt")),
            quote = F, sep = '\t', row.names = F, col.names = T)

})

# total time charged
print(ptime)

