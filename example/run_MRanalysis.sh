#!/usr/bin/env bash
#@auther: chenbin
#@created time: 2023/04/25
#@email: chenbin_6901@163.com/a1030539294@gmail.com
set -uxeo pipefail

vcf_input="Wheat_SNPs.vcf.gz"

## Extract the GT from the VCF file using bcftools
bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%TGT]\n' $vcf_input \
 | sed -r -e '1 s@:GT@@g;1 s@\[[0-9]+\]@@g;1 s@# @@' \
 | bgzip > SNPgeno.GT.vcf.gz

## Get the *.bim file corresponding to the VCF file using plink
plink --vcf $vcf_input --make-bed -out SNPgeno --allow-extra-chr

###
python3 ../scripts/VCFgt2num_converter.py\
 -o SNPgeno_MR.txt\
 -v SNPgeno.GT.vcf.gz\
 -b SNPgeno.bim

###
Rscript ../scripts/MR_analysis.R -G SNPgeno_MR.txt -E MR_expr.txt -P MR_pheno.txt -Q MR_eQTLs.txt -d ./ -p MR_result_test
