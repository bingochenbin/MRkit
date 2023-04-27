#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################################
# Copyright 2023-2033 Bin Chen
#  @author:chenbin
#  @created time:2023/02/26
#  @email: chenbin_6901@163.com/a1030539294@gmail.com
#  @comment: Script to convert genotypes in VCF to numeric codes.
#            convert variates with double minor alleles to 2, variates with double major alleles to 0,
#            variates with mixing alleles to 1, and variates without alleles (eg. ./.) to NA.
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
__author__ = "chenbin"
__email__ = "chenbin_6901@163.com"
__date__ = "2023-02-26"
__version__ = "0.0.1"

import sys
import os
from os.path import join, isdir
from functools import wraps
from resource import getrusage, RUSAGE_SELF, RUSAGE_CHILDREN
import time
from collections import defaultdict
from itertools import islice
import gzip
import argparse
from multiprocessing import Pool
from typing import Generator, List, Dict


def pyversion_checker() -> None:
    """Check if the Python version is lower than v3.8.0"""
    if sys.hexversion < 0x30800f0:
        sys.exit('Error: The Python version being used is too low. You need to update the Python to v3.8.0 or newer.')
    else:
        print('Info: The Python version being used is satisfactory.')


def outdir_checker(outdir) -> str:
    """check the output directory"""
    if not isdir(outdir):
        raise OSError(f"The input outdir {outdir} is not an existing directory! Please check!")
    if not os.listdir(outdir):
        print(f"Warning: The input outdir {outdir} is an empty directory!", flush=True)
    outdir = os.path.join(outdir, '')  # ensure trailing slash

    return outdir


# src:Python-Cookbook-3rd-Edition
class Timer:
    def __init__(self, func=time.perf_counter):
        self.elapsed = 0.0
        self._func = func
        self._start = None

    def start(self):
        if self._start is not None:
            raise RuntimeError('Already started')
        self._start = self._func()

    def stop(self):
        if self._start is None:
            raise RuntimeError('Not started')
        end = self._func()
        self.elapsed += end - self._start
        self._start = None

    def reset(self):
        self.elapsed = 0.0

    @property
    def running(self):
        return self._start is not None

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, *args):
        self.stop()


def stat_resource_usage(func):
    """
    To stat the time and memory usage consumed by the process and its child processes roughly, which executing the
    function
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        # start the timer
        t_perf_counter = Timer(time.perf_counter)
        t_perf_counter.start()
        # the function to be wrapped
        result = func(*args, **kwargs)
        # stop the timer and retrieve resource usage information
        t_perf_counter.stop()
        calling_ru_utime, calling_ru_stime, calling_ru_maxrss, *ign = getrusage(RUSAGE_SELF)
        child_ru_utime, child_ru_stime, child_ru_maxrss, *ign = getrusage(RUSAGE_CHILDREN)
        print('real time: {0}s\nuser time: {1}s\nsys time: {2}s\nmemory usage: {3}K'.format(
            t_perf_counter.elapsed, calling_ru_utime + child_ru_utime,
            calling_ru_stime + child_ru_stime, calling_ru_maxrss + child_ru_maxrss))

        return result
    return wrapper


def parser_bim(genobim) -> Dict:
    """parse *.bim file corresponding to the VCF.
    """
    converter_dict = defaultdict(dict)
    with open(genobim, 'rt') as bf:
        for line in bf:
            _, snpid, _, _, minoral, majoral = line.strip().split('\t')
            # unphased genotypes
            double_minoral = minoral + '/' + minoral
            double_majoral = majoral + '/' + majoral
            mix_major_minoral = majoral + '/' + minoral
            mix_minor_majoral = minoral + '/' + majoral
            # phased genotypes
            pdouble_minoral = minoral + '|' + minoral
            pdouble_majoral = majoral + '|' + majoral
            pmix_major_minoral = majoral + '|' + minoral
            pmix_minor_majoral = minoral + '|' + majoral
            if snpid not in converter_dict:
                converter_dict[snpid][double_minoral] = '2'
                converter_dict[snpid][double_majoral] = '0'
                converter_dict[snpid][mix_minor_majoral] = '1'
                converter_dict[snpid][mix_major_minoral] = '1'
                converter_dict[snpid][pdouble_minoral] = '2'
                converter_dict[snpid][pdouble_majoral] = '0'
                converter_dict[snpid][pmix_minor_majoral] = '1'
                converter_dict[snpid][pmix_major_minoral] = '1'


    return converter_dict


def vcfslicer(genovcf, blocksize) -> Generator:
    """When a dataset(VCF) is loaded, it is sliced in blocks of 1,000 rows (default size)"""
    with gzip.open(genovcf, 'r') as fh:
        try:
            fh.read(2)
            opener = gzip.open
        except gzip.BadGzipFile:
            opener = open
    # slice the dataset into blocks
    with opener(genovcf, 'rt') as vf:
        global HEADER
        HEADER = vf.readline()
        blocks_num = 0
        while True:
            blocks_num += 1
            block = list(islice(vf, blocksize))
            if block:
                yield blocks_num, block
            else:
                break


def vcfconverter(vcfblock) -> List:
    """convert variates with double minor alleles to 2, variates with double major alleles to 0,
    variates with mixing alleles to 1, and variates without alleles (eg. ./.) to NA based on the
    *.bim file corresponding to the VCF."""
    converter_dict = gconverter_dict
    txt_chunk = []
    txt_chunk_append = txt_chunk.append
    for item in vcfblock:
        chrom, pos, snpid, _ = item.split('\t', 3)
        firstmap, secondmap, thirdmap, lastmap, pfirstmap, psecondmap, pthirdmap, plastmap = converter_dict[snpid].items()
        txt = item.replace(*firstmap).replace(*secondmap).replace(*thirdmap).replace(*lastmap).\
            replace(*pfirstmap).replace(*psecondmap).replace(*pthirdmap).replace(*plastmap).replace('./.', 'NA')
        txt_chunk_append(txt)

    return txt_chunk


@stat_resource_usage
def main():
    """Run the script"""
    parser = get_parser()
    args = parser.parse_args()

    # check the python version
    pyversion_checker()
    # check the outdir
    outdir = outdir_checker(args.outdir)
    # parse *.bim file corresponding to the VCF
    global gconverter_dict
    gconverter_dict = parser_bim(args.genobim)
    ##########
    # keep all chunks
    alltxt_chunk = []
    alltxt_chunk_extend = alltxt_chunk.extend
    exc_txt = []  # for adding raised exceptions
    print('Start conversion operation for genotypes ...', flush=True)
    mp = Pool(processes=args.cpu)
    result_ls = []
    for block_num, block in vcfslicer(args.genovcf, args.blocksize):
        result = mp.apply_async(vcfconverter, args=(block,))
        result_ls.append((block_num, result))
        print("Info: Processing block {0} ...".format(block_num), flush=True)
    mp.close()
    mp.join()
    for block_num, res in result_ls:
        try:
            txt_chunk = res.get()
            alltxt_chunk_extend(txt_chunk)
        except Exception:
            exc_type, exc_info, _ = sys.exc_info()
            exc_txt.append('Warning: An unexpected Exception {0} was thrown: {1} when processing block {2}'
                           ''.format(exc_type, exc_info, block_num))
        else:
            print("Info: block {0} done".format(block_num), flush=True)
    if exc_txt:
        print('\n'.join(exc_txt))
    print("Conversion operation for genotypes finished ...", flush=True)
    # Output all txt chunks
    with open(join(outdir, args.outfile), 'wt') as out:
        out.write(HEADER)
        out.write(''.join(alltxt_chunk))


def get_parser():
    """Get options."""
    parser = argparse.ArgumentParser(description='Convert the genotypes in the VCF to numeric codes.',
                                     epilog='Note: Please see the https://github.com/bingochenbin/MRkit for details.')
    parser.add_argument('-d',
                        '--outdir',
                        metavar='DIRECTORY',
                        dest='outdir',
                        default='.',
                        action='store',
                        help='Output directory (default: %(default)s)')
    parser.add_argument('-o',
                        '--outfile',
                        metavar='FILENAME',
                        required=True,
                        action='store',
                        help='Output filename')
    parser.add_argument('-v',
                        '--vcf_input',
                        metavar='FILENAME',
                        dest='genovcf',
                        required=True,
                        action='store',
                        help='The SNP genotyping data file generated by bcftools, can be gzipped')
    parser.add_argument('-b',
                        '--bim_input',
                        metavar='FILENAME',
                        dest='genobim',
                        required=True,
                        action='store',
                        help='*.bim file corresponding to the VCF')
    parser.add_argument('-c',
                        '--blocksize',
                        metavar='INT',
                        dest='blocksize',
                        action='store',
                        default=1000,
                        type=int,
                        help='Integer, read file in slices of blocksize rows (default: %(default)s)')
    parser.add_argument('-n',
                        '--cpu',
                        metavar='INT',
                        dest='cpu',
                        action='store',
                        default=1,
                        type=int,
                        help='The number of cpu cores to use (default: %(default)s)')

    return parser


if __name__ == '__main__':
    main()
