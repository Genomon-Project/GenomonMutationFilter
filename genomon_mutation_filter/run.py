#! /usr/bin/env python

import sys
import os
import math
import argparse
import logging
from . import realignment_filter as rf
from . import indel_filter as idf
from . import breakpoint_filter as brf
from . import simple_repeat_filter as sif
from . import oxog_filter as of
from . import position_filter as pf


#
# Main
#
def run_realignment_filter(arg):

    is_anno = True if arg.print_format == 'anno' else False

    logging.info( 'realignment filter start')
    realignf = rf.Realignment_filter(arg.ref_genome, arg.tumor_min_mismatch, arg.normal_max_mismatch, arg.window_size, arg.score_difference, arg.blat_path, arg.header_flag, arg.max_depth, arg.exclude_sam_flags, arg.thread_num, arg.blat, arg.edlib)
    if is_anno == True:
        realignf.filter(arg.bam1, arg.bam2, arg.output, arg.target_mutation_file)
    else:
        realignf.filter_vcf(arg.bam1, arg.bam2, arg.output, arg.target_mutation_file, arg.sample1, arg.sample2)
    logging.info( 'realignment filter end')


def run_indel_filter(arg):

    is_anno = True if arg.print_format == 'anno' else False

    logging.info( 'indel filter start')
    indelf = idf.Indel_filter(arg.search_length, arg.min_depth, arg.min_mismatch, arg.af_thres, arg.neighbor, arg.header_flag, arg.samtools_path, arg.samtools_params, arg.ref_genome, arg.thread_num)
    if is_anno == True:
        indelf.filter(arg.target_mutation_file, arg.bam2, arg.output)
    else:
        indelf.filter_vcf(arg.target_mutation_file, arg.bam2, arg.output, arg.sample1, arg.sample2)
    logging.info( 'indel filter end')


def run_breakpoint_filter(arg):

    is_anno = True if arg.print_format == 'anno' else False

    logging.info( 'breakpoint filter start')
    bpf = brf.Breakpoint_filter(arg.max_depth, arg.min_clip_size, arg.junc_num_thres, arg.mapq_thres, arg.header_flag, arg.exclude_sam_flags, arg.ref_genome)
    if is_anno == True:
        bpf.filter(arg.target_mutation_file, arg.bam2, arg.output)
    else:
        bpf.filter_vcf(arg.target_mutation_file, arg.bam2, arg.output, arg.sample1, arg.sample2)
    logging.info( 'breakpoint filter end')


def run_simple_repeat_filter(arg):

    is_anno = True if arg.print_format == 'anno' else False

    logging.info( 'simple repeat filter start')
    simplef = sif.Simple_repeat_filter(arg.simple_repeat_db, arg.header_flag)
    if is_anno == True:
        simplef.filter(arg.target_mutation_file, arg.output)
    else:
        simplef.filter_vcf(arg.target_mutation_file, arg.output)
    logging.info( 'simple repeat filter end')


def run_oxog_filter(arg):

    is_tsv = True if arg.print_format == 'tsv' else False

    logging.info( 'oxog filter start')
    oxogf = of.Oxog_filter(arg.ref_genome,arg.samtools_path, arg.mpileup_params)
    if is_tsv == True:
        oxogf.filter(arg.target_mutation_file, arg.bam1, arg.output, arg.bam3)
    else:
        oxogf.filter_vcf(arg.target_mutation_file, arg.bam1, arg.output, arg.bam3)
    logging.info( 'oxog filter end')


def run_position_filter(arg):

    is_tsv = True if arg.print_format == 'tsv' else False

    logging.info( 'position filter start')
    posf = pf.Position_filter(arg.ref_genome, arg.samtools_path, arg.mpileup_params)
    if is_tsv == True:
        posf.filter(arg.target_mutation_file, arg.bam1, arg.output)
    else:
        posf.filter_vcf(arg.target_mutation_file, arg.bam1, arg.output)
    logging.info( 'position filter end')


