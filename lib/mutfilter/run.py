#! /usr/bin/env python

import sys
import os
import math
import argparse
import logging
from realignment_filter import realignment_filter
from indel_filter import indel_filter
from breakpoint_filter import breakpoint_filter
from simple_repeat_filter import simple_repeat_filter


#
# Main
#
def run_realignment_filter(arg):

    logging.info( 'realignment filter start')
    realignf = realignment_filter(arg.ref_genome, arg.tumor_min_mismatch, arg.normal_max_mismatch, arg.window_size, arg.score_difference, arg.blat_path, arg.header_flag, arg.max_depth, arg.exclude_sam_flags)
    realignf.filter(arg.bam1, arg.bam2, arg.output, arg.target_mutation_file)
    logging.info( 'realignment filter end')


def run_indel_filter(arg):

    logging.info( 'indel filter start')
    indelf = indel_filter(arg.search_length, arg.min_depth, arg.min_mismatch, arg.af_thres, arg.neighbor, arg.header_flag, arg.samtools_path, arg.samtools_params)
    indelf.filter(arg.target_mutation_file, arg.bam2, arg.output)
    logging.info( 'indel filter end')


def run_breakpoint_filter(arg):

    logging.info( 'breakpoint filter start')
    bpf = breakpoint_filter(arg.max_depth, arg.min_clip_size, arg.junc_num_thres, arg.mapq_thres, arg.header_flag, arg.exclude_sam_flags)
    bpf.filter(arg.target_mutation_file, arg.bam2, arg.output)
    logging.info( 'breakpoint filter end')


def run_simple_repeat_filter(arg):

    logging.info( 'simple repeat filter start')
    simplef = simple_repeat_filter(arg.simple_repeat_db, arg.header_flag)
    simplef.filter(arg.target_mutation_file, arg.output)
    logging.info( 'simple repeat filter end')


