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
    realignf = realignment_filter(arg.inFasta, arg.tumor_min_mismatch, arg.normal_max_mismatch, arg.window_size)
    realignf.filter(arg.targetTumorBam, arg.targetNormalBam, arg.outputDir, arg.targetMutationFile)
    logging.info( 'realignment filter end')


def run_indel_filter(arg):

    logging.info( 'indel filter start')
    indelf = indel_filter(arg.search_length, arg.min_depth, arg.min_mismatch)
    indelf.filter(arg.targetMutationFile, arg.targetBam, arg.outputDir)
    logging.info( 'indel filter end')


def run_breakpoint_filter(arg):

    logging.info( 'breakpoint filter start')
    bpf = breakpoint_filter(arg.max_depth, arg.min_clip_size, arg.junc_num_thres)
    bpf.filter(arg.targetMutationFile, arg.targetBam, arg.outputDir)
    logging.info( 'breakpoint filter end')


def run_simple_repeat_filter(arg):

    logging.info( 'simple repeat filter start')
    simplef = simple_repeat_filter(arg.simpleRepeatDB)
    simplef.filter(arg.targetMutationFile, arg.outputDir)
    logging.info( 'simple repeat filter end')


