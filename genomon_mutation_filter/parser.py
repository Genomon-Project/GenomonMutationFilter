#! /usr/bin/env python

"""
@author: ken0-1n
"""

import argparse
from .version import __version__
from .run import run_realignment_filter
from .run import run_indel_filter
from .run import run_breakpoint_filter
from .run import run_simple_repeat_filter


def create_parser():
    prog = "mutfilter"
    parser = argparse.ArgumentParser(prog = prog)
    parser.add_argument("--version", action = "version", version = prog + "-" + __version__)
    subparsers = parser.add_subparsers()


    def _create_realignment_parser(subparsers):

        realign_parser = subparsers.add_parser("realignment")
        realign_parser.add_argument( '-t', '--target_mutation_file', help = 'mutation text', type = str, default = None, required = True )
        realign_parser.add_argument( '-1', '--bam1', help = '1st bam file ( tumor )', type = str, default = None, required = True )
        realign_parser.add_argument( '-2', '--bam2', help = '2nd bam file ( control )', type = str, default = None)
        realign_parser.add_argument( '-A', '--sample1', help = '1st sample name ( disease )', type = str, default = None)
        realign_parser.add_argument( '-B', '--sample2', help = '2nd sample name ( control )', type = str, default = None)
        realign_parser.add_argument( '-o', '--output', help = 'Output text file', type = str, default = None, required = True)
        realign_parser.add_argument( '-r', '--ref_genome', help = 'Reference genome', type = str, default = None , required = True)
        realign_parser.add_argument( '--blat', action='store_true', default=False)
        realign_parser.add_argument( '-b', '--blat_path', type = str, default = None)
        realign_parser.add_argument( '-m', '--tumor_min_mismatch', metavar = "tumor_min_mismatch", default='0', type=int)
        realign_parser.add_argument( '-M', '--normal_max_mismatch', metavar = "normal_max_mismatch", default='100000', type=int)
        realign_parser.add_argument( '-s', '--score_difference', metavar = "score_difference", default='5', type=int)
        realign_parser.add_argument( '-w', '--window_size', metavar = "window_size", default='200', type=int)
        realign_parser.add_argument( '-d', '--max_depth', metavar = "max_depth", default='5000', type=int)
        realign_parser.add_argument( '-F', '--exclude_sam_flags', metavar = "exclude_sam_flags", default='3328', type=int)
        realign_parser.add_argument( '-O', '--print_format', choices = ['vcf','anno'], help = 'Print VCF or anno(TSV) format',  default = 'anno' )
        realign_parser.add_argument( '--header', action="store_true", default=False,  dest='header_flag')
        realign_parser.add_argument( '-T', '--thread_num', metavar = "number_of_threads", default='1', type=int)

        return realign_parser
        

    def _create_indel_parser(subparsers):

        indel_parser = subparsers.add_parser("indel")
        indel_parser.add_argument( '-t', '--target_mutation_file', help = 'mutation text', type = str, default = None, required = True )
        indel_parser.add_argument( '-2', '--bam2', help = 'normal bam file', type = str, default = None, required = True)
        indel_parser.add_argument( '-A', '--sample1', help = '1st sample name ( disease )', type = str, default = None)
        indel_parser.add_argument( '-B', '--sample2', help = '2nd sample name ( control )', type = str, default = None)
        indel_parser.add_argument( '-o', '--output', help = 'Output text file', type = str, default = None, required = True)
        indel_parser.add_argument( '-l', '--search_length', metavar = "search_length", default='40', type=int)
        indel_parser.add_argument( '-n', '--neighbor', metavar = "neighbor", default='5', type=int)
        indel_parser.add_argument( '-d', '--min_depth', metavar = "min_depth", default='8', type=int)
        indel_parser.add_argument( '-m', '--min_mismatch', metavar = "min_mismatch", default='100000', type=int)
        indel_parser.add_argument( '-a', '--af_thres', metavar = "allele_frequency_thres", default='1', type=float)
        indel_parser.add_argument( '--header', action="store_true", default=False,  dest='header_flag')
        indel_parser.add_argument( '-s', '--samtools_path', type = str, default = None, required = True)
        indel_parser.add_argument( '-S', '--samtools_params', type = str, default = "-q 20 -BQ0 -d 10000000")
        indel_parser.add_argument( '-O', '--print_format', choices = ['vcf','anno'], help = 'Print VCF or anno(TSV) format',  default = 'anno' )
        indel_parser.add_argument( '-r', '--ref_genome', help = 'Reference genome', type = str, default = None)

        return indel_parser


    def _create_breakpoint_parser(subparsers):
        
        breakpoint_parser = subparsers.add_parser("breakpoint")
        breakpoint_parser.add_argument( '-t', '--target_mutation_file', help = 'mutation text', type = str, default = None, required = True )
        breakpoint_parser.add_argument( '-2', '--bam2', help = 'normal bam file', type = str, default = None, required = True)
        breakpoint_parser.add_argument( '-A', '--sample1', help = '1st sample name ( disease )', type = str, default = None)
        breakpoint_parser.add_argument( '-B', '--sample2', help = '2nd sample name ( control )', type = str, default = None)
        breakpoint_parser.add_argument( '-o', '--output', help = 'Output text file', type = str, default = None, required = True)
        breakpoint_parser.add_argument( '-d', '--max_depth', metavar = "max_depth", default='1000', type=int)
        breakpoint_parser.add_argument( '-c', '--min_clip_size', metavar = "min_clip_size", default='20', type=int)
        breakpoint_parser.add_argument( '-j', '--junc_num_thres', metavar = "junc_num_thres", default='0', type=int)
        breakpoint_parser.add_argument( '-m', '--mapq_thres', metavar = "mapping_quality_thres", default='10', type=int)
        breakpoint_parser.add_argument( '-F', '--exclude_sam_flags', metavar = "exclude_sam_flags", default='3332', type=int)
        breakpoint_parser.add_argument( '--header', action="store_true", default=False,  dest='header_flag')
        breakpoint_parser.add_argument( '-O', '--print_format', choices = ['vcf','anno'], help = 'Print VCF or anno(TSV) format',  default = 'anno' )
        breakpoint_parser.add_argument( '-r', '--ref_genome', help = 'Reference genome', type = str, default = None)

        return breakpoint_parser


    def _create_simplerepeat_parser(subparsers):
        
        simplerepeat_parser = subparsers.add_parser("simplerepeat")
        simplerepeat_parser.add_argument( '-t', '--target_mutation_file', help = 'mutation text', type = str, default = None, required = True )
        simplerepeat_parser.add_argument( '-o', '--output', help = 'Output text file', type = str, default = None, required = True)
        simplerepeat_parser.add_argument( '-S', '--simple_repeat_db', help = 'simple_repeat_database', type = str, default = None, required = True)
        simplerepeat_parser.add_argument('--header', action="store_true", default=False,  dest='header_flag')
        simplerepeat_parser.add_argument( '-O', '--print_format', choices = ['vcf','anno'], help = 'Print VCF or anno(TSV) format',  default = 'anno' )
        
        return simplerepeat_parser


    realign_parser = _create_realignment_parser(subparsers)
    realign_parser.set_defaults(func = run_realignment_filter)
    indel_parser = _create_indel_parser(subparsers)
    indel_parser.set_defaults(func = run_indel_filter)
    breakpoint_parser = _create_breakpoint_parser(subparsers)
    breakpoint_parser.set_defaults(func = run_breakpoint_filter)
    simplerepeat_parser = _create_simplerepeat_parser(subparsers)
    simplerepeat_parser.set_defaults(func = run_simple_repeat_filter)
    return parser
    
