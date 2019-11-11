#! /usr/bin/env python

import sys
import unittest
import os, tempfile, shutil, filecmp
import subprocess
from genomon_mutation_filter import breakpoint_filter as breakf
import genomon_mutation_filter

class TestBreakpoint(unittest.TestCase):

    ######################################
    # Tumor/Normal Pair, Annoformat
    ######################################
    def test1_1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        max_depth = 1000
        min_clip_size = 20
        junc_num_thres = 0
        mapq_thres = 10
        header_flag = True
        exclude_sam_flags = 3332
        ref_genome = cur_dir + "/../data/GRCh37.fa"
        
        bpf = breakf.Breakpoint_filter(max_depth, min_clip_size, junc_num_thres, 
        mapq_thres, header_flag, exclude_sam_flags, ref_genome)
        
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_breakpoint_result_test1_1.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test1.txt"
        bpf.filter(target_mutation_file, bam2, output)
        
        answer_file = cur_dir + "/../data/5929_small_breakpoint_result_answer_test1_1.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))


    def test1_2(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        max_depth = 1000
        min_clip_size = 20
        junc_num_thres = 0
        mapq_thres = 10
        header_flag = True
        exclude_sam_flags = 3332
        ref_genome = cur_dir + "/../data/GRCh37.fa"
        
        bpf = breakf.Breakpoint_filter(max_depth, min_clip_size, junc_num_thres, 
        mapq_thres, header_flag, exclude_sam_flags, ref_genome)
        
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_breakpoint_result_test1_2.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test2.txt"
        bpf.filter(target_mutation_file, bam2, output)
        
        answer_file = cur_dir + "/../data/5929_small_breakpoint_result_answer_test1_2.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        

    def test1_3(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        max_depth = 1000
        min_clip_size = 20
        junc_num_thres = 0
        mapq_thres = 10
        header_flag = True
        exclude_sam_flags = 3332
        ref_genome = cur_dir + "/../data/GRCh37.fa"
        
        bpf = breakf.Breakpoint_filter(max_depth, min_clip_size, junc_num_thres, 
        mapq_thres, header_flag, exclude_sam_flags, ref_genome)
        
        bam2 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_breakpoint_result_test1_3.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test2.txt"
        bpf.filter(target_mutation_file, bam2, output)
        
        answer_file = cur_dir + "/../data/5929_small_breakpoint_result_answer_test1_3.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        
    def test1_4(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        max_depth = 1000
        min_clip_size = 20
        junc_num_thres = 0
        mapq_thres = 10
        header_flag = False
        exclude_sam_flags = 3332
        ref_genome = cur_dir + "/../data/GRCh37.fa"
        
        bpf = breakf.Breakpoint_filter(max_depth, min_clip_size, junc_num_thres, 
        mapq_thres, header_flag, exclude_sam_flags, ref_genome)
        
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_breakpoint_result_test1_4.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test4.txt"
        bpf.filter(target_mutation_file, bam2, output)
        
        answer_file = cur_dir + "/../data/5929_small_breakpoint_result_answer_test1_4.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        

    ######################################
    # Tumor/Normal Pair, VCF format
    ######################################
    def test2_1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        max_depth = 1000
        min_clip_size = 20
        junc_num_thres = 0
        mapq_thres = 10
        header_flag = True
        exclude_sam_flags = 3332
        ref_genome = cur_dir + "/../data/GRCh37.fa"
        
        bpf = breakf.Breakpoint_filter(max_depth, min_clip_size, junc_num_thres, 
        mapq_thres, header_flag, exclude_sam_flags, ref_genome)
        
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_breakpoint_result_test2_1.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test21.txt"
        tumor_sample = "5929_tumor"
        normal_sample = "5929_control"
        bpf.filter_vcf(target_mutation_file, bam2, output, tumor_sample, normal_sample)
        
        answer_file = cur_dir + "/../data/5929_small_breakpoint_result_answer_test2_1.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))


    def test2_2(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        max_depth = 1000
        min_clip_size = 20
        junc_num_thres = 0
        mapq_thres = 10
        header_flag = True
        exclude_sam_flags = 3332
        ref_genome = cur_dir + "/../data/GRCh37.fa"
        
        bpf = breakf.Breakpoint_filter(max_depth, min_clip_size, junc_num_thres, 
        mapq_thres, header_flag, exclude_sam_flags, ref_genome)
        
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_breakpoint_result_test2_2.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test22.txt"
        tumor_sample = "5929_tumor"
        normal_sample = "5929_control"
        bpf.filter_vcf(target_mutation_file, bam2, output, tumor_sample, normal_sample)
        
        answer_file = cur_dir + "/../data/5929_small_breakpoint_result_answer_test2_2.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        

    def test2_3(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        max_depth = 1000
        min_clip_size = 20
        junc_num_thres = 0
        mapq_thres = 10
        header_flag = True
        exclude_sam_flags = 3332
        ref_genome = cur_dir + "/../data/GRCh37.fa"
        
        bpf = breakf.Breakpoint_filter(max_depth, min_clip_size, junc_num_thres, 
        mapq_thres, header_flag, exclude_sam_flags, ref_genome)
        
        bam2 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_breakpoint_result_test2_3.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test22.txt"
        tumor_sample = "5929_tumor"
        normal_sample = "5929_control"
        bpf.filter_vcf(target_mutation_file, bam2, output, tumor_sample, normal_sample)
        
        answer_file = cur_dir + "/../data/5929_small_breakpoint_result_answer_test2_3.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))

    def test2_4(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        max_depth = 1000
        min_clip_size = 20
        junc_num_thres = 0
        mapq_thres = 10
        header_flag = False
        exclude_sam_flags = 3332
        ref_genome = cur_dir + "/../data/GRCh37.fa"
        
        bpf = breakf.Breakpoint_filter(max_depth, min_clip_size, junc_num_thres, 
        mapq_thres, header_flag, exclude_sam_flags, ref_genome)
        
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_breakpoint_result_test2_4.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test24.txt"
        tumor_sample = "5929_tumor"
        normal_sample = "5929_control"

        with self.assertRaises(ValueError) as er:
            bpf.filter_vcf(target_mutation_file, bam2, output, tumor_sample, normal_sample)
        
        
    def test3_1(self):
        
        cmd = ['mutfilter', '--version']
        subprocess.check_call(cmd)
        
        
    def test3_2(self):
        
        cmd = ['mutfilter', 'breakpoint', '--help']
        subprocess.check_call(cmd)


    def test3_3(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        
        max_depth = '1000'
        min_clip_size = '20'
        junc_num_thres = '0'
        mapq_thres = '10'
        exclude_sam_flags = '3332'
        ref_genome = cur_dir + "/../data/GRCh37.fa"
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_breakpoint_result_test3_1.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test1.txt"

        cmd = ['mutfilter', 'breakpoint',
            '-t', target_mutation_file,
            '-2', bam2,
            '-o', output,
            '-r', ref_genome,
            '-d', max_depth,
            '-c', min_clip_size,
            '-j', junc_num_thres,
            '-m', mapq_thres,
            '-F', exclude_sam_flags,
            '--header']
        
        subprocess.check_call(cmd)

        answer_file = cur_dir + "/../data/5929_small_breakpoint_result_answer_test1_1.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))


    def test4_1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        
        max_depth = '1000'
        min_clip_size = '20'
        junc_num_thres = '0'
        mapq_thres = '10'
        exclude_sam_flags = '3332'
        ref_genome = cur_dir + "/../data/GRCh37.fa"
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_breakpoint_result_test4_1.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test1.txt"

        cmd = ['breakpoint',
            '-t', target_mutation_file,
            '-2', bam2,
            '-o', output,
            '-r', ref_genome,
            '-d', max_depth,
            '-c', min_clip_size,
            '-j', junc_num_thres,
            '-m', mapq_thres,
            '-F', exclude_sam_flags,
            '--header']
        
        parser = genomon_mutation_filter.parser.create_parser()
        args = parser.parse_args(cmd)
        args.func(args)

        answer_file = cur_dir + "/../data/5929_small_breakpoint_result_answer_test1_1.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))

if __name__ == "__main__":
    unittest.main()

