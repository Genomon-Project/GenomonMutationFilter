#! /usr/bin/env python

import sys
import unittest
import os, tempfile, shutil, filecmp
import subprocess
# from lib.mutfilter import realignment_filter as rf
from genomon_mutation_filter import realignment_filter as rf
import genomon_mutation_filter

class TestRealignment(unittest.TestCase):

    def setUp(self):
        self.blat = '/home/ubuntu/environment/tools/blat_linux.x86_64.v369/blat'

    ######################################
    # Tumor/Normal Pair, Annoformat
    ######################################
    def test1_1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        ref_genome = cur_dir + "/../data/GRCh37.fa"
        tumor_min_mismatch = 0
        normal_max_mismatch = 100000
        window_size = 500
        score_difference = 5
        blat_path = self.blat
        header_flag = True
        max_depth = 5000
        exclude_sam_flags = 3328
        thread_num = 1
        
        realignf = rf.Realignment_filter(ref_genome, 
        tumor_min_mismatch, normal_max_mismatch, window_size, score_difference, 
        blat_path, header_flag, max_depth, exclude_sam_flags, thread_num)

        bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_realignment_result_test1_1.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test1.txt"
        realignf.filter(bam1, bam2, output, target_mutation_file)
        
        answer_file = cur_dir + "/../data/5929_small_realignment_result_answer_test1_1.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        
    def test1_2(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        ref_genome = cur_dir + "/../data/GRCh37.fa"
        tumor_min_mismatch = 0
        normal_max_mismatch = 100000
        window_size = 200
        score_difference = 5
        blat_path = self.blat
        header_flag = True
        max_depth = 5000
        exclude_sam_flags = 3328
        thread_num = 1
        
        realignf = rf.Realignment_filter(ref_genome, 
        tumor_min_mismatch, normal_max_mismatch, window_size, score_difference, 
        blat_path, header_flag, max_depth, exclude_sam_flags, thread_num)

        bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_realignment_result_test1_2.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test2.txt"
        realignf.filter(bam1, bam2, output, target_mutation_file)
        
        answer_file = cur_dir + "/../data/5929_small_realignment_result_answer_test1_2.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        
    def test1_3(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        ref_genome = cur_dir + "/../data/GRCh37.fa"
        tumor_min_mismatch = 0
        normal_max_mismatch = 100000
        window_size = 200
        score_difference = 5
        blat_path = self.blat
        header_flag = True
        max_depth = 5000
        exclude_sam_flags = 3328
        thread_num = 2
        
        realignf = rf.Realignment_filter(ref_genome, 
        tumor_min_mismatch, normal_max_mismatch, window_size, score_difference, 
        blat_path, header_flag, max_depth, exclude_sam_flags, thread_num)

        bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_realignment_result_test1_3.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test2.txt"
        realignf.filter(bam1, bam2, output, target_mutation_file)
        
        answer_file = cur_dir + "/../data/5929_small_realignment_result_answer_test1_2.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        
    def test1_4(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        ref_genome = cur_dir + "/../data/GRCh37.fa"
        tumor_min_mismatch = 0
        normal_max_mismatch = 100000
        window_size = 200
        score_difference = 5
        blat_path = self.blat
        header_flag = False
        max_depth = 5000
        exclude_sam_flags = 3328
        thread_num = 1
        
        realignf = rf.Realignment_filter(ref_genome, 
        tumor_min_mismatch, normal_max_mismatch, window_size, score_difference, 
        blat_path, header_flag, max_depth, exclude_sam_flags, thread_num)

        bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_realignment_result_test1_4.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test4.txt"
        realignf.filter(bam1, bam2, output, target_mutation_file)
        
        answer_file = cur_dir + "/../data/5929_small_realignment_result_answer_test1_4.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        
    def test1_5(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        ref_genome = cur_dir + "/../data/GRCh37.fa"
        tumor_min_mismatch = 0
        normal_max_mismatch = 100000
        window_size = 200
        score_difference = 5
        blat_path = self.blat
        header_flag = False
        max_depth = 5000
        exclude_sam_flags = 3328
        thread_num = 2
        
        realignf = rf.Realignment_filter(ref_genome, 
        tumor_min_mismatch, normal_max_mismatch, window_size, score_difference, 
        blat_path, header_flag, max_depth, exclude_sam_flags, thread_num)

        bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_realignment_result_test1_5.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test4.txt"
        realignf.filter(bam1, bam2, output, target_mutation_file)
        
        answer_file = cur_dir + "/../data/5929_small_realignment_result_answer_test1_4.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        
        
    def test1_6(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        ref_genome = cur_dir + "/../data/GRCh37.fa"
        tumor_min_mismatch = 0
        normal_max_mismatch = 100000
        window_size = 200
        score_difference = 10
        blat_path = self.blat
        header_flag = True
        max_depth = 5000
        exclude_sam_flags = 3328
        thread_num = 2
        
        realignf = rf.Realignment_filter(ref_genome, 
        tumor_min_mismatch, normal_max_mismatch, window_size, score_difference, 
        blat_path, header_flag, max_depth, exclude_sam_flags, thread_num)

        bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_realignment_result_test1_6.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test2.txt"
        realignf.filter(bam1, bam2, output, target_mutation_file)
        
        answer_file = cur_dir + "/../data/5929_small_realignment_result_answer_test1_6.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))


    ######################################
    # Tumor Single, Annoformat
    ######################################
    def test2_1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        ref_genome = cur_dir + "/../data/GRCh37.fa"
        tumor_min_mismatch = 0
        normal_max_mismatch = 100000
        window_size = 500
        score_difference = 5
        blat_path = self.blat
        header_flag = True
        max_depth = 5000
        exclude_sam_flags = 3328
        thread_num = 1
        
        realignf = rf.Realignment_filter(ref_genome, 
        tumor_min_mismatch, normal_max_mismatch, window_size, score_difference, 
        blat_path, header_flag, max_depth, exclude_sam_flags, thread_num)

        bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        bam2 = None
        output = cur_dir + "/../data/5929_small_realignment_result_test2_1.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test11.txt"
        realignf.filter(bam1, bam2, output, target_mutation_file)
        
        answer_file = cur_dir + "/../data/5929_small_realignment_result_answer_test2_1.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        
    def test2_2(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        ref_genome = cur_dir + "/../data/GRCh37.fa"
        tumor_min_mismatch = 0
        normal_max_mismatch = 100000
        window_size = 200
        score_difference = 5
        blat_path = self.blat
        header_flag = True
        max_depth = 5000
        exclude_sam_flags = 3328
        thread_num = 1
        
        realignf = rf.Realignment_filter(ref_genome, 
        tumor_min_mismatch, normal_max_mismatch, window_size, score_difference, 
        blat_path, header_flag, max_depth, exclude_sam_flags, thread_num)

        bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        bam2 = None
        output = cur_dir + "/../data/5929_small_realignment_result_test2_2.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test12.txt"
        realignf.filter(bam1, bam2, output, target_mutation_file)
        
        answer_file = cur_dir + "/../data/5929_small_realignment_result_answer_test2_2.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        
    def test2_3(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        ref_genome = cur_dir + "/../data/GRCh37.fa"
        tumor_min_mismatch = 0
        normal_max_mismatch = 100000
        window_size = 200
        score_difference = 5
        blat_path = self.blat
        header_flag = True
        max_depth = 5000
        exclude_sam_flags = 3328
        thread_num = 2
        
        realignf = rf.Realignment_filter(ref_genome, 
        tumor_min_mismatch, normal_max_mismatch, window_size, score_difference, 
        blat_path, header_flag, max_depth, exclude_sam_flags, thread_num)

        bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        bam2 = None
        output = cur_dir + "/../data/5929_small_realignment_result_test2_3.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test12.txt"
        realignf.filter(bam1, bam2, output, target_mutation_file)
        
        answer_file = cur_dir + "/../data/5929_small_realignment_result_answer_test2_2.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))

    def test2_4(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        ref_genome = cur_dir + "/../data/GRCh37.fa"
        tumor_min_mismatch = 0
        normal_max_mismatch = 100000
        window_size = 200
        score_difference = 5
        blat_path = self.blat
        header_flag = False
        max_depth = 5000
        exclude_sam_flags = 3328
        thread_num = 1
        
        realignf = rf.Realignment_filter(ref_genome, 
        tumor_min_mismatch, normal_max_mismatch, window_size, score_difference, 
        blat_path, header_flag, max_depth, exclude_sam_flags, thread_num)

        bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        bam2 = None
        output = cur_dir + "/../data/5929_small_realignment_result_test2_4.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test14.txt"
        realignf.filter(bam1, bam2, output, target_mutation_file)
        
        answer_file = cur_dir + "/../data/5929_small_realignment_result_answer_test2_4.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        
    def test2_5(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        ref_genome = cur_dir + "/../data/GRCh37.fa"
        tumor_min_mismatch = 0
        normal_max_mismatch = 100000
        window_size = 200
        score_difference = 5
        blat_path = self.blat
        header_flag = False
        max_depth = 5000
        exclude_sam_flags = 3328
        thread_num = 2
        
        realignf = rf.Realignment_filter(ref_genome, 
        tumor_min_mismatch, normal_max_mismatch, window_size, score_difference, 
        blat_path, header_flag, max_depth, exclude_sam_flags, thread_num)

        bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        bam2 = None
        output = cur_dir + "/../data/5929_small_realignment_result_test2_5.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test14.txt"
        realignf.filter(bam1, bam2, output, target_mutation_file)
        
        answer_file = cur_dir + "/../data/5929_small_realignment_result_answer_test2_4.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))

    def test2_6(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        ref_genome = cur_dir + "/../data/GRCh37.fa"
        tumor_min_mismatch = 0
        normal_max_mismatch = 100000
        window_size = 200
        score_difference = 10
        blat_path = self.blat
        header_flag = True
        max_depth = 5000
        exclude_sam_flags = 3328
        thread_num = 2
        
        realignf = rf.Realignment_filter(ref_genome, 
        tumor_min_mismatch, normal_max_mismatch, window_size, score_difference, 
        blat_path, header_flag, max_depth, exclude_sam_flags, thread_num)

        bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        bam2 = None
        output = cur_dir + "/../data/5929_small_realignment_result_test2_6.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test12.txt"
        realignf.filter(bam1, bam2, output, target_mutation_file)
        
        answer_file = cur_dir + "/../data/5929_small_realignment_result_answer_test2_6.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))

    ######################################
    # Tumor/Normal Pair, VCF 
    ######################################
    def test3_1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        ref_genome = cur_dir + "/../data/GRCh37.fa"
        tumor_min_mismatch = 0
        normal_max_mismatch = 100000
        window_size = 500
        score_difference = 5
        blat_path = self.blat
        header_flag = True
        max_depth = 5000
        exclude_sam_flags = 3328
        thread_num = 1
        
        realignf = rf.Realignment_filter(ref_genome, 
        tumor_min_mismatch, normal_max_mismatch, window_size, score_difference, 
        blat_path, header_flag, max_depth, exclude_sam_flags, thread_num)

        bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_realignment_result_test3_1.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test21.txt"
        tumor_sample = "5929_tumor"
        normal_sample = "5929_control"
        realignf.filter_vcf(bam1, bam2, output, target_mutation_file, tumor_sample, normal_sample)
        
        answer_file = cur_dir + "/../data/5929_small_realignment_result_answer_test3_1.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))

    def test3_2(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        ref_genome = cur_dir + "/../data/GRCh37.fa"
        tumor_min_mismatch = 0
        normal_max_mismatch = 100000
        window_size = 500
        score_difference = 5
        blat_path = self.blat
        header_flag = True
        max_depth = 5000
        exclude_sam_flags = 3328
        thread_num = 1
        
        realignf = rf.Realignment_filter(ref_genome, 
        tumor_min_mismatch, normal_max_mismatch, window_size, score_difference, 
        blat_path, header_flag, max_depth, exclude_sam_flags, thread_num)

        bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_realignment_result_test3_2.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test22.txt"
        tumor_sample = "5929_tumor"
        normal_sample = "5929_control"
        realignf.filter_vcf(bam1, bam2, output, target_mutation_file, tumor_sample, normal_sample)
        
        answer_file = cur_dir + "/../data/5929_small_realignment_result_answer_test3_2.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        
    def test3_3(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        ref_genome = cur_dir + "/../data/GRCh37.fa"
        tumor_min_mismatch = 0
        normal_max_mismatch = 100000
        window_size = 500
        score_difference = 5
        blat_path = self.blat
        header_flag = True
        max_depth = 5000
        exclude_sam_flags = 3328
        thread_num = 2
        
        realignf = rf.Realignment_filter(ref_genome, 
        tumor_min_mismatch, normal_max_mismatch, window_size, score_difference, 
        blat_path, header_flag, max_depth, exclude_sam_flags, thread_num)

        bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_realignment_result_test3_3.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test22.txt"
        tumor_sample = "5929_tumor"
        normal_sample = "5929_control"
        realignf.filter_vcf(bam1, bam2, output, target_mutation_file, tumor_sample, normal_sample)
        
        answer_file = cur_dir + "/../data/5929_small_realignment_result_answer_test3_2.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))

    def test3_4(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        ref_genome = cur_dir + "/../data/GRCh37.fa"
        tumor_min_mismatch = 0
        normal_max_mismatch = 100000
        window_size = 500
        score_difference = 5
        blat_path = self.blat
        header_flag = False
        max_depth = 5000
        exclude_sam_flags = 3328
        thread_num = 1
        
        realignf = rf.Realignment_filter(ref_genome, 
        tumor_min_mismatch, normal_max_mismatch, window_size, score_difference, 
        blat_path, header_flag, max_depth, exclude_sam_flags, thread_num)

        bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_realignment_result_test3_4.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test24.txt"
        tumor_sample = "5929_tumor"
        normal_sample = "5929_control"
        
        with self.assertRaises(ValueError) as er:
            realignf.filter_vcf(bam1, bam2, output, target_mutation_file, tumor_sample, normal_sample)
        # self.assertTrue(er.exception.message.endswith('is not in list'))

    def test3_5(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        ref_genome = cur_dir + "/../data/GRCh37.fa"
        tumor_min_mismatch = 0
        normal_max_mismatch = 100000
        window_size = 500
        score_difference = 5
        blat_path = self.blat
        header_flag = False
        max_depth = 5000
        exclude_sam_flags = 3328
        thread_num = 2
        
        realignf = rf.Realignment_filter(ref_genome, 
        tumor_min_mismatch, normal_max_mismatch, window_size, score_difference, 
        blat_path, header_flag, max_depth, exclude_sam_flags, thread_num)

        bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_realignment_result_test3_5.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test24.txt"
        tumor_sample = "5929_tumor"
        normal_sample = "5929_control"
        
        with self.assertRaises(RuntimeError) as er:
            realignf.filter_vcf(bam1, bam2, output, target_mutation_file, tumor_sample, normal_sample)
        # self.assertTrue(er.exception.message.endswith('There was an error!'))


    def test3_6(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        ref_genome = cur_dir + "/../data/GRCh37.fa"
        tumor_min_mismatch = 0
        normal_max_mismatch = 100000
        window_size = 200
        score_difference = 10
        blat_path = self.blat
        header_flag = True
        max_depth = 5000
        exclude_sam_flags = 3328
        thread_num = 2
        
        realignf = rf.Realignment_filter(ref_genome, 
        tumor_min_mismatch, normal_max_mismatch, window_size, score_difference, 
        blat_path, header_flag, max_depth, exclude_sam_flags, thread_num)

        bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_realignment_result_test3_6.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test22.txt"
        tumor_sample = "5929_tumor"
        normal_sample = "5929_control"
        realignf.filter_vcf(bam1, bam2, output, target_mutation_file, tumor_sample, normal_sample)

        answer_file = cur_dir + "/../data/5929_small_realignment_result_answer_test3_6.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))

    ######################################
    # Tumor Single, VCF
    ######################################
    def test4_1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        ref_genome = cur_dir + "/../data/GRCh37.fa"
        tumor_min_mismatch = 0
        normal_max_mismatch = 100000
        window_size = 500
        score_difference = 5
        blat_path = self.blat
        header_flag = True
        max_depth = 5000
        exclude_sam_flags = 3328
        thread_num = 1
        
        realignf = rf.Realignment_filter(ref_genome, 
        tumor_min_mismatch, normal_max_mismatch, window_size, score_difference, 
        blat_path, header_flag, max_depth, exclude_sam_flags, thread_num)

        bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        bam2 = None
        output = cur_dir + "/../data/5929_small_realignment_result_test4_1.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test31.txt"
        tumor_sample = "5929_tumor"
        normal_sample = None
        realignf.filter_vcf(bam1, bam2, output, target_mutation_file, tumor_sample, normal_sample)
        
        answer_file = cur_dir + "/../data/5929_small_realignment_result_answer_test4_1.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))

    def test4_2(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        ref_genome = cur_dir + "/../data/GRCh37.fa"
        tumor_min_mismatch = 0
        normal_max_mismatch = 100000
        window_size = 200
        score_difference = 5
        blat_path = self.blat
        header_flag = True
        max_depth = 5000
        exclude_sam_flags = 3328
        thread_num = 1
        
        realignf = rf.Realignment_filter(ref_genome, 
        tumor_min_mismatch, normal_max_mismatch, window_size, score_difference, 
        blat_path, header_flag, max_depth, exclude_sam_flags, thread_num)

        bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        bam2 = None
        output = cur_dir + "/../data/5929_small_realignment_result_test4_2.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test32.txt"
        tumor_sample = "5929_tumor"
        normal_sample = None
        realignf.filter_vcf(bam1, bam2, output, target_mutation_file, tumor_sample, normal_sample)
        
        answer_file = cur_dir + "/../data/5929_small_realignment_result_answer_test4_2.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        
    def test4_3(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        ref_genome = cur_dir + "/../data/GRCh37.fa"
        tumor_min_mismatch = 0
        normal_max_mismatch = 100000
        window_size = 200
        score_difference = 5
        blat_path = self.blat
        header_flag = True
        max_depth = 5000
        exclude_sam_flags = 3328
        thread_num = 2
        
        realignf = rf.Realignment_filter(ref_genome, 
        tumor_min_mismatch, normal_max_mismatch, window_size, score_difference, 
        blat_path, header_flag, max_depth, exclude_sam_flags, thread_num)

        bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        bam2 = None
        output = cur_dir + "/../data/5929_small_realignment_result_test4_3.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test32.txt"
        tumor_sample = "5929_tumor"
        normal_sample = None
        realignf.filter_vcf(bam1, bam2, output, target_mutation_file, tumor_sample, normal_sample)
        
        answer_file = cur_dir + "/../data/5929_small_realignment_result_answer_test4_2.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        
    def test4_4(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        ref_genome = cur_dir + "/../data/GRCh37.fa"
        tumor_min_mismatch = 0
        normal_max_mismatch = 100000
        window_size = 200
        score_difference = 5
        blat_path = self.blat
        header_flag = False
        max_depth = 5000
        exclude_sam_flags = 3328
        thread_num = 1
        
        realignf = rf.Realignment_filter(ref_genome, 
        tumor_min_mismatch, normal_max_mismatch, window_size, score_difference, 
        blat_path, header_flag, max_depth, exclude_sam_flags, thread_num)

        bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        bam2 = None
        output = cur_dir + "/../data/5929_small_realignment_result_test4_4.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test34.txt"
        tumor_sample = "5929_tumor"
        normal_sample = None

        with self.assertRaises(ValueError) as er:
            realignf.filter_vcf(bam1, bam2, output, target_mutation_file, tumor_sample, normal_sample)
        # self.assertTrue(er.exception.message.endswith('is not in list'))
        
        
    def test4_5(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        ref_genome = cur_dir + "/../data/GRCh37.fa"
        tumor_min_mismatch = 0
        normal_max_mismatch = 100000
        window_size = 200
        score_difference = 5
        blat_path = self.blat
        header_flag = False
        max_depth = 5000
        exclude_sam_flags = 3328
        thread_num = 2
        
        realignf = rf.Realignment_filter(ref_genome, 
        tumor_min_mismatch, normal_max_mismatch, window_size, score_difference, 
        blat_path, header_flag, max_depth, exclude_sam_flags, thread_num)

        bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        bam2 = None
        output = cur_dir + "/../data/5929_small_realignment_result_test4_5.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test34.txt"
        tumor_sample = "5929_tumor"
        normal_sample = None
        
        with self.assertRaises(RuntimeError) as er:
            realignf.filter_vcf(bam1, bam2, output, target_mutation_file, tumor_sample, normal_sample)
        # self.assertTrue(er.exception.message.endswith('There was an error!'))
        
    def test4_6(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        ref_genome = cur_dir + "/../data/GRCh37.fa"
        tumor_min_mismatch = 0
        normal_max_mismatch = 100000
        window_size = 200
        score_difference = 10
        blat_path = self.blat
        header_flag = True
        max_depth = 5000
        exclude_sam_flags = 3328
        thread_num = 2
        
        realignf = rf.Realignment_filter(ref_genome, 
        tumor_min_mismatch, normal_max_mismatch, window_size, score_difference, 
        blat_path, header_flag, max_depth, exclude_sam_flags, thread_num)

        bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        bam2 = None
        output = cur_dir + "/../data/5929_small_realignment_result_test4_6.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test32.txt"
        tumor_sample = "5929_tumor"
        normal_sample = None
        realignf.filter_vcf(bam1, bam2, output, target_mutation_file, tumor_sample, normal_sample)
        
        answer_file = cur_dir + "/../data/5929_small_realignment_result_answer_test4_6.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        
    
    def test5_1(self):
        
        cmd = ['mutfilter', '--version']
        subprocess.check_call(cmd)
        
        
    def test5_2(self):
        
        cmd = ['mutfilter', 'realignment', '--help']
        subprocess.check_call(cmd)


    def test5_3(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        
        ref_genome = cur_dir + "/../data/GRCh37.fa"
        tumor_min_mismatch = '0'
        normal_max_mismatch = '100000'
        window_size = '500'
        score_difference = '5'
        blat_path = self.blat
        max_depth = '5000'
        exclude_sam_flags = '3328'
        thread_num = '1'
        bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_realignment_result_test3_1.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test1.txt"

        cmd = ['mutfilter', 'realignment',
            '-t', target_mutation_file,
            '-1', bam1,
            '-2', bam2,
            '-o', output,
            '-r', ref_genome,
            '-b', self.blat,
            '-m', tumor_min_mismatch,
            '-M', normal_max_mismatch,
            '-w', window_size,
            '-s', score_difference,
            '-d', max_depth,
            '-F', exclude_sam_flags,
            '-T', thread_num,
            '--header']
        
        subprocess.check_call(cmd)

        answer_file = cur_dir + "/../data/5929_small_realignment_result_answer_test1_1.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))

    ######################################
    # Tumor/Normal Pair, Execute
    ######################################
    
    def test6_1(self):
        
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        ref_genome = cur_dir + "/../data/GRCh37.fa"
        tumor_min_mismatch = 0
        normal_max_mismatch = 100000
        window_size = 500
        score_difference = 5
        blat_path = self.blat
        header_flag = True
        max_depth = 5000
        exclude_sam_flags = 3328
        thread_num = 1
        
        bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_realignment_result_test6_1.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test1.txt"
        answer_file = cur_dir + "/../data/5929_small_realignment_result_answer_test1_1.txt"
        
        parser = genomon_mutation_filter.parser.create_parser()
        args = parser.parse_args(["realignment", "-1", bam1, "-2", bam2, "-t", target_mutation_file, "-o", output, "-r", ref_genome, "-b", blat_path, "-m", str(tumor_min_mismatch),  "-M", str(normal_max_mismatch), "-s", str(score_difference), "-w", str(window_size), "-d", str(max_depth), "-F", str(exclude_sam_flags), "--header"])
        args.func(args)
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))


if __name__ == "__main__":
    unittest.main()

