#! /usr/bin/env python

import sys
import unittest
import os, tempfile, shutil, filecmp
from lib.mutfilter import indel_filter as idf

class TestIndel(unittest.TestCase):

    def setUp(self):
        if sys.version_info.major == 2:
            samtools = "/home/ubuntu/miniconda2/bin/samtools"
        else:
            samtools = "/home/ubuntu/miniconda3/bin/samtools"
        self.samtools = samtools 

    ######################################
    # Tumor/Normal Pair, Annoformat
    ######################################
    def test1_1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        search_length = 40
        min_depth = 8
        min_mismatch = 100000
        af_thres = 1
        neighbor = 5
        header_flag = True
        samtools_path = self.samtools
        samtools_params = '-q 20 -BQ0 -d 10000000'
        ref_genome = cur_dir + "/../data/GRCh37.fa"
        
        indelf = idf.Indel_filter(search_length, min_depth, 
        min_mismatch, af_thres, neighbor, header_flag, 
        samtools_path, samtools_params, ref_genome)
        
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_indel_result_test1_1.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test1.txt"
        indelf.filter(target_mutation_file, bam2, output)
        
        answer_file = cur_dir + "/../data/5929_small_indel_result_answer_test1_1.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))


    def test1_2(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        search_length = 40
        min_depth = 8
        min_mismatch = 100000
        af_thres = 1
        neighbor = 5
        header_flag = True
        samtools_path = self.samtools
        samtools_params = '-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP'
        ref_genome = cur_dir + "/../data/GRCh37.fa"
        
        indelf = idf.Indel_filter(search_length, min_depth, 
        min_mismatch, af_thres, neighbor, header_flag, 
        samtools_path, samtools_params, ref_genome)
        
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_indel_result_test1_2.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test2.txt"
        indelf.filter(target_mutation_file, bam2, output)
        
        answer_file = cur_dir + "/../data/5929_small_indel_result_answer_test1_2.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        

    def test1_3(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        search_length = 40
        min_depth = 8
        min_mismatch = 100000
        af_thres = 1
        neighbor = 5
        header_flag = True
        samtools_path = self.samtools
        samtools_params = '-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP'
        ref_genome = cur_dir + "/../data/GRCh37.fa"
        
        indelf = idf.Indel_filter(search_length, min_depth, 
        min_mismatch, af_thres, neighbor, header_flag, 
        samtools_path, samtools_params, ref_genome)
        
        bam2 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_indel_result_test1_3.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test2.txt"
        indelf.filter(target_mutation_file, bam2, output)
        
        answer_file = cur_dir + "/../data/5929_small_indel_result_answer_test1_3.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        
    def test1_4(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        search_length = 40
        min_depth = 8
        min_mismatch = 100000
        af_thres = 1
        neighbor = 5
        header_flag = False
        samtools_path = self.samtools
        samtools_params = '-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP'
        ref_genome = cur_dir + "/../data/GRCh37.fa"
        
        indelf = idf.Indel_filter(search_length, min_depth, 
        min_mismatch, af_thres, neighbor, header_flag, 
        samtools_path, samtools_params, ref_genome)
        
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_indel_result_test1_4.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test4.txt"
        indelf.filter(target_mutation_file, bam2, output)
        
        answer_file = cur_dir + "/../data/5929_small_indel_result_answer_test1_4.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        

    ######################################
    # Tumor/Normal Pair, VCF format
    ######################################
    def test2_1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        search_length = 40
        min_depth = 8
        min_mismatch = 100000
        af_thres = 1
        neighbor = 5
        header_flag = True
        samtools_path = self.samtools
        samtools_params = '-q 20 -BQ0 -d 10000000'
        ref_genome = cur_dir + "/../data/GRCh37.fa"
        
        indelf = idf.Indel_filter(search_length, min_depth, 
        min_mismatch, af_thres, neighbor, header_flag, 
        samtools_path, samtools_params, ref_genome)
        
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_indel_result_test2_1.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test21.txt"
        tumor_sample = "5929_tumor"
        normal_sample = "5929_control"
        indelf.filter_vcf(target_mutation_file, bam2, output, tumor_sample, normal_sample)
        
        answer_file = cur_dir + "/../data/5929_small_indel_result_answer_test2_1.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))


    def test2_2(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        search_length = 40
        min_depth = 8
        min_mismatch = 100000
        af_thres = 1
        neighbor = 5
        header_flag = True
        samtools_path = self.samtools
        samtools_params = '-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP'
        ref_genome = cur_dir + "/../data/GRCh37.fa"
        
        indelf = idf.Indel_filter(search_length, min_depth, 
        min_mismatch, af_thres, neighbor, header_flag, 
        samtools_path, samtools_params, ref_genome)
        
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_indel_result_test2_2.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test22.txt"
        tumor_sample = "5929_tumor"
        normal_sample = "5929_control"
        indelf.filter_vcf(target_mutation_file, bam2, output, tumor_sample, normal_sample)
        
        answer_file = cur_dir + "/../data/5929_small_indel_result_answer_test2_2.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        

    def test2_3(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        search_length = 40
        min_depth = 8
        min_mismatch = 100000
        af_thres = 1
        neighbor = 5
        header_flag = True
        samtools_path = self.samtools
        samtools_params = '-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP'
        ref_genome = cur_dir + "/../data/GRCh37.fa"
        
        indelf = idf.Indel_filter(search_length, min_depth, 
        min_mismatch, af_thres, neighbor, header_flag, 
        samtools_path, samtools_params, ref_genome)
        
        bam2 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_indel_result_test2_3.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test22.txt"
        tumor_sample = "5929_tumor"
        normal_sample = "5929_control"
        indelf.filter_vcf(target_mutation_file, bam2, output, tumor_sample, normal_sample)
        
        answer_file = cur_dir + "/../data/5929_small_indel_result_answer_test2_3.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))


    def test2_4(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        search_length = 40
        min_depth = 8
        min_mismatch = 100000
        af_thres = 1
        neighbor = 5
        header_flag = False
        samtools_path = self.samtools
        samtools_params = '-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP'
        ref_genome = cur_dir + "/../data/GRCh37.fa"
        
        indelf = idf.Indel_filter(search_length, min_depth, 
        min_mismatch, af_thres, neighbor, header_flag, 
        samtools_path, samtools_params, ref_genome)
        
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_indel_result_test2_4.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test24.txt"
        tumor_sample = "5929_tumor"
        normal_sample = "5929_control"
        with self.assertRaises(ValueError) as er:
            indelf.filter_vcf(target_mutation_file, bam2, output, tumor_sample, normal_sample)
        

if __name__ == "__main__":
    unittest.main()

