#! /usr/bin/env python

import sys
import unittest
import os, tempfile, shutil, filecmp
import subprocess
from genomon_mutation_filter import indel_filter as idf
import genomon_mutation_filter

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
        thread_num = 1
        
        indelf = idf.Indel_filter(search_length, min_depth, 
        min_mismatch, af_thres, neighbor, header_flag, 
        samtools_path, samtools_params, ref_genome, thread_num)
        
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
        thread_num = 1
        
        indelf = idf.Indel_filter(search_length, min_depth, 
        min_mismatch, af_thres, neighbor, header_flag, 
        samtools_path, samtools_params, ref_genome, thread_num)
        
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
        thread_num = 1
        
        indelf = idf.Indel_filter(search_length, min_depth, 
        min_mismatch, af_thres, neighbor, header_flag, 
        samtools_path, samtools_params, ref_genome, thread_num)
        
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
        thread_num = 1
        
        indelf = idf.Indel_filter(search_length, min_depth, 
        min_mismatch, af_thres, neighbor, header_flag, 
        samtools_path, samtools_params, ref_genome, thread_num)
        
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_indel_result_test1_4.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test4.txt"
        indelf.filter(target_mutation_file, bam2, output)
        
        answer_file = cur_dir + "/../data/5929_small_indel_result_answer_test1_4.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))

        
    def test1_5(self):
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
        thread_num = 2
        
        indelf = idf.Indel_filter(search_length, min_depth, 
        min_mismatch, af_thres, neighbor, header_flag, 
        samtools_path, samtools_params, ref_genome, thread_num)
        
        bam2 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_indel_result_test1_3.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test2.txt"
        indelf.filter(target_mutation_file, bam2, output)
        
        answer_file = cur_dir + "/../data/5929_small_indel_result_answer_test1_3.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))

        
    def test1_6(self):
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
        thread_num = 2
        
        indelf = idf.Indel_filter(search_length, min_depth, 
        min_mismatch, af_thres, neighbor, header_flag, 
        samtools_path, samtools_params, ref_genome, thread_num)
        
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
        thread_num = 1
        
        indelf = idf.Indel_filter(search_length, min_depth, 
        min_mismatch, af_thres, neighbor, header_flag, 
        samtools_path, samtools_params, ref_genome, thread_num)
        
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
        thread_num = 1
                
        indelf = idf.Indel_filter(search_length, min_depth, 
        min_mismatch, af_thres, neighbor, header_flag, 
        samtools_path, samtools_params, ref_genome, thread_num)
        
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
        thread_num = 1
        
        indelf = idf.Indel_filter(search_length, min_depth, 
        min_mismatch, af_thres, neighbor, header_flag, 
        samtools_path, samtools_params, ref_genome, thread_num)
        
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
        thread_num = 1
        
        indelf = idf.Indel_filter(search_length, min_depth, 
        min_mismatch, af_thres, neighbor, header_flag, 
        samtools_path, samtools_params, ref_genome, thread_num)
        
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_indel_result_test2_4.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test24.txt"
        tumor_sample = "5929_tumor"
        normal_sample = "5929_control"
        with self.assertRaises(ValueError) as er:
            indelf.filter_vcf(target_mutation_file, bam2, output, tumor_sample, normal_sample)
        

    def test2_5(self):
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
        thread_num = 2
        
        indelf = idf.Indel_filter(search_length, min_depth, 
        min_mismatch, af_thres, neighbor, header_flag, 
        samtools_path, samtools_params, ref_genome, thread_num)
        
        bam2 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_indel_result_test2_3.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test22.txt"
        tumor_sample = "5929_tumor"
        normal_sample = "5929_control"
        indelf.filter_vcf(target_mutation_file, bam2, output, tumor_sample, normal_sample)
        
        answer_file = cur_dir + "/../data/5929_small_indel_result_answer_test2_3.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))


    def test2_6(self):
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
        thread_num = 2
        
        indelf = idf.Indel_filter(search_length, min_depth, 
        min_mismatch, af_thres, neighbor, header_flag, 
        samtools_path, samtools_params, ref_genome, thread_num)
        
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_indel_result_test2_4.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test24.txt"
        tumor_sample = "5929_tumor"
        normal_sample = "5929_control"
        with self.assertRaises(RuntimeError) as er:
            indelf.filter_vcf(target_mutation_file, bam2, output, tumor_sample, normal_sample)
            
            
    def test3_1(self):
        
        cmd = ['mutfilter', '--version']
        subprocess.check_call(cmd)
        
        
    def test3_2(self):
        
        cmd = ['mutfilter', 'indel', '--help']
        subprocess.check_call(cmd)


    def test3_3(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        
        search_length = '40'
        min_depth = '8'
        min_mismatch = '100000'
        af_thres = '1'
        neighbor = '5'
        header_flag = True
        samtools_path = self.samtools
        samtools_params = '-q 20 -BQ0 -d 10000000'
        ref_genome = cur_dir + "/../data/GRCh37.fa"
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_indel_result_test3_1.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test1.txt"

        cmd = ['mutfilter', 'indel',
            '-t', target_mutation_file,
            '-2', bam2,
            '-o', output,
            '-r', ref_genome,
            '-l', search_length,
            '-d', min_depth,
            '-m', min_mismatch,
            '-a', af_thres,
            '-n', neighbor,
            '-s', samtools_path,
            '-S', samtools_params,
            '--header']
        
        subprocess.check_call(cmd)

        answer_file = cur_dir + "/../data/5929_small_indel_result_answer_test1_1.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))


    ######################################
    # Tumor/Normal Pair, Execute
    ######################################
    
    def test4_1(self):
        
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
        bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        output = cur_dir + "/../data/5929_small_indel_result_test4_1.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test1.txt"
        answer_file = cur_dir + "/../data/5929_small_indel_result_answer_test1_1.txt"
        
        parser = genomon_mutation_filter.parser.create_parser()
        args = parser.parse_args(["indel", "-2", bam2, "-t", target_mutation_file, "-o", output, "-r", ref_genome, "-l", str(search_length), "-n", str(neighbor),  "-d", str(min_depth), "-m", str(min_mismatch), "-a", str(af_thres), "-s", samtools_path, "-S", samtools_params, "--header"])
        args.func(args)
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))



    ######################################
    # Tumor/Normal Pair, Annoformat CRAM
    ######################################
    def test5_1(self):
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
        thread_num = 1
        
        indelf = idf.Indel_filter(search_length, min_depth, 
        min_mismatch, af_thres, neighbor, header_flag, 
        samtools_path, samtools_params, ref_genome, thread_num)
        
        bam2 = cur_dir + "/../data/5929_control_small.markdup.cram"
        output = cur_dir + "/../data/5929_small_indel_result_test1_1.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test1.txt"
        indelf.filter(target_mutation_file, bam2, output)
        
        answer_file = cur_dir + "/../data/5929_small_indel_result_answer_test1_1.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))


    def test5_2(self):
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
        thread_num = 1
        
        indelf = idf.Indel_filter(search_length, min_depth, 
        min_mismatch, af_thres, neighbor, header_flag, 
        samtools_path, samtools_params, ref_genome, thread_num)
        
        bam2 = cur_dir + "/../data/5929_control_small.markdup.cram"
        output = cur_dir + "/../data/5929_small_indel_result_test1_2.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test2.txt"
        indelf.filter(target_mutation_file, bam2, output)
        
        answer_file = cur_dir + "/../data/5929_small_indel_result_answer_test1_2.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        

    def test5_3(self):
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
        thread_num = 1
        
        indelf = idf.Indel_filter(search_length, min_depth, 
        min_mismatch, af_thres, neighbor, header_flag, 
        samtools_path, samtools_params, ref_genome, thread_num)
        
        bam2 = cur_dir + "/../data/5929_tumor_small.markdup.cram"
        output = cur_dir + "/../data/5929_small_indel_result_test1_3.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test2.txt"
        indelf.filter(target_mutation_file, bam2, output)
        
        answer_file = cur_dir + "/../data/5929_small_indel_result_answer_test1_3.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))

        
    def test5_4(self):
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
        thread_num = 1
        
        indelf = idf.Indel_filter(search_length, min_depth, 
        min_mismatch, af_thres, neighbor, header_flag, 
        samtools_path, samtools_params, ref_genome, thread_num)
        
        bam2 = cur_dir + "/../data/5929_control_small.markdup.cram"
        output = cur_dir + "/../data/5929_small_indel_result_test1_4.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test4.txt"
        indelf.filter(target_mutation_file, bam2, output)
        
        answer_file = cur_dir + "/../data/5929_small_indel_result_answer_test1_4.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))

        
    def test5_5(self):
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
        thread_num = 2
        
        indelf = idf.Indel_filter(search_length, min_depth, 
        min_mismatch, af_thres, neighbor, header_flag, 
        samtools_path, samtools_params, ref_genome, thread_num)
        
        bam2 = cur_dir + "/../data/5929_tumor_small.markdup.cram"
        output = cur_dir + "/../data/5929_small_indel_result_test1_3.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test2.txt"
        indelf.filter(target_mutation_file, bam2, output)
        
        answer_file = cur_dir + "/../data/5929_small_indel_result_answer_test1_3.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))

        
    def test5_6(self):
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
        thread_num = 2
        
        indelf = idf.Indel_filter(search_length, min_depth, 
        min_mismatch, af_thres, neighbor, header_flag, 
        samtools_path, samtools_params, ref_genome, thread_num)
        
        bam2 = cur_dir + "/../data/5929_control_small.markdup.cram"
        output = cur_dir + "/../data/5929_small_indel_result_test1_4.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test4.txt"
        indelf.filter(target_mutation_file, bam2, output)
        
        answer_file = cur_dir + "/../data/5929_small_indel_result_answer_test1_4.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))




if __name__ == "__main__":
    unittest.main()

