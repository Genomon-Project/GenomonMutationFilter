#! /usr/bin/env python

import sys
import unittest
import os, tempfile, shutil, filecmp
import subprocess
from genomon_mutation_filter import simple_repeat_filter as sif
import genomon_mutation_filter

class TestSimpleRepeat(unittest.TestCase):


    ######################################
    # Tumor/Normal Pair, Annoformat
    ######################################
    def test1_1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        header_flag = True
        simple_repeat_db = cur_dir + "/../data/simpleRepeat.bed.gz"
        simplef = sif.Simple_repeat_filter(simple_repeat_db, header_flag)

        output = cur_dir + "/../data/5929_small_simple_result_test1_1.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test2.txt"
        simplef.filter(target_mutation_file, output)
        
        answer_file = cur_dir + "/../data/5929_small_simple_result_answer_test1_1.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))


    def test1_2(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        header_flag = False
        simple_repeat_db = cur_dir + "/../data/simpleRepeat.bed.gz"
        simplef = sif.Simple_repeat_filter(simple_repeat_db, header_flag)

        output = cur_dir + "/../data/5929_small_simple_result_test1_2.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test4.txt"
        simplef.filter(target_mutation_file, output)
        
        answer_file = cur_dir + "/../data/5929_small_simple_result_answer_test1_2.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        

    ######################################
    # Tumor/Normal Pair, VCF format
    ######################################
    def test2_1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        header_flag = True
        simple_repeat_db = cur_dir + "/../data/simpleRepeat.bed.gz"
        simplef = sif.Simple_repeat_filter(simple_repeat_db, header_flag)

        output = cur_dir + "/../data/5929_small_simple_result_test2_1.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test22.txt"
        simplef.filter_vcf(target_mutation_file, output)
        
        answer_file = cur_dir + "/../data/5929_small_simple_result_answer_test2_1.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))


    def test2_2(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        header_flag = False
        simple_repeat_db = cur_dir + "/../data/simpleRepeat.bed.gz"
        simplef = sif.Simple_repeat_filter(simple_repeat_db, header_flag)

        output = cur_dir + "/../data/5929_small_simple_result_test2_2.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test24.txt"
        simplef.filter_vcf(target_mutation_file, output)
        
        answer_file = cur_dir + "/../data/5929_small_simple_result_answer_test2_2.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))
        


    def test3_1(self):
        
        cmd = ['mutfilter', '--version']
        subprocess.check_call(cmd)
        
        
    def test3_2(self):
        
        cmd = ['mutfilter', 'simplerepeat', '--help']
        subprocess.check_call(cmd)
        
        
    def test3_3(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        
        simple_repeat_db = cur_dir + "/../data/simpleRepeat.bed.gz"
        output = cur_dir + "/../data/5929_small_simple_result_test3_1.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test2.txt"

        cmd = ['mutfilter', 'simplerepeat',
            '-t', target_mutation_file,
            '-o', output,
            '-S', simple_repeat_db,
            '--header']
        
        subprocess.check_call(cmd)

        answer_file = cur_dir + "/../data/5929_small_simple_result_answer_test1_1.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))


    def test4_1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        
        simple_repeat_db = cur_dir + "/../data/simpleRepeat.bed.gz"
        output = cur_dir + "/../data/5929_small_simple_result_test4_1.txt"
        target_mutation_file = cur_dir + "/../data/5929_small_mutation_result_test2.txt"

        cmd = ['simplerepeat',
            '-t', target_mutation_file,
            '-o', output,
            '-S', simple_repeat_db,
            '--header']
        
        parser = genomon_mutation_filter.parser.create_parser()
        args = parser.parse_args(cmd)
        args.func(args)

        answer_file = cur_dir + "/../data/5929_small_simple_result_answer_test1_1.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))

if __name__ == "__main__":
    unittest.main()

