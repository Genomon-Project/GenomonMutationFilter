#! /usr/bin/env python

import sys
import unittest
import os, shutil, filecmp
import subprocess
import genomon_mutation_filter.oxog_filter as of

class TestOxog(unittest.TestCase):

    def test1_1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter(None,None,None)

        bases = "A"
        quals = "F"
        flags = "64"
        var2nm_1st, var2nm_2nd = oxof.parse_bases(bases, quals, flags)

        self.assertTrue(var2nm_1st == {"A":1})
        self.assertTrue(var2nm_2nd == {"A":0})

    def test1_2(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter(None,None,None)

        bases = "AA"
        quals = "FF"
        flags = "64,99"
        var2nm_1st, var2nm_2nd = oxof.parse_bases(bases, quals, flags)

        self.assertTrue(var2nm_1st == {"A":2})
        self.assertTrue(var2nm_2nd == {"A":0})

    def test1_3(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter(None,None,None)

        bases = "tt"
        quals = "FF"
        flags = "64,75"
        var2nm_1st, var2nm_2nd = oxof.parse_bases(bases, quals, flags)

        self.assertTrue(var2nm_1st == {"t":2})
        self.assertTrue(var2nm_2nd == {"t":0})

    def test1_4(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter(None,None,None)

        bases = "A$a"
        quals = "FF"
        flags = "64,71"
        var2nm_1st, var2nm_2nd = oxof.parse_bases(bases, quals, flags)

        self.assertTrue(var2nm_1st == {"A":1,"a":1})
        self.assertTrue(var2nm_2nd == {"A":0,"a":0})

    def test1_5(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter(None,None,None)

        bases = "AT-1A"
        quals = "FF"
        flags = "64,66"
        var2nm_1st, var2nm_2nd = oxof.parse_bases(bases, quals, flags)

        self.assertTrue(var2nm_1st == {'A': 1, 'T': 1})
        self.assertTrue(var2nm_2nd == {"A":0 , "T":0})

    def test1_6(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter(None,None,None)

        bases = "At+1a"
        quals = "FF"
        flags = "65,64"
        var2nm_1st, var2nm_2nd = oxof.parse_bases(bases, quals, flags)

        self.assertTrue(var2nm_1st == {'A': 1, 't': 1})
        self.assertTrue(var2nm_2nd == {"A":0, "t":0})

    def test1_7(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter(None,None,None)

        bases = "*N"
        quals = "*F"
        flags = "63,64"
        var2nm_1st, var2nm_2nd = oxof.parse_bases(bases, quals, flags)

        self.assertTrue(var2nm_1st == {'N': 1})
        self.assertTrue(var2nm_2nd == {"N":0})

    def test1_8(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter(None,None,None)

        bases = "A"
        quals = "F"
        flags = "128"
        var2nm_1st, var2nm_2nd = oxof.parse_bases(bases, quals, flags)

        self.assertTrue(var2nm_1st == {"A":0})
        self.assertTrue(var2nm_2nd == {"A":1})

    def test1_9(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter(None,None,None)

        bases = "AA"
        quals = "FF"
        flags = "128,163"
        var2nm_1st, var2nm_2nd = oxof.parse_bases(bases, quals, flags)

        self.assertTrue(var2nm_1st == {"A":0})
        self.assertTrue(var2nm_2nd == {"A":2})

    def test1_10(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter(None,None,None)

        bases = "tt"
        quals = "FF"
        flags = "64,128"
        var2nm_1st, var2nm_2nd = oxof.parse_bases(bases, quals, flags)

        self.assertTrue(var2nm_1st == {"t":1})
        self.assertTrue(var2nm_2nd == {"t":1})

    def test1_11(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter(None,None,None)

        bases = "A$a"
        quals = "FF"
        flags = "2211,64"
        var2nm_1st, var2nm_2nd = oxof.parse_bases(bases, quals, flags)

        self.assertTrue(var2nm_1st == {"A":0,"a":1})
        self.assertTrue(var2nm_2nd == {"A":1,"a":0})

    def test1_12(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter(None,None,None)

        bases = "t+1aA"
        quals = "FF"
        flags = "128,128"
        var2nm_1st, var2nm_2nd = oxof.parse_bases(bases, quals, flags)

        self.assertTrue(var2nm_1st == {'A': 0, 't': 0})
        self.assertTrue(var2nm_2nd == {"A":1, "t":1})

    # Error is collect because base length and quals length are not same.
    # def test1_13(self):
    #     cur_dir = os.path.dirname(os.path.abspath(__file__))

    #     oxof = of.Oxog_filter(None,None,None)

    #     bases = "AA"
    #     quals = "F"
    #     flags = "64,64"
    #     var2nm_1st, var2nm_2nd = oxof.parse_bases(bases, quals, flags)

    def test2_1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter(None,None,None)

        ref = "C"
        alt = "A"
        f1r2 = 2
        f2r1 = 17 
        f_oxog = oxof.flag_oxog(ref, alt, f1r2, f2r1)
        self.assertTrue(f_oxog == 0)

    def test2_2(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter(None,None,None)

        ref = "C"
        alt = "A"
        f1r2 = 2
        f2r1 = 18 
        f_oxog = oxof.flag_oxog(ref, alt, f1r2, f2r1)
        self.assertTrue(f_oxog == 1)

    def test2_3(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter(None,None,None)

        ref = "C"
        alt = "A"
        f1r2 = 1
        f2r1 = 8
        f_oxog = oxof.flag_oxog(ref, alt, f1r2, f2r1)
        self.assertTrue(f_oxog == 1)

    def test2_4(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter(None,None,None)

        ref = "G"
        alt = "T"
        f1r2 = 17
        f2r1 = 2 
        f_oxog = oxof.flag_oxog(ref, alt, f1r2, f2r1)
        self.assertTrue(f_oxog == 0)

    def test2_5(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter(None,None,None)

        ref = "G"
        alt = "T"
        f1r2 = 18 
        f2r1 = 2 
        f_oxog = oxof.flag_oxog(ref, alt, f1r2, f2r1)
        self.assertTrue(f_oxog == 1)

    def test2_6(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter(None,None,None)

        ref = "G"
        alt = "T"
        f1r2 = 8
        f2r1 = 1
        f_oxog = oxof.flag_oxog(ref, alt, f1r2, f2r1)
        self.assertTrue(f_oxog == 1)

    def test3_1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter("samtools","-q 20 -B -Q15 -d 10000000 --output-extra FLAG -x",1)
        bam = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        reg = "chr1:12345678-12345678"

        with open(os.devnull, 'w') as FNULL:
            d1, d2 = oxof.call_mpileup(reg, bam, FNULL)

        self.assertTrue(True)


    def test4_1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter("samtools","-q 20 -B -Q15 -d 10000000 --output-extra FLAG -x",1)
        bam = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_vcf = cur_dir + "/../data/5929_small_mutation_result_test21.txt"
        output = cur_dir + "/../data/5929_small_mutation_result_test21_oxogout.txt"

        with open(os.devnull, 'w') as FNULL:
            oxof.filter_main_vcf(in_vcf, bam, output, "5929_tumor", "5929_control")

        self.assertTrue(True)


