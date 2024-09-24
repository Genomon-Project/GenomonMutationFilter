#! /usr/bin/env python

import sys
import unittest
import os, shutil, filecmp
import subprocess
import genomon_mutation_filter.oxog_filter as of

class TestOxog(unittest.TestCase):

    def test1_1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter(None,None)

        bases = "A"
        quals = "F"
        var2nm = oxof.parse_bases(bases, quals)

        self.assertTrue(var2nm == {"A":1})

    def test1_2(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter(None,None)

        bases = "AA"
        quals = "FF"
        var2nm = oxof.parse_bases(bases, quals)

        self.assertTrue(var2nm == {"A":2})

    def test1_3(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter(None,None)

        bases = "tt"
        quals = "FF"
        var2nm = oxof.parse_bases(bases, quals)

        self.assertTrue(var2nm == {"t":2})

    def test1_4(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter(None,None)

        bases = "A$a"
        quals = "FF"
        var2nm = oxof.parse_bases(bases, quals)

        self.assertTrue(var2nm == {"A":1,"a":1})

    def test1_5(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter(None,None)

        bases = "AT-1A"
        quals = "FF"
        var2nm = oxof.parse_bases(bases, quals)

        self.assertTrue(var2nm == {'A': 1, 'T': 1})

    def test1_6(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter(None,None)

        bases = "At+1a"
        quals = "FF"
        var2nm = oxof.parse_bases(bases, quals)

        self.assertTrue(var2nm == {'A': 1, 't': 1})

    def test1_7(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter(None,None)

        bases = "*N"
        quals = "*F"
        var2nm = oxof.parse_bases(bases, quals)

        self.assertTrue(var2nm == {'N': 1})

    # Error is collect because base length and quals length are not same.
    # def test1_8(self):
    #     cur_dir = os.path.dirname(os.path.abspath(__file__))

    #     oxof = of.Oxog_filter(None,None)

    #     bases = "AA"
    #     quals = "F"
    #     var2nm = oxof.parse_bases(bases, quals)

    def test2_1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter(None,None)

        ref = "C"
        alt = "A"
        f1r2 = 2
        f2r1 = 17 
        f_oxog = oxof.flag_oxog(ref, alt, f1r2, f2r1)
        self.assertTrue(f_oxog == 0)

    def test2_2(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter(None,None)

        ref = "C"
        alt = "A"
        f1r2 = 2
        f2r1 = 18 
        f_oxog = oxof.flag_oxog(ref, alt, f1r2, f2r1)
        self.assertTrue(f_oxog == 1)

    def test2_3(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter(None,None)

        ref = "C"
        alt = "A"
        f1r2 = 1
        f2r1 = 8
        f_oxog = oxof.flag_oxog(ref, alt, f1r2, f2r1)
        self.assertTrue(f_oxog == 1)

    def test2_4(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter(None,None)

        ref = "G"
        alt = "T"
        f1r2 = 17
        f2r1 = 2 
        f_oxog = oxof.flag_oxog(ref, alt, f1r2, f2r1)
        self.assertTrue(f_oxog == 0)

    def test2_5(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter(None,None)

        ref = "G"
        alt = "T"
        f1r2 = 18 
        f2r1 = 2 
        f_oxog = oxof.flag_oxog(ref, alt, f1r2, f2r1)
        self.assertTrue(f_oxog == 1)

    def test2_6(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter(None,None)

        ref = "G"
        alt = "T"
        f1r2 = 8
        f2r1 = 1
        f_oxog = oxof.flag_oxog(ref, alt, f1r2, f2r1)
        self.assertTrue(f_oxog == 1)

    def test3_1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        oxof = of.Oxog_filter("samtools","-q 20 -B -Q15 -d 10000000")
        bam = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        reg = "chr1:12345678-12345678"

        with open(os.devnull, 'w') as FNULL:
            d1, d2 = oxof.call_mpileup(reg, bam, FNULL)

        self.assertTrue(True)



