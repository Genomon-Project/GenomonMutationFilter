#! /usr/bin/env python

import sys
import unittest
import os, tempfile, shutil, filecmp
import subprocess
import genomon_mutation_filter.position_filter as pf

class TestPosition(unittest.TestCase):

    def test1_1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        posf = pf.Position_filter(None,None,None,None)

        bases = "A"
        positions = "50"
        qnames = "name1"
        flags = "64"
        ref = "T"
        var2nm, var2pos, var2qname, var2flag = posf.parse_bases(bases, positions, qnames, flags, ref)

        self.assertTrue(var2nm == {"A":1})
        self.assertTrue(var2pos== {"A":['50']})
        self.assertTrue(var2qname== {"A":['name1']})
        self.assertTrue(var2flag== {"A":['64']})

    def test1_2(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        posf = pf.Position_filter(None,None,None,None)

        bases = "AA"
        positions = "50,49"
        qnames = "name1,name2"
        flags = "64,128"
        ref = "T"
        var2nm, var2pos, var2qname, var2flag = posf.parse_bases(bases, positions, qnames, flags, ref)

        self.assertTrue(var2nm == {"A":2})
        self.assertTrue(var2pos== {"A":['50','49']})
        self.assertTrue(var2qname== {"A":['name1','name2']})
        self.assertTrue(var2flag== {"A":['64','128']})

    def test1_3(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        posf = pf.Position_filter(None,None,None,None)

        bases = ",."
        positions = "50,49"
        qnames = "name1,name2"
        flags = "64,128"
        ref = "T"
        var2nm, var2pos, var2qname, var2flag = posf.parse_bases(bases, positions, qnames, flags, ref)

        self.assertTrue(var2nm == {"T":2})
        self.assertTrue(var2pos== {"T":['50','49']})
        self.assertTrue(var2qname== {"T":['name1','name2']})
        self.assertTrue(var2flag== {"T":['64','128']})

    def test1_4(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        posf = pf.Position_filter(None,None,None,None)

        bases = "A$a"
        positions = "50,49"
        qnames = "name1,name2"
        flags = "64,128"
        ref = "T"
        var2nm, var2pos, var2qname, var2flag = posf.parse_bases(bases, positions, qnames, flags, ref)

        self.assertTrue(var2nm == {"A":2})
        self.assertTrue(var2pos== {"A":['50','49']})
        self.assertTrue(var2qname== {"A":['name1','name2']})
        self.assertTrue(var2flag== {"A":['64','128']})

    def test1_5(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        posf = pf.Position_filter(None,None,None,None)

        bases = "A.-1A"
        positions = "50,49"
        qnames = "name1,name2"
        flags = "64,128"
        ref = "T"
        var2nm, var2pos, var2qname, var2flag = posf.parse_bases(bases, positions, qnames, flags, ref)

        self.assertTrue(var2nm == {'A': 1, 'T': 1, '-A': 1})
        self.assertTrue(var2pos == {'A': ['50'], 'T': ['49'], '-A': ['49']})
        self.assertTrue(var2qname == {'A': ['name1'], 'T': ['name2'], '-A': ['name2']})
        self.assertTrue(var2flag == {'A': ['64'], 'T': ['128'], '-A': ['128']})

    def test1_6(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        posf = pf.Position_filter(None,None,None,None)

        bases = "A,+1a"
        positions = "50,49"
        qnames = "name1,name2"
        flags = "64,128"
        ref = "T"
        var2nm, var2pos, var2qname, var2flag = posf.parse_bases(bases, positions, qnames, flags, ref)

        self.assertTrue(var2nm == {'A': 1, 'T': 1, '+A': 1})
        self.assertTrue(var2pos == {'A': ['50'], 'T': ['49'], '+A': ['49']})
        self.assertTrue(var2qname == {'A': ['name1'], 'T': ['name2'], '+A': ['name2']})
        self.assertTrue(var2flag == {'A': ['64'], 'T': ['128'], '+A': ['128']})

    def test1_7(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        posf = pf.Position_filter(None,None,None,None)

        bases = "*N"
        positions = "48,49"
        qnames = "name1,name2"
        flags = "64,128"
        ref = "T"
        var2nm, var2pos, var2qname, var2flag = posf.parse_bases(bases, positions, qnames, flags, ref)

        self.assertTrue(var2nm == {'N': 1})
        self.assertTrue(var2pos == {'N': ['49']})
        self.assertTrue(var2qname == {'N': ['name2']})
        self.assertTrue(var2flag == {'N': ['128']})

    def test2_1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        ref_genome = cur_dir + "/../database/GRCh37/GRCh37.fa"
        posf = pf.Position_filter(ref_genome,"samtools","-q 20 -B -Q15 -d 10000000 --output-BP --output-QNAME --output-extra FLAG", 1)
        bam = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        reg = "chr1:12345678-12345678"

        with open(os.devnull, 'w') as FNULL:
            l_ret = posf.call_mpileup(reg, bam, FNULL)

        self.assertTrue(l_ret == None)


    def test3_1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        posf = pf.Position_filter(None,None,None,None)

        ref = "T"
        alt = "A"
        ret = posf.get_alt_pileup_key(ref, alt)

        self.assertTrue(ret == "A")

    def test3_2(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        posf = pf.Position_filter(None,None,None,None)

        ref = "T"
        alt = "TA"
        ret = posf.get_alt_pileup_key(ref, alt)

        self.assertTrue(ret == "+A")

    def test3_3(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        posf = pf.Position_filter(None,None,None,None)

        ref = "TA"
        alt = "T"
        ret = posf.get_alt_pileup_key(ref, alt)

        self.assertTrue(ret == "-A")

    def test3_4(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        posf = pf.Position_filter(None,None,None,None)

        ref = "T"
        alt = "TAA"
        ret = posf.get_alt_pileup_key(ref, alt)

        self.assertTrue(ret == "+AA")

    def test3_5(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        posf = pf.Position_filter(None,None,None,None)

        ref = "TAA"
        alt = "T"
        ret = posf.get_alt_pileup_key(ref, alt)

        self.assertTrue(ret == "-AA")

    def test3_6(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        posf = pf.Position_filter(None,None,None,None)

        ref = "TT"
        alt = "AA"
        ret = posf.get_alt_pileup_key(ref, alt)

        self.assertTrue(ret == "A")

    def test3_7(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        posf = pf.Position_filter(None,None,None,None)

        ref = "TTT"
        alt = "AAA"
        ret = posf.get_alt_pileup_key(ref, alt)

        self.assertTrue(ret == "A")

    def test3_8(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        posf = pf.Position_filter(None,None,None,None)

        ref = "TT"
        alt = "AAA"
        ret = posf.get_alt_pileup_key(ref, alt)

        self.assertTrue(ret == None)

    def test4_1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        posf = pf.Position_filter(None,None,None,None)

        cigar = [(4,10),(0,40)]
        l_pos, r_pos = posf.get_cigar_size(cigar)

        self.assertTrue(l_pos == 10)
        self.assertTrue(r_pos == 0)

    def test4_2(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        posf = pf.Position_filter(None,None,None,None)

        cigar = [(4,20),(0,40),(4,10)]
        l_pos, r_pos = posf.get_cigar_size(cigar)

        self.assertTrue(l_pos == 20)
        self.assertTrue(r_pos == 10)

    def test4_3(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        posf = pf.Position_filter(None,None,None,None)

        cigar = [(0,40),(5,10)]
        l_pos, r_pos = posf.get_cigar_size(cigar)

        self.assertTrue(l_pos == 0)
        self.assertTrue(r_pos == 10)

    def test4_4(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        posf = pf.Position_filter(None,None,None,None)

        cigar = [(0,40)]
        l_pos, r_pos = posf.get_cigar_size(cigar)

        self.assertTrue(l_pos == 0)
        self.assertTrue(r_pos == 0)

    def test5_1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        ref_genome = cur_dir + "/../database/GRCh37/GRCh37.fa"
        posf = pf.Position_filter(ref_genome,"samtools","-q 20 -B -Q15 -d 10000000 --output-BP --output-QNAME --output-extra FLAG", 1)
        bam = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_vcf = cur_dir + "/../data/5929_small_mutation_result_test21.txt"
        output = cur_dir + "/../data/5929_small_mutation_result_test21_posout.txt"

        with open(os.devnull, 'w') as FNULL:
            posf.filter_main_vcf(in_vcf, bam, output, "5929_tumor", "5929_control")

        self.assertTrue(True)


