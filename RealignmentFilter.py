import sys
import os
import re
import pysam
import scipy.special
import argparse
import logging
import subprocess
import math
from AutoVivification import AutoVivification

#
# Class definitions
#
class RealignmentFilter:

    def __init__(self, mutationFileInfo):
        self.mfi = mutationFileInfo
        self.max_depth = int(1000)
        self.max_mismatch = int(5)

        self.reference_genome = '/home/ogawaprj/ngs/ref/GRCh37-lite_PCAWG_bwa-0.7.10/GRCh37-lite_PCAWG.fa'
        self.split_refernece_thres = int(200)
        self.window = int(200)

        self.tumor_min_mismatch = int(3)
        self.normal_max_mismatch = int(1)

        self.complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        self.blat_cmds = ['/usr/local/bin/blat', '-stepSize=5', '-repMatch=2253']
        
        
    ############################################################
    def set_filter_info(self, mutation_line_info, tumor_ref,tumor_alt,tumor_other, normal_ref, normal_alt, normal_other):
            
        if(tumor_alt >= self.tumor_min_mismatch and normal_alt <= self.normal_max_mismatch):
            mutation_line_info[2].append(0)
        else:
            mutation_line_info[2].append(1)
           
        mutation_line_info[3].append(str(tumor_ref)  +"\t"+ str(tumor_alt)  +"\t"+ str(tumor_other) +"\t" +
                                     str(normal_ref) +"\t"+ str(normal_alt) +"\t"+ str(normal_other))


    ############################################################
    def makeTwoReference(self, chr,start,end,ref,alt, outputFilePath):

        hOUT = open(outputFilePath, 'w')
        
        seq = ""
        label = ','.join([chr, str(start), str(end), ref, alt])
        range = chr + ":" + str(int(start) - self.window + 1) +"-"+ str(int(end) + self.window)
        for item in pysam.faidx(self.reference_genome, range):
            if item[0] == ">": continue
            seq = seq + item.rstrip('\n').upper()

        print >> hOUT, '>' + label + "_ref"
        print >> hOUT, seq
        # print '>' + label + "_ref"
        # print seq

        # for insertion
        if ref == "-":   seq = seq[0:(self.window + 1)] + alt + seq[-self.window:]
        # for deletion
        elif alt == "-": seq = seq[0:self.window] + seq[-self.window:]
         # for SNV
        else:            seq = seq[0:self.window] + alt + seq[-self.window:]

        print >> hOUT, '>' + label + "_alt"
        print >> hOUT, seq
        # print '>' + label  + "_alt"
        # print seq

        hOUT.close()


    ############################################################
    def extractRead(self, bamfile, chr,start,end, outputFilePath):

        hOUT = open(outputFilePath, 'w')

        for read in bamfile.fetch(chr,start,end):
          
            # get the flag information
            flags = format(int(read.flag), "#014b")[:1:-1]

            # skip unmapped read 
            if flags[2] == "1" or flags[3] == "1": continue 

            # skip supplementary alignment
            if flags[8] == "1" or flags[11] == "1": continue

            # skip duplicated reads
            if flags[10] == "1": continue

            # chr_current = bamfile.getrname(read.tid)
            # pos_current = int(read.pos + 1)
            # dir_current = ("-" if flags[4] == "1" else "+")
            # chr_pair = bamfile.getrname(read.rnext)
            # pos_pair = int(read.pnext + 1)
            # dir_pair = ("-" if flags[5] == "1" else "+")
            
            chr_pair = bamfile.getrname(read.rnext)
            print chr_pair

            tempSeq = ""
            if flags[4] == "1":
                tempSeq = "".join(self.complement.get(base) for base in reversed(str(read.seq)))
            else:
                tempSeq = read.seq

            # the first read
            if flags[6] == "1":
                print >> hOUT, '>' + read.qname + '/1'
                print >> hOUT, tempSeq
                # print '>' + read.qname + '/1'
                # print tempSeq
            else:
                print >> hOUT, '>' + read.qname + '/2'
                print >> hOUT, tempSeq
                # print '>' + read.qname + '/2'
                # print tempSeq

        hOUT.close()


    ############################################################
    def checkScore(self, align):

        tempScore = 100
        for i in range(0, len(align)):
            tempScore = min(tempScore, align[i][0])

        return(tempScore)


    ############################################################
    def summarizeRefAlt(self, inputPsl):
        
        hIN = open(inputPsl, 'r')

        numOther = 0
        numAlt = 0
        numRef = 0

        tempID = ""
        tempAlt = []
        tempRef = []

        for line in hIN:
            F = line.rstrip('\n').split('\t')
            if F[0].isdigit() == False: continue

            # remove the read pair num info ("/1", "/2") 
            F[9] = F[9][0:-2]
            if tempID != F[9]:
                if tempID != "":
                
                    ####
                    tempAltNM = self.checkScore(tempAlt)
                    tempRefNM = self.checkScore(tempRef)
                    if tempAltNM >= self.max_mismatch and tempRefNM >= self.max_mismatch: numOther += 1
                    elif tempAltNM <  tempRefNM: numAlt += 1
                    elif tempRefNM <= tempAltNM: numRef += 1

                tempID = F[9]
                tempAlt = []
                tempRef = []

            tNM = int(F[10]) - int(F[0]) + int(F[5]) + int(F[7])
            tpos = int(F[15])
            tdir = F[8]

            if F[13][-3:] == "alt":
                tempAlt.append((tNM, tpos, tdir))
            elif F[13][-3:] == "ref":
                tempRef.append((tNM, tpos, tdir))

        ####
        tempAltNM = self.checkScore(tempAlt)
        tempRefNM = self.checkScore(tempRef)
        if tempAltNM >= self.max_mismatch and tempRefNM >= self.max_mismatch: numOther += 1
        elif tempAltNM <  tempRefNM: numAlt += 1
        elif tempRefNM <= tempAltNM: numRef += 1

        hIN.close()

        return([numRef, numAlt, numOther])
               

    ############################################################
    def filter(self, in_tumor_bam, in_normal_bam, outputFilePath):
        logging.info( 'filter start')

        tumor_samfile = pysam.Samfile(in_tumor_bam, "rb")
        normal_samfile = pysam.Samfile(in_normal_bam, "rb")
        
        for mutation_line_info in self.mfi['line_info']:
            chr = mutation_line_info[1]['chr']
            start = mutation_line_info[1]['start']
            end = mutation_line_info[1]['end']
            ref = mutation_line_info[1]['ref']
            alt = mutation_line_info[1]['alt']

            if tumor_samfile.count(chr, start, (start+1)) >= self.max_depth:
                self.set_filter_info(mutation_line_info,0,0,0,0,0,0)
                continue 

            if normal_samfile.count(chr, start, (start+1)) >= self.max_depth:
                self.set_filter_info(mutation_line_info,0,0,0,0,0,0)
                continue 

            self.makeTwoReference(chr,start,end,ref,alt,outputFilePath + "/tmp.refalt.fa")
            
            # extract short reads from tumor sequence data around the candidate
            self.extractRead(tumor_samfile,chr,start,end,outputFilePath + "/tmp.fa")

            # alignment tumor short reads to the reference and alternative sequences
            FNULL = open(os.devnull, 'w')
            subprocess.call(self.blat_cmds + [outputFilePath + "/tmp.refalt.fa", outputFilePath + "/tmp.fa", outputFilePath + "/tmp.psl"], 
                            stdout = FNULL, stderr = subprocess.STDOUT)
            FNULL.close()

            # summarize alignment results
            tumor_ref, tumor_alt, tumor_other = self.summarizeRefAlt(outputFilePath + "/tmp.psl")

            # extract short reads from normal sequence data around the candidate
            self.extractRead(normal_samfile,chr,start,end,outputFilePath + "/tmp.fa")

            # alignment normal short reads to the reference and alternative sequences
            FNULL = open(os.devnull, 'w')
            subprocess.call(self.blat_cmds + [outputFilePath + "/tmp.refalt.fa", outputFilePath + "/tmp.fa", outputFilePath + "/tmp.psl"], 
                            stdout = FNULL, stderr = subprocess.STDOUT)
            FNULL.close()

            # summarize alignment results
            normal_ref, normal_alt, normal_other = self.summarizeRefAlt(outputFilePath + "/tmp.psl")
    
            self.set_filter_info(mutation_line_info, tumor_ref,tumor_alt,tumor_other, normal_ref, normal_alt, normal_other)
            

    logging.info( 'filter end')

