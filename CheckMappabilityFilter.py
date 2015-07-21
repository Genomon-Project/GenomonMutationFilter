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
class CheckMappabilityFilter:

    def __init__(self, mutationFileInfo):
        self.mfi = mutationFileInfo

        self.reference_genome = '/home/ogawaprj/ngs/ref/GRCh37-lite_PCAWG_bwa-0.7.10/GRCh37-lite_PCAWG.fa'
        self.window = int(50)
        self.blat_cmds = ['/home/ogawaprj/ngs/bin/blat_x86_64/blat', '-fine']
        

    ############################################################
    def getScore(self, align):
        if len(align) >= 1:
            align.sort(key=lambda x:x[0],reverse=True)
            return (align[0][0], align[0][1])
        return (0,0)


    ############################################################
    def makePosFasta(self, outputFilePath):

        hOUT = open(outputFilePath, 'w')

        for mutation_line_info in self.mfi['line_info']:
            chr = mutation_line_info[1]['chr']
            start = mutation_line_info[1]['start']
            end = mutation_line_info[1]['end']
            ref = mutation_line_info[1]['ref']
            alt = mutation_line_info[1]['alt']

            seq = ""
            label = ','.join([chr, str(start), str(end), ref, alt])
            range = chr + ":" + str(int(start) - self.window + 1) +"-"+ str(int(end) + self.window)
            for item in pysam.faidx(self.reference_genome, range):
                if item[0] == ">": continue
                seq = seq + item.rstrip('\n').upper()

            # for insertion
            if ref == "-":   seq = seq[0:(self.window + 1)] + alt + seq[-self.window:]
            # for deletion
            elif alt == "-": seq = seq[0:self.window] + seq[-self.window:]
            # for SNV
            else:            seq = seq[0:self.window] + alt + seq[-self.window:]

            print >> hOUT, '>' + label
            print >> hOUT, seq

        hOUT.close()


    ############################################################
    def summarizeAlt(self, inputPsl):
        
        hIN = open(inputPsl, 'r')

        resultDict = {}
        
        tempID = ""
        tempMyScore = int(0)
        tempElseScore = []
        ####
        for line in hIN:
            F = line.rstrip('\n').split('\t')
            if F[0].isdigit() == False: continue

            # remove the read pair num info ("/1", "/2") 
            if tempID != F[9]:
                if tempID != "":
                
                    ####
                    tempElseMaxScore, tempElseMaxRegion = self.getScore(tempElseScore)
                    resultDict[tempID] = str(tempMyScore) +'\t'+ str(tempElseMaxScore) +'\t'+ str(tempElseMaxRegion) 
        
                tempID = F[9]
                tempMyScore = int(0)
                tempElseScore = []
            
            tchr = F[13]
            tstart = F[15]
            tend = F[16]
            FF = tempID.split(','); 

            # print tempID +"\t"+ F[0] +"\t"+tchr + "\t" + str(tstart) + "\t" + str(tend)
            # print tempID +"\t"+ F[0] +"\t"+FF[0] + "\t" + str(int(FF[1]) - self.window) + "\t" + str(int(FF[2]) + self.window)

            if (str(FF[0]) == str(tchr) and (int(FF[1]) - self.window) == int(tstart) and (int(FF[2]) + self.window) == int(tend)):
                tempMyScore = F[0];
            else:
                tempElseScore = []
                tempElseScore.append((F[0], tchr +':'+ tstart +'-'+ tend))

        hIN.close()

        tempElseMaxScore, tempElseMaxRegion = self.getScore(tempElseScore)
        resultDict[tempID] = str(tempMyScore) +'\t'+ str(tempElseMaxScore) +'\t'+ tempElseMaxRegion 

        # for key in resultDict:
        #  print key + "\t" + resultDict[key]

        return resultDict


    ############################################################
    def set_map_ability_info(self, mapAbilityDict):
        
        # for key in mapAbilityDict:
        #   print key + "\t" + mapAbilityDict[key]

        for mutation_line_info in self.mfi['line_info']:
            chr = mutation_line_info[1]['chr']
            start = mutation_line_info[1]['start']
            end = mutation_line_info[1]['end']
            ref = mutation_line_info[1]['ref']
            alt = mutation_line_info[1]['alt']

            score = ''
            label = ','.join([chr, str(start), str(end), ref, alt])
            if label in mapAbilityDict: 
                scores = mapAbilityDict[label]
                tempMyScore, tempElseMaxScore, tempElseMaxRegion = scores.split('\t')

                if int(tempMyScore) == 0:
                    score = 'unaligned' +'\t'+ '---'
                    
                elif int(tempElseMaxScore) == 0:
                    score = '---' +'\t'+ '---'

                else:
                    score = str(tempElseMaxRegion) +'\t'+ str(int(tempMyScore) - int(tempElseMaxScore))

            else:
                score = 'unaligned' +'\t'+ '----'

            mutation_line_info[3].append(score)


    ############################################################
    def filter(self,  outputFilePath):
        logging.info( 'filter start')

        self.makePosFasta(outputFilePath + "/tmp.posalt.fa")

        # alignment alternative sequences to the human reference genome
        FNULL = open(os.devnull, 'w')
        subprocess.call(self.blat_cmds + [self.reference_genome, outputFilePath + "/tmp.posalt.fa", outputFilePath + "/tmp.psl"], 
                        stdout = FNULL, stderr = subprocess.STDOUT)
        FNULL.close()

        mapAbilityDict = self.summarizeAlt(outputFilePath + "/tmp.psl")
            
        self.set_map_ability_info(mapAbilityDict)

    logging.info( 'filter end')

