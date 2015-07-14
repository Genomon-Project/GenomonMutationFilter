import sys
import os
import re
import pysam
import scipy.special
import argparse
import logging
import math
from AutoVivification import AutoVivification

#
# Class definitions
#
class IndelFilter:

    def __init__(self, mutationFileInfo):
        self.mfi = mutationFileInfo
        self.search_length = int(10)
        self.min_depth = int(10)
        self.max_depth = int(1000)
        self.min_mismatch = int(3)
        self.target = re.compile( '([\+\-])([0-9]+)([ACGTNacgtn]+)' )
        self.remove_chr = re.compile( '\^.' )
        

    ############################################################
    def set_filter_info(self, mutation_line_info, max_mismatch_count, max_mismatch_rate):
        if(max_mismatch_count <= 0):
            mutation_line_info[2].append(0)
            mutation_line_info[3].append("0\t0")
            
        else:
            if(max_mismatch_count >= self.min_mismatch):
                mutation_line_info[2].append(1)
            else:
                mutation_line_info[2].append(0)

            mutation_line_info[3].append(str(max_mismatch_count)+"\t"+str(max_mismatch_rate))


    ############################################################
    def filter(self, in_bam):
        logging.info( 'filter start')

        samfile = pysam.Samfile(in_bam, "rb")
        
        for mutation_line_info in self.mfi['line_info']:
            chr = mutation_line_info[1]['chr']
            start = mutation_line_info[1]['start']
            end = mutation_line_info[1]['end']
        
            max_mismatch_count = 0
            max_mismatch_rate = 0
            
            if (mutation_line_info[1]['ref'] == '-' or mutation_line_info[1]['alt'] == '-'):
                self.set_filter_info(mutation_line_info, max_mismatch_count, max_mismatch_rate)
                continue

            if samfile.count(chr, start, (start+1)) >= self.max_depth:
                self.set_filter_info(mutation_line_info, max_mismatch_count, max_mismatch_rate)
                continue
            
            region = chr +":"+str(max(0, (start - self.search_length))) +"-"+ str(end + self.search_length)

            ####
            # print region
            for mpileup in pysam.mpileup( '-BQ', '0', '-d', '10000000', "-r", region, in_bam ):
                
                #
                # Prepare mpileup data
                #
                mp_list = str( mpileup.translate( None, '\n' ) ).split( '\t' )
                mp_list_len = len( mp_list )
                coordinate = mp_list[ 0:3 ]

                #
                # skip if depth is 0
                #
                if mp_list[ 3 ] == '0' or ( mp_list_len > 6 and mp_list[ 6 ] == '0' ):
                    continue
                
                #
                # skip if depth < min_depth
                #
                if int(mp_list[ 3 ]) < self.min_depth:
                    continue

                #
                depth = mp_list[ 3 ]
                read_bases = mp_list[ 4 ]
                qual_list = mp_list[ 5 ]

                #
                # position id,
                # mpileup output 4th row(number of read covering the site),
                # 5th row(read bases),
                # 6th row(base quality)
                #
                indel = AutoVivification()

                #
                # Look for deletion/insertion and save info in 'indel' dictionary
                #
                #   ([\+\-])[0-9]+[ACGTNacgtn]+
                #
                # m.group(1): + or - (deletion/insertion)
                # m.group(2): number of deletion/insertion
                # m.group(3): nucleotides
                #
                deleted = 0
                iter = self.target.finditer( read_bases )
                for m in iter:
                    site = m.start()
                    type = m.group( 1 )
                    num = m.group( 2 )
                    bases = m.group( 3 )[ 0:int( num ) ]
                    if bases.islower():
                        strand = ( '-', '+' )
                    else:
                        strand = ( '+', '-' )

                    key = '\t'.join( coordinate + [ bases.upper() ] )
                    if type in indel and key in indel[ type ]:
                        indel[ type ][ key ][ strand[ 0 ] ] += 1
                    else:
                        indel[ type ][ key ][ strand[ 0 ] ] = 1
                        indel[ type ][ key ][ strand[ 1 ] ] = 0

                    read_bases = read_bases[ 0:site - deleted ] + read_bases[ site + int(num) + len( num ) + 1 - deleted: ]
                    deleted += 1 + len( num ) + int( num )

                #
                # Remove '^.' and '$'
                #
                read_bases = self.remove_chr.sub( '', read_bases )
                read_bases = read_bases.translate( None, '$' ) 

                #
                # Error check
                #
                if len( read_bases ) != len( qual_list ):
                    logging.error( "mpileup data is not good: {0}, {1}".format( mpileup, read_bases ) )
                    exit(1)

                #
                for type in ( '+', '-' ):
                    if type in indel:
                        for key in indel[ type ].keys():

                            mismatch_count = ( indel[ type ][ key ][ '-' ] + indel[ type ][ key ][ '+' ])
                            mismatch_rate = (mismatch_count / int(depth))

                            if mismatch_count >= max_mismatch_rate:
                                max_mismatch_count = mismatch_count
                                max_mismatch_rate  = mismatch_rate

            ####
            self.set_filter_info(mutation_line_info, max_mismatch_count, max_mismatch_rate)


        logging.info( 'filter end')

