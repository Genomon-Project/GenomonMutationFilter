import sys
import os
import re
import pysam
import argparse
import logging
from auto_vivification import auto_vivification

#
# Class definitions
#
class indel_filter:

    def __init__(self, search_length, min_depth, min_mismatch, af_thres, base_qual_thres, neighbor):
        self.search_length = search_length
        self.min_depth = min_depth
        self.min_mismatch = min_mismatch
        self.af_thres = af_thres
        self.base_qual_thres = str(base_qual_thres)
        self.neighbor = neighbor
        self.target = re.compile( '([\+\-])([0-9]+)([ACGTNacgtn]+)' )
        self.remove_chr = re.compile( '\^.' )
    

    def write_result_file(self, line, file_handle, max_mismatch_count, max_mismatch_rate):
        print >> file_handle, (line +"\t"+ str(max_mismatch_count) +"\t"+ str(max_mismatch_rate)) 


    def filter(self, in_mutation_file, in_bam, output):

        samfile = pysam.Samfile(in_bam, "rb")
        
        chrIndex = 0
        srcfile = open(in_mutation_file,'r')
        header = srcfile.readline()  
        headerlist = header.split('\t')
        for colname in headerlist:
            if (colname == 'Chr'): 
                break
            chrIndex += 1
        
        hResult = open(output,'w')
        
        # print header
        print >> hResult, (header.rstrip('\n')
                       + "\tmismatch_count\tmismatch_rate")

        for line in srcfile:
            line = line.rstrip()
            itemlist = line.split('\t')
    
            # input file is annovar format (not zero-based number)
            chr = itemlist[chrIndex]
            start = (int(itemlist[chrIndex + 1]) - 1)
            end = int(itemlist[chrIndex + 2])
            ref = itemlist[chrIndex + 3]
            alt = itemlist[chrIndex + 4]
            
            chr = chr.replace('chr', '') 
            
            max_mismatch_count = 0
            max_mismatch_rate = 0
            
            # if (ref == '-' or alt == '-'):
            #     self.write_result_file(line, hResult, '---', '---')
            #     continue
            region = chr +":"+str(max(0, (start - self.search_length + 1))) +"-"+ str(end + self.search_length)

            ####
            # print region
            for mpileup in pysam.mpileup( '-BQ', '0', '-d', '10000000', "-q", self.base_qual_thres, "-r", region, in_bam ):
                # print mpileup.rstrip()

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
                indel = auto_vivification()

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
                            start_pos = mp_list[ 1 ]
                            
                            mismatch_count = ( indel[ type ][ key ][ '-' ] + indel[ type ][ key ][ '+' ])
                            mismatch_rate = (float(mismatch_count) / float(depth))

                            if mismatch_rate >= max_mismatch_rate:
                                start_pos = int(start_pos)
                                end_pos   = int(start_pos)

                                if (type == '-'):
                                    start_pos = int(start_pos) + 1
                                    end_pos = int(start_pos) + len((key.split('\t'))[3]) - 1 

                                # print "m: " + str(start_pos) +"-"+ str(end_pos)
                                # print "o: " + str(start) +"-"+ str(end)
                                # print mismatch_count
                                # print mismatch_rate

                                if ((start_pos - self.neighbor <= int(start) + 1 and int(start) + 1 <= self.neighbor + end_pos) 
                                  or(start_pos - self.neighbor <= int(end)      and  int(end)       <= self.neighbor + end_pos)): 

                                    max_mismatch_count = mismatch_count
                                    max_mismatch_rate  = mismatch_rate

            ####
            # print "mmc: " + str(max_mismatch_count)
            # print "mm:  " + str(self.min_mismatch)
            if(max_mismatch_count <= self.min_mismatch or max_mismatch_rate <= self.af_thres):
                self.write_result_file(line, hResult, max_mismatch_count, max_mismatch_rate)

        hResult.close()


