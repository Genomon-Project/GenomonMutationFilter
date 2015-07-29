import tabix
import sys
import os
import re
import logging

#
# Class definitions
#
class simple_repeat_filter:


    def __init__(self, simple_repeat_db):
        self.simple_repeat_db = simple_repeat_db


    def write_result_file(self, line, file_handle, db_pos, db_seq):
        print >> file_handle, (line + "\t" +db_pos+ "\t" +db_seq)
        

    ############################################################
    def filter(self, in_mutation_file, output):
    
        tb = tabix.open(self.simple_repeat_db)

        # tabix open
        chrIndex = 0
        srcfile = open(in_mutation_file,'r')
        header = srcfile.readline()  
        headerlist = header.split('\t')
        for colname in headerlist:
            if (colname == 'Chr'): 
                break
            chrIndex += 1
        
        hResult = open(output,'w')

        print >> hResult, (header.rstrip('\n')
                       + "\tsimple_repeat_pos\tsimple_repeat_seq")

        for line in srcfile:
            line = line.rstrip()
            itemlist = line.split('\t')
    
            # input file is annovar format (not zero-based number)
            chr = itemlist[chrIndex]
            start = (int(itemlist[chrIndex + 1]) - 1)
            end = int(itemlist[chrIndex + 2])

            # chr = chr.replace('chr', '') 

            # tabix databese is a zero-based number 
            position_array = []
            sequence_array = []
            try:
                records = tb.query(chr,start,end)
                for record in records:
                    position_array.append(record[0] +":"+ record[1] +"-"+ record[2])
                    sequence_array.append(record[15])
                    
            except tabix.TabixError:
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                logging.error( ("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) ) )
                self.write_result_file(line, hResult, '---', '---')
                continue
          
            ####
            db_pos = ";".join(map(str,position_array[0:5]))
            db_seq = ";".join(map(str,sequence_array[0:5]))
            if db_pos == "": db_pos = "---"
            if db_seq == "": db_seq = "---"
            self.write_result_file(line, hResult, db_pos, db_seq)

        hResult.close()


