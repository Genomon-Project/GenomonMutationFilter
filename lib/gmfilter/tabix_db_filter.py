import tabix
import sys
import os
import re
import logging

#
# Class definitions
#
class tabix_db_filter:


    def __init__(self, tabix_db_file):
        self.tabix_db_file = tabix_db_file


    ############################################################
    def filter(in_mutation_file, output):
    
        tb = tabix.open(self.tabix_db_file)

        # tabix open
        chrIndex = 0
        srcfile = open(in_mutation_file,'r')
        header = srcfile.readline()  
        headerlist = header.split('\t')
        for colname in headerlist:
            if (colname == 'Chr'): 
                break
            chrIndex += 1
        
        for line in srcfile:
            line = line.rstrip()
            itemlist = line.split('\t')
    
            # input file is annovar format (not zero-based number)
            chr = itemlist[chrIndex]
            start = (int(itemlist[chrIndex + 1]) - 1)
            end = int(itemlist[chrIndex + 2])
            
            # tabix databese is a zero-based number 
            try:
                records = tb.query(chr,start,end)
            except tabix.TabixError:
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                logging.error( ("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) ) )
                continue
          
            ####
            print line + "\tobj"   

