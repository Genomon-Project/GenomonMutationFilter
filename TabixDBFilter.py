
import tabix
import sys
import os
import re
import logging
# from mongo_utility import mongo_utility 

class TabixDBFilter:


    def __init__(self):
        logging.info( 'create TabixDBFilter')


    def filter(self, mutationFileInfo, dbTabixFile):
    
        logging.info( 'database='+dbTabixFile)
        tb = tabix.open(dbTabixFile)

        # tabix open
        for mutation_line_info in mutationFileInfo['line_info']:

            # input file is annovar format (not zero-based number)
            chr = mutation_line_info[1]['chr']
            start = mutation_line_info[1]['start']
            end = mutation_line_info[1]['end']
            
            # tabix databese is a zero-based number 
            try:
                records = tb.query(chr,start,end)
            except tabix.TabixError:
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                logging.error( ("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) ) )
                raise
               #  continue
          
            for record in records:
                if record:
                    mutation_line_info[2].append(0)
                else:
                    mutation_line_info[2].append(1)
                break
                
