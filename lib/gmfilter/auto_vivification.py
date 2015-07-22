import os
import re
import scipy.special

#
# Class definitions
#

############################################################
class auto_vivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)

        except KeyError:
            value = self[item] = type(self)()

        return value

