"""Classes for lastRunFile (output file from last run of gasic_batch)"""

import os
import sys
import pandas as pd

class lastRunFile(object):
    """lastRunFile: tab-delimited output (no header)
    produced by gasic batch script"""

    def __init__(self, inFile):
        self.fileName = inFile            
        self.tbl = pd.read_csv(inFile, sep='\t', header=None)

    def mgID_iter(self):
        """Iterator for all metagenome ID in table"""
        for i,row in self.tbl.iterrows():
            yield row[0]

    def get_mgIDs(self):
        """Get list of all metagenome IDs in table"""
        return list(self.tbl[0])

    def mgID_exists(self, mgIDs):
        """Does metagenome ID(s) exist in table?

        Args:
        mgIDs -- list of metagenome IDs to check for.
        
        Return:
        boolean of success
        """
        return self.tbl[0].isin(mgIDs)


    def mgEntries(self, mgIDs):
        """selecting rows corresponding to metagenomes listed in mgIDs.

        Args:
        mgIDs -- list of metagenome IDs to select rows

        Return:
        pandas DF of rows conresponding to metagenome
        """
        return self.tbl.loc[self.tbl[0].isin(mgIDs)]


        
        
