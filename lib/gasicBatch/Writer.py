import sys
import os
import pandas as pd


class OutputWriter(object):
    """Writing functions for gasic batch"""

    def __init__(self, metagenome_ID, nCol=8, sep='\t'):
        """
        Args:
        mgID -- MGRAST metagenome ID
        nCol -- number of columns in table to write
        sep -- value separator
        """
        self.mgID = metagenome_ID
        self.nCol = nCol - 2
        self.sep = sep

    def writeValues(self, outvals):
        """Standard writing of row for output table.
        Args:
        outvals -- dict of output values
        TODO:
        make more flexible
        """        
        # metagenome_id, refSequences, n-mapped, corrected-abundacne, error, pval, mg_platform
        print '{mgID}\t{ref}\t{total}\t{mapped}\t{corr}\t{error}\t{pval}\t{mg_platform}'.format(**outvals)
        
    def lastRun(self, df):
        """If metagenome in last run output, write old output
        Args:
        df -- pandas dataframe of last run output
        """
        msg =  '  Metagenome "{}" in last-run file. Writing old output; moving to next metagenome\n\n'
        sys.stderr.write(msg.format(self.mgID))
        df.to_csv(sys.stdout, sep=self.sep, header=None, index=None)

    def noReadFile(self):
        """If the read file could not be downloaded"""
        msg = 'No read file downloaded for Metagenome {0}. Skipping metagenome\n\n'
        sys.stderr.write(msg.format(self.mgID))
        line = [self.mgID, 'ERROR:no_metagenome_read_file'] + ['NA'] * self.nCol
        print self.sep.join(line) 

    def platformUserSkip(self, platform):
        """If read from the sequencing platform should be skipped"""
        msg = ' The platform "{}" is in the --platform list. Skipping metagenome.\n\n'
        sys.stderr.write(msg.format(platform))
        line = [self.mgID, 'ERROR:undesired_platform'] + ['NA'] * self.nCol
        print self.sep.join(line)

    def platformUnknown(self, platform):
        """If the sequencing platform could not be determined"""
        msg = ' The platform "{}" could not be determined. Skipping metagenome.\n\n'
        sys.stderr.write(msg.format(platform))
        line = [self.mgID, 'ERROR:undetermined_platform'] + ['NA'] * self.nCol
        print self.sep.join(line)
        
    def simReadError(self):
        """If error(s) during read simulation"""
        msg = '\n  WARNING: Read simulation error for metagenome "{}". Skipping metagenome.\n\n'
        sys.stderr.write(msg.format(self.mgID))
        line = [self.mgID, 'ERROR:read_simulation_error'] + ['NA'] * self.nCol
        print self.sep.join(line)

    def readFileFormatConversionError(self):
        """If error(s) during conversion of read file format (e.g. fastq to fasta)"""
        msg = '\n  WARNING: error during read file format conversion for metagenome "{}". Skipping metagenome.\n\n'
        sys.stderr.write(msg.format(self.mgID))
        line = [self.mgID, 'ERROR:read_file_format_conversion_error'] + ['NA'] * self.nCol
        print self.sep.join(line)
        
        