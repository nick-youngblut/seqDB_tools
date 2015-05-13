"""Classes for nameFile"""

import os
import sys
from collections import defaultdict

class NameFile(object):
    """nameFile = 1 or 2 tab-delimited columns
    1st column = refFile (refernece sequence)
    2nd column = indexFile (read mapper index)
    """

    def __init__(self, nameFile):
        """init
        Args:
        nameFile -- name of nameFile (string)
        """
        self.read_nameFile(nameFile)


    def read_nameFile(self, nameFile):
        """Reading in nameFile
        Args:
        nameFile -- name of nameFile        
        """
        with open(nameFile, 'rb') as fh:
            self.names = []
            lineCount = 0
            for l in fh.readlines():
                row = l.rstrip().split('\t')
                try: 
                    row[1]
                except: 
                    # assuming that the refFasta is the same as the index file
                    row.append(row[0]) 
                name = Name(row[0], row[1], lineCount)
                self.names.append(name)

                lineCount += 1


    def iter_names(self):
        """Interate through all names
        
        Return:
        Object of Names class
        """
        for Name in self.names:
            yield Name


    def set_names_keypair(self, i, k, v):
        """Adding a key-pair to dict in names (array of dicts)
        
        Args:
        i = row index
        k = key
        v = value
        """
        try:
            self.names[i]
        except IndexError as err:
            raise type(err)(err.message + '. Index: "{0}" not found'.format(i))

        self.names[i][k] = v
        return 1        


    def get_name(self, i):
        """
        Args:
        i -- index in names list
        """
        return self.names[i]

    def get_names(self):
        return self.names
        
    def len(self):
        return len(self.names)



class Name(object):
    """Class for individual reference file metadata"""

    def __init__(self, fastaFile, indexFile, rowIndex):
        """init
        Args:
        fastaFile -- fastaFile name
        indexFile -- indexFile name
        rowIndex -- rowIndex (depreciated?)
        """
        self.fastaFile = fastaFile
        self.indexFile = indexFile
        self.rowIndex = rowIndex


    def __repr__(self):
        tmp = '\n'.join( ['{0} : {1}'.format(k,v) for k,v in sorted(self.__dict__.items())] )
        return tmp        

    # getters         
    def get_fastaFile(self):
        return self.fastaFile

    def get_indexFile(self):
        return self.indexFile        

    def get_simReadsFile(self):
        return self.simReadsFile

    def get_simReadsFileType(self):
        return self.simReadsFileType

    def get_simReadsCount(self):
        return self.simReadsCount
        
    def get_refSamFile(self):
        return self.refSamFile

    # setters
    def set_refSamFile(self, refSamFile):
        self.refSamFile = refSamFile

    def set_simReadsFile(self, simReadsFile):
        self.simReadsFile = simReadsFile

    def set_simReadsFileType(self, fileType):
        self.simReadsFileType = fileType

    def set_simReadsCount(self, count):
        self.simReadsCount = count
        
    def add_simSamFile(self, i, j, samSimFile):
        """setting simulation sam files (pairwise mappings to subject (meta)genomes)

        Args:
        i = index i
        j = index j
        samSimFile = sam file name resulting from mapping
        """
        if not hasattr(self, 'samSimFiles'):
            self.samSimFiles = defaultdict(dict)

        self.samSimFiles[i][j] = samSimFile
