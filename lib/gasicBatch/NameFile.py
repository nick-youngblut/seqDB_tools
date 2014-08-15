"""Classes for nameFile"""

import os
import sys


class NameFile(object):
    """nameFile = 1 or 2 tab-delimited columns
    1st column = refFile (refernece sequence)
    2nd column = indexFile (read mapper index)
    """

    def __init__(self, nameFile):
        self.read_nameFile(nameFile)


    def read_nameFile(self, nameFile):
        """Reading in nameFile"""
        with open(nameFile, 'rb') as fh:
            self.names = []
            lineCount = 0
            for l in fh.readlines():
                row = l.rstrip().split('\t')
                try: 
                    row[1]
                except: 
                    row.append(None)
                self.names.append(dict(refFile=row[0], 
                                   indexFile=row[1], 
                                   rowIndex=lineCount))
                lineCount += 1


    def iter_names(self):
        """Interate through all names
        
        Output:
        [{refFile=refFile, indexFile=indexFile},]
        """
        for row in self.names:
            yield row


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


    def get_names_row(self, i):
        return self.names[i]

    def len(self):
        return len(self.names)
