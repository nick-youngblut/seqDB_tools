#!/usr/bin/env python

#--- Option parsing ---#
"""
gasic_seqDB_batch.py: run GASiC with multiple query metagenomes

indexFile -- bowtie2 index file 
metaFile -- metadata file with mg-rast or sra metagenome sequence links

Usage:
  gasic_seqDB_batch.py <indexFile> <metaFile>
  gasic_seqDB_batch.py -h | --help
  gasic_seqDB_batch.py --version

Options:
  -h --help     Show this screen.
  --version     Show version.
  --speed=<kn>  Speed in knots [default: 10].
  --moored      Moored (anchored) mine.
  --drifting    Drifting mine.
"""

from docopt import docopt

if __name__ == '__main__':
    args = docopt(__doc__, version='0.1')


#--- Package import ---#
import os
import sys

scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)
libDir = os.path.join(scriptDir, '../lib/gasic-r16/')
sys.path.append(libDir)

#import core.tools

import gasicBatch.objects as gobj



#--- Option error testing ---#

# readable files
def fileExists(fileName):
    if not os.path.isfile(fileName):
        raise IOError('"{0}" does not exist'.format( fileName ) )

fileExists(args['<metaFile>'])

# bowtie2 in path?



#--- Main ---#


# reading 
## metaFile
metaFileO = gobj.metaFile(fileName=args['<metaFile>'], seqDB='MGRAST', sep='\t')


# each metagenome
for ID in metaFileO.iterByRow():
    # downloading     
    metaFileO.download( ID )

    sys.exit()
