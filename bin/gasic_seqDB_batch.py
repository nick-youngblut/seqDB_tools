#!/usr/bin/env python

#--- Option parsing ---#
"""
gasic_seqDB_batch.py: run GASiC with multiple query metagenomes

indexFile -- bowtie2 index file 
metaFile -- metadata file with mg-rast or sra metagenome sequence links

Usage:
  gasic_seqDB_batch.py [options] <indexFile> <metaFile>
  gasic_seqDB_batch.py -h | --help
  gasic_seqDB_batch.py --version

Options:
  --seqDB=<seqDB>       Sequence database name ('MGRAST' or 'SRA'). [default: MGRAST]
  --stages=<stages>...  MG-RAST processing stages to attempt download of. 
  Each in list will be tried until successful download. [default: 200,150]
  --version     Show version.
  -h --help     Show this screen.
"""

from docopt import docopt
import re

if __name__ == '__main__':
    args = docopt(__doc__, version='0.1')
    args['--seqDB'] = args['--seqDB'].upper()
    args['--stages'] = re.split(r'\s*,\s*', args['--stages'])


#--- Package import ---#
import os
import sys
import tempfile
import shutil

scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)
libDir = os.path.join(scriptDir, '../lib/gasic-r16/')
sys.path.append(libDir)

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
if args['--seqDB'] == 'MGRAST':
    metaFileO = gobj.MetaFile_MGRAST(fileName=args['<metaFile>'], 
                                     stages=args['--stages'], 
                                     sep='\t')
elif args['--seqDB'] == 'SRA':
    metaFileO = gobj.MetaFile_SRA(fileName=args['<metaFile>'],
                                  sep='\t')


# each metagenome
for row in metaFileO.iterByRow():

    # making temp directory and chdir
    tmpdir = tempfile.mkdtemp()
    
    try:
        # downloading     
        metaFileO.download( ID=row['ID'] )
         
        # determine read stats
        metaFileO.getReadStats( fileName=metaFileO.outFile, fileFormat='fasta')

        # read mapping with gasic
        

        # read simulation and mapping

        # creating similarity matrix

        # similarity correction

    finally:
        # clean up tempdir
        shutil.rmtree(tmpdir)

    # debug
    sys.exit()
