#!/usr/bin/env python

#--- Option parsing ---#
"""
gasic_seqDB_batch.py: run GASiC with multiple query metagenomes

refFile -- fasta of reference sequence (used to create bowtie2 index file). Use '%s' to match multiple files.
indexFile -- bowtie2 index file. Use '%s' to match multiple files.
metaFile -- metadata file with mg-rast or sra metagenome sequence links

Usage:
  gasic_seqDB_batch.py [options] <nameFile> <refFile> <indexFile> <metaFile>
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

import gasicBatch.MetaFile as MetaFile
from gasicBatch.ReadMapper import ReadMapper
from gasicBatch.ReadSimulator import ReadSimulator

#--- Option error testing ---#

# readable files
def fileExists(fileName):
    if not os.path.isfile(fileName):
        raise IOError('"{0}" does not exist'.format( fileName ) )

fileExists(args['<metaFile>'])


#--- Main ---#
# metaFile loading
if args['--seqDB'] == 'MGRAST':
    metaF = MetaFile.MetaFile_MGRAST(fileName=args['<metaFile>'], 
                                     stages=args['--stages'], 
                                     sep='\t')
elif args['--seqDB'] == 'SRA':
    metaF = MetaFile.MetaFile_SRA(fileName=args['<metaFile>'],
                                  sep='\t')


# each metagenome
for row in metaF.iterByRow():

    # making temp directory and chdir
    #tmpdir = tempfile.mkdtemp()
    
    # downloading     
    ret = metaF.download( ID=row['ID'] )
    if ret is None: 
        sys.stderr.write('No read file downloaded for Metagenome {ID}. Moving to next metagenome'.format(**row))
        continue
            
    # determine read stats
    ## skipping if no downloaded file found
    metaF.getReadStats(fileFormat='fasta')

    # read mapping
    ## creating object for specific mapper
    mapper = ReadMapper.getMapper('bowtie2')    # factory class
    ## calling mapper
    indexFile = args['<indexFile>']
    readFile = metaF.get_readFile()    
    samFile = mapper.run_mapper(indexFile, readFile, fileType='fasta')   # f = bowtie2 flag for file type


    # similarity existimation
    ## generate reads for every reference genome

    ## find out how many reads were generated
    ### Attention: Here we assume that all files contain the same number of read and are stored in fastq format
    #num_reads = len( [ True for i in SeqIO.parse(sim_files[0],'fastq') ] )

    ## map the reads of every reference to all references

    ## parse SAM files to create numpy array
    
    # read simulation 
    #simulator = ReadSimulator.getSimulator('mason')
    #outDir = os.path.abspath(os.path.curdir)
    #simulator.run_simulator(args['<refFile>'], outDir=outDir)

    sys.exit()


    # creating similarity matrix

    # similarity correction

    # clean up tempdir
    #shutil.rmtree(tmpdir)

    # debug
    sys.exit()
