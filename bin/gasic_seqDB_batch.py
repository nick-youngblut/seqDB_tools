#!/usr/bin/env python

#--- Option parsing ---#
"""
gasic_seqDB_batch.py: run GASiC with multiple query metagenomes

nameFile -- tab-delim table (refFile<tab>indexFile).
            refFile = fasta of reference sequences used to make bowtie2 index
            indexFile = read mapper index file
metaFile -- metadata file with mg-rast or sra metagenome sequence links

Usage:
  gasic_seqDB_batch.py [options] <nameFile> <metaFile>
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
from collections import defaultdict

from Bio import SeqIO
import pysam
import numpy as np

scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)
libDir = os.path.join(scriptDir, '../lib/gasic-r16/')
sys.path.append(libDir)

import gasicBatch.MetaFile as MetaFile
import gasicBatch.NameFile as NameFile
from gasicBatch.ReadMapper import ReadMapper
from gasicBatch.ReadSimulator import ReadSimulator
from gasicBatch.CorrectAbundances import CorrectAbundances


#--- Option error testing ---#

# readable files?
def fileExists(fileName):
    if not os.path.isfile(fileName):
        raise IOError('"{0}" does not exist'.format( fileName ) )

fileExists(args['<metaFile>'])
fileExists(args['<nameFile>'])



#--- Main ---#
# nameFile loading
nameF = NameFile.NameFile(args['<nameFile>'])


# metaFile loading
if args['--seqDB'] == 'MGRAST':
    metaF = MetaFile.MetaFile_MGRAST(fileName=args['<metaFile>'], 
                                     stages=args['--stages'], 
                                     sep='\t')
elif args['--seqDB'] == 'SRA':
    metaF = MetaFile.MetaFile_SRA(fileName=args['<metaFile>'],
                                  sep='\t')



# each metagenome (getting from certain seqDB)
for row in metaF.iterByRow():

    # making temp directory and chdir
    #tmpdir = tempfile.mkdtemp()

    # unpack
    metagenomeID = row['ID']
    
    # downloading     
    ret = metaF.download( ID=metagenomeID )
    if ret is None: 
        sys.stderr.write('No read file downloaded for Metagenome {0}. Moving to next metagenome'.format(metagenomeID))
        continue
            
    # determine read stats
    ## skipping if no downloaded file found
    metaF.getReadStats(fileFormat='fasta')


    # read mapping
    ## creating object for specific mapper
    mapper = ReadMapper.getMapper('bowtie2')    # factory class

    ## calling mapper for each index file
    for name in nameF.iter_names():
        # unpack
        indexFile = name.get_indexFile()
        readFile = metaF.get_readFile()    
        # mapping
        samFile = mapper.run_mapper(indexFile, readFile, fileType='fasta')   # f = bowtie2 flag for file type
        name.set_refSamFile(samFile)
    

    # similarity estimation
    ## select simulator
    simulator = ReadSimulator.getSimulator('mason')

    ## foreach refFile: simulate reads 
    simReadsFiles = dict()
    for name in nameF.iter_names():
        # unpack
        fastaFile = name.get_fastaFile()
        indexFile = name.get_indexFile()
        readFile = metaF.get_readFile()    

        # read simulation 
        outDir = os.path.abspath(os.path.curdir)
        (simReadsFile,fileType) = simulator.run_simulator(fastaFile, outDir=outDir)
        
        # find out how many reads were generated
        ### Attention: Here we assume that all files contain the same number of read and are stored in fastq format
        num_reads = len( [ True for i in SeqIO.parse(simReadsFile, fileType) ] )

        ## convert readFile to fasta if needed
        if fileType.lower() == 'fastq':
            fastaFile = os.path.splitext(simReadsFile)[0] + '.fna'
            SeqIO.convert(simReadsFile, 'fastq', fastaFile, 'fasta')
            simReadsFile = fastaFile
        
        # saving reads file names in names class
        name.set_simReadsFile(simReadsFile)


    # pairwise mapping of the simulated reads from each ref to all references 
    ## resulting SAM file names saves as numpy array
    n_refs = nameF.len()    
    simSamFiles = np.array([['' for i in range(n_refs)] for j in range(n_refs)], dtype=object)
    for i in range(n_refs):
        # getting simulated reads from first query reference taxon 
        nameQuery = nameF.get_name(i)        
        simReadsFile = nameQuery.get_simReadsFile()

        for j in range(n_refs):
            # getting index file of subject for mapping to subject
            nameSubject = nameF.get_name(j)
            indexFile = nameSubject.get_indexFile()

            # mapping
            simSamFile = mapper.run_mapper(indexFile, simReadsFile, fileType='fasta')   # f = bowtie2 flag for file type

            # adding samSimFile to names file
            #nameQuery.add_simSamFile(i, j, simSamFile)            
            simSamFiles[i,j] = simSamFile
    
    # parse SAM files to create numpy array
    mappedReads = np.zeros((n_refs, n_refs, num_reads))
    for i in range(n_refs):
        for j in range(n_refs):            
            # count the reads in i mapping to subject j
            samfh = pysam.Samfile(simSamFiles[i,j], "r")
            mappedReads[i,j,:] = np.array( [int(not read.is_unmapped) for read in samfh] )


    # save the similarity matrix
    matrixOutFile = metagenomeID + '_simMtx'
    np.save(matrixOutFile, mappedReads)
    matrixOutFile += '.npy'
    sys.stderr.write('Wrote similarity matrix: {}'.format(matrixOutFile))
            

    # similarity correction 
    ## input: matrix & original reads -> ref sam file
    refSamFiles = [name.get_refSamFile() for name in nameF.iter_names()]
    CorAbund = CorrectAbundances()
    
    result = CorAbund.similarityCorrection(refSamFiles, matrixOutFile, 3)
    

    print result; sys.exit()
    
    
    # clean up tempdir
    #shutil.rmtree(tmpdir)

    # debug
    sys.exit()
