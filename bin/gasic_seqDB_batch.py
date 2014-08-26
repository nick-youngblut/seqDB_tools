#!/usr/bin/env python

#--- Option parsing ---#
"""
gasic_seqDB_batch.py: run GASiC with multiple query metagenomes

nameFile -- tab-delim table (refFile<tab>indexFile).
            refFile = fasta of reference sequences used to make bowtie2 index
            indexFile = read mapper index file
metaFile -- metadata file with mg-rast or sra metagenome sequence links

Usage:
  gasic_seqDB_batch.py [options] <nameFile> <metaFile> [<stages>...]
  gasic_seqDB_batch.py -h | --help
  gasic_seqDB_batch.py --version

Options:
  <nameFile>         Name file (See Description).
  <metaFile>         Tab-delim table of metadata for the sequence database of interest.
  <stages>...        MG-RAST processing stages to attempt download of. 
                     Each in list will be tried until successful download. [default: 200,150]
  --seqDB=<seqDB>    Sequence database name ('MGRAST' or 'SRA'). [default: MGRAST]
  --version          Show version.
  -h --help          Show this screen.
  --debug            Debug mode

Description:
  Run the GASiC pipeline on >= reference sequences (>=1 nucleotid fasta files),
  with the query reads being those selected from a sequence database (eg., MGRAST).
  Selecting which (meta)genome reads are used for the GASiC pipeline is determined
  from the user-provide metaFile.
  
  <nameFile> -- Format: 1 or 2 columns; 1st column: 'reference_fasta';
                optional 2nd column: 'reference_index'
                The index is the mapper index file (eg., bowtie2 index),
                which is only needed if the mappers requires an index file.
"""

from docopt import docopt
import os, sys

if __name__ == '__main__':
    args = docopt(__doc__, version='0.1')
    args['--seqDB'] = args['--seqDB'].upper()
    if len(args['<stages>']) == 0:
        args['<stages>'] = [200, 150]

#--- Package import ---#
import os
import sys
import tempfile
import shutil
from collections import defaultdict
import multiprocessing as mp

from Bio import SeqIO
import pysam
import numpy as np
import parmap

scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)
libDir = os.path.join(scriptDir, '../lib/gasic-r16/')
sys.path.append(libDir)

import gasicBatch.MetaFile as MetaFile
import gasicBatch.NameFile as NameFile
from gasicBatch.ReadMapper import ReadMapper, PairwiseMapper
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
                                     stages=args['<stages>'], 
                                     sep='\t')
elif args['--seqDB'] == 'SRA':
    metaF = MetaFile.MetaFile_SRA(fileName=args['<metaFile>'],
                                  sep='\t')

# current working directory
origWorkDir = os.path.abspath(os.curdir)


# each metagenome (getting from certain seqDB)
for mg in metaF.iterByRow():

    ## unpack
    mgID = mg.get_ID()
    
    # making temp directory and chdir; all processes for metagenome done in temp dir
    if args['--debug'] == False:
        tmpdir = tempfile.mkdtemp()
        os.chdir(tmpdir)

    # downloading metagenome reads
    dw_success = mg.download( )
    if not dw_success: 
        sys.stderr.write('No read file downloaded for the metagenome "{0}".\n\tMoving to next metagenome\n'.format(mgID))
        continue
        
    # determine read stats
    ## skipping if no downloaded file found
    mg.get_readStats(fileFormat='fasta')

    
    #-- mapping metagenome reads to reference fasta(s) --#
    ## creating object for specific mapper
    mapper = ReadMapper.getMapper('bowtie2')    # factory class
   # mapper.set_paramsByReadStats(mg)
            
    ### getting iterable and args
    indexFiles = [name.get_indexFile() for name in nameF.iter_names()]
    ### calling mapper via parmap
    readFile = mg.get_readFile()
    ret = parmap.map(mapper, indexFiles, readFile, tmpdir,  processes=4)
    ### assigning sam file names to name instances
    ret = dict(ret)
    for name in nameF.iter_names():
        indexFile = name.get_indexFile()
        name.set_refSamFile(ret[indexFile])
    
    
    #for name in nameF.iter_names():
        # unpack
#        indexFile = name.get_indexFile()
#        readFile = mg.get_readFile()    
#        # mapping
#        samFile = mapper.run_mapper(indexFile, readFile)   # f = bowtie2 flag for file type (fasta)
#        name.set_refSamFile(samFile)
    

    #-- similarity estimation by pairwise mapping simulated reads --#
    ## select simulator
    simulator = ReadSimulator.getSimulator('mason')
    ## setting params based on metagenome read stats & platform
    platform, simParams = simulator.get_paramsByReadStats(mg)
    
    ## foreach refFile: simulate reads 
    #simReadsFiles = dict()
    #num_reads = None
    
#    simulator.run_simulatorMP(nameF, platform=platform, params=simParams)

#    sys.exit()
    
    # for name in nameF.iter_names():
    #     # unpack
    #     fastaFile = name.get_fastaFile()
    #     indexFile = name.get_indexFile()
    #     #readFile = mg.get_readFile()    

    #     # read simulation 
    #     outDir = os.path.abspath(os.path.curdir)
    #     (simReadsFile,fileType) = simulator.run_simulator(fastaFile, outDir=outDir,
    #                                                       platform=platform,
    #                                                       params={platform : simParams})
        
    #     # find out how many reads were generated
    #     ### Attention: Here we assume that all files contain the same number of read and are stored in fastq format
    #     num_reads = len( [ True for i in SeqIO.parse(simReadsFile, fileType) ] )

    #     ## convert readFile to fasta if needed
    #     if fileType.lower() == 'fastq':
    #         fastaFile = os.path.splitext(simReadsFile)[0] + '.fna'
    #         SeqIO.convert(simReadsFile, 'fastq', fastaFile, 'fasta')
    #         simReadsFile = fastaFile
        
    #     # saving reads file names in names class
    #     name.set_simReadsFile(simReadsFile)

        
    # pairwise mapping of the simulated reads from each ref to all references 
    pwm = PairwiseMapper(nameF, mapper)
    simSamFiles = pwm.pairwiseMap()

    
    # parse sam files to create a numpy array of reads mapped
    (I,J) = simSamFiles.shape
    mappedReads = np.zeros((I, J, num_reads))    
    for i in range(I):
        for j in range(J):
            #print simSamFiles[i,j]
            
            # count the reads in i mapping to subject j
            samfile = pysam.Samfile(simSamFiles[i,j], "r")
            mappedReads[i,j,:] = np.array( [int(not read.is_unmapped) for read in samfile] )
            samfile.close()
               
            
    # save the similarity matrix
    matrixOutFile = mgID + '_simMtx'
    np.save(matrixOutFile, mappedReads)
    matrixOutFile += '.npy'
    sys.stderr.write('Wrote similarity matrix: {}\n'.format(matrixOutFile))
            

    # similarity correction 
    ## input: matrix & original reads -> ref sam file
    refSamFiles = [name.get_refSamFile() for name in nameF.iter_names()]
    CorAbund = CorrectAbundances()

    nBootstrap = 3
    result = CorAbund.similarityCorrection(refSamFiles, matrixOutFile, nBootstrap)
    

    # writing output
    for i,name in enumerate(nameF.get_names()):
        total = result['total']
        outvals = dict(
            ref = name.get_fastaFile(),
            mapped = result['num_reads'][i],
            corr = result['corr'][i] * total,
            error = result['err'][i] * total,
            pval = result['p'][i],
            mgID = mgID  # metagenome containing the reads used
            )

        print '{mgID}\t{ref}\t{mapped}\t{corr}\t{error}\t{pval}'.format(**outvals)

        
    # moving back to original working directory
    if args['--debug'] == False:
        os.chdir(origWorkDir)
        
    # debug
    #sys.exit()
