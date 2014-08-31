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
  --ncores-map=<nm>  Number of parallel read mapping calls. [default: 1]  
  --ncores-sim=<ns>  Number of parallel read simulations. [default: 1]
  --ncores-3rd=<nt   Number of cores used by 3rd party software. [default: 1] 
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

import pysam
import numpy as np

scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)
libDir = os.path.join(scriptDir, '../lib/gasic-r16/')
sys.path.append(libDir)

import gasicBatch.MetaFile as MetaFile
import gasicBatch.NameFile as NameFile
from gasicBatch.ReadMapper import ReadMapper #, PairwiseMapper
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
# unpack args
ncores_map = int(args['--ncores-map'])
ncores_sim = int(args['--ncores-sim'])
ncores_3rd = int(args['--ncores-3rd'])

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

    # unpack
    mgID = mg.get_ID()
    
    # making temp directory and chdir
    if args['--debug'] == False:
        tmpdir = tempfile.mkdtemp()
        os.chdir(tmpdir)
    else:
        tmpdir = os.curdir

    # downloading metagenome reads
    ret = mg.download( )
    if ret is None: 
        sys.stderr.write('No read file downloaded for Metagenome {0}. Moving to next metagenome'.format(mgID))
        continue
            
    #-- determine read stats --#
    ## skipping if no downloaded file found
    mg.get_ReadStats(fileFormat='fasta')


    #-- read mapping --#
    ## creating object for specific mapper
    mapper = ReadMapper.getMapper('bowtie2')    # factory class
   # mapper.set_paramsByReadStats(mg)
    
    ## calling mapper for each index file
    mapper.parallel(nameF, mg, nprocs=ncores_map, params={'-f':'', '-p':ncores_3rd})


    #-- similarity estimation by pairwise mapping simulated reads --#
    ## select simulator
    simulator = ReadSimulator.getSimulator('mason')
    ## setting params based on metagenome read stats & platform
    platform, simParams = simulator.get_paramsByReadStats(mg)
    ## calling simulator using process pool
    simulator.parallel(nameF, nprocs=ncores_sim, outDir=tmpdir, platform=platform, params=simParams)
    # finding out how many reads were generated
    num_reads = [name.get_simReadsCount() for name in nameF.iter_names()]
    if len(set(num_reads)) > 1:
        sys.stderr.write('Warning: differing numbers of reads generated by simulator')
    num_reads = num_reads[0]

        
    #-- pairwise mapping of the simulated reads from each ref to all references --#
    # making list of tuples for all pairwise comparisons
    pairwiseComps = []
    n_refs = nameF.len()
    for i in range(n_refs):
        # getting simulated reads from first query reference taxon
        nameQuery = nameF.get_name(i)
        simReadsFile = nameQuery.get_simReadsFile()
        
        for j in range(n_refs):
        # getting index file of subject for mapping to subject
            nameSubject = nameF.get_name(j)
            indexFile = nameSubject.get_indexFile()
            # append values
            pairwiseComps.append((i,j,indexFile,simReadsFile,))

    # pairwise mapping
    pairwiseComps = mapper.pairwise(pairwiseComps, nprocs=ncores_map,
                                    tmpFile=True, params={'-f':'', '-p':ncores_3rd})

    # convert to numpy array
    simSamFiles = np.array([['' for i in range(n_refs)] for j in range(n_refs)], dtype=object)
    for row in pairwiseComps:
        i = row[0]
        j = row[1]
        simSamFile = row[4]
        simSamFiles[i,j] = simSamFile
    
    
    # parse SAM files to create numpy array of number reads mapped    
    mappedReads = np.zeros((n_refs, n_refs, num_reads))
    for i in range(n_refs):
         for j in range(n_refs):
             # count the reads in i mapping to subject j
             samfh = pysam.Samfile(simSamFiles[i,j], "r")
             mappedReads[i,j,:] = np.array( [int(not read.is_unmapped) for read in samfh] )
             samfh.close()

    print mappedReads
    sys.exit()
             
    
    ## resulting SAM file names saves as numpy array
    # n_refs = nameF.len()
    # simSamFiles = np.array([['' for i in range(n_refs)] for j in range(n_refs)], dtype=object)
    # for i in range(n_refs):
    #     # getting simulated reads from first query reference taxon
    #     nameQuery = nameF.get_name(i)
    #     simReadsFile = nameQuery.get_simReadsFile()
        
    #     for j in range(n_refs):
    #     # getting index file of subject for mapping to subject
    #         nameSubject = nameF.get_name(j)
    #         indexFile = nameSubject.get_indexFile()
            
    #         # mapping
    #         simSamFile = mapper.run_mapper(indexFile, simReadsFile, fileType='fasta')   # f = bowtie2 flag for file type
            
    #         # adding samSimFile to names file
    #         #nameQuery.add_simSamFile(i, j, simSamFile)
    #         simSamFiles[i,j] = simSamFile
            
    # # parse SAM files to create numpy array of number reads mapped
    # mappedReads = np.zeros((n_refs, n_refs, num_reads))
    # for i in range(n_refs):
    #     for j in range(n_refs):
    #         # count the reads in i mapping to subject j
    #         samfh = pysam.Samfile(simSamFiles[i,j], "r")
    #         mappedReads[i,j,:] = np.array( [int(not read.is_unmapped) for read in samfh] )

            
    # # parse sam files to create a numpy array of reads mapped
    # (I,J) = simSamFiles.shape
    # mappedReads = np.zeros((I, J, num_reads))    
    # for i in range(I):
    #     for j in range(J):
    #         #print simSamFiles[i,j]
            
    #         # count the reads in i mapping to subject j
    #         samfile = pysam.Samfile(simSamFiles[i,j], "r")
    #         mappedReads[i,j,:] = np.array( [int(not read.is_unmapped) for read in samfile] )
    #         samfile.close()
               



            
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
