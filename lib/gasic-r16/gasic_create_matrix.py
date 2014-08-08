#!/usr/bin/python
"""
Copyright (c) 2012, Martin S. Lindner, LindnerM@rki.de, 
Robert Koch-Institut, Germany,
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * The name of the author may not be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL MARTIN S. LINDNER BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import sys
import os
import glob
import optparse

import numpy as np
from Bio import SeqIO
import pysam

from core import tools


def similarity_matrix_raw(names, ref_pattern, index_pattern, temp_dir, simulator, mapper):
    """
    Perform the read generation and mapping step for the similarity matrix
    calculation. The raw output is used to bootstrap a similarity matrix.
    
    INPUT:
    names:             array of genome names
    ref_pattern:       pattern pointing to the reference genome files used by the read simulator
    index_pattern:     pattern pointing to the mapping index files used by the read mapper
    temp_dir:          directory where temporary files are stored (simulated reads, SAM files).
                       Attention: make sure that there is enough space available!

    OUTPUT:
    mapped_reads:      Mapping information about mapped reads readily usable by
                       bootstrap_similarity_matrix().
    """

    if not mapper in tools.run_mapper:
        raise Exception('Aborting. Mapper "%s" not found in tools.py'%mapper)
    if not simulator in tools.run_simulator:
        raise Exception('Aborting. Simulator "%s" not found in tools.py'%simulator)
    if not os.path.isdir(temp_dir):
        os.makedirs(temp_dir)
    
    n_seq = len(names)
    rng = range(n_seq)

    # construct arrays with real file names from pattern
    ref_files = [ ref_pattern%nm for nm in names ] # filenames of reference sequences
    index_files = [ index_pattern%nm for nm in names ]   # filenames of mapper index files
    sim_files = [ temp_dir+'/'+nm+'.fastq' for nm in names ] # filenames of simulated read files
    
    # generate reads for every reference genome
    for i in rng:
        tools.run_simulator[simulator](ref_files[i], sim_files[i])

    # find out how many reads were generated
    # Attention: Here we assume that all files contain the same number of reads
    # and are stored in fastq format
    num_reads = len( [ True for i in SeqIO.parse(sim_files[0],'fastq') ] )
    
    # map the reads of every reference to all references
    for i in rng:
        for j in rng:
            samfile = temp_dir+'/'+names[i]+'-'+names[j]+'.sam'
            tools.run_mapper[mapper](index_files[j], sim_files[i], samfile)

    # parse SAM files
    mapped_reads = np.zeros((n_seq,n_seq, num_reads))
    for i in rng:
        for j in rng:
            # count the reads in i mapping to refernce j
            samfile = temp_dir+'/'+names[i]+'-'+names[j]+'.sam'
            samhandle = pysam.Samfile(samfile, "r")
            mapped_reads[i,j,:] = np.array( [int(not read.is_unmapped) for read in samhandle] )

    return mapped_reads





if __name__=="__main__":
    usage = """%prog [options] NAMES

Calculate the similarity matrix.

First, a set of reads is simulated for every reference genome using a read
simulator from core/tools.py specified via -s.
Second, the simulated reads of each species are mapped against all reference
genomes using the mapper specified with -m.
Third, the resulting SAM-files are analyzed to calculate the similarity
matrix. The similarity matrix is stored as a numpy file (-o).

Input:
NAMES:  Filename of the names file; the plain text names file should
        contain one name per line. The name is used as identifier in
        the whole algorithm.

See the provided LICENSE file or source code for license information.
"""

    # configure the parser
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-s', '--simulator', type='string', dest='simulator', default=None, help='Identifier of read simulator defined in core/tools.py [default: %default]')
    parser.add_option('-r', '--reference', type='string', dest='ref', default='./ref/%s.fasta', help='Reference sequence file pattern for the read simulator. Placeholder for the name is "%s". [default: %default]')    
    parser.add_option('-m', '--mapper', type='string', dest='mapper', default=None, help='Identifier of mapper defined in core/tools.py [default: %default]')
    parser.add_option('-i', '--index', type='string', dest='index', default='./ref/%s.fasta', help='Reference index files for the read mapper. Placeholder for the name is "%s". [default: %default]')    
    parser.add_option('-t', '--temp', type='string', dest='temp', default='./temp', help='Directory to store temporary simulated datasets and SAM files. [default: %default]')
    parser.add_option('-o', '--output', type='string', dest='out', default='./similarity_matrix.npy', help='Output similarity matrix file. [default: %default]')
    # parse arguments
    opt, args = parser.parse_args()
    
    numArgs = len(args)
    if numArgs == 1:
        # read the Names file
        names_file = args[0]
        names = tools.read_names(names_file)

        # construct the similarity matrix
        smat = similarity_matrix_raw(names, opt.ref, opt.index, opt.temp, opt.simulator, opt.mapper)

        # save the similarity matrix
        np.save(opt.out, smat)
        print "Wrote similarity matrix to",opt.out
    else:
        parser.print_help()
        sys.exit(1)
