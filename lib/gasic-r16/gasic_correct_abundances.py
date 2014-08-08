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

import numpy as np
import pysam
import sys
import glob
import optparse

from core import gasic
from core import tools


def similarity_correction(names, smat_raw, sam_pattern, bootstrap_samples):
    """
    Perform similarity correction step. The similarity matrix and mapping
    results must be available.
    
    INPUT:
    names:             array of genome names
    smat_raw:          mapping information for similarity matrix with same ordering of genomes as in 'names'
    sam_pattern:       pattern pointing to the filenames of the SAM-files to analyze
    bootstrap_samples: number of bootstrap samples, use 1 to disable bootstrapping

    OUTPUT:
    total:             total number of reads in the dataset
    num_reads:         number of reads mapped to each genome (array)
    corr:              abundance of each genome after similarity correction
    err:               estimated standard error
    p:                 p-value for the confidence, that the true abundance is above some threshold
    """

    # find out the number of reads
    total = len( [1 for read in pysam.Samfile(sam_pattern%(names[0]), "r")] )
    print ". found %i reads"%total

    # initialize some arrays
    #   mapping information; mapped[i,j]=1 if read j was successfully mapped to i.
    mapped = np.zeros( (len(names), total) )

    #   total number of successfully mapped reads per reference
    num_reads = np.zeros( (len(names),) )

    # analyze the SAM files
    for n_ind,nm in enumerate(names):
        print ".. analyzing SAM-File %i of %i"%(n_ind+1,len(names))
        # samfile filename
        sf = pysam.Samfile(sam_pattern%nm, "r")

        # go through reads in samfile and check if it was successfully mapped
        mapped[n_ind,:] = np.array([int(not rd.is_unmapped) for rd in sf])
        num_reads[n_ind] = sum(mapped[n_ind,:])

    # run similarity correction step
    p,corr,var = gasic.bootstrap(mapped, smat_raw, bootstrap_samples)
    err = np.sqrt(var)
    return total,num_reads,corr,err,p
    

def unique(names, sam_pattern):
    """ Determine the number of unique reads for every species based on the read names.
    INPUT:
    names:             array of genome names
    sam_pattern:       pattern pointing to the filenames of the SAM-files to analyze

    OUTPUT:
    unique:            number of unique reads per species.
    """
    # one set for the names of mapped reads for each species
    mapped_read_names = [set() for nm in names]
    
    for n,nm in enumerate(names):
        # parse the samfile
        sf = pysam.Samfile(sam_pattern%nm, "r")
        for read in sf:
            # add the hash of read name to the set, if read was mapped
            if not read.is_unmapped:
                mapped_read_names[n].add(hash(read.qname))

    unique_read_names = [set() for nm in names]
    for n in range(len(names)):
        others = set()
        for m in range(len(names)):
            if n!=m:
                others |= mapped_read_names[m]
        unique_read_names[n] = mapped_read_names[n] - others

    return np.array([len(unq) for unq in unique_read_names])

        
if __name__=="__main__":
    usage = """%prog NAMES

Run the similarity correction step.


Note: Although it is possible to run the read mappers by hand or to create the
similarity matrix manually, we strongly recommend to use the provided Python
scripts 'run_mappers.py' and 'create_similarity_matrix.py'.

Input:
NAMES:  Filename of the names file; the plain text names file should
        contain one name per line. The name is used as identifier in
        the whole algorithm.

See the provided LICENSE file or source code for license information.
"""

    # configure the parser
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-m', '--similarity-matrix', type='string', dest='smat', default='./similarity_matrix.npy', help='Path to similarity matrix file. The similarity matrix must be created with the same NAMES file. [default: %default]')
    parser.add_option('-s', '--samfiles', type='string', dest='sam', default='./SAM/%s.sam', help='Pattern pointing to the SAM files created by the mapper. Placeholder for the name is "%s". [default: %default]')

    parser.add_option('-b', '--bootstrap-samples', type='int', dest='boot', default=100, help='Set the number of bootstrap samples. Use 1 to disable bootstrapping [default: %default]')
    parser.add_option('-o', '--output', type='string', dest='out', default='./results.txt', help='Plain text output file containing the results. [default: %default]')
    # parse arguments
    opt, args = parser.parse_args()
    
    numArgs = len(args)
    if numArgs == 1:
        # read the Names file
        names_file = args[0]
        names = tools.read_names(names_file)

        # load the similarity matrix
        smat = np.load(opt.smat)

        # start similarity correction
        total,num_reads,corr,err,p = similarity_correction(names, smat, opt.sam, opt.boot)

        # write results into tab separated file.
        ofile = open(opt.out,'w')
        ofile.write("#genome name\tmapped reads\testimated reads\testimated error\tp-value\n")
        for n_ind,nm in enumerate(names):
            # Name, mapped reads, corrected reads, estimated error, p-value
            out = "{name}\t{mapped}\t{corr}\t{error}\t{pval}\n"
            ofile.write(out.format(name=nm,mapped=num_reads[n_ind],corr=corr[n_ind]*total,error=err[n_ind]*total,pval=p[n_ind]))
        ofile.close()
        print ". wrote results to", opt.out

    else:
        parser.print_help()
        sys.exit(1)
