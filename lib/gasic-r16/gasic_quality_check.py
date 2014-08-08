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

import optparse
import glob
import os
import sys
try:
    import matplotlib.pyplot as plt
except ImportError:
    raise Exception('Import Error: Failed to load matplotlib module. Please install matplotlib properly (http://matplotlib.sourceforge.net/).')

try:
    import pysam
except ImportError:
    raise Exception('Import Error: Failed to load pysam module. Please install pysam properly (http://code.google.com/p/pysam/).')

import numpy as np
import scipy.stats as st
from core import tools

if __name__=="__main__":
    usage = """%prog [Options] SAMFILE

Perform a sanity check on the mapping results (SAM file).

This tool analyzes the output of the read mapper and provides useful
information to the user to decide whether the set of reference sequences
and the mapper settings are feasible for the given dataset.

Input:
SAMFILE:  Name of the SAM file to analyze
        
"""

    # configure the parser
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-o', '--outpath', type='string', dest='out', default='./', help='Path to the directory for the analysis output. [default: %default]')
    # parse arguments
    options, args = parser.parse_args()

    numArgs = len(args)
    if numArgs >= 1:
        sam_names = args

        # create output directory if necessary
        if not os.path.exists(os.path.dirname(options.out)):
            os.makedirs(os.path.dirname(options.out))
    
        # analyze all files and gather statistics
        for i,nm in enumerate(sam_names):
            dataset_name = os.path.splitext(os.path.split(nm)[1])[0]
            sf = pysam.Samfile(nm,'r')
            cov = np.zeros( (sum(sf.lengths),) )
            start_pos = np.cumsum(sf.lengths)-sf.lengths[0]
            total_reads = 0
            mapped_reads = 0
            for read in sf:
                total_reads += 1
                if not read.is_unmapped:
                    r_start = start_pos[read.tid] + read.pos
                    r_end = start_pos[read.tid] + read.pos + read.qlen
                    cov[r_start:r_end] += 1
                    mapped_reads += 1

            # calculate coverage
            max_cov = np.max(cov)
            mean_cov = np.mean(cov)

            # plot the coverage
            plt.figure()
            h = plt.hist(cov, bins=min(max_cov,100), range=(0,max_cov))
            plt.xlabel('coverage')
            plt.ylabel('# occurrences')
            plt.title('Coverage histogram for '+dataset_name)
            plt.savefig(options.out+'/'+dataset_name+'_coverage.png', dpi=300)

            # write out user information
            f = open(options.out+'/'+dataset_name+'_info.txt','w')
            f.write(' Name:\t'+nm+'\n')
            f.write('Total Genome Length:\t%i\n'%(sum(sf.lengths)))
            f.write('Num. Contigs:\t%i\n\n'%(len(sf.lengths)))
            f.write('Mapping Results:\n')
            f.write('Total Reads:\t%i\n'%(total_reads))
            f.write('Mapped Reads:\t%i\n'%(mapped_reads))
            f.write('Average Coverage:\t%f\n'%(mean_cov))
            zero_cov = h[0][0]/float(sum(h[0]))
            f.write('Fraction of bases with 0 Coverage:\t%f\n\n'%(zero_cov) )

            f.write('Comments:\n')
            if mean_cov < 1.:
                f.write('Low overall coverage. ')
                if mapped_reads/float(total_reads)>0.01:
                    f.write('Consider using a larger dataset.\n')
                else:
                    f.write('Abundance of this or related sequences is below 0.01.\n')
            
            p = st.poisson(mean_cov)
            if zero_cov > 2*p.pmf(0):
                f.write('Unnaturally many bases with zero coverage. Consider that this species is not present in the dataset.\n')
            
            f.close()
    else:
        parser.print_help()
        sys.exit(1)
    
