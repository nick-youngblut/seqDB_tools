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

from core import tools

if __name__=="__main__":
    usage = """%prog NAMES READS -i REF -o OUT -m MAPPER

Run a read mapper to map reads to reference genomes.

The names in the NAMES file will be inserted in the provided string
patterns. Each pattern must contain exactly one "%s" placeholder
(python string formatting). 

Input:
NAMES:  Filename of the names file; the plain text names file should
        contain one name per line. The name is used as identifier in
        the whole algorithm.
        
READS:  File containing the reads to be mapped.

"""

    # configure the parser
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-m', '--mapper', type='string', dest='mapper', default=None, help='Identifier of mapper defined in core/tools.py [default: %default]')
    parser.add_option('-i', '--index', type='string', dest='ref', default='./ref/%s.fasta', help='Pattern, that points to the reference sequences/indices when used with a name. Placeholder for the name is "%s". [default: %default]')
    parser.add_option('-o', '--output', type='string', dest='out', default='./SAM/%s.sam', help='Pattern, that points to the output SAM file, when used with a name. Placeholder for the name is "%s". [default: %default]')
    # parse arguments
    options, args = parser.parse_args()

    numArgs = len(args)
    if numArgs == 2:
        reads = args[1]
        names = tools.read_names(args[0])
        out = [options.out%nm for nm in names]
        ref = [options.ref%nm for nm in names]

        if not options.mapper:
            parser.print_help()
            raise Exception('Aborting. No mapper specified.')
        
        # create output directory if necessary
        if not os.path.exists(os.path.dirname(options.out)):
            os.makedirs(os.path.dirname(options.out))
    
        for i in range(len(names)):
            print ". mapping reads with %s to %s"%(options.mapper, ref[i]) 
            tools.run_mapper[options.mapper](ref[i], reads, out[i])
    else:
        parser.print_help()
        sys.exit(1)
    
