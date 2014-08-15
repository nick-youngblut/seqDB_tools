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

import os

# Names File Reader
#
# Reads a names file (plain text, one name per line) and returns an array of strings
def read_names(filename):
    f = open(filename,'r')
    names = [nm.rstrip('\n\r') for nm in f]
    f.close()
    return names


# Read Mappers
#
# Define caller functions for the read mappers here. These fuctions call the
# read mappers on the command line with the correspondig parameters. Output
# is a SAM file at the specified location. Each function takes 3 mandatory
# string arguments:
#   index: name of the reference sequence or reference index
#   reads: name of the file containing the reads
#   out:   name of the output SAM file
# Additional mapper parameters may be speciefied:
#   param: string with parameters for the mapper [default: ""]
#
# Note: the mapper should be configured in that way, that it only reports ONE
# match for each read (e.g. the best match).

def run_bowtie2(index, reads, out, param=""):
    command = "bowtie2 -U {reads} -x {index} -S {samfile} {param} --local -M 0 ".format(reads=reads, index=index, samfile=out, param=param)
    print "Executing:",command
    os.system(command)
    return 1

def run_bowtie(index, reads, out, param=""):
    command = "bowtie -S -p 2 -q -3 30 -v 2 {param} {index} {reads} > {samfile}".format(index=index, reads=reads, samfile=out, param=param)
    print "Executing:",command
    os.system(command)
    return 1

def run_bwa(index, reads, out, param=""):
    command = "bwa aln {param} {index} {reads} > /tmp/res.sai".format(index=index, reads=reads, param=param)
    print "Executing:",command
    os.system(command)
    # convert the bwa output to SAM and remove temporary file
    command_sam = "bwa samse %s /tmp/res.sai %s > %s && rm /tmp/res.sai"%(index,reads,out)
    os.system(command_sam)
    return 1

def run_bwasw(index, reads, out, param=""):
    command = "bwa bwasw {param} {index} {reads} > {samfile}".format(index=index, reads=reads, samfile=out, param=param)
    print "Executing:",command
    os.system(command)
    return 1


run_mapper = dict( bowtie=run_bowtie,
                   bowtie2=run_bowtie2,
                   bwa = run_bwa,
                   bwasw = run_bwasw,)

"""
How to add your custom mapper

1. Create a caller function
   - Create a copy of one of the existing functions, e.g. run_bowtie
   - Rename it, customize it, but DO NOT TOUCH the interface!

2. Add the caller function to the run_mapper dict
   - The dict entry should have the format: [name] = [caller function]

3. Now you can use your mapper:
 >>> import tools_lib
 >>> tools_lib.run_mapper["name"]
"""

                    

# Read Simulators
#
# Define caller functions for the read simulators here. These fuctions call
# the read simulators on the command line with the correspondig parameters.
# Output is a SAM file at the specified location. Each function takes 3
# mandatory string arguments:
#   ref: name of the reference sequence file
#   out: name of the output reads file
#
# Note: since simulators tend to have many tuning parameters, we encourage
# to create a seperate caller function for every scenario.

def run_mason_illumina(ref, out):
    command = "mason illumina -N 10000 -hi 0 -hs 0 -n 72  -sq -o %s %s"%(out,ref)
    print "Executing:",command
    os.system(command)
    # remove the needless SAM file
    command = "rm %s.sam"%out
    os.system(command)
    return 1

def run_dwgsim(ref, out):
    command = "dwgsim -c 2 -1 80 -2 0 -r 0 -y 0 -e 0.002 -N 100000 -f TACG %s %s"%(ref, out)
    print "Executing:",command
    os.system(command)
    # remove all additional files and rename reads file
    command = "mv {out}.bfast.fastq {out} && rm {out}.bwa.read1.fastq && rm {out}.bwa.read2.fastq && rm {out}.mutations.txt".format(out=out)
    os.system(command)
    return 1


run_simulator = dict( mason_illumina=run_mason_illumina,
                      dwgsim=run_dwgsim,)
                    
"""
How to add your custom read simulator

1. Create a caller function
   - Create a copy of one of the existing functions, e.g. run_dwgsim
   - Rename it, customize it, but DO NOT TOUCH the interface!

2. Add the caller function to the run_simulator dict
   - The dict entry should have the format: [name] = [caller function]

3. Now you can use your simulator:
 >>> import tools_lib
 >>> tools_lib.run_simulator["name"]
"""
