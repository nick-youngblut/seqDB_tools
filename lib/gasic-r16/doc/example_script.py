"""
File: example_script.py

Description:
    This file contains commented analysis script demonstrating how GASiC
    can be used via the scripting interface. The contained code only has
    exemplary character and is not meant to be executed directly.

Author: Martin Lindner, LindnerM@rki.de
"""


"""
    Example 1:
    ----------

    Simple GASiC workflow. Reproduces the command line tool
    functionalities.

    Assumptions:
    - Data: FASTQ reads in './data/reads.fastq'
    - Reference genomes: FASTA files in './ref/'
    - Mapper index: Index files in './ref/index/'
    - bowtie mapper
    - Mason simulator
    - Names File in 'names.txt'
"""

# load module for system operations
import os


# load GASiC modules
from core import tools
import create_matrix
import correct_abundances


# set all paths and patterns
p_data  = './data/reads.fastq'
p_ref   = './ref/%s.fasta'
p_index = './ref/index/%s'
p_SAM   = './SAM/%s.sam'
p_temp  = './distance'


# create directories for SAM files and temporary files
os.makedirs(p_SAM)
os.makedirs(p_temp)


# read the Names File
names = tools.read_names('names.txt')


# Step 1: map the reads to every reference genome
for name in names:
    # fill the current name into the pattern of reference index and SAM path
    index_name = p_index%name
    SAM_name = p_SAM%name

    # run the read mapper
    tools.run_mapper['bowtie'](index_name, p_data, SAM_name)


# Step 2: calculate the similarity matrix
sim_matrix = create_matrix.similarity_matrix_raw(names, p_ref, p_index, p_temp, 'Mason', 'bowtie')


# Step 3: estimate the true abundances
total,mapped,est,err,p = correct_abundances.similarity_correction(names, sim_matrix, p_SAM, 100)

print "Finished."




"""
    Example 2:
    ----------

    Apply GASiC to 5 datasets, e.g. containing repetitions of the
    same experiment

    Assumptions:
    - Data: FASTQ reads in './data/' with names 'rep0.fastq' - 'rep5.fastq'
    - Reference genomes: FASTA files in './ref/'
    - Mapper index: Index files in './ref/index/'
    - a custom mapper
    - Mason simulator
    - Names File in 'names.txt'
"""

# load module for system operations
import os


# load GASiC modules
from core import tools
import create_matrix
import correct_abundances


# set all paths and patterns
p_data  = './data/rep%i.fastq'
p_ref   = './ref/%s.fasta'
p_index = './ref/index/%s'
p_temp  = './distance'


# create directories for SAM files and temporary files
os.makedirs('./SAM')
os.makedirs(p_temp)


# read the Names File
names = tools.read_names('names.txt')


# define a wrapper function for 'my_mapper'
def run_my_mapper(index, reads, out, param=""):
    command = "my_mapper -d some_default {param} {index} {reads} > {samfile}".format(index=index, reads=reads, samfile=out, param=param)
    print "Executing:",command
    os.system(command)
    return 1


# register wrapper function temporarily as 'my_mapper'
tools.run_mapper['my_mapper'] = run_my_mapper


# calculate the similarity matrix only once
sim_matrix = create_matrix.similarity_matrix_raw(names, p_ref, p_index, p_temp, 'Mason', 'my_mapper')


# set up numpy arrays to store all results
num_genomes = len(names)
import numpy as np
total = np.zeros((5,))
mapped= np.zeros((5,num_genomes))
est   = np.zeros((5,num_genomes))
err   = np.zeros((5,num_genomes))
p     = np.zeros((5,num_genomes))
unique= np.zeros((5,num_genomes))

# iterate over all datasets
for it in range(5):
    dataset = p_data%it

    # map the reads to every reference genome
    for name in names:
        # fill the current name into the pattern of reference index and SAM path
        index_name = p_index%name
        SAM_name = './SAM/rep%i_%s.sam'%(it,name)
        
        # run the read mapper
        tools.run_mapper['bowtie'](index_name, dataset, SAM_name)

        # estimate the true abundances and count unique reads
        SAM_pattern = './SAM/rep{num}_%s.sam'.format(num=it)
        total[it],mapped[it,:],est[it,:],err[it,:],p[it,:] = correct_abundances.similarity_correction(names, sim_matrix, SAM_pattern, 100)
        unique[it,:] = correct_abundances.unique(names, SAM_pattern)

# GASiC calculation finished. Now we can post-process the results
# for example, calculate mean and variance for each genome over all datasets
mean_est = np.mean(est,axis=0)
var_est = np.var(est,axis=0)

print mean_est
print var_est

print "Finished."
