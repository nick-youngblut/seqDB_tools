"""Read mapping classes. Using strategy pattern."""

import sys
import os
from distutils.spawn import find_executable
import multiprocessing as mp
from functools import partial
import parmap
import numpy as np
import uuid



def randomString(string_length=10):
    """Returns a random string"""
    rs = str(uuid.uuid4())
    rs = rs.replace('-','')
    return rs[0:string_length]

    

class ReadMapper(object):

    @staticmethod
    def getMapper(mapper=None):
        """factory designating subclass to use for mapping

        Args:
        mapper -- designates mapper subclass

        Mappers:
        bowtie2
        """
        
        mappers = dict(bowtie2=MapperBowtie2)

        if mapper in mappers:
            return mappers[mapper]()
        else:
            raise TypeError('mapper: "{0}" not implemented')


    def exeExists(self, exe):
        """checking to see if executable is in path"""
        if find_executable(exe):
            return 1
        else:
            raise IOError('"{0}" is not in your $PATH'.format(exe))

            

class MapperBowtie2(ReadMapper):

    def __init__(self, executable=['bowtie2-build', 'bowtie2']):
        # checking for bowtie2 in path
        [self.exeExists(x) for x in executable]
        # attr
        self.exe = executable
                

    def __call__(self, indexFile, readFile, outFile=None,
                   subproc=None, samFile=None, params={'-f': ''}):
        """Calling bowtie2 for mapping

        Args:
        indexFile -- bowtie2 index file
        readFile -- read file provided to bowtie2
        outFile -- sam output file. If None: using indexFile basename.
        subproc -- subprocess?
        samFile -- output SAM file name. Default to edited indexFile name.
        params -- bowtie2 parameters. Value = '' if boolean parameter
        """
        # outFile name
        if outFile is None:
            (basename, ext) = os.path.splitext(indexFile)
            outFile = basename + '.sam'
        if samFile is not None:
            outFile = samFile

        # setting params if any exist
        params = ' '.join( ['{0} {1}'.format(k,v) for k,v in params.items()] )
        
        # calling bowtie2
        cmd = 'bowtie2 -U {reads} -x {index} -S {samfile} --local {params}'
        cmd = cmd.format(reads=readFile, index=indexFile, samfile=outFile, params=params)
        sys.stderr.write( 'Executing: "{0}"\n'.format(cmd) )
        os.system(cmd)

        # return
        if subproc is not None:
            subproc.send(outFile)
            subproc.close()
        else:
            return outFile
        
            
    def parallel(self, names, mg, nprocs=1, **kwargs):
        """Calling mapper using multiple processors.

        Args:
        names -- NameFile instance with iter_names() method
        mg -- MetaFile instance with readFile attrib
        nprocs -- number of parallel calls
        kwargs -- passed to mapper call
        
        Return:
        refSamFile attrib set for each name in names 
        """

        # making list of tuples (indexFile, readFile)
        lt = [(name.get_indexFile(), mg.get_readFile()) for name in names.iter_names()]

        # altering function kwargs
        new_mapper = partial(self, **kwargs)

        # calling mapper
        samFiles = parmap.starmap(new_mapper, lt, processes=nprocs)

        # adding samFile attrib to name instances
        for i,name in enumerate(names.iter_names()):
            name.set_refSamFile(samFiles[i])
                                
            
    def make_index(self,subjectFile, outFile=None, **kwargs):
        """Making index file for subject fasta file

        Args:
        subjectFile -- subject sequence to call bowtie2 (sequence being mapped to)
        outFile -- output file name
        kwargs -- passed to bowtie2-build
        """
        # outFile
        if outFile is None:
            (basename, ext) = os.path.splitext(indexFile)
            outFile = basename + '.sam'            

        # setting params if any exist
        params = ' '.join( ['-{0} {1}'.format(k,v) for k,v in kwargs.items()] )

        # call command
        cmd = 'bowtie2-build {subject} {outFile}'
        cmd = cmd.format(subject=subjectFile, outFile=outFile)        
        sys.stderr.write( 'Executing: "{0}"\n'.format(cmd) )
        os.system(cmd)
        return outFile

    

        


class PairwiseMapper_OLD(object):
    """Class for pairwise read mapping"""

    def __init__(self, name, mapper):
        """Providing nameFile object and mapper to determine
        what to map and how to do it.
        
        Args:
        name -- nameFile object
        mapper -- ReadMapper object

        Return:
        2d numpy array of SAM files produced by mapping
        """
        self.name = name
        self.n_refs = name.len()
        self.mapper = mapper    # readmapper class
            
            
    def pairwiseMap(self):
        """ resulting SAM file names saves as numpy array"""

        # setting up multiprocessing
        parentConns = []
        procs = []
        
        for i in range(self.n_refs):
            # getting simulated reads from first query reference taxon
            nameQuery = self.name.get_name(i)
            simReadsFile = nameQuery.get_simReadsFile()
        
            for j in range(self.n_refs):
                # getting index file of subject for mapping to subject
                nameSubject = self.name.get_name(j)
                indexFile = nameSubject.get_indexFile()

                # unique names for sam file
                samFile = randomString() + '.sam'
                
                # setting proc (pipe)
                parentConn, childConn = mp.Pipe()
                parentConns.append([i,j,parentConn])
                p = mp.Process(target=self.mapper.run_mapper,
                               args=[indexFile, simReadsFile],
                               kwargs=dict(fileType='fasta',subproc=childConn,
                                       samFile = samFile))
                procs.append(p)
                p.start()
                
        # getting data from children and adding to numpy array
        simSamFiles = np.array([['' for i in range(self.n_refs)] for j in range(self.n_refs)], dtype=object)
        for row in parentConns:
            i = row[0]
            j = row[1]
            parentConn = row[2]
            simSamFiles[i,j] = parentConn.recv()

        # waiting for all processes to finish
        for p in procs:
            p.join()

        # return array
        return simSamFiles
                
        
        
            