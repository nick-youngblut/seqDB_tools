"""Read mapping classes. Using strategy pattern."""

import sys
import os
from distutils.spawn import find_executable

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
                

    def run_mapper(self, indexFile, readFile, outFile=None, fileType='fastq', **kwargs):
        """Calling bowtie2 for mapping

        Args:
        indexFile -- bowtie2 index file
        readFile -- read file provided to bowtie2
        outFile -- sam output file. If None: using indexFile basename.
        """
        # outFile
        if outFile is None:
            (basename, ext) = os.path.splitext(indexFile)
            outFile = basename + '.sam'            

        # setting params if any exist
        params = ' '.join( ['-{0} {1}'.format(k,v) for k,v in kwargs.items()] )
        ## fileType
        if fileType.lower() == 'fasta':
            params += ' -f'
        
        # calling bowtie2
        cmd = 'bowtie2 -U {reads} -x {index} -S {samfile} --local {params}'
        cmd = cmd.format(reads=readFile, index=indexFile, samfile=outFile, params=params)
        sys.stderr.write( 'Executing: "{0}"\n'.format(cmd) )
        os.system(cmd)
        return outFile


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


        
