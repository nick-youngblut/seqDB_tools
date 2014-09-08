"""Classes for read simulation"""

import sys
import os

from distutils.spawn import find_executable
import multiprocessing as mp
import parmap
from functools import partial
from Bio import SeqIO


class ReadSimulator(object):
    """Factory class for read simulator"""
    
    @staticmethod
    def getSimulator(simulator):
        """Factor method for selecting simulator class"""
        
        simulators = dict(mason=mason, griner=grinder)

        simulator = simulator.lower()
        if simulator in simulators:
            return simulators[simulator]()
        else:
            raise TypeError('Simulator: "{0}" not supported\n'.format(simulator))

    def exeExists(self, exe):
        """checking to see if executable is in path"""
        if find_executable(exe):
            return 1
        else:
            raise IOError('"{0}" is not in your $PATH'.format(exe))


class grinder(ReadSimulator):
    """Class for calling Grinder simulator"""

    def run_simulator(refFile, outFile=None, params=''):
        """Calling grinder simulator"""
        print 'TODO'
        pass


class mason(ReadSimulator):
    """Class for calling mason simulator"""

    def __init__(self, executable='mason'):
        # in PATH?
        self.exeExists(executable)
        # attr
        self.exe = executable


    def get_paramsByReadStats(self, mg, params=dict()):
        """Getting simulator params based on read stats (e.g., read lengths &
        sequencing platform

        Args:
        mg -- MetaFile row class
        params -- dict that sets initial params
        """
        # params for platform
        platform = None
        if mg.platform == 'illumina':
            platform = 'illumina'
        elif mg.platform == '454' or mg.platform.startswith('pyro'):
            platform = '454'
        elif mg.platform == 'sanger':
            platform = 'sanger'
        else:
            raise KeyError('"{}" platform is not supported!\n'.format(mg.platform))
            
        # stats
        if hasattr(mg, 'readStats'):
            if mg.platform == 'illumina':
                params['--read-length'] =  int(mg.readStats['median'])
            elif mg.platform == '454' or mg.platform == 'sanger':
                params['--read-length-mean'] = mg.readStats['mean']
                params['--read-length-error'] = mg.readStats['stdev']
        return platform, params
            

    def __call__(self, refFile, outFile=None, outDir=None,
                      platform='illumina', params=None):
        """Calling mason; default for creating illumina reads

        Default keyword params for mason (override any of them with params):

        ADD PARAMS HERE
        
        Args:
        refFile -- fasta file using for generating reads
        outFile -- output
        outDir -- directory to write output. Default: same as refFile
        platform -- sequencing platform
        params -- parameters passed to mason. {param : value}
                  value = '' if param is boolean
        
        Return:
        string with file name written by mason

        """
        # setting params
        ## defaults
        defaultParams = {'illumina' : {
            '--num-reads' : 10000,
            '--haplotype-indel-rate' : 0,
            '--haplotype-snp-rate' : 0,
            '--read-length' : 100 },
            '454' : {
            '--num-reads' : 10000,
            '--haplotype-indel-rate' : 0,
            '--haplotype-snp-rate' : 0 },
            'sanger' : {
            '--num-reads' : 10000,
            '--haplotype-indel-rate' : 0,
            '--haplotype-snp-rate' : 0 }}
        ## changing defaults
        #for k,v in defaultParams.items():
        #    if params.has_key(k):
        defaultParams[platform].update(params)

        
        # output filetype
        fileType = 'fasta'
        for k,v in defaultParams.items():
            if '-sq' in v.keys() or '--simulate-qualities' in v.keys():
                fileType = 'fastq'        

                
        # output file: if None: creating from refFile
        if outFile is None:
            (basename, ext) = os.path.splitext(refFile)
            if fileType == 'fasta':
                ext = '.fa'
            else:
                ext = '.fq'
            outFile = basename + '_simReads' + ext
        if outDir is not None:
            outFile = os.path.join(outDir, os.path.basename(outFile))
        logFile = outFile + '.log'

        
        ## params as string
        cmdParams = ' '.join(['{} {}'.format(k,v) for k,v in defaultParams[platform].items()])
                    
        # calling
        cmd = 'mason {platform} {params} -o {outFile} {refFile} > {logFile}'
        cmd = cmd.format(platform=platform, params=cmdParams, refFile=refFile, outFile=outFile, logFile=logFile)
        
        sys.stderr.write("Executing: {0}\n".format(cmd))
        os.system(cmd)


        # return 
        return dict(simReadsFile=outFile, simReadsFileType=fileType)

        
    def parallel(self, names, fileType='fasta', nprocs=1, **kwargs):
        """Running simulator using apply_async

        Args:
        names -- NameFile class with iter_names() method
        nprocs -- max number of parallel simulation calls
        kwargs -- passed to simulator

        Attribs added to each name instance in names:
        simReadsFile -- file name of simulated reads
        simReadsFileType -- file type (eg., 'fasta' or 'fastq')
        simReadsFileCount -- number of simulated reads

        Return:
        boolean of run success/fail
        """
        # making list of fasta file to provide simulator call
        fastaFiles = [name.get_fastaFile() for name in names.iter_names()]

        # settig kwargs
        new_simulator = partial(self, **kwargs)

        # calling simulator
        res = parmap.map(new_simulator, fastaFiles, processes=nprocs)

        # checking that simulated reads were created for all references
        for row in res:
            if not os.path.isfile(row['simReadsFile']):
                return False
            elif os.stat(row['simReadsFile'])[0] == 0:
                return False
        
        # converting reads to fasta if needed
        if fileType.lower() == 'fasta':
            for result in res:
                simFile = result['simReadsFile']
                fileType = result['simReadsFileType'].lower()
                if fileType != 'fasta':
                    fastaFile = os.path.splitext(simFile)[0] + '.fna'
                    SeqIO.convert(simFile, fileType, fastaFile, 'fasta')
                    result['simReadsFile'] = fastaFile
                    result['simReadsFileType'] = 'fasta'
                    
        # setting attribs in name instances                    
        for i,name in enumerate(names.iter_names()):
            # read file
            simReadsFile = res[i]['simReadsFile']
            name.set_simReadsFile(simReadsFile)
            # file type
            fileType = res[i]['simReadsFileType'].lower()
            name.set_simReadsFileType(fileType)
            # number of simulated reads            
            num_reads = len([True for i in SeqIO.parse(simReadsFile, fileType)])
            name.set_simReadsCount(num_reads)
            
        return True