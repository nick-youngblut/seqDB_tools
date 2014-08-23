"""Classes for read simulation"""

import sys
import os

from distutils.spawn import find_executable


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


    def get_paramsByReadStats(self, mg):
        """Getting simulator params based on read stats (e.g., read lengths &
        sequencing platform

        Args:
        mg -- MetaFile row class
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
        params = dict()
        if hasattr(mg, 'readStats'):
            if mg.platform == 'illumina':
                params['--read-length'] =  mg.readStats['median']
            elif mg.platform == '454' or mg.platform == 'sanger':
                params['--read-length-mean'] = mg.readStats['mean']
                params['--read-length-error'] = mg.readStats['stdev']
        return platform, params
            

    def run_simulator(self, refFile, outFile=None, outDir=None,
                      platform='illumina', params=None):
        """Calling mason; default for creating illumina reads

        Default keyword params for mason (override any of them with params):

        ADD HERE
        
        Args:
        refFile -- fasta file using for generating reads
        outFile -- output
        outDir -- directory to write output. Default: same as refFile
        platform -- sequencing platform
        params -- parameters passed to mason. dict of dict: {platform : {param : value}}.
                  value = '' if param is boolean
        
        Return:
        string with file name written by mason
        """

        # output file: if None: creating from refFile
        if outFile is None:
            (basename, ext) = os.path.splitext(refFile)
            outFile = basename + '_simReads.fq'        
        if outDir is not None:
            outFile = os.path.join(outDir, os.path.basename(outFile))
        logFile = outFile + '.log'

        # setting params
        ## defaults
        defaultParams = {'illumina' : {
            '-N' : 10000,
            '-hi' : 0,
            '-hs' : 0,
            '-n' : 100,
            '-sq' : ''
        },'454' : {
            '-N' : 10000,
            '-hi' : 0,
            '-hs' : 0,
            '-sq' : ''
        },'sanger' : {
            '-N' : 10000,
            '-hi' : 0,
            '-hs' : 0,
            '-o' : 'fastq'
        }}
        ## changing defaults
        for k,v in defaultParams.items():
            if params.has_key(k):
                defaultParams[k].update(params[k])
        
        ## params as string
        cmdParams = ' '.join(['{} {}'.format(k,v) for k,v in defaultParams[platform].items()])
                    
        # calling
        cmd = 'mason {platform} {params} -o {outFile} {refFile} > {logFile}'
        cmd = cmd.format(platform=platform, params=cmdParams, refFile=refFile, outFile=outFile, logFile=logFile)
        
        sys.stderr.write("Executing: {0}\n".format(cmd))
        os.system(cmd)
        return outFile, 'fastq'


