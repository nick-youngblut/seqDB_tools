"""Classes for read simulation"""

import sys
import os


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
        

    def run_simulator(self, refFile, outFile=None, outDir=None,
                      params='illumina -N 10000 -hi 0 -hs 0 -n 72 -sq'):
        """Calling mason; default for creating illumina reads

        Args:
        refFile: fasta file using for generating reads
        outFile: output

        Return:
        outFile
        """

        # output file: if None: creating from refFile
        if outFile is None:
            (basename, ext) = os.path.splitext(refFile)
            outFile = basename + '_simReads.fq'        
        if outDir is not None:
            outFile = os.path.join(outDir, os.path.basename(outFile))
        logFile = outFile + '.log'

        # calling
        cmd = 'mason {params} -o {outFile} {refFile} > {logFile}'
        cmd = cmd.format(params=params, refFile=refFile, outFile=outFile, logFile=logFile)
        sys.stderr.write("Executing: {0}\n".format(cmd))
        os.system(cmd)
        return outFile
