import sys
import os
import re
import gzip
import zlib
from StringIO import StringIO

from Bio import SeqIO
import scipy
import requests
import pandas as pd



class MetaFile(object):
    """metadata table file object class
    """
    
    def __init__(self, fileName=None, fileObject=None, 
                 stages=None, **pandas_kwargs):
        """Loading file as pandas dataFrame. Checking for necessary columns
        
        Kwargs:
        seqDB -- datafile type: 'MGRAST' or 'SRA'
        stages -- mg-rast processing stage for downloading files. If multiple stages provided,
          will try each in succession until an entry is returned.
        """

        # attributes
        self.fileName = fileName
        self.fileObject = fileObject
        try:
            self.stages = [str(x) for x in stages]
        except TypeError:
            self.stages = None
        self.pd_kwargs = pandas_kwargs

        # loading MetaFile file
        self.loadFile()        


    def loadFile(self):        
        """Loading file as pandas DataFrame """
        if self.fileName is not None:
            self._tbl = pd.read_csv(self.fileName, **self.pd_kwargs)
        elif self.fileObject is not None:            
            self._tbl = pd.read_csv(self.fileObject, **self.pd_kwargs)
        else:
            raise IOError( 'Provide either fileName or fileObject' )

            
    @classmethod
    def gunzip(cls, inFile, outFile=None):
        """gunzip of a provided file"""
        if outFile is None:
            outFile, ext = os.path.splitext(inFile)
        print outFile

        
    #--- setters/getters ---#
    # readFile attr = downloaded read file
    def set_readFile(self, fileName): self.readFile = fileName        
    def get_readFile(self): return self.readFile
    

    
class MetaFile_MGRAST(MetaFile):
    """Subclass of MetaFile class for importing mgrast metafile """

    def __init__(self, fileName=None, fileObject=None, 
                 stages=[150, 100], **pandas_kwargs):
        MetaFile.__init__(self, fileName=fileName, fileObject=fileObject, 
                          stages=stages, **pandas_kwargs)
        self.checkFile()
    

    def checkFile(self):
        """Checking for needed columns in mg-rast metadata file"""
        try:
            self._tbl['id']
        except:
            raise IOError( 'MGRAST metaFile does not have "id" column' )


    def iterByRow(self):
        """Return iterator for processing table by row for mg-rast metadata table
        requires metaFile in pandas dataframe class
        
        Return:
        metagenome_id
        """
        reC = re.compile(r'((mgm)*\d\d\d\d\d\d\d\.\d)')
        
        for rowIndex, row in self._tbl.iterrows():
            metagenome_id = row['id']
            
            # check that metagenome ID is in correct format
            m = reC.match(metagenome_id)
            if not m:
                raise ValueError('id: "{0}" is not in correct format'.format(metagenome_id))
            else:
                metagenome_id = m.group(1)

            # create metaFile row class
            mfr = MetaFile_MGRAST_row(row, self.stages, rowIndex=rowIndex)
            
            # yield
            yield mfr


            
class MetaFile_MGRAST_row(object):
    """Just single row of metafile table"""

    def __init__(self, row, stages, rowIndex=None):
        """
        Attribs:
        row -- pandas dataframe row
        stages -- list of MGRAST processing stages
        """
        if type(row) is not pd.core.series.Series:
            raise TypeError('row must be pd.core.series.Series\n')
            
        self.ID = row['id']
        self.row = row
        self.rowIndex = rowIndex
        self.stages = stages
        
        # determining sequencing platform
        self.platform = self.getPlatform()

        
    def getPlatform(self):
        """Determining sequencing platform from metadata file.
        NoneType returned if the platform cannot be determined.

        Column(s) to check:
        seq_method
        
        Args:

        Return:
        string identifying sequencing platform.
        """
        
        # check for column
        try: 
            self.row['seq_method']  = str(self.row['seq_method'])
        except KeyError as err:
            raise type(err)(err.message + '. "seq_meth" column not in metadata file!')
            
        # determine sequencing platform
        platform = None
        if re.search(r'454|pyro', self.row['seq_method'], flags=re.I) is not None:
            platform = '454'
        elif re.search( r'illumina', self.row['seq_method'], flags=re.I) is not None:
            platform = 'illumina'
        elif re.search( r'sanger', self.row['seq_method'], flags=re.I) is not None:
            platform = 'sanger'
        # return
        return platform


    def download(self):
        """MG-RAST API request for downloading metagenome reads.
        Reads are downloaded as gzipped files.

        Attrib edit:
        readFile -- string with downloaded file name
        readFileFormat -- sequence file format
        """        

        # input check
        if self.ID is None:
            raise TypeError('ID cannot be None type')
        
        # trying each stage
        for stage in self.stages:
            ret = self._dlStage(stage, self.ID)            
            if ret:  # success!
                return ret
            else:
                sys.stderr.write(' Requested content was empty! Stage{} did not return anything\n'.format(stage))

        # if fail: set readFile as NoneType
        self.set_readFile(None)        
        return False
            

    def _dlStage(self, stage, ID, iteration=1, lastIter=9):
        """Downloading the fasta file based on a particular stage and iteration.
        Will recursively try next iteration until a sequence file is downloaded.

        Args:
        stage -- mgrast pipeline stage
        iteration -- iteration through the MGRAST pipeline
        lastIter -- the last iteration that will be tried

        Return:
        string with downloaded file name; Nonetype if no file
        """
        # lastIter
        if iteration > lastIter:
            sys.stderr.write(' Exceeded iterations. Giving up\n')
            return False
        
        # initialize url
        url = 'http://api.metagenomics.anl.gov/1/download/'
        url = url + ID + '?file={0}.{1}'.format(stage, iteration)    
        
        # send request to mgrast
        sys.stderr.write( 'For ID: "{0}", trying stage: "{1}"\n'.format(ID, stage) )
        sys.stderr.write( 'Sending request: "{0}"\n'.format(url) )
        req = requests.get(url)
        sys.stderr.write(' Request status: {0}\n'.format(str(req.status_code)))
        if req.status_code != 200:   # request failed; try next stage 
            sys.stderr.write('  Request != 200. Stage{} did not return a valid file!\n'.format(stage))
            return False
                
        # writing content to file
        outFile = ID + '_stage' + stage + '.fasta'
        d = zlib.decompressobj(16+zlib.MAX_WBITS)
        
        # trying to write
        with open(outFile, 'wb') as fd:
            for chunk in req.iter_content(10000):
                try:
                    fd.write( d.decompress(chunk) )
                except zlib.error:
                    fd.write( chunk )

        # determine sequence file format
        self._set_seqFormat(outFile)
                    
        # check that downloaded file is non-empty; if not trying next pipeline iteration
        if os.stat(outFile)[6] == 0:
            sys.stderr.write(' File empty. Trying next pipeline iteration!\n')            
            self._dlStage(stage, ID, iteration=iteration + 1)  # recursive
        elif not self.get_readFileFormat():
            sys.stderr.write(' File not a sequence file. Trying next pipeline iteration!\n')            
            self._dlStage(stage, ID, iteration=iteration + 1)  # recursive            
        else:
            sys.stderr.write(' File written: {0}\n'.format(outFile))            
            
        # setting readFile (outFile)
        self.set_readFile(outFile)
        return True


    def _set_seqFormat(self, inFile, nlines=100):
        """Determining the format of the seuqence file.

        Args:
        inFile -- file name
        nlines -- number of lines in file to check (starting from top)

        Attrib set:
        readFileForamt -- set to fasta or fastq or NoneType
        """
        with open(inFile, 'r') as fd:
            head = ''.join([fd.readline() for x in xrange(nlines)])

        # format?
        nseqs_fasta = len( [seq_rec.id for seq_rec in SeqIO.parse(StringIO(head), 'fasta')] )
        nseqs_fastq = len( [seq_rec.id for seq_rec in SeqIO.parse(StringIO(head), 'fastq')] )
        
        if nseqs_fasta > 0 and nseqs_fastq > 0:
            if nseqs_fasta > nseqs_fastq:
                self.set_readFileFormat('fasta')
            elif nseqs_fasta < nseqs_fastq:
                self.set_readFileFormat('fastq')
            else:
                raise IOError('  The file appears to be both fasta and fastq\n')                
        elif nseqs_fasta > 0:
            self.set_readFileFormat('fasta')
        elif nseqs_fastq > 0:
            self.set_readFileFormat('fastq')
        else:
            self.set_readFileFormat(None)
            
            
    def is_readFileEmpty(self):
        """Determine whether read file is empty.

        Return:
        True if empty or no read file"""
        if self.get_readFile() is None:
            return True
        elif os.stat(self.get_readFile())[6] == 6:
            return True
        else:
            return False
            
            
    def to_fasta(self, rmFile=False):
        """Converting from fastq to fasta.

        Args:
        rmFile -- remove old version of file?

        Attrib edit:
        readFile name set to new file (*.fasta)
        readFileFormat set to fasta format        
        """
        # unpack
        readFile = self.get_readFile()
        readFileFormat = self.get_readFileFormat()

        # new file name
        basename,ext = os.path.splitext(readFile)
        ## rename file to prevent overwrite
        if ext == '.fasta':
            os.rename(readFile, basename + '.tmp')
            readFile = basename + '.tmp'        
        newFile = basename + '.fasta'

        # convert
        SeqIO.convert(readFile, readFileFormat, newFile, 'fasta')

        # remove
        if rmFile:
            os.remove(readFile)

        # setting attributes
        self.set_readFile(newFile)
        self.set_readFileFormat('fasta')

        return True
        
            
    def get_ReadStats(self, readFile=None, fileFormat='fasta'):
        """Loading read file and getting stats:
        Length distribution: min, max, mean, median, stdev
        
        Args:
        readFile -- if None using self.readFile

        Output:
        readStats attribute -- dict of stats (eg., 'mean' : mean_value); returns None if no fileName provided
        """

        # args
        if readFile is None:
            try:
                readFile = self.get_readFile()
            except AttributeError:
                sys.stderr.write('No readFile found! Returning 0\n')
                return False
                
        # get seq lengths & stats
        seqLens = []
        for seq_record in SeqIO.parse(readFile, fileFormat):
            seqLens.append( len(seq_record) )

        if len(seqLens) == 0:
            return False
            
        # stats on read lengths
        self.readStats = {
            'n' : len(seqLens),
            'min' : min(seqLens),
            'mean' : scipy.mean(seqLens),
            'median' : scipy.median(seqLens),
            'max' : max(seqLens),
            'stdev' : scipy.std(seqLens)
            }
        return True


    #-- getters/setters --#
    def get_platform(self):
        return self.platform
    def get_ID(self):
        return self.ID

    # readFile attr = downloaded read file
    def set_readFile(self, fileName):
        self.readFile = fileName        
    def get_readFile(self):
        return self.readFile
    def set_readFileFormat(self, fileFormat):
        self.readFileFormat = fileFormat
    def get_readFileFormat(self):
        return self.readFileFormat
    def get_readCount(self):
        try:
            return self.readStats['n']
        except AttributeError or KeyError:
            return None
