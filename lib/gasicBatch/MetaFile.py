import sys
import os
import re
import gzip
import zlib

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
        outFile -- string with downloaded file name
        """        

        # input check
        if self.ID is None:
            raise TypeError('ID cannot be None type')
        
        # trying each stage
        for stage in self.stages:
            outFile = self._dlStage(stage, self.ID)            
            if outFile:
                self.set_readFile(outFile)
                return 1
            else:
                sys.stderr.write(' Requested content was empty! Stage{} did not return anything\n'.format(stage))
        return 0
            

    def _dlStage(self, stage, ID, iteration=1, lastIter=9):
        """Downloading the fasta file based on a particular stage and iteration.

        Args:
        stage -- mgrast pipeline stage
        iteration -- iteration through the MGRAST pipeline
        lastIter -- the last iteration that will be tried
        """
        # lastIter
        if iteration > lastIter:
            sys.stderr.write(' Exceeded iterations. Giving up\n')
            return 0
        
        # initialize url
        url = 'http://api.metagenomics.anl.gov/1/download/'
        url = url + ID + '?file={0}.{1}'.format(stage, iteration)    
        
        # send request to mgrast
        sys.stderr.write( 'For ID: "{0}", trying stage: "{1}"\n'.format(ID, stage) )
        sys.stderr.write( 'Sending request: "{0}"\n'.format(url) )
        req = requests.get(url)
        sys.stderr.write(' Request status: {0}\n'.format(str(req.status_code)))
        if req.status_code != 200:   # try next stage 
            sys.stderr.write('  Request != 200. Stage{} did not return a valid file!\n'.format(stage))
            return 0
                
        # writing content to file
        outFile = ID + '_stage' + stage + '.fasta'
        d = zlib.decompressobj(16+zlib.MAX_WBITS)
        
        # trying to write
        try:
            with open(outFile, 'wb') as fd:
                for chunk in req.iter_content(1000):
                    fd.write( d.decompress( chunk) )
            sys.stderr.write(' File written: {0}\n'.format(outFile))
        except:  # trying next iteration
            sys.stderr.write(' File was not a compress sequence file! Trying next pipeline iteration!\n')
            self._dlStage(stage, ID, iteration=iteration + 1)

        # checking that content was written
        ## else: try next stage
        if os.stat(outFile)[6] != 0:
            #self.set_readFile(outFile)
            return outFile
        else:
            return 0

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
                return 0
                
        # get seq lengths & stats
        seqLens = []
        for seq_record in SeqIO.parse(readFile, fileFormat):
            seqLens.append( len(seq_record) )
        
        # stats on read lengths
        self.readStats = { 
            'min' : min(seqLens),
            'mean' : scipy.mean(seqLens),
            'median' : scipy.median(seqLens),
            'max' : max(seqLens),
            'stdev' : scipy.std(seqLens)
            }
        return 1


        

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

        