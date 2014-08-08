import sys
import os
import re
#import urllib
#import urllib2
import requests

import pandas as pd


class metaFile(object):
    """metaFile object class
    """
    
    def __init__(self, fileName=None, fileObject=None, 
                 seqDB='MGRAST', stages=[150, 100], **pandas_kwargs):
        """Loading file as pandas dataFrame. Checking for necessary columns
        
        Kwargs:
        seqDB -- datafile type: 'MGRAST' or 'SRA'
        stages -- mg-rast processing stage for downloading files. If multiple stages provided,
          will try each in succession until an entry is returned.
        """

        # attributes
        self.fileName = fileName
        self.fileObject = fileObject
        self.seqDB = seqDB
        self.stages = [str(x) for x in stages]
        

        # checking seqDB
        if self.seqDB is not 'MGRAST' and self.seqDB is not 'SRA':
            raise IOError( 'seqDB must be "MGRAST" or "SRA"' )

        # loading file
        sep = '\t'
        if fileName is not None:
            self._tbl = pd.read_csv(fileName, **pandas_kwargs)
        elif fileObject is not None:            
            self._tbl = pd.read_csv(fileObject, **pandas_kwargs)
        else:
            raise IOError( 'Provide either fileName or fileObject' )

        # check for needed columns
        ## MGRAST
        ### columns needed: id
        if self.seqDB is 'MGRAST':
            try:
                self._tbl['id']
            except:
                raise IOError( 'MGRAST metaFile does not have "id" column' )

        ## SRA
        ### columns needed: ftp
        elif self.seqDB is 'SRA':
            try:
                self._tbl['ftp']
            except:
                raise IOError( 'SRA metaFile does not have "ftp" column' )

        ## seqDB not recognized
        else:
            raise IOError( 'seqDB "{0}" is not recognized'.format(self.seqDB) )


    def iterByRow(self):
        """Return iterator for processing table by row
        
        Return:
        metagenome_id
        """
        reC = re.compile(r'((mgm)*\d\d\d\d\d\d\d\.\d)')
        
        for count, row in self._tbl.iterrows():
            metagenome_id = row['id']

            # check that metagenome ID is in correct format
            m = reC.match(metagenome_id)
            if not m:
                raise ValueError('id: "{0}" is not in correct format'.format(metagenome_id))
            else:
                metagenome_id = m.group(1)

            # yield
            yield metagenome_id
            

    def download(self, ID):
        """MG-RAST API request for downloading metagenome

        Kwargs:
        ID -- metagenome ID
        """        
        
        # trying each stage
        for stage in self.stages:
            # initialize url
            url = 'http://api.metagenomics.anl.gov/1/download/'
            url = url + ID
            url = url + '?stage={0}'.format(stage)            

            # send request to mgrast
            sys.stderr.write( 'Sending request: "{0}"\n'.format(url) )
            req = requests.get(url)
            
            ## need to use curl get request to download as compressed file

            print(req.status_code)
            #urllib.urlretrieve(url, ID)
            #response = urllib2.urlopen(url)
            #print(response)
            #html = response.read()
            #print(url)
        
