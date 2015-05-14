#!/usr/bin/env python

"""
BaseSpace_get.py: get files from BaseSpace

Usage:
  BaseSpace_get.py project <Id> [options]
  BaseSpace_get.py sample <Id> [options]
  BaseSpace_get.py -h | --help
  BaseSpace_get.py --version

Options:
  <Id>            Project/Sample Id
  --user=<u>      BaseSpace user name
  --config=<cfg>  Config file path.
                  [default: ~/.basespace.cfg]
  --profile=<pf>  Profile to use from the config file.
                  [default: DEFAULT]
  --version       Show version.
  -h --help       Show this screen.
"""

from docopt import docopt
import sys,os
import configobj
from StringIO import StringIO
from configobj import ConfigObj
from validate import Validator
from BaseSpacePy.api.BaseSpaceAPI import BaseSpaceAPI
from BaseSpacePy.api.BaseSpaceException import CredentialsException


def get_configspec():
    configspec = """
    [__many__]
    clientKey    = string()
    clientSecret = string()     
    appSessionId = string()
    apiServer    = string(default=None)
    """
    return ConfigObj(StringIO(configspec))


def file_download(myAPI, sampleID):
    files = myAPI.getSampleFilesById(sampleID)
    for f in files:
        msg = 'Downloading: {}\n'
        sys.stderr.write(msg.format(f.Name))
        myAPI.fileDownload(f.Id, '.')
    
    

def main(uargs):    
    uargs['--config'] = os.path.abspath(os.path.expanduser(uargs['--config']))
    conf = ConfigObj(uargs['--config'], configspec=get_configspec())    

    try:
        conf_args = conf[uargs['--profile']]
    except KeyError:
        msg = 'Profile "{}" not found in config file'
        raise KeyError, msg.format(uargs['--profile'])
    
 
    # selecting config args
    keys = ['clientKey', 'clientSecret', 'apiServer', 'apiVersion', 'appSessionId', 'accessToken']
    conf_args = [conf_args[x] for x in keys]
    myAPI = BaseSpaceAPI(*conf_args)    

    # setting user
    if uargs['--user'] is not None:
        user = myAPI.getUserById(uargs['--user'])
    else:
        user = myAPI.getUserById('current')

    # all samples for project
    if uargs['project']:
        samples = myAPI.getSamplesByProject(uargs['<Id>'])
        for sample in samples:
            file_download(myAPI, sample.Id)

    if uargs['sample']:
        file_download(myAPI, uargs['<Id>'])
            
    

if __name__ == '__main__':
    uargs = docopt(__doc__, version='0.1')
    main(uargs)

