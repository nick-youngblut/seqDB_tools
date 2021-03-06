#!/usr/bin/env python

"""
BaseSpace_list.py: list user information

Usage:
  BaseSpace_list.py [options]
  BaseSpace_list.py -h | --help
  BaseSpace_list.py --version

Options:
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

    # user projects
    projects = myAPI.getProjectByUser()
    print "## The projects for this user are:" 
    project_headers = ['ProjectID','Name','UserOwnedBy','DateCreated']
    sample_headers = ['SampleID','NumReadsRaw','NumReadsPF','IsPairedEnd','Status']
    print '\t'.join(project_headers + sample_headers)
    for project in projects:
        project_data = [project.Id,
                         project.Name,
                         project.UserOwnedBy,
                         project.DateCreated]

        samples = myAPI.getSamplesByProject(project.Id)
        for sample in samples:
            sample_data = [sample.Id,
                           sample.NumReadsRaw,
                           sample.NumReadsPF,
                           sample.IsPairedEnd,
                           sample.Status]
            print '\t'.join([str(x) for x in project_data + sample_data])
    
            #files = myAPI.getSampleFilesById(sample.Id)
            #for f in files:
                #fo = myAPI.getFileById(f.Id)
            #    myAPI.fileDownload(f.Id, '.')
                

    # user runs
    print "## The runs for this user are:"
    print '\t'.join(['RunID','ExperimentName','UserOwnedBy',
                     'DateCreated', 'Status'])
    runs = user.getRuns(myAPI)
    for run in runs:
        data = [run.Id,
                run.ExperimentName,
                run.UserOwnedBy,
                run.DateCreated,
                run.Status]
        print '\t'.join([str(x) for x in data])

    

if __name__ == '__main__':
    uargs = docopt(__doc__, version='0.1')
    main(uargs)

