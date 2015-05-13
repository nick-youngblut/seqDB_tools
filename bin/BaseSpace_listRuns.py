#!/usr/bin/env python

"""
BaseSpace_listRuns.py: list runs for a User

Usage:
  BaseSpace_listRuns.py [options]
  BaseSpace_listRuns.py -h | --help
  BaseSpace_listRuns.py --version

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
    myProjects = myAPI.getProjectByUser()
    print "The projects for this user are:" 
    print '\n'.join([str(x) for x in myProjects])

    # user runs
    runs = user.getRuns(myAPI)
    print "The runs for this user are:"
    print '\n'.join([str(x) for x in runs])


if __name__ == '__main__':
    uargs = docopt(__doc__, version='0.1')
    main(uargs)

