#!/usr/bin/env python

"""
BaseSpace_summary_format.py: formatting copied & pasted run summary table from BaseSpace

Usage:
  BaseSpace_summary_format.py [options] <table_files>...
  BaseSpace_summary_format.py -h | --help
  BaseSpace_summary_format.py --version

Options:
  <table_files>   >=1 summary table file.
  -o=<o>          Basename of output files.
                  [Default: BaseSpace_summary]
  --version       Show version.
  -h --help       Show this screen.
"""

from docopt import docopt
import sys,os
import re
import pandas as pd

def parse_file_name(file_name):
    x = file_name.split('__',2)
    try:
        len(x[1])
    except IndexError:
        x.append('NA')
    assert len(x)==2, 'File name parsing error: "{}"'.format(file_name)
    return x


def read_file(inFile):
    data = []

    if inFile == '-':
        iFH = sys.stdin
    else:
        iFH = open(inFile, 'rb')    
    for line in iFH:
        line = line.rstrip()
        if re.match(r'^\s*$', line):
            continue
        line = [x.lstrip().rstrip() for x in line.split('\t')]
        line = [re.sub(' +', ' ', x) for x in line]
        line = [re.sub('[()]', '', x) for x in line]

        data.append(line)        

    iFH.close()
    return data


def parse_data(data):
    # parsing tables
    tbl1 = []
    tbl2 = []
    end = 0
    level = ''
    for line in data:
        if end == 0:
            tbl1.append(line)
        elif end == 1:
            if line[0] == 'Lane':
                line[0] = 'Level'
            if len(line) == 1 and line[0] != 'Level':
                level = line[0]
                continue
            else:
                line[0] = level
            tbl2.append(line)
        if line[0] == 'Total':
            end = 1
    return tbl1, tbl2


def write_table(fh, tbl, file_name, runID, date, header=False):    
    # header
    if header:
        x = '\t'.join(['file_name','runID','run_date'] + tbl[0])
        fh.write(x + '\n')
    # body
    for line in tbl[1:]:
        x = '\t'.join([file_name, runID, date] + line)
        fh.write(x + '\n')


def main(uargs):
    outbase = uargs['-o']
    outFile1 = outbase + '_tbl1.txt'
    outFile2 = outbase + '_tbl2.txt'
    oFH1 = open(outFile1, 'wb')
    oFH2 = open(outFile2, 'wb')
    header = True
    for file_name in uargs['<table_files>']:
        msg = 'Cannot find: "{}"'
        assert os.path.isfile(file_name), msg.format(file_name)

        runID,date = parse_file_name(file_name)
        data = read_file(file_name)
        tbl1,tbl2 = parse_data(data)

        write_table(oFH1, tbl1, file_name, runID, date, header=header)
        write_table(oFH2, tbl2, file_name, runID, date, header=header)
        header = False
        
    oFH1.close()
    oFH2.close()

    sys.stderr.write('File written: {}\n'.format(outFile1))
    sys.stderr.write('File written: {}\n'.format(outFile2))
    

if __name__ == '__main__':
    uargs = docopt(__doc__, version='0.1')
    main(uargs)

