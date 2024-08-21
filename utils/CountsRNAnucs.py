#! /usr/bin/env python3

import os, sys, argparse
import pandas as pd
from datetime import datetime

def argParse():
    des='Get Nucleotide composition from sRNAs. Starting T with 3rd G etc.'
    parser = argparse.ArgumentParser(description=des,
                                     allow_abbrev=False,
                                     add_help=False,
                                     usage=argparse.SUPPRESS,
                                     formatter_class=argparse.RawTextHelpFormatter
                                     )

    req=parser.add_argument_group('Required')

    req.add_argument('-d','--directory', action='store', help='Directory containing the FastQ (ungzipped) files.', required=True)

    opt=parser.add_argument_group('Option')

    opt.add_argument('-s','--sizes', nargs='+', type=int, help='Sizes you wish to have checked. Default all FastQs in directory.')
    opt.add_argument('-h','--help', action='help', help='Display this message and exit.')

    args = parser.parse_args()

    return args, parser


def get_Files(path, sizes=None):
    if sizes:
        dir = [file for file in os.listdir(path) if (file.endswith('fq') or file.endswith('fastq')) and (int(file.split('bp')[0].split('_')[-1]) in sizes)]
    else:
        dir = [file for file in os.listdir(path) if (file.endswith('fq') or file.endswith('fastq'))]
    if not dir:
        print('ERROR! No FastQ files found! Make sure you have the correct path or the correct format.')
        sys.exit()
    return dir

def getSizes(path):
    sizes=sorted([int(i.split('bp')[0].split('_')[-1]) for i in get_Files(path)])
    return sizes

def getSeqs(file):
    infile=open(file)#'MACIESmapped_25bp_sRNAs.fq')

    seqs=[]
    header=infile.readline()
    while header:
        seqs.append(infile.readline().strip())
        PLUS=infile.readline()
        QUAL=infile.readline()
        header = infile.readline()
    return seqs

def getNucDist(seqs):
    dct={'T':0,'A':0,'G':0,'C':0, 'N':0}
    for i in seqs:
        dct[i[0]]+=1
    dct_norm={f'{k} %':v/sum(dct.values()) for k,v in dct.items()}
    return {**dct, **dct_norm}

def getTGperc(dct, seqs):
    T=0
    nT=0
    cntTG=0
    cntTnG=0
    cntnTG=0
    cntnTnG=0
    for i in seqs:
        if i.startswith('T'):
            T+=1
            if i[2].upper() == 'G':
                cntTG+=1
            else:
                cntTnG+=1
        else:
            nT+=1
            if i[2].upper() == 'G':
                cntnTG+=1
            else:
                cntnTnG+=1

    cntTG_perc=cntTG/len(seqs) #T#
    cntTnG_perc=cntTnG/len(seqs) #T#
    cntnTnG_perc=cntnTnG/len(seqs) #nT#
    cntnTG_perc=cntnTG/len(seqs) #nT#
    return {'T/G':cntTG_perc,'T/noG':cntTnG_perc,'noT/noG':cntnTnG_perc,'noT/G':cntnTG_perc}

def updateDataFrame(df,dct, dct_TGperc,  size):
    for k,v in dct.items():
        df.loc[size, k]=v
    for k,perc in dct_TGperc.items():
        df.loc[size, k]=perc

def masterFREQprogram(path, sizes):
    df=pd.DataFrame(index=sizes, columns=['T','A','G','C','N','T %','A %','G %','C %','N %','T/G','T/noG','noT/noG','noT/G'])

    dir=get_Files(path, sizes)
    for file in dir:
        size=int(file.split('bp')[0].split('_')[-1])
        seqs=getSeqs(file)
        dct=getNucDist(seqs)
        dct_TGperc=getTGperc(dct,seqs)
        updateDataFrame(df,dct,dct_TGperc,size)

    df.to_csv(f'{str(datetime.now()).split()[0].replace("-","")}_NucCompII.tsv', sep='\t')


if __name__ == '__main__':
    print('STARTING')
    args, parser = argParse()
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit()

    path=args.directory
    if args.sizes:
        sizes=args.sizes
    else:
        sizes=getSizes(path)

    masterFREQprogram(path, sizes)
    print('DONE')

