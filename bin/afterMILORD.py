#! /usr/bin/env python3
## Author: Sebastian Bechara
## Last Modified: 07.2024
## Mail: ask Mariusz

import os, sys, argparse
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from Bio import SeqIO
from datetime import datetime

def argParse():
    '''
    '''
    des = '\
This script will generate the cryptic and alternative excision plots when supplied\n\
with a GFF3 file produced by the ParTIES submodule MILORD. \n\b\
Usage:\n\t\
afterMILORD.py -g <inGFF3>'
    parser = argparse.ArgumentParser(description=des,
                                allow_abbrev=False,
                                add_help=False,
                                usage=argparse.SUPPRESS,
                                formatter_class=argparse.RawTextHelpFormatter)
    req = parser.add_argument_group('Required')
    req.add_argument('-g','--GFF',
                     action='store',
                     metavar='',
                     type=str,
                     help='Specify the GFF-file containing the excision info')
    
    opt = parser.add_argument_group('Optional')
    opt.add_argument('-c','--cutOff',
                     action='store',
                     metavar='',
                     type=float,
                     help='Specify cut-off to use for regarding an event as "true".\nCut-off is calculated by forming the ratio of support reads for\nthe variant and the reference. Default: 0.1')
    opt.add_argument('-o', '--outPath',
                     action='store',
                     metavar='',
                     type=str,
                     help='Specify the path where to save the output. If it does not\nexist it will be created. Default is current directory.')
    opt.add_argument('-f', '--format',
                     action='store',
                     metavar='',
                     type=str,
                     help='specify which format to save the plot in. \nAvailable: SVG (default), PNG, JPG')
    opt.add_argument('-t', '--TSV',
                     action='store_true',
                     help='Wether to write the specific location of excision events\nto file or not. It will not by default.')
    opt.add_argument('-h','--help',
                     action='help',
                     help='Print this message and exit.')
    # opt.add_argument('-d','--double',
    #                  action='store_true',
    #                  help='Wether to put naive and cut-off data on the same figure side by side or not')

    return parser

def checkArgs(parser, args):
    ARGS = {'infile':None,
            'outPath':'.',
            'cutOff':0.1,
            'F':'svg',
            'TSV':False} 
    ## required
    if not args.GFF:
        print('ERROR!! No GFF specified!')
        parser.print_help(sys.stderr)
        sys.exit()
    elif args.GFF.split('.')[-1] not in ['GFF3','GFF', 'gff', 'gff3']:
        print('ERROR!! Input file not supported! You need a GFF3 file!')
        parser.print_help(sys.stderr)
        sys.exit()
    else: ARGS['infile'] = args.GFF
    ## optional
    if args.outPath:
        if not os.path.isdir(args.outPath):
            os.mkdir(args.outPath)
            ARGS['outPath'] = args.outPath
        else: ARGS['outPath'] = args.outPath 
    if args.cutOff:
        if type(args.cutOff) != float:
            print('ERROR!! If you specify a custom cut-off, please specify a float!')
            parser.print_help(sys.stderr)
            sys.exit()
        else: ARGS['cutOff'] = args.cutOff
    if args.format:
        if args.format.upper() not in ['SVG','PNG','JPG']:
            print('ERROR!! Format not supported!!')
            parser.print_help(sys.stderr)
            sys.exit()
        else: ARGS['F'] = args.format.lower()
    if args.TSV:
        ARGS['TSV'] = True
    
    return ARGS
      

def mainLoop(infile, co=0.1, TSV=False, outPath='.', F='svg'):
    inFile=[i for i in open(infile) if not i.startswith('#')]
    cols=['dodgerblue','forestgreen','orange']

    fig, axs = plt.subplots(3,2,figsize=(9,6), layout='constrained')
    for i in range(2):
        types = ['IES', 'ALTERNATIVE', 'CRYPTIC']
        for index in range(3):
            inin = [i for i in inFile if i.split('type=')[-1].split(';')[0]== types[index]]
            l = [int(i.split('\t')[4])-int(i.split('\t')[3]) for i in inFile if i.split('type=')[-1].split(';')[0]== types[index]]
            sVar = [int(i.split('variant=')[-1].split(';')[0]) for i in inFile if i.split('type=')[-1].split(';')[0]== types[index]]
            sRef = [int(i.split('_ref=')[-1].split(';')[0]) for i in inFile if i.split('type=')[-1].split(';')[0]== types[index]]
            sRat = np.array(sVar) / (np.array(sVar)+np.array(sRef)) 
            if i == 0:     
                axs[index,i].hist(l , np.arange(5,150,1), edgecolor='white', label=f'{types[index]}', color=cols[index], snap=False)
                axs[index,i].set_ylabel(types[index])
                if index == 0:
                    axs[index,i].set_title('w/o cut-off')   
            else:
                axs[index,i].hist(np.array(l)[sRat>co] , np.arange(5,150,1), edgecolor='white', label=f'{types[index]}', color=cols[index], snap=False)
                if index == 0:
                    axs[index,i].set_title(f'w/ ratio cut-off > {co}')   
    ## write outRes
            if TSV:
                if i== 0:
                    with open(f'{outPath}/{types[index]}_ResLoc.tsv', 'w') as w:
                        out=['\t'.join([i.split('\t')[0], i.split('\t')[3], i.split('\t')[4]]) for i in inFile if i.split('type=')[-1].split(';')[0]== types[index]]
                        w.write('\n'.join(out)+'\n')
                else:
                    with open(f'{outPath}/{types[index]}_ResLocCutOff.tsv', 'w') as w:
                        ll = np.array(inin)[sRat>0.1]
                        out=['\t'.join([i.split('\t')[0], i.split('\t')[3], i.split('\t')[4]]) for i in ll if i.split('type=')[-1].split(';')[0]== types[index]]
                        w.write('\n'.join(out)+'\n')
    fig.supxlabel('Size of exision event', fontsize=16)
    fig.supylabel('Occurrence', fontsize=16)
    plt.savefig(f'{outPath}/CrypticPlot.{F}')

def getTime():
    return str(datetime.now()).split('.')[0].replace('-','.').replace(' ','|')

def main(args=None):
    parser= argParse()
    num=1
    if args:
        num=2
    if not args:
        args = parser.parse_args()
    if len(sys.argv) == num:
        parser.print_help(sys.stderr)
        sys.exit()

    ARGS = checkArgs(parser, args)
    infile = ARGS['infile']
    outPath = ARGS['outPath']
    CO = ARGS['cutOff']
    F = ARGS['F']
    TSV = ARGS['TSV']
    print(f'[{getTime()}] - Running . . . ')    
    mainLoop(infile, CO, TSV, outPath, F)
    print(f'[{getTime()}] - Done!')

if __name__ == '__main__':
    main()
    