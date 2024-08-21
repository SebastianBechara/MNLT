#! /usr/bin/env python3
# Author: Sebastian Bechara
# Last modified: 04.2024
import sys , os , argparse
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime

def argParse():
    '''
    Generates Parser using argparser
    Returns the arguments specified by comman-line call and the parser itself
    '''
    des=f'\n{"#"*90}\n\nScript to map presorted sRNAs and to generate SIMPLE (MDS, IES, Other, Unknown)\n\
or COMPLEX (Vector, MDS, IES, OES, Mitochondria, Klebsiella, rRNA, Unknown) maping distribution.\n\
Reads need to be split into corresponding size, their size in the file name (e.g. _15bp.fastq) and NOT compressed.\n\n\
Dependencies: Hisat2, Python-libraries (matplotlib; numpy)\n\
\tUsage: sRNAHistSimple.py -d <directory> -min <min_size> -max <max_size>\n\
\t       -o <output_directory> -t <plot_title> -p <simple|complex|both> -s <separator>\n      '

    parser = argparse.ArgumentParser(description=des,
                        allow_abbrev=False,
                        add_help=False,
                        usage=argparse.SUPPRESS,
                        formatter_class=argparse.RawTextHelpFormatter)
    par = parser.add_argument_group('Arguments')
    par.add_argument('-d','--inDir',
                        metavar='',#argparse.SUPPRESS,#'\b',
                        action='store',
                        type=str,
                        help='[Required] Directory containing sRNA sequencing data (not compressed)')
    par.add_argument('-min','--MinSize',
                        metavar='',#argparse.SUPPRESS,#'\b',
                        action='store',
                        type=int,
                        default=15,
                        help='[Optional] Minimum Insert size to use (default = 15)')
    par.add_argument('-max','--MaxSize',
                        metavar='',#argparse.SUPPRESS,#'\b',
                        action='store',
                        type=int,
                        default=35,
                        help='[Optional] Maximum Insert size to use (default = 35)')
    par.add_argument('-o','--outDir',
                        metavar='',#argparse.SUPPRESS,#'\b',
                        action='store',
                        type=str,
                        help=f'[Optional] Name for the directory that is going to be generated.\n{" "*11}If "." then working directory is used as out directory (default)')
    par.add_argument('-t','--title',
                        metavar='',#argparse.SUPPRESS,#'\b',
                        action='store',
                        type=str,
                        help=f'[Optional] Title of the plot.\n{" "*11}If not set, base filename will be displayed')
    par.add_argument('-p','--plotComp',
                        #metavar='\b',
                        metavar='',#argparse.SUPPRESS,#'\b',
                        action='store',
                        default='simple',
                        type=str,
                        help=f'[Optional] Whether to have the complex or simple plot, or both.\n{" "*11}Specify by one of these: "simple", "complex", "both" (default: simple)')
    par.add_argument('-s', '--sep',
                        metavar='',#argparse.SUPPRESS,#'\b',
                        action='store',
                        type=str,
                        default='bp.fastq',
                        help=f'[Optional] File endings after the size number, used for separating\n{" "*11}and filtering (Default: bp.fastq, so bp is the effective separator)')
    par.add_argument('-th', '--threads',
                        metavar='',#argparse.SUPPRESS,#'\b',
                        action='store',
                        type=int,
                        default=10,
                        help='[Optional] Specify the amount of threads to use while mapping data (default = 10)')
    par.add_argument('-k','--keep_tmp',
                        #metavar='\b',
                        action='store_true',
                        help=f'[Optional] Specify whether to keep tmp directory.\n{" "*11}Contains intermediate sam files (removed by default)')
    par.add_argument('-h', '--help',
                        #metavar='\b',
                        action='help',
                        help='Display this message and exit')
    return  parser

def checkArgs(args, parser):
    '''
    Checks the arguments given to the parser and prepares arguments w/o a default value.
    Returns specified input and output directories and the title of the plot to generate.
    '''
    if not args.inDir:# or not args.MinSize or not args.MaxSize:
        print('Error! There are required arguments missing\n################\n\nUsage:\n\n')
        parser.print_help(sys.stderr)
        sys.exit()
    if len(args.inDir.split('/')) == 1:
        inDir = f'{os.getcwd()}/{args.inDir}'
    else: inDir = args.inDir
    if args.outDir:
        outDir = args.outDir
    else:
        outDir = os.getcwd()
    if args.title:
        title = args.title
    else:
        title = inDir.split('/')[-1]
    return inDir, outDir, title

def checkFolders(outDir=None):
    '''
    Creates Directories including the tmp directory and an output directory if specified.
    '''
    if os.path.isdir('tmp') != True:
        os.mkdir('tmp')
    if os.path.isdir("tmp/mapped") != True:
        os.mkdir("tmp/mapped")
    if os.path.isdir("tmp/unmapped") != True:
        os.mkdir("tmp/unmapped")
    if outDir and outDir != os.getcwd():
        if not os.path.isdir(outDir):
            os.mkdir(outDir)

def getSeqFiles(parser, direc, minS, maxS, sep):
    '''
    Fetches sequencing files present in the given input directory.
    Takes the directory as well as minimum and maximum size of sRNAs to use as input.
    Also take a separator as input, in order to be able to properly bin the files,
    this needs to be something like _15bp.fastq; meaning '_' then the size then the separator.
    Separator is also used as file ending, so it essentially needs to be xxx.fastq.
    Returns the sequencing files to use as list.
    '''
    seqFiles=[]
    # print(f'getSeqFiles {direc}')
    for file in os.listdir(direc):
        # print(file, sep)
        # if sep:
        if file.endswith(sep):
            ss=sep.split('.')[0]
            if int(minS) <= int(file.split('_')[-1].split(ss)[0]) <= int(maxS):
                seqFiles.append(file)
    return seqFiles

def mapSeq(seqFiles, directory, threads):
    '''
    This maps the reads to the references and saves mapped and unmapped read for each reference
    in the tmp/mapped|unmapped directory.
    System call to Hisat2. By default it uses 10 threads.  
    '''
    tags=['vec', 'mac', 'ies', 'oes', 'mito', 'kleb', 'rRNA']
    refs=['L4440/vector','para_MAC/mac','para_MIC/ies','new_para_MIC/oes', 'mito/mito', 'klebs/klebs', 'para_rRNA/rRNA']
    files = seqFiles
    p = '/'.join(os.path.abspath(__file__).split('/')[:-2])
    for file in files:
        for i in range(len(tags)):
            if i==0:
                cmd=f'hisat2 -p {threads} --un tmp/unmapped/{file}.{tags[i]}_unmapped.fastq --al tmp/mapped/{file}.{tags[i]}_mapped.fastq -x {p}/res/ref/hisat2/{refs[i]} -U {directory}/{file} > tmp/shite{tags[i]}.sam'
            else:
                cmd=f'hisat2 -p {threads} --un tmp/unmapped/{file}.{tags[i]}_unmapped.fastq --al tmp/mapped/{file}.{tags[i]}_mapped.fastq -x {p}/res/ref/hisat2/{refs[i]} -U tmp/unmapped/{file}.{tags[i-1]}_unmapped.fastq > tmp/shite{tags[i]}.sam'
            os.system(cmd)
    
def countReads(file, sep):
    '''
    Counts reads of specified file and return a key count pair for the following dictionary.
    '''
    s = sep.split('.')[0]
    with open(file, 'r') as infile:
        key = file.split(s)[0].split('_')[-1]
        cnt=0
        for line in infile:
            cnt+=1
        if cnt%4 != 0:
            print('ERROR! There was a FastQ file that was not divisible by 4. File in question: {file}')
            sys.exit()
        else:
            return int(key), cnt/4

def calcFreq(minS, maxS, sep):
    '''
    This calculates the actual frequencies of reads mapping to each reference.
    Takes threshhold sizes and the separator as input and returns a dictionary containing the ferequencies.
    '''
    tags=['vec', 'mac', 'ies', 'oes', 'mito', 'kleb', 'rRNA', 'unmapped']
    dcts = {tag: {n:0 for n in range(minS,maxS+1)} for tag in tags }
    for file in os.listdir('tmp/mapped'):
        tag = file.split('_mapped')[0].split('.')[-1]
        if tag in tags:
            key_cnt = countReads('tmp/mapped/'+file, sep)
            dcts[tag][key_cnt[0]] += key_cnt[1]
    
    for file in os.listdir('tmp/unmapped'):
        if 'rRNA_unmapped' in file:
            key_cnt = countReads('tmp/unmapped/'+file, sep)
            dcts['unmapped'][key_cnt[0]] += key_cnt[1]
    
    totReads = sum( [ sum(dct.values()) for tag, dct in dcts.items() ] )
    for tag , dct in dcts.items():
        for key, cnt in dct.items():
            dcts[tag][key] = cnt / totReads
    return dcts

def writeTSV(dcts,title, indir, outdir):
    '''
    Writes file containing the frequencies of each size and the corresponding reference.
    Will have the same name as the title for the plots and a .log ending, but is essentially a text file.
    '''
    tab = '\t'
    d = {i:[] for t, dct in dcts.items() for i in dct.keys()}
    for t, dct in dcts.items(): # t==tag; d==dct
        for b, c in dct.items(): # b==bp; c==cnt 
            d[b].append(c)

    with open(outdir+'/'+'_'.join(title.split())+'.log', 'w') as w:
        w.write(f'# Used directory: {indir}\n#\n#\n\n')
        w.write('sRNA size\tVector\tMDS\tIES\tOES\tMitochondria\tKlebsiella\trRNA\tUnmapped\n')
        for bp, lst in d.items():
            w.write(f'{bp}\t{tab.join([str(i) for i in lst])}\n')

def extYvals(dct):
    '''
    Extracts the y values for plotting from nested dictionary.
    '''
    ll = []
    for k,d in dct.items():
        l=[]
        for kk, v in d.items():
            l.append(v)
        ll.append(np.array(l))
    # print('extract y vals\n', ll)
    return ll

def getSimpleData(dcts):
    '''
    Generates the y values for simple plots.
    '''
    # print('SIMPLE DATA', '\n', dcts)
    cols=['forestgreen','red','grey','lightgrey']
    lgd=['MDS','IES','Other','Unmapped']
    dctSimple={'mac':None, 
                'ies':None,
                'other':{i:0 for i in dcts['vec'].keys()},
                'unmapped':{i:0 for i in dcts['vec'].keys()}}
    for tag, d in dcts.items():
        if tag in dctSimple.keys():
            dctSimple[tag] = d
        elif tag != 'unmapped':
            for bp , cnt in d.items():
                dctSimple['other'][bp]+=cnt
        else:
            dctSimple['unmapped']=d    
    y_val = extYvals(dctSimple) #[np.array(list(vals)) for t,d in dctSimple.items() for vals in d.values()]
    return y_val

def makePlots(x_val, y_val, cols, title, outdir,lgd, comp):
    '''
    This generates the plots and saves them either in the CWD or the specified out directory.
    The y-axis limit is set by default to (0, 0.8)
    '''
    plt.figure(figsize=(9,6))
    for i in range(len(y_val)):
        if i==0:
            plt.bar(x_val, y_val[i],color=cols[i])
        else:
            plt.bar(x_val, y_val[i],color=cols[i], bottom=np.sum(y_val[i-1::-1], axis=0))
    
    plt.xlabel('sRNA length [bp]', fontsize=14)
    plt.xlabel('Reads (normalised to total reads)', fontsize=14)
    plt.ylim(0,0.8)
    plt.legend(lgd, loc='best')
    # if title:
    #     ttl = title
    # else:
    #     ttl = dir.split('/')[-1]
    plt.title(title, fontsize=16)
    # plt.show()
    if comp == 'simple':
        plt.savefig(f'{outdir}/{"_".join(title.split())}_simple.svg', bbox_inches='tight')
    if comp == 'complex':
        plt.savefig(f'{outdir}/{"_".join(title.split())}_complex.svg', bbox_inches='tight')

def makeRes(dcts, minS, masS, outdir, indir, title, comp):
    '''
    Generates all the results, so log file and plots.
    '''
    writeTSV(dcts,title, indir, outdir)

    x_val = np.array(list(dcts['vec'].keys()))
    if comp=='simple':
        cols=['forestgreen','red','grey','lightgrey']
        lgd=['MDS','IES','Other','Unmapped']
        y_val = getSimpleData(dcts)
        makePlots(x_val, y_val, cols, title, outdir, lgd, comp) 
    elif comp=='complex':
        lgd=['Vector', 'MDS','IES','OES','Mitochondria','Klebsiella','rRNA','Unmapped']
        cols=['darkviolet', 'forestgreen','red','dodgerblue','orange','black','gold','grey']
        y_val = extYvals(dcts)#[np.array(list(vals)) for t,d in dcts.items() for v in d.values()]
        makePlots(x_val, y_val, cols, title, outdir, comp)
    elif comp=='both':
        y_valSimple = getSimpleData(dcts)
        colsSimple=['forestgreen','red','grey','lightgrey']
        lgdSimple=['MDS','IES','Other','Unmapped']
        makePlots(x_val, y_valSimple, colsSimple, title, outdir, lgdSimple, 'simple')
        lgdComp=['Vector', 'MDS','IES','OES','Mitochondria','Klebsiella','rRNA','Unmapped']
        colsComp=['darkviolet', 'forestgreen','red','dodgerblue','orange','black','gold','grey']
        y_valComp = extYvals(dcts)#[np.array(list(vals)) for t,d in dcts.items() for v in d.values()]
        makePlots(x_val, y_valComp, colsComp, title, outdir, lgdComp, 'complex')

def main(args=None):
    parser = argParse()
    num=1
    if args:
        num=2
    if not args:
        args = parser.parse_args()
    if len(sys.argv) == num:
        parser.print_help(sys.stderr)
        sys.exit()
    
    def getTime():
        return str(datetime.now()).split('.')[0].replace('-','.').replace(' ','|')
        
    print(f'\n\n[{getTime()}] - Starting Analysis.\n{parser.description}')
    print(f'[{getTime()}] - Preparing arguments and generating directories')
    inDir, outDir, title = checkArgs(args, parser)
    # print(args)
    checkFolders(outDir)
    print(f'[{getTime()}] - Fetching sequencing files')
    seqFiles = getSeqFiles(parser, inDir, args.MinSize, args.MaxSize, args.sep)  #(parser, direc, minS, maxS, sep)
    #print(seqFiles)
    print(f'[{getTime()}] - Mapping reads')
    mapSeq(seqFiles, inDir, args.threads) ### change function to work
    print(f'[{getTime()}] - Calculating frequencies')
    dcts = calcFreq(args.MinSize, args.MaxSize, args.sep)  ### change: remove default sep = None
    # print(dcts)
    print(f'[{getTime()}] - Generating results')
    makeRes(dcts, args.MinSize, args.MaxSize,outDir, inDir, title, args.plotComp)  # (dcts, minS, masS, outdir, indir, title, comp)
    if args.keep_tmp:
        print(f'[{getTime()}] - Keeping tmp directory')
        sys.exit()
    else: 
        print(f'[{getTime()}] - Removing tmp directory')
        os.system('rm -r tmp')
        sys.exit()

if __name__ == '__main__':
    main()
    