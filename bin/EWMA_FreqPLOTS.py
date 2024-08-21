#! /usr/bin/env python3
## Author: Sebastian Bechara
## Last Modified: 06.2024

import sys, os, argparse, shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

def argParse():
    '''
    Generates Parser using argparser
    Returns the arguments specified by comman-line call and the parser itself
    '''
    des='\
This script will generate the rolling mean IRS plots and the base frequency plots\n\
as found in Swart et al. 2017.\n\
You may generate the plots for multiple KDs, but be aware that you will get quite a few\n\
plots. Also the rolling mean over the IRSs may get overly cluttered if you put to many\n\
experiments in them.'
    parser = argparse.ArgumentParser(description = des,
                                     allow_abbrev = False,
                                     add_help = False,
                                     usage=argparse.SUPPRESS,
                                     formatter_class=argparse.RawTextHelpFormatter)
    req = parser.add_argument_group('Required')
    req.add_argument('-s','--switch',
                     action = 'store',
                     metavar = '',
                     type = str,
                     help = 'Specify what to do. Either "EWMA", "FREQ" or "BOTH"')
    req.add_argument('-i', '--IRS',
                     action = 'store',
                     metavar = '',
                     type = str,
                     help = 'File containing all calculated IRSs for wished KDs. Typically a tsv file')
    req.add_argument('-k', '--KDs',
                     nargs = '+',
                     metavar = '',
                     help = 'Specify KDs for which to generate the plots. If multiple KDs are wished \nprovide them separated by a whitespace. They need to be supplied exactly as they are \nin the IRS file (case-sensitive)')

    opt = parser.add_argument_group('Optional')
    opt.add_argument('-a','--IES_annotaion',
                     action='store',
                     metavar='',
                     type=str,
                     help='Annotation file containing info on IESs. If not specified will use the one in \nthe resource directory')
    opt.add_argument('-p','--peaks',
                     nargs='+',
                     metavar='',
                     help='[Freq] Specify IES-peaks to plot. Provide them separated by a whitespace,\nordered and multiples of 2. Peaks will be infered by input, e.g. if input = "15 25 98\n1110 2000 4000" you will get three peaks with (15,25), (98,1110),(2000,4000).\nDefault are (26,31),(44,52)&(54,62)')
    opt.add_argument('-o','--outDir',
                    action='store',
                    metavar='',
                    type=str,
                    help='Specify directory & path in which to save the plots. If it does not exist it will\nbe created. Default is current work directory')
    # opt.add_argument('-l','--labels',
    #                  actrion='store_true',
    #                  metavar='',
    #                  help='[Freq] Specify wether to denote the corresponding nucleotide in the plot')
    opt.add_argument('-h','--help',
                     action='help',
                     help='Display this message and exit'
                     )

    return parser

def checkArgs(args, parser):
    '''
    Checks the arguments given to the parser 
    '''
    if not args.switch or (args.switch not in ['BOTH', 'EWMA', 'FREQ']):
        print(args.switch, len(args.switch))
        print('ERROR! You need to specify what to do!\n\n####################\n\n')
        parser.print_help(sys.stderr)
        sys.exit()
    if not args.IRS:
        print('ERROR! No IRS file specified!\n\n####################\n\n')
        parser.print_help(sys.stderr)
        sys.exit()
    if not args.KDs:
        print('ERROR! You need to specify at least one KD!\nOtherwise your plots are going to be quite empty.\n\n####################\n\n')
        parser.print_help(sys.stderr)
        sys.exit()
    if args.peaks:
        if len(args.peaks) % 2 != 0:
            print('ERROR! You need to specify start AND end of peaks to concider!\n\n####################\n\n')
            parser.print_help(sys.stderr)
            sys.exit()

def splitPeaks(peaks): 
    '''
    Splits the list provided by the user into "peaks"; meaning tuples of two
    '''
    l=[]
    t=()
    for i in range(len(peaks)):
        if (i+1) %2 == 0:
            t += (int(peaks[i]),)
            l.append(t)
            t = ()
            continue
        t += (int(peaks[i]),)
    return l

def filterDF(df, KDs, ctrls):
    '''
    Filters the provided DataFrame by the controls present in the IRS file
    '''
    for spl in KDs:
        for i in df.index:
            cFlag=False
            for ctrl in ctrls:
                if df.loc[i,ctrl] > 0.05:
                    cFlag=True
            if cFlag:
                df.loc[i, spl] = np.nan
            # else:
            #     continue
    return df

def prepData(IRSfile, KDs):
    '''
    Creates DataFrame that is used in later functions
    '''
    df = pd.read_table(IRSfile, sep='\t', header = 0)
    ctrls = [i for i in df.columns if 'Ctrl' in i or 'ctrl' in i]
    df = df[['ID','Length']+KDs+ctrls]
    df = filterDF(df, KDs, ctrls)
    return df

def prep1(df,lst):
    '''
    Additional preparation step needed for the EWMA calculation
    '''
    DF = pd.DataFrame()
    for spl in lst:
        dct={k:np.nanmean(g[spl].tolist()) for k,g in df.groupby('Length')}
        DF.index=dct.keys()
        DF[spl]=dct.values()
    return DF

def calcEWMA(df, KDs):
    '''
    Calculates EWMAs for the provided KDs 
    '''
    for spl in KDs:
        df[f'{spl} 5EWMA'] = df[spl].ewm(span=5).mean()
        df[f'{spl} 50EWMA'] = df[spl].ewm(span=50).mean()
    return df

def generateEWMA(df, lst, outPath='.'):
    '''
    Generates the plots for the EWMA calculation
    '''
    df = calcEWMA( prep1(df, lst) ,lst)
    colours=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    for i in [5, 50]:
        fig = plt.figure()
        ax = fig.add_axes([0.1,0.1,0.8,0.8])
        for cnt, spl in enumerate(lst):
            # plt.scatter(x=df2_1['Length'], y=df2_1[spl], color=colours[cnt], s=1, alpha=0.005)
            plt.scatter(x=df.index, y=df[spl], color=colours[cnt], s=1, alpha=0.2)
            plt.plot(df.index, df[f'{spl} {i}EWMA'], linewidth=1, color=colours[cnt], label=spl) # dff
        if i == 5:
            plt.xlim(22,150)
            ax.set_xticks([25]+[i for i in range(40,150, 20)]+[150])
            ax.set_xticklabels([25]+[i for i in range(40,150, 20)]+[150])
        else:
            plt.xlim(22,1000)
            ax.set_xticks([i for i in range(20,1001, 100)]+[1000])
            ax.set_xticklabels([i for i in range(20,1001, 100)]+[1000])
        plt.xlabel('IES Length [bps]', fontsize=16)
        plt.ylabel('IES Retention Scores', fontsize=16)
        leg=plt.legend(ncol=3)#, loc='upper center')
        #plt.show()
        # if outPath:
        plt.savefig(f'{outPath.rstrip("/")}/EWMA_{i}bpWindow.svg', bbox_inches='tight')

def baseFreq(df, KDs, IESAnnoFile, peaks=[(26,31),(44,52),(54,62)], outPath='.'):
    '''
    Calculates and plots basefrequency for the first three bases that follow after the TA
    The code is a bit confusing and not really nicely readable. If you have a idea to tweak it feel free 
    '''
    DF = pd.DataFrame(index=[round(i,2) for i in np.arange(0,1.01,0.01)])
    DF['IRS'] = DF.index
    IESAnno = {i.split('ID=')[-1].split(';')[0] : i.split('sequence=')[-1].strip() for i in open(IESAnnoFile) if not i.startswith('#')}

    for spl in KDs:
        for MIN , MAX in peaks:
            m = np.mean(df[spl])
            sd = np.std(df[spl])
            m_min_sd = m-2*sd
            m_plus_sd = m+2*sd
            dct = { round(i,2):[] for i in np.arange(0,1.01,0.01) } ## dct w/ k==IRS
            for ind , i in enumerate(df[spl]):
                if (MIN <= df.loc[ind, 'Length'] <= MAX) and (not np.isnan(i)):
                    dct[round(i,2)].append(df.loc[ind, 'Length']) ## values are the actual IRSs
            DF[f'{spl} Length Means'] = [np.mean(v) for v in dct.values()]
            fig, ax = plt.subplots(4, figsize=(8,7), sharex=True)
            ax[3].scatter(DF[f'{spl} Length Means'].index, y=DF[f'{spl} Length Means'], marker='d', s=5)
            ax[3].axvline(x=m, color='grey', linestyle='dashed')
            ax[3].axvline(x=m_min_sd, color='grey', linestyle='dotted')
            ax[3].axvline(x=m_plus_sd, color='grey', linestyle='dotted')
            ax[3].set_ylim(MIN-5,MAX+5)
            ax[3].set_ylabel('IES Length', fontsize=16)
            for pos in range(3):
                dct = { round(i,2):[] for i in np.arange(0,1.01,0.01) }
                for ind, i in enumerate(df[spl]):
                    if (MIN <= df.loc[ind, 'Length'] <= MAX) and (not np.isnan(i)):
                        dct[round(i,2)].append(IESAnno[df.loc[ind, 'ID']][pos+2])
                dct_freq = { k : {'A':0, 'T':0, 'G':0, 'C':0} for k in dct.keys()}
                for k,v in dct.items():
                    for nuc in v:
                        if nuc in dct_freq[k]:
                            dct_freq[k][nuc]+=1
                dct_freq_norm = { irs : { nuc : amount/len(dct[irs]) for nuc, amount in dct_freq[irs].items() if dct[irs]} for irs in dct.keys() }
                a=pd.Series(dct_freq_norm).reset_index().set_axis(['IRS','Nucs'],axis=1)
                b=a['Nucs'].apply(pd.Series).set_axis(['A','T','G','C'],axis=1)
                b.index=[round(i,2) for i in np.arange(0,1.01,0.01)]
                dd = pd.concat([DF,b], axis=1)
                for i in ['A','T','G','C']:
                    dd[f'{i} 10EWMA']=dd[i].ewm(span=10).mean()
                if pos==0:
                    ax[pos].plot(dd.loc[m_min_sd:m_plus_sd,'IRS'], dd.loc[m_min_sd:m_plus_sd,'C 10EWMA'], color='blue')
                    ax[pos].plot(dd.loc[m_min_sd:m_plus_sd,'IRS'], dd.loc[m_min_sd:m_plus_sd,'T 10EWMA'], color='red')
                    ax[pos].scatter(dd.index, dd['C'], color='blue', s=5)
                    ax[pos].scatter(dd.index, dd['T'], color='red', s=5)
                    ax[pos].axvline(x=m, color='grey', linestyle='dashed')
                    ax[pos].axvline(x=m_min_sd, color='grey', linestyle='dotted')
                    ax[pos].axvline(x=m_plus_sd, color='grey', linestyle='dotted')
                    ax[pos].set_title(spl, fontsize=16)
                    ax[pos].set_ylim(0,1)
                    ax[pos].set_ylabel(f'% {pos+1}. base', fontsize=16)
                if pos==1:
                    ax[pos].plot(dd.loc[m_min_sd:m_plus_sd,'IRS'], dd.loc[m_min_sd:m_plus_sd,'A 10EWMA'], color='green')
                    ax[pos].plot(dd.loc[m_min_sd:m_plus_sd,'IRS'], dd.loc[m_min_sd:m_plus_sd,'T 10EWMA'], color='red')
                    ax[pos].scatter(dd.index, dd['A'], color='green', s=5)
                    ax[pos].scatter(dd.index, dd['T'], color='red', s=5)
                    ax[pos].axvline(x=m, color='grey', linestyle='dashed')
                    ax[pos].axvline(x=m_min_sd, color='grey', linestyle='dotted')
                    ax[pos].axvline(x=m_plus_sd, color='grey', linestyle='dotted')
                    ax[pos].set_ylim(0,1)
                    ax[pos].set_ylabel(f'% {pos+1}. base', fontsize=16)
                if pos==2:
                    ax[pos].plot(dd.loc[m_min_sd:m_plus_sd,'IRS'], dd.loc[m_min_sd:m_plus_sd,'G 10EWMA'], color='orange')
                    ax[pos].scatter(dd.index, dd['G'], color='orange', s=5)
                    ax[pos].axvline(x=m, color='grey', linestyle='dashed')
                    ax[pos].axvline(x=m_min_sd, color='grey', linestyle='dotted')
                    ax[pos].axvline(x=m_plus_sd, color='grey', linestyle='dotted')
                    ax[pos].set_ylim(0,1)
                    ax[pos].set_ylabel(f'% {pos+1}. base', fontsize=16)
            # fig.suptitle(splt, fontsize=16)
            # plt.title(splt, fontsize=16)
            plt.xlim(0,1)
            plt.xlabel('IES Retension Score', fontsize=16)
            outSpl = '_'.join(spl.replace('/','_').split())
            plt.savefig(f'{outPath.rstrip("/")}/{outSpl}_baseFreqPlot_peak_{MIN}_{MAX}.svg', bbox_inches='tight')
    return True

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
    #args = parser.parse_args()
    
    def getTime():
        return str(datetime.now()).split('.')[0].replace('-','.').replace(' ','|')

    checkArgs(args,parser)
    KDs = args.KDs
    IRS = args.IRS
    if args.IES_annotaion:
        IESAnnoFile = args.IES_annotaion
    else:
        p = '/'.join(os.path.abspath(__file__).split('/')[:-2])
        IESAnnoFile = f'{p}/res/anno/internal_eliminated_sequence_PGM_ParTIES.pt_51_with_ies.gff3'
    peaks = None
    outPath = None
    if args.outDir:
        outPath = args.outDir
        if os.path.exists(outPath):
            shutil.rmtree(outPath)
        os.mkdir(outPath)
    if args.peaks:
        peaks = args.peaks
        peaks = splitPeaks(peaks)

    print(f'\n\n[{getTime()}] - Starting Analysis\n')
    print(f'\n\n[{getTime()}] - Preparing Data\n')
    df = prepData(IRS, KDs) 
    print(f'\n\n[{getTime()}] - Starting Analysis\n')
    if args.switch == 'EWMA':
        print(f'\n\n[{getTime()}] - Calculating EWMA and generating plots\n')
        if outPath:
            generateEWMA(df, KDs, outPath)
        else: 
            generateEWMA(df, KDs, outPath='.')
        print(f'\n\n[{getTime()}] - Done!\n')
    
    if args.switch == 'FREQ':
        print(f'\n\n[{getTime()}] - Calculating and ploting base frequencies\n')
        if outPath and peaks:
            baseFreq(df, KDs, IESAnnoFile, peaks=peaks, outPath=outPath)
        elif outPath: 
            baseFreq(df, KDs, IESAnnoFile, outPath=outPath)
        elif peaks:
            baseFreq(df, KDs, IESAnnoFile, peaks=peaks)
        else:
            baseFreq(df, KDs, IESAnnoFile)
        print(f'\n\n[{getTime()}] - Done!\n')

    if args.switch == 'BOTH':
        print(f'\n\n[{getTime()}] - Calculating EWMA and generating plots\n')
        if outPath:
            generateEWMA(df, KDs, outPath)
        else: 
            generateEWMA(df, KDs, outPath='.')
        print(f'\n\n[{getTime()}] - Done!\n')

        print(f'\n\n[{getTime()}] - Calculating and ploting base frequencies\n')
        if outPath and peaks:
            baseFreq(df, KDs, IESAnnoFile, peaks=peaks, outPath=outPath)
        elif outPath: 
            baseFreq(df, KDs, IESAnnoFile, outPath=outPath)
        elif peaks:
            baseFreq(df, KDs, IESAnnoFile, peaks=peaks)
        else:
            baseFreq(df, KDs, IESAnnoFile)
        print(f'\n\n[{getTime()}] - Done!\n')


if __name__ == '__main__':
    main()

    
