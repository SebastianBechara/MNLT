#! /usr/bin/env python3
import os, sys, argparse
import regex as re
# from pychara.annotparse.annoTool import parseAnnoIES
from lib.pychara.annotparse.annoTool import parseAnnoIES

def argParse():
    des='This script will return a FastQ of reads spanning MDS|IES junctions. It takes a sam as input and only outputs reads that map perfectly (25M in CIGAR).'
    parser=argparse.ArgumentParser(description=des,
                            allow_abbrev=False,
                            add_help=False,
                            usage=argparse.SUPPRESS)

    req=parser.add_argument_group('Required')
    req.add_argument('-i','--inFile', action='store', help='Sam file to infer overlapping reads of')

    opt=parser.add_argument_group('Optional')
    opt.add_argument('-g','--GFF',action='store', help='GFF file specifying IESs. Default is for my machine')
    opt.add_argument('-g2','--GFF_2', action='store', help='GFF file specifying genes. Default is for my machine')
    opt.add_argument('-o','--outReads', action='store', help='Name of FastQ file of the boundary reads')
    opt.add_argument('-o2','--outTSV', action='store', help='Name of TSV file containing info on the overlapped IESs')
    opt.add_argument('-h','--help', action='help', help='Display this message and exit')
    args=parser.parse_args()
    return args, parser

def importIESanno(path='/home/bechara/Bioinformatics/Annotation_files/internal_eliminated_sequence_PGM_ParTIES.pt_51_with_ies.gff3',
                  path2='/home/bechara/Bioinformatics/Annotation_files/ptetraurelia_mac_51_annotation_v2.0.gff3'):
    print(' . . . generating IES annotation dictionary')
    IESs, IESdct = parseAnnoIES(path, path2, '/home/bechara/Bioinformatics/Annotation_files/IESRetentionScores/All_IRSs_20240123.tsv')
    return IESdct

def evalCigar(s):
    if len(re.findall(r'\d+\w',s))==1:
        return True
    else:
        return False

def masterBOUNDARYProgram(insam,IES,out=None,out2=None):
    checked=[]
    if out:
        out=out
    else:
        out=insam.split('.')[0]+'MDSIESBoundaryReads.fq'
    if out2:
        out2=out2
    else:
        out2=insam.split('.')[0]+'MDSIESBoundaryIESs.tsv'
    overlappedIESs={ies:0 for scaf,lst in IES.items() for ies in lst}
    with open(out,'w')as w:
        print(' . . . writing MDS|IES boundary Reads to file')
        for line in open(insam):
            scaf=line.split('\t')[2]
            if (not line.startswith('@')) and (scaf != '*') and (scaf in IES):
                start=int(line.split('\t')[3])
                #evalCig=evalCigar(line.split('\t')[5]
                if evalCigar(line.split('\t')[5]):#)==list:
                    stop=start+25
                else:
                    continue
                if line.split('\t')[0] in checked:
                    continue
                for ies in IES[scaf]:
                    if (start+5 < ies.start < stop-5) or (start+5<ies.stop<stop-5):
                        overlappedIESs[ies]+=1
                        id=line.split('\t')[0]
                        seq=line.split('\t')[9]
                        qual=line.split('\t')[10]
                        checked.append(id)
                        w.write('@'+id+'\n'+seq+'\n+\n'+qual+'\n')
    with open(out2, 'w') as w:
        for ies,cnt in overlappedIESs.items():
            if cnt!=0:
                w.write('\t'.join([str(ies.id),str(ies.scaf),str(ies.start),str(ies.stop),str(ies.length),str(ies.GC),str(cnt)])+'\n')

if __name__ == '__main__':
    args, parser=argParse()
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit()
    if not args.inFile:
        print('\n\nERROR! You need to specify an input SAM file! \n\n')
        print('#'*90+'\n\n')
        parser.print_help(sys.stderr)
        sys.exit()
    if args.GFF:
        if args.GFF_2:
            IES = importIESanno(args.GFF, args.GFF_2)
        else:
            print('\n\nERROR! You need to specify a GFF for the genes!\nDefault: /home/bechara/Bioinformatics/Annotation_files/ptetraurelia_mac_51_annotation_v2.0.gff3\n\n')
            print('#'*90+'\n\n')
            parser.print_help(sys.stderr)
            sys.exit()
    if args.GFF_2:
        if not args.GFF:
            print('\n\nERROR! You have a GFF for the genes but not for the IESs!\nDefault:/home/bechara/Bioinformatics/Annotation_files/internal_eliminated_sequence_PGM_ParTIES.pt_51_with_ies.gff3 \n\n')
            print('#'*90+'\n\n')
            parser.print_help(sys.stderr)
            sys.exit()
    if (not args.GFF) and (not args.GFF_2):
        IES = importIESanno()
    if args.outReads:
        if args.outTSV:
            masterBOUNDARYProgram(args.inFile,IES, args.outReads, args.outTSV)
        else:
            masterBOUNDARYProgram(args.inFile,IES, args.outReads)
    else:
        masterBOUNDARYProgram(args.inFile,IES)
