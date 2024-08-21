#! /usr/bin/env python3

import argparse, os, sys
from Bio import SeqIO
from Bio import Seq

def argParse():
    des = 'This script outputs ORF from an input fasta containing transcripts\n#######\n\
finished:\n-return transcripts with lo to no orfs'
    parser = argparse.ArgumentParser(description=des,
                                      allow_abbrev=False,
                                      add_help=False,
                                      formatter_class=argparse.RawTextHelpFormatter,
                                      usage=argparse.SUPPRESS)
    req = parser.add_argument_group('Required')
    req.add_argument('-i','--input', action='store', help='Input FastA file')

    opt = parser.add_argument_group('Optional')
    opt.add_argument('-n','--ncRNA', action ='store_true', help='Additionally outputs transcripts with no or only small ORFs (<100)')
    opt.add_argument('-os', '--ORFsize', action='store', help='specifies the size of the ORFs to filter. Default 100 aa')
    opt.add_argument('-o', '--output', action='store', help='Filebase to be used for output. (default=inputfilebase+ORF.fa)')
    opt.add_argument('-v','--verbose', action='store_true', help='Prints a neat overview and outputs a tiny log file [fastaBase.log]')
    opt.add_argument('-h','--help', action='help', help='Display this message and exit')

    args = parser.parse_args()
    return args, parser

def parseFasta(FILE):
    return SeqIO.parse(FILE, 'fasta')

def get6Frames(rec):
    # seqs=[]
    seqsfwd={i:[] for i in range(1,4)}
    seqsrev={i:[] for i in range(1,4)}
    for i in range(3):
        seqsfwd[i+1].append(Seq.translate(rec.seq[i:], table=6, to_stop=False, stop_symbol='*'))
        seqsrev[i+1].append(Seq.translate(rec.seq.reverse_complement()[i:], table=6, to_stop=False, stop_symbol='*'))
        # seqs.append(Seq.translate(rec.seq[i:], table=6, to_stop=False, stop_symbol='*'))
        # seqs.append(Seq.translate(rec.seq.reverse_complement()[i:],table=6, to_stop=False, stop_symbol='*'))
    return seqsfwd, seqsrev

# def getMpos(seq):
#     for i in seq:
#         if i.upper() == 'M':
#             return seq.index(i)

# for rec in infasta:
#     seqs = get6Frames(rec)
#     for i in seqs:
#         ORF=''
#         for aa in i:
#             if aa == '*':
#                 break
#             if aa.upper() == 'M':
#                 ORF+=aa
#
#             if ORF:
#                 ORF+=aa


def get_ORFs(seqsfwd):#, seqsrev):
    FWD_ORF={i:[] for i in range(1,4)}
    # REV_ORF=[i:[] for i in range(1,4]}
    for frame,seqs in seqsfwd.items():
        for seq in seqs:
            for subseq in seq.split('*'):
                if subseq.find('M') < 0:
                    continue
                start=subseq.find('M')
                stop=len(subseq[subseq.find('M'):])+start
                FWD_ORF[frame].append((start, stop, subseq[subseq.find('M'):]))# len(subseq[subseq.find('M'):]), subseq[subseq.find('M'):]))      #########  add start and stop in sequence
    # for frame,seqs in seqsrev.items():
    #     for seq in seqs:
    #         for subseq in seq.split('*'):
    #             if subseq.find('M') < 0:
    #                 continue
    #             FWD_REV[frame].append((subseq.find('M'), len(subseq[subseq.find('M'):]), subseq[subseq.find('M'):]))

    # for cnt,i in enumerate(seqs):
    #     # ORF_over200=[]
    #     # ORF_under200=[]
    #     for subseq in i.split('*'):
    #         # print(subseq.find('M'),type(subseq.find('M')))
    #         if subseq.find('M') <0:
    #             continue
    #         if cnt <3:
    #             FWD_ORF.append(subseq[subseq.find('M'):])
    #         else:
    #             REV_ORF.append(subseq[subseq.find('M'):])
            # if subseq[subseq.find('M'):] >= 200:
            #     ORF_over200.append(subseq[subseq.find('M'):])
            # elif subseq[subseq.find('M'):] < 200:
            #     ORF_under200.append(subseq[subseq.find('M'):])
            # else:
            #     ORF.append(subseq[subseq.find('M'):])
    return FWD_ORF#, REV_ORF

def masterORFprogram(inFile):
    dct={}
    for rec in inFile:#parseFasta(inPath):
        fwd, rev = get6Frames(rec)
        ORFfwd = get_ORFs(fwd)
        ORFrev = get_ORFs(rev)
        dct[rec.id]=[ORFfwd,ORFrev]   ###  { ID : [ FWDORF{frame:[(start,stop,aaSeq)]} , REVORF{frame:[(start,stop,aaSeq)]} ]}
        #dct[rec.id]= [i for i in get_ORFs(get6Frames(rec))]
    return dct

def writeOutFile(args, res, out=None, orfsize=100):
    infasta={i.id:str(i.seq) for i in SeqIO.parse(args.input,'fasta')}
    if out:
        outfile=out
    else:
        outfile=args.input.split('.fa')[0]
    if args.ncRNA:
        with open(outfile+'_NUCmin100ORFs.fa','w') as w:
            for id,dcts in res.items():
                cnt=0
                idCnt=len( [i for frame, lst in dcts[0].items() for i in lst] ) + len( [i for frame, lst in dcts[0].items() for i in lst] )
                for dct in dcts:
                    for frame , lst in dct.items():
                        for tpl in lst:
                            if cnt == idCnt:
                                w.write('>'+id+'\n'+infasta[id]+'\n') ##### link nuc size to aa size
                            if len(tpl[2])*3 <int(orfsize):#100:
                                cnt+=1


    with open(outfile+'_6RF_ORFs.fa','w') as w:
        for id,dcts in res.items():
            for dct in dcts:
                for frame, lst in dct.items():
                    for tpl in lst:
                        fid= id+ '_RF'+str(frame)+'_'+str(tpl[0])+'..'+str(tpl[1])
                        w.write('>'+fid+'\n'+str(tpl[2])+'\n')
                # w.write('>'+str(k)+'_RF'+'\n'+str(v)+'\n')

if __name__ == '__main__':

    args, parser = argParse()
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit()
    if not args.input:
        print('ERROR! No input file specified')
        parser.print_help(sys.stderr)
        sys.exit()
    if not args.input.endswith('.fa') and not args.input.endswith('.fasta'):
        print('ERROR! Wrong input format. Provide a FastA file')
        parser.print_help(sys.stderr)
        sys.exit()

    infasta = parseFasta(args.input)
    res=masterORFprogram(infasta)

    if args.output:
        if args.ORFsize:
            writeOutFile(args, res, args.output, args.ORFsize)
        else:
            writeOutFile(args, res, args.output)
    else:
        writeOutFile(args, res)
    # infasta = parseFasta(args.input)
    # get6Frames(infasta)
