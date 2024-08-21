#!/usr/bin/env python3

# Author: Xyrus Maurer-Alcala
# Contact: maurerax@gmail.com or xyrus.maurer-alcala@izb.unibe.ch
# Modifications: Sebastian Bechara
# Contact: ask Mariusz
# Last Modified: 2024.07
# usage: python nlab.IRS.py -fwd FwdReads -rev RevReads -k KD

'''A dirty pythonic re-write of the MIRET module of the parTIES program
(Wilkes, Arnaiz, and Sperling 2016 Bioinformatics).

Author: Xyrus Maurer-Alcala
Modifications: Sebastian Bechara 

Note that the intent was to solely run the MIRET module, but with a dramatic
reduction in computational time (e.g. ~16 hours to ~1 hour). Best practice
is to use a random sample of ~10-15 million reads (most commonly used amount).
This can easily be done with seqtk (recommended).

Will perform the read mapping and IES retention score calulations for you.


Dependencies:
    Hisat2, samtools'''

import argparse
import subprocess
import sys, os

from argparse import RawTextHelpFormatter,SUPPRESS

def check_args():
	parser = argparse.ArgumentParser(description =
        '\nRuns "streamlined" version of ParTIES (Wilkes, Arnaiz, and Sperling ' \
        '2016 Bioinformatics).\nAdjust the PATH to your germline/somatic ' \
        'Hisat2 databases (lines 66/69) and gff3 of IES positions (line 173) ' \
        'as needed.', usage=SUPPRESS, formatter_class=RawTextHelpFormatter)

	required_arg_group = parser.add_argument_group("Required Arguments")
	required_arg_group.add_argument("--kd", "-kd", action = "store",
		help = " Knockdown name.\n")

	optional_arg_group = parser.add_argument_group("Optional Arguments")
	optional_arg_group.add_argument("--soma_sam", "-s", action = "store",
        help = " SamFile for the reads mapped against the somatic genome.\n\n")
	optional_arg_group.add_argument("--germ_sam", "-g", action = "store",
        help = " SamFile for the reads mapped against the germline genome.\n\n")
	optional_arg_group.add_argument("--fwd", "-fwd", action = "store",
        help = " 'Forward' paired end reads.\n\n")
	optional_arg_group.add_argument("--rev", "-rev", action = "store",
        help = " 'Reverse' paired end reads.\n\n")
	optional_arg_group.add_argument("--single", "-single", action = "store",
        help = " Single end reads.\n\n")
	optional_arg_group.add_argument("--na", "-na", action = "store_true",
        help = " Record IESs with low support as 'NA' rather than '0'.\n\n")
	optional_arg_group.add_argument('-t', '--threads', action='store', 
        help = 'Specify number of threads')

	return parser 


# Map the reads with "default" parameters: no spliced alignments (use MILORD
# from parties if interested).
def map_reads(reads, KD, soma = True, threads = 8):
    p = '/'.join(os.path.abspath(__file__).split('/')[:-2])
    if soma:
        db = f"{p}/res/ref/hisat2/para_MAC/mac"
        out_sam = f'{KD}.SomaMap.sam'
    else:
        db = f"{p}/res/ref/hisat2/para_MIC/ies"
        out_sam = f'{KD}.GermMap.sam'
    if len(reads) == 2:
        hisat_cmd = f'hisat2 -x {db} --no-spliced-alignment -1 {reads[0]} ' \
                    f'-2 {reads[1]} -p {threads} | samtools view -SF 4 -h - '\
                    f' > {out_sam}'
    elif len(reads) == 1:
        hisat_cmd = f'hisat2 -x {db} --no-spliced-alignment -U {reads[0]} ' \
                    f' -p {threads} | samtools view -SF 4 -h - > {out_sam}'
    else:
        print('Somthing went wrong...')
        sys.exit()

    hisat_call = subprocess.call(hisat_cmd, shell=True)#,
        #stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

    return out_sam


# Coordinates in the somatic and germline genomes of the IES boundaries.
def parseGFF(gff):
    iesBySCF = {}
    iesDict = {}

    for i in open(gff).readlines():
        if i.startswith("scaffold51_"):
            soma = (int(i.split(";")[0].split(".")[-1]),)
            germ = tuple(map(int,i.split("\t")[3:5]),)+(i.split(";")[0].split("=")[1],)
            iesBySCF.setdefault(i.split("_with_IES")[0],[]).append(soma+germ)
            iesDict[germ[-1]] = [0,0]

    return iesBySCF, iesDict


# Parse through sam-formatted files and evaluate coordinates that span
# MDS-IES boundaries as either IES- or IES+ (follows MIRET module).
def parseSam(samfile, iesDict, iesBySCF, soma=True):
    linesProcessed = 0

    for i in open(samfile):
        if not i.startswith("@"):
            scf = i.split("\t")[2].split("_with")[0]
            print(f'\rLooked at {linesProcessed} lines', end='')
            if scf in iesBySCF.keys():
                linesProcessed += 1

                if len(i.split("\t")[5]) == 4:
                    fpos = int(i.split("\t")[3])
                    tpos = fpos+len(i.split("\t")[9])

                    for j in iesBySCF[scf]:
                        if soma and tpos > j[0]:
                            if fpos <= j[0]-15 and tpos >= j[0]+15:
                                iesDict[j[-1]][0] += 1
                        else:
                            if tpos < j[1]:
                                break
                            else:
                                boundary = 0
                                if fpos <= j[1]-15 and tpos >= j[1]+15:
                                    boundary += 1
                                if fpos <= j[2]-15 and tpos >= j[2]+15:
                                    boundary += 1
                                if boundary:
                                    iesDict[j[-1]][1] += 1

    print(f'\rLooked at {linesProcessed} lines')   #### iesDict = { ID : [soma count, germline count]}


# Calculates the IES Retention Scores (IRS) for each IES.
def calcIRS(iesDict, KD, nan = False):
    irs_scores = {}
    poor_ies = []
    for k, v in iesDict.items():
        if sum(v) < 5:
            poor_ies.append(k)
            # If "nan" is set to true, then IESs with poor support/information
            # are recorded as "NA", otherwise, recoreded as "0".
            if nan:
                irs_scores[k] = "NA"
            else:
                irs_scores[k] = 0
        else:
            irs_scores[k] = v[1]/sum(v)
    if not nan:
        print(f'There were {len(poor_ies)} IESs with too few reads (<5) to evalute '\
            f'its retention score.\nTheir IRSs have been recorded as "0".\n')
    else:
        print(f'There were {len(poor_ies)} IESs with too few reads (<5) to evalute '\
            f'its retention score.\nTheir IRSs have been recorded as "NA".\n')

    with open(f'{KD}.IRS.tsv','w+') as w:
        w.write(f'IES\t{KD}\n')
        for k, v in irs_scores.items():
            if v != 'NA':
                w.write(f'{k}\t{format(float(v),".4f")}\n')
            else:
                w.write(f'{k}\t{v}\n')

    with open(f'{KD}_cov.tsv' , 'w') as w:
        w.write(f'IES ID\tIESMin\tIESPlus\n')
        for k,v in iesDict.items():
            w.write(f'{k}\t{v[0]}\t{v[1]}\n')

def main(args=None):
    parser = check_args()
    num=1
    if args:
        num=2
    if not args:
        args = parser.parse_args()
    if len(sys.argv) == num:
        parser.print_help(sys.stderr)
        sys.exit()

    if not args.soma_sam and not args.germ_sam:
        if args.threads:
            T = args.threads
        else: T = 8
        if args.fwd and args.rev:
            reads = [args.fwd, args.rev]
        elif args.single:
            reads = [args.single]
        else:
            print("[Error] Invalid combination of reads files. Accepts either " \
                "paired-end or single-end reads.")
            sys.exit(1)
        print('\n#------------------------------------------#\n# Mapping reads against SOMA (MAC) \
#\n#------------------------------------------#\n')
        args.soma_sam = map_reads(reads, args.kd, True, T)
        print('\n#------------------------------------------#\n# Mapping reads against GERMLINE (MIC) \
#\n#------------------------------------------#\n')
        args.germ_sam = map_reads(reads, args.kd, False, T)
    p = '/'.join(os.path.abspath(__file__).split('/')[:-2])
    ies_gff = f"{p}/res/anno/internal_eliminated_sequence_PGM_ParTIES.pt_51_with_ies.gff3"
    print('\n#------------------------------------------#\n# Reading IES annotations \
#\n#------------------------------------------#\n')
    iesBySCF, iesDict = parseGFF(ies_gff)

    print("\n#------------------------------------------#\n# Parsing Soma-Mapped Reads \
#\n#------------------------------------------#\n")
    parseSam(args.soma_sam, iesDict, iesBySCF)

    print("\n#------------------------------------------#\n# Parsing Germ-Mapped Reads \
#\n#------------------------------------------#\n")
    parseSam(args.germ_sam, iesDict, iesBySCF, soma=False)

    calcIRS(iesDict, args.kd, args.na)

if __name__ == "__main__":
    main()
