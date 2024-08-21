#! /usr/bin/env python3

import sys, os, argparse
from Bio import SeqIO
from BCBio import GFF

def argParse():
    des='Generate a genbank file out of a genomic fasta and the respective gff3 file'

    parser = argparse.ArgumentParser(description=des,
                                     add_help = False,
                                     allow_abbrev=False,
                                     usage=argparse.SUPPRESS)

    req = parser.add_argument_group('Required')
    req.add_argument('-f','--fasta',action='store', help='FastA file containing the genome sequence')
    req.add_argument('-g', '--gff', action='store', help='GFF3 file containing the coresponding annotations to the genome file')

    opt = parser.add_argument_group('Optional')
    opt.add_argument('-t', '--moleculeType', action='store', help='Molecule type of gff. If genomic then DNA (also the default).')

    args = parser.parse_args()

    return args, parser

def checkMoleculeType(iterator, molType='DNA'):
    """
    Set the molecule type to DNA.
    """
    for rec in iterator:
        if 'molecule_type' not in rec.annotations:
            rec.annotations['molecule_type']=molType
        yield rec

def convertSubfeatur2Feature(iterator):
    """
    Converts subfeatures (CDS,UTRs,etc.) of genes to regular features.
    Seems that the genbank file will only contain gene annotations otherwise.
                       !!!!!NOT DONE YET!!!!!
    """
    for rec in iterator:
        newRec=[]
        for feat in rec.features:
            newRec.append(rec)
            while len(newRec) > 0:



def main(args):
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit()
    if not args.fasta and not args.gff:
        print('ERROR!! Missing input!')
        parser.print_help(sys.stderr)
        sys.exit()
    if not args.fasta.endswith('fa'):
        if not args.fasta.endswith('fasta'):
            print('ERROR!! Wrong Fromat!')
            parser.print_help(sys.stderr)
            sys.exit()

    inFasta = SeqIO.to_dict(SeqIO.parse(args.fasta,'fasta'))

    gff_iter = GFF.parse(args.gff, inFasta)

    out = args.gff.split('.fa')[0]+'Scrp.genbank'
    SeqIO.write(checkMoleculeType(gff_iter), out, 'genbank')


if __name__ == '__main__':
    args, parser = argParse()
    main(args)
