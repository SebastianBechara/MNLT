#! /usr/bin/env python3

# input protein sequence for just a optimisation. If you want to recodonise input the nucleotide sequence

import sys
from Bio import SeqIO


codon_usage = {
                'M' : ['ATG'],
                'E' : ['GAA', 'GAG'],
                'S' : ['TCG', 'AGC',  'TCC', 'TCT', 'AGT', 'TCA'],  ## 'TCN',
                'N' : ['AAT', 'AAC'],
                'K' : ['AAA', 'AAG'],
                'D' : ['GAT', 'GAC'],
                'L' : ['CTT',  'TTA', 'TTG', 'CTA', 'CTC', 'CTG'],   ## 'CTN',
                'G' : ['GGT',  'GGA', 'GGC', 'GGG'],  ## 'GGN',
                'I' : ['ATC', 'ATT', 'ATA'],
                'T' : ['ACA', 'ACT', 'ACC', 'ACG'],
                'C' : ['TGT', 'TGC'],
                'F' : ['TTC', 'TTT'],
                'R' : ['CGG', 'CGA', 'CGC', 'AGA', 'CGT', 'AGG'],
                'Y' : ['TAT', 'TAC'],
                'W' : ['TGG'],
                'Q' : ['TAG', 'CAA', 'TAA', 'CAG'],
                'A' : [ 'GCC', 'GCT', 'GCG', 'GCA'],  ##  'GCN',
                'H' : ['CAC', 'CAT'],
                'V' : ['GTT', 'GTC',  'GTG', 'GTA'],   ##  'GTN',
                'P' : ['CCT', 'CCC', 'CCA', 'CCG'],
                '*' : ['TGA']
                }

codon_usage_high = {
                'M' : 'ATG',
                'E' : 'GAA', #, 'GAG'],
                'S' : 'TCA', #G', 'AGC',  'TCC', 'TCT', 'AGT', 'TCA'],  ## 'TCN',
                'N' : 'AAT', # 'AAC'],
                'K' : 'AAA', # 'AAG'],
                'D' : 'GAT', #, 'GAC'],
                'L' : 'TTA', #'CTT',  'TTA', 'TTG', 'CTA', 'CTC', 'CTG'],   ## 'CTN',
                'G' : 'GGA', #T',  'GGA', 'GGC', 'GGG'],  ## 'GGN',
                'I' : 'ATT', #C', 'ATT', 'ATA'],
                'T' : 'ACT', #A', 'ACT', 'ACC', 'ACG'],
                'C' : 'TGT', #, 'TGC'],
                'F' : 'TTT', #C', 'TTT'],
                'R' : 'AGA', #CGG', 'CGA', 'CGC', 'AGA', 'CGT', 'AGG'],
                'Y' : 'TAT', # 'TAC'],
                'W' : 'TGG',
                'Q' : 'TAA', #G', 'CAA', 'TAA', 'CAG'],
                'A' : 'GCT', #'GCC', 'GCT', 'GCG', 'GCA'],  ##  'GCN',
                'H' : 'CAT', #C', 'CAT'],
                'V' : 'GTT', # 'GTC',  'GTG', 'GTA'],   ##  'GTN',
                'P' : 'CCA', #T', 'CCC', 'CCA', 'CCG'],
                '*' : 'TGA'
                }


def find_new_cod(cod, codons):
    for k,v in codons.items():
        if cod.upper() in v:
            new_cod = v[rd.randint(0, len(v)-1)]
    return new_cod

def get_aminoacid(cod, codons):
    for k,v in codons.items():
        if cod in v:
            return k

def get_new_seq(seq, codons):
    new_seq=[]
    ## iteration=[]
    for i in range(0, len(seq)+1,3):
        print(i)
        cod = seq[i:i+3]
        new_cod = cod
        print(cod, new_cod)
        if cod != '' or new_cod != '' and len(cod) == 3:
            if len(list(set([i[0] for i in codons[get_aminoacid(cod, codons)]]))) != 1:
                while cod[0] == new_cod[0]:
                    new_cod = find_new_cod(cod, codons)#, codon_usage)
                else:
                    new_seq.append(new_cod)
                    ## iteration.append(new_seq)
            elif len(list(set([i[1] for i in codons[get_aminoacid(cod, codons)]]))) != 1:
                while cod[1] == new_cod[1]:
                    new_cod = find_new_cod(cod, codons)#, codon_usage)
                else:
                    new_seq.append(new_cod)
            elif len(list(set([i[2] for i in codons[get_aminoacid(cod, codons)]]))) != 1:
                while cod[2] == new_cod[2]:
                    new_cod = find_new_cod(cod, codons)#, codon_usage)
                else:
                    new_seq.append(new_cod)
            else:
                new_seq.append(cod)
    return ''.join(new_seq)

def check_GC_content(seq):
    GC_only = seq.replace('A','').replace('T','')
    GC_perc = len(GC_only) / len(seq)
    return GC_perc

def write_outseq(seq):
    with open('temp/new_fancy_sequence.fasta', 'w') as w:
        w.write('>new fancy codon optimised sequence\n'+str(seq)+'\n')

def run_sil_recod(seq, codons):
    ##
    # seqs = {i.description:str(i.seq) for i in SeqIO.parse(str(seq),'fasta')}
    ##
    og_gc = check_GC_content(seq)
    new_gc = og_gc+.6
    iteration = []
    # new_seqs=[]
       ############ Why The Fuck is new_seq referenced before assignment ????????
    ###
    # for k,v in seqs.items():
    ###
    if int(new_gc*100000)  in range(int((og_gc-.05)*100000) , int((og_gc+.05)*100000)):
        seq = new_seq
    else:
        while int(new_gc*100000)  not in range(int((og_gc-.05)*100000) , int((og_gc+.05)*100000)):
            new_seq = get_new_seq(seq, codons)
            new_gc= check_GC_content(new_seq)
            iteration.append(new_seq)
            print(new_seq)
            if len(iteration) == 100:
                new_seqs_gc = {i:int(check_GC_content(i)*100000) for i in iteration}
                min_gc = min(list(new_seqs_gc.values()))
                for k1,v in new_seqs_gc.items():
                    if min_gc == v:
                        key = k1
                # key = [k  for k,v in new_seqs_gc.items() if min_gc in v  ]
                print(type (key))
                new_seq = new_seqs_gc[  key  ]# [k  for k,v in new_seqs_gc.items() if min(list(new_seqs_gc.values())) in v  ][0]   ]
                # write_outseq(new_seq))
                print(new_seq)
                # new_seqs.append('>'+str(k)+'_recodonised\n'+str(key)+'\n')
                break
            # print(new_seq)
        else:
            # new_seqs.append('>'+str(k)+'_recodonised\n'+str(new_seq)+'\n')
            print('\nDONE !\n\n', new_seq)
        # return new_seqs

def just_recod(seq, codons):
    # seqs = {i.description:str(i.seq) for i in SeqIO.parse(str(seq),'fasta')}
    # new_seqs = []
    # for k,v in seqs.items():
    #     new_seq=[codons[i] for i in v]
    #     new_seqs.append('>'+str(k)+'_recodonised\n'+str(''.join(new_seq))+'\n')

    new_seq=[codons[i.upper()] for i in seq]
                                                                    # for i in seq:
                                                                    #     # print(i)
                                                                    #     new_seq.append(codons[i])
    # return new_seqs
    return ''.join(new_seq)

# def write_file(filename, seqs):
#     with open(str(filename).split('.fasta')[0]+'_recodonised.fasta') as w:
#         w.write('\n'.join(seqs))

def main():
    if len(sys.argv[1:]) != 2:
        print('too many or missing arguments ! \n\n\t>>>Usage\n\tcodon_opt.py <recodonise | optimise> <plain sequence in correct RF>')
        sys.exit()

    seq = sys.argv[2]
    flag = sys.argv[1]

    # seqs = {i.description:str(i.seq) for i in SeqIO.parse(str(file),'fasta')}
    # for k,seq in seqs.items():
    if flag == 'recodonise':
        run_sil_recod(seq, codon_usage)
        # new = run_sil_recod(seq, codon_usage)
        # write_file(seq, new)

    elif flag == 'optimise':
        print('\nDONE !\n\nHere is you new sequence:\n',just_recod(seq, codon_usage_high))
        
