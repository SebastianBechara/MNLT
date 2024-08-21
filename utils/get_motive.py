import os ,sys
import regex as re
from Bio import SeqIO

def calc(trans_fasta):

    prot_fasta = {str(i.description):str(i.seq) for i in SeqIO.parse(trans_fasta, 'fasta')}
    #print(prot_fasta)
    #pattern = r'([DE]{1}YK[DE]{4}K){e<=1}'   ##   for flag related stuff
    pattern = r'[PVILMFYW]{1}K[RHKDESTNQCUGPAVILMFYW]{1}[DE]{1}'   ##   for sumoylation
    #seq_pos = []
    res = {str(i.description):[] for i in SeqIO.parse(trans_fasta, 'fasta')}
    #print(res)
    #cnt=0
    print('\n\n####\n\nSearching pattern . . .\n\n####\n\n\n')
    for k,v in prot_fasta.items():
        print('\t##\tSearching pattern for: '+str(k), end='\t\t\r')
        #try:
        #    re.compile(pattern).search(v).string   ### key : [(no_match|containing_string)]
        #except:
        #    res[k].append('no_match')   ### key : [(no_match|containing_string)]
        #else:
        #    res[k] += [re.compile(pattern).search(v).string]
        iterator = re.compile(pattern).finditer(v)
        if all(False for _ in iterator) == True:
            res[k] += ['no_match', r'None', r'None', r'None', r'None']
        else:
            cnt=1
            for i in re.compile(pattern).finditer(v):
                if len(res[k]) >= 5:
                    #res[k+'_'+str(cnt)] += res[k] + [i.start(), i.end(), i.group()] #v[i.start():i.end()]]
                    key = k+'_'+str(cnt)
                    res[key] = [str(v), i.start(), i.end(), i.group()] #v[i.start():i.end()]]
                    cnt+=1
                    try:
                        res[key] += [i.fuzzy_counts]   ### key : [(no_match|containing_string), start_pos, end_pos, hit_substring, fuzzy_count(substitutions, insertions, deletions)]
                   # .join([j for j in i.fuzzy_counts])
                    except:
                        res[key] += [r'None']
                else:
                    res[k] += [str(v), i.start(), i.end(), i.group()] #v[i.start():i.end()]]   ### key : [(no_match|containing_string), start_pos, end_pos, hit_substring]
                    #seq_pos.append((i.start(), i.end()))
                    try:
                        res[k] += [i.fuzzy_counts]   ### key : [(no_match|containing_string), start_pos, end_pos, hit_substring, fuzzy_count(substitutions, insertions, deletions)]
                    except:
                        res[k] += [r'None']
    print('\t##\tSearching pattern for: '+str(k)+'\n')
    return prot_fasta, res

def calc_mass(prot_fasta, res):

    AA_sizes = {    ## electircally charged AA (+ or neu)
                'R':0.1742,
                'H':0.1552,
                'K':0.1462,
                    ## electrically charged AA (-)
                'D':0.1331,
                'E':0.1471,
                    ## polar 'charged' AA, hydrophilic
                'S':0.1051,
                'T':0.1191,
                'N':0.1321,
                'Q':0.1462,
                    ## special AA
                'C':0.1212,
                'U':0.168,
                'G':0.0751,
                'P':0.1151,
                    ## apolar AA, hydrophobic
                'A':0.0891,
                'V':0.1171,
                'I':0.1312,
                'L':0.1312,
                'M':0.1492,
                'F':0.1652,
                'Y':0.1812,
                'W':0.2042,
                #'X':sum(list(AA_sizes.values()))/len(list(AA_sizes.keys()))
                }

    print('\n\n####\n\nCalculating Mass . . .\n\n####\n\n')
    for k,v in res.items():
        print('\t##\tCalculating mass for: '+str(k), end='\t\t\r')
        mass = 0
        #length = 0
        if res[k][0] != 'no_match':
            if '_' in k:
                for aa in prot_fasta[k.split('_')[0]].strip('*'):
                    if aa == 'X':
                        mass += sum(list(AA_sizes.values()))/len(list(AA_sizes.keys()))
                    else:
                        mass += AA_sizes[aa]
            else:
                for aa in prot_fasta[k].strip('*'):
                    if aa == 'X':
                        mass += sum(list(AA_sizes.values()))/len(list(AA_sizes.keys()))
                    else:
                        mass += AA_sizes[aa]
                    # res[k] += [len(v), mass]   ### key : [(no_match|containing_string), start_pos, end_pos, hit_substring, length_of_seq, kDa]
        else:
            res[k] += [r'None', r'None']
            continue
        if '_' in k:
            res[k] += [len(prot_fasta[k.split('_')[0]]), mass]
        else:
            res[k] += [len(prot_fasta[k]), mass]
    print('\t##\tCalculating mass for: '+str(k)+'\n')
    return res

#def h_dist(str, str2):
    ##d
#    for bla in bla:
#        print()

def write_to_file(results):
    output = 'Flag_search_RM.tsv'
    with open(str(output),'w') as outfile:
        outfile.write('Assecion Number\tThe whole sequence or a no match flag\tstart position of the matching sequence\tend position\
         of matching sequence\tthe actual matching string\t1.Substitutions 2.Insertions 3.Deletions\tlength of input sequence\tmass in kDa of input\n')
        for k,v in results.items():
            if v[0] != 'no_match':
                outfile.write(str(k)+'\t'+str(v[0])+'\t'+str(v[1])+'\t'+str(v[2])+'\t'+str(v[3])+'\t'+str(v[4])+'\t'+str(v[5])+'\t'+str(v[6])+'\n')
def main():

    file = sys.argv[1]

    prot, diction = calc(file)
    result = calc_mass(prot, diction)
    write_to_file(result)
    #cnt=0
    #for k,v in diction.items():
        #if diction[k][0] != 'no_match':
            #cnt+=1
            #print('number of containg seqs: '+str(cnt))
            #print(diction[k])
            #print(diction[k][4])
    #print(diction)
main()
