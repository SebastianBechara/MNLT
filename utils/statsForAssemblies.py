#!/usr/bin/env python3
from Bio import SeqIO
from Bio.Seq import Seq
import regex as re
import os,sys, argparse

def argParse():
    descr = '\n****\n\nScript for analysing several assembly stats like, N50, Telomeres etc.\nUsage: Name.py [Option] -A AssemblyFile -D DirectoryWithAssemblies -g GenomeEstimate -T TelomereSequence\n\n****'
    parser = argparse.ArgumentParser(prog='STATSforASSEMBLIES',description=descr, allow_abbrev=False)
    parser.add_argument('-A','--assembly', metavar='', action='store', default='', type=str,help='The assembly file to analyse; FastA format')
    parser.add_argument('-D','--directory', metavar='', action='store', default='', type=str,help='A directory containing one or multiple assembly file(s) to analyse; FastA format')
    parser.add_argument('-g', '--genome_estimate', metavar='', action='store',type=int, help='Specify for NG50 calculation')
    parser.add_argument('-B', '--onlyProtRec', metavar='', action='store_true', default=False, help='Specify this flag if you only want the protein recovery analysis')
    parser.add_argument('-b', '--BLAST_', metavar='',action='store_true', default=False, help='Add this flag if protein recovery rate using BLAST is wanted - off by default')
    parser.add_argument('-f', '--fasta_for_BLAST', metavar='', action='store',type=str, help='Requires -b, The FastA file containing protein sequences to use for protein recovery')
    parser.add_argument('-e', '--E_val', metavar='', action='store', type=float, default=0.0001, help='Specify E-value cut-off for BLAST recovery [0.0001]')
    parser.add_argument('-o', '--outpuFile', metavar='', action='store_true', default=False, help="Location and name of output file; directories that don't exist will be generated")
    parser.add_argument('-T','--telomere', metavar='', default=['CCCAAA','CCCCAA'], help='Specify the telemoric sequence to search for, you may specify multiple ones; default = CCCAAA & CCCCAA', nargs='+')

    if len(sys.argv[1:]) == 0:
        print('\n*****\n\n')
        print(parser.description)
        sys.exit()

    args=parser.parse_args()

    return args



def prepFolders(flag=False):
    if os.path.isdir('temp')!=True:
        os.mkdir('temp')
    if flag:
        os.mkdir('temp/BLASTDB')

def getFolders(inpath):
    files = os.path.listdir()

def loadGenome(infasta):
    return {i.description:str(i.seq) for i in SeqIO.parse(str(infasta),'fasta')}

def getLens(dct):
    return sorted([len(i) for i in dct.values()], reverse=True)

def calc_N50(lens):
    print('calculating N50 . . .')
    total = sum(lens)
    seq=0
    for i in lens:
        seq+=i
        if (seq/total)*1000000000 >= 0.5*1000000000:
            return i

def calc_NG50(lens, genome_estimate):
    print('calculating NG50 . . .')
    total = int(genome_estimate)
    seq=0
    for i in lens:
        seq+=i
        if (seq/total)*1000000000 >= 0.5*1000000000:
            return i

def max_minScaff(lens):
    return max(lens), min(lens)

def total_bp(lens):#scaffs):
    # tot = 0
    # for v in list(scaffs.values()):
    #     tot += len(v)
    # return tot
    tot = sum(lens)
    print('total size of assembly is:',tot/1000000,'mbp')
    return tot


def checkGaps(assembly):
    Ns={k:[] for k in assembly.keys()}
    for k,v in assembly.items():
        for i in re.compile(r'N{5,}').finditer(v):
            Ns[k]+=[[i.start(), i.end()]]
    N={k:v for k,v in Ns.items() if len(v) >= 1}
    gaps=[]
    for k,v in N.items():
        for i in v:
            gaps.append(i[1]-i[0])
    return N, gaps

def telomerePrep():#teloSeq):
    # telos=[]
    # telosRev=[]
    # for seq in teloSep:
    #     seqRev=str(Seq(seq).reverse_complement())
    #     if seq not in telos or seqRev not in telos:
    #         telos.append(seq)
    #         telosRev.append(seqRev)
    C4A2='CCCCAA'
    C3A3='CCCAAA'
    T2G4=str(Seq(C4A2).reverse_complement())
    T3G3=str(Seq(C3A3).reverse_complement())
    return  C4A2, C3A3, T2G4, T3G3 #telos, telosRev

def checkTelomeres(assembly):
    # telos, telosRev = telomerePrep(teloSeq)
    # telomeresFWD=[]
    # telomeresREV=[]
    #
    # for telo in telos:
    #     for k,v in assembly.items():
    #         if len(list(re.compile(r'('+str(telo)+'){2,}').finditer(v[0:251]))) > 0:
    #             telomeresFWD.append(k)
    # for telo in telosRev:
    #     for k,v in assembly.items():
    #         if len(list(re.compile(r'('+str(telo)+'){2,}').finditer(v[0:251]))) > 0:
    #             telomeresREV.append(k)
    # TelFwd=set(telomeresFWD)
    # TelRev=set(telomeresREV)
    # both=list(TelFwd & TelRev)
    # single=list(TelFwd ^ TelRev)
    # return both, single
    c4a2Telos=[]
    c3a3Telos=[]
    t2g4Telos=[]
    t3g3Telos=[]
    C4A2, C3A3, T2G4, T3G3 = telomerePrep()
    for k,v in assembly.items():
        if len(list(re.compile(r'('+C4A2+'){2,}').finditer(v[0:251]))) > 0:
            c4a2Telos.append(k)
        if len(list(re.compile(r'('+C3A3+'){2,}').finditer(v[0:251]))) > 0:
            c3a3Telos.append(k)
        if len(list(re.compile(r'('+T2G4+'){2,}').finditer(v[-251:-1]))) > 0:
            t2g4Telos.append(k)
        if len(list(re.compile(r'('+T3G3+'){2,}').finditer(v[-251:-1]))) > 0:
            t3g3Telos.append(k)
    c4a2Telos=list(set(c4a2Telos))
    c3a3Telos=list(set(c3a3Telos))
    t2g4Telos=list(set(t2g4Telos))
    t3g3Telos=list(set(t3g3Telos))
    fwdBoth=set(c4a2Telos+c3a3Telos)
    revBoth=set(t2g4Telos+t3g3Telos)
    bothTelos=list(fwdBoth & revBoth)
    singleTelos=list(fwdBoth ^ revBoth)
    return bothTelos, singleTelos

# def checkTelomeres(assembly):
#     c4a2Telos=[]
#     c3a3Telos=[]
#     t2g4Telos=[]
#     t3g3Telos=[]
#     C4A2, C3A3, T2G4, T3G3 = telomerePrep()
#     for k,v in assembly.items():
#         if len(list(re.compile(r'('+C4A2+'){2,}').finditer(v[251:-251]))) > 0:
#             c4a2Telos.append(k)
#         if len(list(re.compile(r'('+C3A3+'){2,}').finditer(v[251:-251]))) > 0:
#             c3a3Telos.append(k)
#         if len(list(re.compile(r'('+T2G4+'){2,}').finditer(v[251:-251]))) > 0:
#             t2g4Telos.append(k)
#         if len(list(re.compile(r'('+T3G3+'){2,}').finditer(v[251:-251]))) > 0:
#             t3g3Telos.append(k)
#     c4a2Telos=list(set(c4a2Telos))
#     print(len(c4a2Telos))
#     c3a3Telos=list(set(c3a3Telos))
#     print(len(c3a3Telos))
#     t2g4Telos=list(set(t2g4Telos))
#     print(len(t2g4Telos))
#     t3g3Telos=list(set(t3g3Telos))
#     print(len(t3g3Telos))
#     fwdBoth=set(c4a2Telos+c3a3Telos)
#     print(len(fwdBoth))
#     revBoth=set(t2g4Telos+t3g3Telos)
#     print(len(revBoth))
#     bothTelos=list(fwdBoth & revBoth)
#     print(len(bothTelos))
#     singleTelos=list(fwdBoth ^ revBoth)
#     print(len(singleTelos))
#     return bothTelos, singleTelos
def makeBlastDB(assembly):
    out = 'temp/BLASTDB/'+str(assembly).split('/')[-1]+'BLASTDB'
    os.system('makeblastdb -in '+str(assembly)+' -dbtype nucl -parse_seqids -out '+str(out))
    return out

def tblastn(inprot, assembly_db, evalue):
    out = 'temp/'+str(inprot).split('/')[-1]+'_v_'+str(assembly_db).split('/')[-1]
    os.system('tblastn -query '+str(inprot)+' -db '+str(assembly_db)+' -evalue '+str(evalue)+' -max_target_seqs 1 -outfmt 6 -out '+out)
    return out

def filterBlast(blrep, pattern): ## NEED TO WRITE NEW!!!
    lt=[]
    for i in blrep: #open(str(blrep)).readlines():
        if i.startswith(str(pattern)):
            if i not in lt:
                lt.append(i)
    lst=list(set([i.split('\t')[0] for i in lt]))
    return lst

def getMappedPerc(blastRep,ref_fasta): ## NEED TO WRITE NEW!!!
    CiliateProtDB = {i.description:str(i.seq) for i in SeqIO.parse(str(ref_fasta),'fasta')}
    cilProtPTET ={k:v for k,v in CiliateProtDB.items() if 'PTET' in k}
    cilProtPCAU={k:v for k,v in CiliateProtDB.items() if 'PCAU' in k}
    cilProtICH={k:v for k,v in CiliateProtDB.items() if 'IMG' in k}
    cilProtTTH={k:v for k,v in CiliateProtDB.items() if 'TTHERM' in k}
    blastRp = open(str(blastRep)).readlines()
    ptet=filterBlast(blastRp, 'PTET')
    pcau=filterBlast(blastRp, 'PCAU')
    ttherm=filterBlast(blastRp, 'TTHER')
    ich=filterBlast(blastRp, 'IMG')
    perc_ptet=len(ptet)/len(cilProtPTET)
    perc_pcau=len(pcau)/len(cilProtPCAU)
    perc_tthm=len(ttherm)/len(cilProtTTH)
    perc_ich=len(ich)/len(cilProtICH)
    return perc_ptet, perc_pcau, perc_tthm, perc_ich

def runDefaults(assembly, genome_estimate):
    gen = loadGenome(assembly)
    lens = getLens(gen)
    N50 = calc_N50(lens)
    NG50 = calc_NG50(lens, genome_estimate)
    Tot = total_bp(lens)
    Max, Min = max_minScaff(lens)
    N,gaps = checkGaps(gen)
    both, single = checkTelomeres(gen)
    return N50, NG50, Tot, Max, Min, N, Gaps, both, single

def runBLASTProtRec(inprot, assembly, e_value):
    path = makeBlastDB(assembly)
    BlastRep = tblastn(inprot, path, e_value)
    getMappedPerc(BlastRep)#, pattern)



def writeOutRepDef(N50, NG50, Tot, Max, Min, N, Gaps, both, single, perc):
    with open('BLA')

def writeOutRepProtRec(perc):
    with open('BLA')

if __name__ == "__main__":
    args = argParse()
    if args.directory != '':
        res={}
        cnt=0
        dirCont = os.listdir(args.directory)
        for i in dirCont:
            if cnt == len(dirCont):
                print('ERROR: No FastA found!!\n\n make sure the files end with .fa or .fasta\n\n')
                print(parser.description)
                sys.exit()
            if i.endswith('.fa') or i.endswith('.fasta'):
                assembly, N50, NG50, Tot, Max, Min, N, Gaps, both, single = runDefaults(i)
                res[i.split('/')[-1]]=[N50, NG50, Tot, Max, Min, N, Gaps, both, single]
                if args.BLAST_:
                    perc = runBLASTProtRec(args.fasta_for_BLAST, assembly, args.E_val)
                    res[i.split('/')[-1]]=[N50, NG50, Tot, Max, Min, N, Gaps, both, single, perc]
            else:
                cnt+=1
    elif args.assembly !='':
        if args.assembly.endswith('.fa') or args.assembly.endswith('.fasta'):
            N50, NG50, Tot, Max, Min, N, Gaps, both, single = runDefaults(args.assembly)
        else:
            print('ERROR: No FastA found!!\n\n make sure the file ends with .fa or .fasta\n\n')
            print(parser.description)
            sys.exit()
    else:
        print('ERROR:No Input Specified!!\n\nPlease Specify a FastA file or a directory containing FastA files\n\n')
        print(parser.description)
        sys.exit()

    if args.onlyProtRec:
        assemblyDB = makeBlastDB(args.fasta_for_BLAST)
        blastRep = tblastn
