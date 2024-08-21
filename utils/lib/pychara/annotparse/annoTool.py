#! /usr/bin/env python3
## to import:
##      from pychara.annotparse.annoTool import *
from datetime import datetime
from time import time
# import os,sys, argparse
from Bio import SeqIO
from Bio.SeqUtils import GC
import pandas as pd

class TE:
    def __init__(self, line):
        self.line=line.strip()
        self.id= self.line.split('ID=')[-1].split(';')[0]
        self.scaf= self.line.split('\t')[0]
        self.start=min( int(self.line.split('\t')[3]), int(self.line.split('\t')[4]))
        self.stop=max(int(self.line.split('\t')[3]), int(self.line.split('\t')[4]))
        self.orientation=self.line.split('\t')[6]
        self.length=self.stop-self.start
        self.sequence=None
        self.GC=None
    def setSeqGC(self, genome):
        self.sequence = mic[self.scaf][self.start-1:self.stop]
        self.GC=GC(self.sequence)

def parseTE(TEanno,MIC=None):
    TEs = []
    TEdct = {}
    if MIC:
        mic = {i.id:str(i.seq) for i in SeqIO.parse(MIC,'fasta')}
    for line in open(TEanno):
        if MIC:
            te=TE(line).setSeqGC(mic)
        else:
            te=TE(line)
        TEs.append(te)
        if te.scaf not in TEdct:
            TEdct[te.scaf]=[te]
        else:
            TEdct[te.scaf].append(te)
    return TEs, TEdct

class IES:
    def __init__(self,id, scaff, start, stop, seq):
        self.id = id
        self.TAScarPos = int(id.split('.')[-1])
        self.scaf = scaff
        self.start = start
        self.stop = stop
        self.length = stop-start
        self.sequence = seq
        self.GC = GC(self.sequence)
        self.inGene = None
        self.allMods = {'A':0, 'T':0, 'G':0, 'C':0, 'N':0}
        self.allModsPos = {'A':[], 'T':[], 'G':[], 'C':[]}
        self.allModsQv = {'A':[], 'T':[], 'G':[], 'C':[], 'N':[]}
        self.mods = {'m6A':0, 'm4C':0, 'm5C':0}
        self.modsPos = {'+':
                            {'m6A':[],'m4C':[],'m5C':[]},
                        '-':
                            {'m6A':[],'m4C':[],'m5C':[]}
                        }
        self.modsIdenQv = {'+':
                                {'m6A':[],'m4C':[],'m5C':[]},
                           '-':
                                {'m6A':[],'m4C':[],'m5C':[]}
                            }
        self.IRSs=None
    def countMods(self):
        allMods, allIdenMods = sum(self.allMods.values()) , sum(self.mods.values())
        return allMods, allIdenMods
    def addIRS(self, IRS):
        # IRS=pd.read_table(file)
        d = {k:v[0] for k,v in IRS[IRS['ID']==self.id].iloc[:,2:].to_dict(orient='list').items()}
        self.IRSs=d
def prepGENEAnnoforIES(anno):
    anno=[i for i in open(anno) if not i.startswith('#')]
    annoDct = {i.split('\t')[0]:[] for i in anno}
    for i in anno:
        if i.split('\t')[2] == 'gene':
            annoDct[i.split('\t')[0]].append(i)
    return annoDct

def checkInGene(anno, ies):
            if ies.scaf.split('_with')[0] in anno:
                for i in anno[ies.scaf.split('_with')[0]]: #== i.split('\t')[0] and 'gene' == i.split('\t')[2]:
                    startPos=min(int(i.split('\t')[3]), int(i.split('\t')[4]))
                    stopPos=max(int(i.split('\t')[3]), int(i.split('\t')[4]))
                    if startPos <= ies.TAScarPos <= stopPos:
                        ies.inGene = True

def parseAnnoIES(annoIES,annoforGENE,IRSfile):
    IRS=pd.read_table(IRSfile)
    start = time()
    annoGENE=prepGENEAnnoforIES(annoforGENE)
    IESs=[IES(i.split('ID=')[-1].split(';')[0],i.split('\t')[0],int(i.split('\t')[3]), int(i.split('\t')[4]), i.split('quence=')[-1].strip()) for i in open(annoIES) if i.startswith('#') != True]
    IESdct={}
    for i in IESs:
        checkInGene(annoGENE, i)
        if i.id in IRS['ID'].to_list():
            i.addIRS(IRS)
        if i.scaf not in IESdct:
            IESdct[i.scaf]=[i]
        else:
            IESdct[i.scaf].append(i)
    stop=time()
    print('Done ! Elapsed time: '+str(stop-start))
    return IESs, IESdct



class GENE:
    def __init__(self,ID, scaf):
        self.id=ID
        self.sid=ID[-7:]
        self.scaf=scaf
        self.mRNApos=None
        self.CDSs=[] ## start/stop positions
        self.CDScnt=None  ## number of CDSs
        self.intronCNT=None  ## number of Introns
        self.hasIntron=None   ## True if it has introns
        self.introns=[]   ## start/stop introns
        self.orientation=None
        self.spliceJunc=[]  ## splicejunction, donor/acceptor positions
        self.multUTR5=False
        self.UTR5=[]#None
        # if self.UTR5 == False:
        #     if elf.orientation == '-':
        #         self.UTR5 = [max(sorted(self.CDSs, key=lambda x: x[0],reverse=True)[0])+500,min(sorted(self.CDSs, key=lambda x: x[0],reverse=True)[0])]
        #     else:
        #         self.UTR5 = [min(sorted(self.CDSs, key=lambda x: x[0])[0])-500,max(sorted(self.CDSs, key=lambda x: x[0],reverse=True)[0])]
        self.UTR3=[]#None
        self.multUTR3=False
        self.allMods={'A':0, 'T':0, 'G':0, 'C':0}
        self.allModsPos={'A':[], 'T':[], 'G':[], 'C':[]}
        self.mods={'m6A':0,'m4C':0,'m5C':0}  ## all modifications in gene
        self.modsPos={'m6A':[],'m4C':[], 'm5C':[]}  ## position of modifications
        self.modsIdenQv={'m6A':[],'m4C':[], 'm5C':[]}
        self.modFeature={'CDS':0,'intron':0,"3'UTR":0,"5'UTR":0}
        self.modsIdenCNT={
        'm6A':{'CDS':0,'intron':0,"3'UTR":0,"5'UTR":0},
        'm4C':{'CDS':0,'intron':0,"3'UTR":0,"5'UTR":0},
        'm5C':{'CDS':0,'intron':0,"3'UTR":0,"5'UTR":0}
        }
        #################
        #  self.modsIdenCNTNorm not implemented yet
        #################
        self.modsIdenCNTNorm={
        'm6A':{'CDS':[],'intron':[],"3'UTR":[],"5'UTR":[]},
        'm4C':{'CDS':[],'intron':[],"3'UTR":[],"5'UTR":[]},
        'm5C':{'CDS':[],'intron':[],"3'UTR":[],"5'UTR":[]}
        }   ######## normalised against feature length; if SINGLE FEATURE has 2 mods, it will
            ######## be normalised against that. So you get different floats that get summed up
        # if self.UTR5 == False:
        #     if elf.orientation == '-':
        #         self.UTR5 = [max(sorted(self.CDSs, key=lambda x: x[0],reverse=True)[0])+500,min(sorted(self.CDSs, key=lambda x: x[0],reverse=True)[0])]
        #     else:
        #         self.UTR5 = [min(sorted(self.CDSs, key=lambda x: x[0])[0])-500,max(sorted(self.CDSs, key=lambda x: x[0],reverse=True)[0])]
        self.ATG=None
        self.stop=None
        self.coding=None
        self.done=False
    def check(self,extUTRs=False):
        if extUTRs:
            checkUTRs(self)
        else:
            checkCoding(self)
        checkCDSsIntronsNC(self)
        # distModCNTForNorm(self)
        # normCNTs(self)
    def checkQuant(self):
        initNormCNTdct(self)
        distModCNTForNorm(self)
        normCNTs(self)
    # def distMods(self):
    #     for modType,pos in self.modsPos.items():
    #         for mod in pos:
    #             if self.CDScnt > 1:
    #                 for cds in self.CDSs:
    #                     if min(cds) <= mod <= max(cds):
    #                         self.modsCNT[modType]['CDS']+=1
    #             else:
    #                 if min(self.CDSs[0]) <= mod <= max(self.CDSs[0]):
    #                     self.modsCNT[modType]['CDS']+=1
    #             if self.intronCNT > 1:
    #                 for intron in self.introns:
    #                     if min(intron) <= mod <= max(intron):
    #                         self.modsCNT[modType]['intron']+=1
    #             else:
    #                 if min(self.introns[0]) <= mod <= max(self.introns[0]):
    #                     self.modsCNT[modType]['intron']+=1
    #             if self.multUTR3:
    #                 for utr3 in self.UTR3:
    #                     if min(utr3) <= mod <= max(utr3):
    #                         self.modsCNT[modType]["3'UTR"]+=1
    #             else:
    #                 if min(self.UTR3[0]) <= mod <= max(self.UTR3[0]):
    #                     self.modsCNT[modType]["3'UTR"]+=1
    #             if self.multUTR5:
    #                 for utr5 in self.UTR5:
    #                     if min(utr5) <= mod <= max(utr5):
    #                         self.modsCNT[modType]["5'UTR"]+=1
    #             else:
    #                 if min(self.UTR5[0]) <= mod <= max(self.UTR5[0]):
    #                     self.modsCNT[modType]["5'UTR"]+=1


def checkCoding(self):
    if self.done:
        if self.coding == None:
            self.coding = True

def checkUTRs(self):
    if self.done:
        if self.coding == None:
            self.coding = True
        # print('CHECKUTRS init: ', self.UTR5, self.UTR3, self.orientation, self.coding)
        if self.UTR5 == list() and self.coding:
            # print(self.UTR5, self.coding)
            if self.orientation == '-':
                # print('CDSs:', self.CDSs, 'scaf & id', self.scaf, self.id)
                self.UTR5.append( (max(sorted(self.CDSs, key=lambda x: x[0],reverse=True)[0])+251,max(sorted(self.CDSs, key=lambda x: x[0],reverse=True)[0])+1)  )
                # print('CHECKUTRS after: ', self.UTR5)
            else:
                self.UTR5.append( (min(sorted(self.CDSs, key=lambda x: x[0])[0])-251,min(sorted(self.CDSs, key=lambda x: x[0])[0])-1) )
                # print('CHECKUTRS after: ', self.UTR5)
        if self.UTR3 == list() and self.coding:
            # print('CHECKUTRS init: ', self.UTR3)
            if self.orientation == '-':
                self.UTR3.append( (min(sorted(self.CDSs, key=lambda x: x[0])[0])-101,min(sorted(self.CDSs, key=lambda x: x[0])[0])-1) )
                # print('CHECKUTRS after: ', self.UTR3)
            else:
                self.UTR3.append( (max(sorted(self.CDSs, key=lambda x: x[0],reverse=True)[0])+101,max(sorted(self.CDSs, key=lambda x: x[0],reverse=True)[0])+1) )
                # print('CHECKUTRS after: ', self.UTR3)


def checkCDSsIntronsNC(self):
    if self.done:
        self.CDScnt=len(self.CDSs)
        # print(self.CDScnt, self.CDScnt-1)
        # self.intronCNT=self.CDScnt-1
        # print(self.intronCNT)
        # if self.intronCNT > 0:
        #     self.hasIntron=True
        # else:
        #     self.hasIntron=False
        for i in range(len(self.CDSs)):
            if i == len(self.CDSs)-1:
                pass
            else:
                self.introns.append((self.CDSs[i][1]+1, self.CDSs[i+1][0]-1))
                # print('CDSs:',self.CDSs, 'Introns:',self.introns)
        # if self.coding == None:
        #     self.coding = True
        if len(self.UTR5) > 1:
            self.multUTR5 = True
            for i in range(len(self.UTR5)):
                if i == len(self.UTR5)-1:
                    pass
                else:
                    self.introns.append( (self.UTR5[i][1]+1, self.UTR5[i+1][0]-1 ))
        if len(self.UTR3) > 1:
            self.multUTR3 = True
            for i in range(len(self.UTR3)):
                if i == len(self.UTR3)-1:
                    pass
                else:
                    self.introns.append( (self.UTR3[i][1]+1, self.UTR3[i+1][0]-1 ))
        self.intronCNT=len(self.introns)
        if self.intronCNT > 0:
            self.hasIntron=True
        else:
            self.hasIntron=False


def parseAnnoGENE(anno, extUTRs=False):
    start= time()
    print('['+str(datetime.now()).replace(' ' , '|').split('.')[0]+'] Parsing gene annotations . . . ')
    genes = {}
    anno = [i for i in open(str(anno)) if not i.startswith('#')]
    id = anno[0].split('ID=')[-1].split(';')[0]
    orient = anno[0].split('\t')[6]
    scaf = anno[0].split('\t')[0]
    gene = GENE(id,scaf)
    pos = (min(int(anno[0].split('\t')[3]) ,int(anno[0].split('\t')[4])), max(int(anno[0].split('\t')[3]) ,int(anno[0].split('\t')[4]))   )
    gene.orientation=orient
    gene.mRNApos = pos
    # print(gene.id, gene.sid, gene.orientation,gene.scaf)
    if 'ncRNA' in anno[0]:
        gene.coding=False
        # print(gene.coding)
    for i in anno:
        # print(i)
        # gene.orientation=i.split('\t')[6]
        if gene.sid not in i:
            ##### finish current gene object, add it to genes list and generate new gene object
            gene.done = True
            gene.check(extUTRs)
            # print('3UTR:',gene.UTR3,'5UTR:',gene.UTR5, )
            # if extUTRs:
            #     gene.checkUTRs()
            # gene.checkCDSs()
            if gene.scaf in genes:
                genes[gene.scaf].append(gene)
            else:
                genes[gene.scaf]=[gene]
            # genes.append(gene)
            ##### load new gene object with ID and orientation
            id = i.split('ID=')[-1].split(';')[0]
            scaf = i.split('\t')[0]
            gene = GENE(id,scaf)
            pos = (min(int(i.split('\t')[3]) ,int(i.split('\t')[4])), max(int(i.split('\t')[3]) ,int(i.split('\t')[4]))   )
            gene.mRNApos = pos
            gene.orientation=i.split('\t')[6]
            # print('new gene:',gene.id, gene.sid, gene.scaf)
            continue
        if 'ncRNA' in i:
            gene.coding = False
        if 'CDS' in i and gene.id == id:
            gene.CDSs.append((int(i.split('\t')[3]), int(i.split('\t')[4])))
            # print('CDSs:',gene.CDSs)
        if 'three_prime_UTR' in i and gene.id == id:
            startUTR = min( int(i.split('\t')[4]),int(i.split('\t')[3]))
            stopUTR = max( int(i.split('\t')[4]),int(i.split('\t')[3]))
            if stopUTR - startUTR == 0:
                continue
            gene.UTR3.append( (int(i.split('\t')[3]), int(i.split('\t')[4])) )
            # print('3UTR:',gene.UTR3)
        if 'five_prime_UTR' in i and gene.id == id:
            startUTR = min( int(i.split('\t')[4]),int(i.split('\t')[3]))
            stopUTR = max( int(i.split('\t')[4]),int(i.split('\t')[3]))
            if stopUTR - startUTR == 0:
                continue
            gene.UTR5.append( (int(i.split('\t')[3]), int(i.split('\t')[4])) )
            # print('5UTR:',gene.UTR5)
    stop=time()
    print('['+str(datetime.now()).replace(' ' , '|').split('.')[0]+'] Parsing gene annotations . . . Done! Elapsed time: '+str(stop-start))
    return genes
