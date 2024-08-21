#! /usr/bin/env python3
import numpy as np
import pandas as pd
import decimal
from decimal import Decimal
decimal.getcontext().rounding='ROUND_HALF_UP'
# theoretically it's not done, can be that I'll finish it in the OG script.
# I'll try to remember to also dump it here
class IES:
    def __init__(self,id, scaff, start, stop):
        self.id = id
        self.TAScarPos = int(id.split('.')[-1])
        self.scaf = scaff
        self.start = start
        self.stop = stop
        self.length = stop-start
        self.IRS = None
        self.modpos = {'fwd':[],'rev':[]}
        self.relModPosFlat = {'fwd':[],'rev':[]}
        self.relmodpos = {orient:np.array([0 for _ in np.arange(-3,3,0.01)]) for orient in ['fwd','rev']}
        self.modProbs={'fwd':[],'rev':[]}
        self.modConsSeq = {'naive': {'fwd':[], 'rev':[]},
                           'relative':{'fwd':[], 'rev':[]}}
        self.plmn50Seq=None
        self.seq=None
        self.Acnt=None
    def addIRS(self, IRS):
        #IRS=pd.read_table(IRSfile)
        d = {k:v[0] for k,v in IRS[IRS['ID']==self.id].iloc[:,2:].to_dict(orient='list').items()}
        self.IRSs=d
    def addSeqs(self, ref):
        self.seq=ref[self.scaf][self.start-1:self.stop]
        self.plmn50Seq=ref[self.scaf][self.start-51:self.stop+50]
    def addModSeq(self, ref, Rel=False):
        if Rel:
            for orient, lst in self.relModPosFlat.items():
                # print(lst)
                for pos in lst:
                    if orient == 'fwd':
                        self.modConsSeq['relative'][orient].append(ref[self.scaf][pos-11:pos+10])
                    else:
                        self.modConsSeq['relative'][orient].append(Seq(ref[self.scaf][pos-11:pos+10]).reverse_complement())
        else:
            for orient, lst in self.modpos.items():
                for pos in lst:
                    if orient == 'fwd':
                        self.modConsSeq['naive'][orient].append(ref[self.scaf][pos-11:pos+10])
                    else:
                        self.modConsSeq['naive'][orient].append(Seq(ref[self.scaf][pos-11:pos+10]).reverse_complement())
    def countAs(self):
        cnt=0
        for nuc in self.plmn50Seq:
            if nuc.upper() == 'A':
                cnt+=1
        self.Acnt=cnt
    def normRelCount(self):
        from decimal import Decimal
        for orient, array in self.relmodpos.items():
                self.relmodpos[orient] = np.array([round(Decimal(i),5) for i in array/self.Acnt])
    def checkHemiMethyl(self):
        pass
