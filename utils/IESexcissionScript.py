#!/usr/bin/env python
from Bio import SeqIO
import os, sys
import regex as re
import argparse as ap
import time
########################
##### for output put in all IRS scores or dcl23, dcl235, dcl5, nowa12, iswi, ezl1
########################
## hardcoded dictionary; contains first nucleotide of first IES and last nucleotide of last IES
pointers={'scaffold51_57_with_IES': [74, 10493],
          'scaffold51_78_with_IES': [55, 6378],
          'scaffold51_150_with_IES': [125, 4609],
          'scaffold51_187_with_IES': [50, 3683]
         }

## IES and Read classes
class IES:
    def __init__(self,id, scaff, start, stop):
        self.id = id
        self.scaff = scaff
        self.start = start
        self.stop = stop
        self.position = 0 ## which ies it is in the fragment -8e.g. first == 0, second == 1)
        self.cntNaive = 0
        self.cntStrc = 0
    def add_cntNaive(self):
        self.cntNaive += 1
    def add_cntStrc(self):
        self.cntStrc += 1
    def setPosition(self, nmb):
        self.position += int(nmb)

class READ:
    def __init__(self, id, length, scaffold):
        self.id = id
        self.length = length
        self.scaff = scaffold
        self.structure = []
        self.parseCIGAR = []
    def defineStructure(self, list):
        self.structure += list
    # def defineStructureV2(self, bool):
    #     self.structure.append(bool)
    def addCIGARinfo(self, lst):
        self.parseCIGAR += lst

def parseIESs(annoFile):
    '''
    gets IES infos and puts them in a list containing IES object
    takes IES annotation file as argument
    '''
    print('. . . generating IES objects')
    anno=open(str(annoFile)).readlines()
    IESs=[IES(i.split('ID=')[-1].split(';')[0],i.split('\t')[0],int(i.split('\t')[3]), int(i.split('\t')[4])) for i in anno if i.startswith('#') != True]
    IESsDict={i.split('\t')[0]:[] for i in anno if i.startswith('#') != True}
    for i in anno:
        if i.startswith('#') != True:
            IESsDict[i.split('\t')[0]].append(IES(i.split('ID=')[-1].split(';')[0],i.split('\t')[0], int(i.split('\t')[3]), int(i.split('\t')[4])))
    return IESs , IESsDict

# def parseIESsdict(annoFile):
#     '''
#     gets IES info and puts them in a dictionary with the corresponding scaffold as key
#     needed later for the structure and to get the IES position in the scaffold
#     '''
#     anno=open(str(annoFile)).readlines()
#     IESs={i.split('\t')[0]:[] for i in anno}
#     for i in anno:
#         IESs[i.split('\t')[0]].append(IES(i.split('ID=')[-1].split(';')[0],i.split('\t')[0], int(i.split('\t')[3]), int(i.split('\t')[4])))
#     return IESs

def getIESPositions(IESdct, IESs):
    '''
    get the IESs position in the scaffold/fragment. E.g. if it's the first on the fragment it has 0, if the second it's 1
    '''
    print('. . . getting IES positions on PCR-fragment')
    for ies in IESs:
        for scaf,iesList in IESdct.items():
            if scaf == ies.scaff:
                for Ies in iesList:
                    if ies.id == Ies.id:
                        ies.setPosition(iesList.index(Ies))

def parseCIGAR(line):
    '''
    gets Cigar info (deletions and their positions) only when the read spans the locus from first to last IES
    takes line of a sam file aas argument
    '''
    cntRead=0
    cntScaf=int(line.split('\t')[3])
    deletions=[]
    if cntScaf <= min(pointers[line.split('\t')[2].split('_loc')[0]]):
        for i in re.findall(r'(\d+)(\w)', line.split('\t')[5]):
            if i[1] == 'N' or i[1]=='D':
                if int(i[0]) > 10:
                    cntScaf+=int(i[0])
                    deletions.append((i, cntScaf)) ### maybe put this before adding the deletions
                    continue ## to see if it adds up
                cntScaf+=int(i[0]) ## same as above
            if i[1] == 'M':
                cntRead+=int(i[0])
                cntScaf+=int(i[0])
            if i[1] == 'S' or i[1] == 'I':
                cntRead+=int(i[0])
    if cntScaf >= max(pointers[line.split('\t')[2].split('_loc')[0]]):
        return line.split('\t')[2].split('_loc')[0], cntRead, cntScaf , deletions
    else:
        return '','','',''

def processParseCIGAR_cntNaive(cigOut, IESs):
    '''
    processes the output of parseCIGAR(), adds the naive counts per IES
    takes ParseCIGAR() output and IES objects list as arguments
    '''
    for i in cigOut:
        if i[0] == '':
            continue  ## i[0]==scaffold, i[1]==readLength, i[2]==coveredScaffold, i[3]==listWithDeletions(CIGs)[(int,D),int]
        for ies in IESs: ## IESs == list of IES objects
            if i[0] == ies.scaff:
                for j in i[3]:
                    print(int(j[1]),'-', int(j[0][0]), '=',int(j[1])- int(j[0][0]),'\t', int(j[1]))
                    if ies.start-4 <= int(j[1])-int(j[0][0]) <= ies.start+4 and ies.stop-4 <= int(j[1]) <= ies.stop+4:# or ies.start-5 <= int(i[3][1])-int(i[3][0][0]) <= ies.stop-5 or ies.start+5 <= int(i[3][1])-int(i[3][0][0]) <= ies.stop+5:
                        ies.add_cntNaive()
                    else:
                        continue

def processParseCIGAR_cntNaive_V2(cigOut, ies):  #### NOT DONE YET!
    '''
    processes the output of parseCIGAR(), adds the naive counts per IES
    takes ParseCIGAR() output and IES objects list as arguments
    '''
    # for i in cigOut:
        # if i[0] == '':
        #     continue  ## i[0]==scaffold, i[1]==readLength, i[2]==coveredScaffold, i[3]==listWithDeletions(CIGs)[(int,D),int]
        # for ies in IESs: ## IESs == list of IES objects
        #     if i[0] == ies.scaff:
    for j in cigOut:
        # print(int(j[1]),'-', int(j[0][0]), '=',int(j[1])- int(j[0][0]),'\t', int(j[1]))
        if ies.start-5 <= int(j[1])-int(j[0][0]) <= ies.start+5 and ies.stop-5 <= int(j[1]) <= ies.stop+5:# or ies.start-5 <= int(i[3][1])-int(i[3][0][0]) <= ies.stop-5 or ies.start+5 <= int(i[3][1])-int(i[3][0][0]) <= ies.stop+5:
            ies.add_cntNaive()
        else:
            continue

def getStructure(line, IESs, cigar):  #### NOT DONE YET !!! WORKS SOMEWHAT, have to check for overly fuzzy deletions(keep them appart and leave them)
    '''
    gets structure of the read by going through the sam line
    takes line of a sam file, IESs object list and split CIGAR as arguments
    '''
    struct=[]
    skippedDeletions=[]
    checkedCig=[]
    weird=[]
    cnt=0
    # print('CIGAR!!: ', cigar)
    for scaf, iess in IESs.items():
        if scaf == line.split('\t')[2].split('_lo')[0]:
            for ies in iess:#range(0,len(iess)-1
                for j in cigar:
                    if j not in checkedCig:# or not in skippedDeletions:
                        # if iess[i] != iess[i-cnt]:
                        if ies.stop+4 < int(j[1])-int(j[0][0]):
                            struct.append(False)
                            # print('if 1: ',ies.id,iess.index(ies), len(struct), struct)
                            break
                        elif ies.start-5 <= int(j[1])-int(j[0][0]) <= ies.start+5 and ies.stop-5 <= int(j[1]) <= ies.stop+5:
                            struct.append(True)
                            checkedCig.append(j)
                            # print('elif 1: ', ies.id,iess.index(ies),len(struct), struct)
                            break
                        elif ies.start > int(j[1]):
                            struct.append(False)
                            checkedCig.append(j)
                            skippedDeletions.append(j)
                            # print('elif 2: ', ies.id,iess.index(ies),len(struct), struct)
                            break
                        else:
                            struct.append(False)
                            checkedCig.append(j)
                            skippedDeletions.append(j)
                            weird.append(j)
                            # print('else: ', ies.id,iess.index(ies),len(struct), struct)
                            break
                if cigar != []:#] or cigar != []:
                    if ies.start > cigar[-1][-1]: # and cigar != []:
                        struct.append(False)
                        # print('if 2: ', ies.id,iess.index(ies),len(struct), struct)
    return struct, checkedCig, skippedDeletions, weird

def getSamInfo(sam, IESs, IESsDict):  #### NOT DONE YET !!! WORKS also w the defineSturcture function
    '''
    goes through samfile and collects data as read structure, CIGAR info and (hopefully) converts them to counts
    takes sam file and IES objects list as arguments
    '''
    print('. . . gathering intel from sam')
    Reads={'scaffold51_57_with_IES': [],
           'scaffold51_78_with_IES': [],
           'scaffold51_150_with_IES': [],
           'scaffold51_187_with_IES': []
           }
    print('\tloading sam')
    insam=open(str(sam)).readlines()
    for line in insam:
        print('\tprocessing line '+str(insam.index(line)+1)+' of '+str(len(insam)), end='\t\t\r')
        #parseCig=[]
        if line.startswith('@')!=True and (line.split('\t')[2].split('_loc')[0] in pointers) == True:
            parseCig=parseCIGAR(line)
            if parseCig[0] == '':
                continue
            # print(parseCig,'\n', parseCig[3], 'len: ', len(parseCig[3]), 'type: ', type(parseCig[3]), '\n', parseCig[0])
            Reads[line.split('\t')[2].split('_lo')[0]].append(READ(line.split('\t')[0], len(line.split('\t')[9]), line.split('\t')[2])) #.addCIGARinfo(parseCig).defineStructure(getStructure(line, IESs, parseCig[3]))
            Reads[line.split('\t')[2].split('_lo')[0]][-1].addCIGARinfo(parseCig[3])
            Reads[line.split('\t')[2].split('_lo')[0]][-1].defineStructure(getStructure(line, IESsDict, parseCig[3])[0])
    print('\tprocessing line'+str(insam.index(line)+1)+' of '+str(len(insam)))
    return Reads

##  MAYBE OBSOLETE?? MAYBE NOT
def strucCnt(structList, structure, ies):
    '''
    gets the counts for the structural stuff ## maybe have the structure in 1,0
    format and just multiply the by the number of reads that have this excission
    '''
    if structure not in structList:
        for i in range(0, len(structure)-1):
            if ies.position == i and structure[i] == True:
                ies.add_cntStrc()
                return True
            # else:
            #     return False

# def getCounts(Reads, IESs):    #### NOT DONE YET!!!
#     '''
#     goes through the read objects that hold the infos and generates the counts
#     '''
#     strucCnts=[]
#     checkedStrucCnts=[]
#     # checked=[]
#     checkedReads=[]
#     total_reads = {i:0 for i in Reads.keys()}
#     total_readsStruc = {i:0 for i in Reads.keys()}
#     for key in Reads.keys():
#         if Reads[key] != []:
#             # checkedReads=[]
#             for read in Reads[key]:
#                 # checked=[]
#                 if read.parseCIGAR != []: #read not in checkedReads and read.parseCIGAR != []:
#                     if
#                     for ies in IESs:
#                         # checked=[]
#                         if ies.scaff==key:#ies not in checked and ies.scaff == key:
#                             # for read in Reads[key]:
#                             #     # print(read.id, read.parseCIGAR)
#                             #     # if read.parseCIGAR[0] != '':# and read not in checkedReads:
#                             #     if read not in checkedReads:
#                             processParseCIGAR_cntNaive_V2(read.parseCIGAR, ies)
#                             # if strucCnt(strucCnts, read.structure, ies):
#                             if strucCnt(strucCnts, read.structure, ies):
#                                 checkedStrucCnts.append(True)
#                             if len(checkedStrucCnts) == len(read.structure):
#                                 total_readsStruc[key] +=1
#                                 # strucCnts.append(read.structure)
#                                 # checked.append(ies)
#                                 # total_readsStruc[key] +=1
#                             # total_readsStruc[key] +=1
#                             # checked.append(ies)
#                             strucCnts.append(read.structure)
#                 # checkedReads.append(read)   #### MAYBE HAVE THE CHECKED READS IN THE LOOP?#
#                 total_reads[key]+=1
#     return total_reads, total_readsStruc
# ###################################################################################
#  alternative for getCounts() that ACTUALLY WORKS !!!!
###################################################################################
#
# checked_structures=[]
# totStrc=0
# totReads=0
# for key in Reads.keys():
#     if Reads[key] != []:
#         for read in Reads[key]:
#             if read.parseCIGAR != []:
#                 # if read.structure not in checked_structures:
#                 for ies in IESs:
#                     if ies.scaff == key:
#                         processParseCIGAR_cntNaive_V2(read.parseCIGAR, ies)
#                         strucCnt(checked_structures, read.structure, ies)
#                 checked_structures.append(read.structure)
#             else:
#                 totReads+=1
#                 continue
#             totStrc+=1
#             totReads+=1

##########
# v2 as above  WORKS BETTER THAN ABOVE;;; MAYBE TAKE THIS ONE !!!!
##########

def strucCnt_V2(checked_structures, structure, ies):
    '''
    gets the counts for the structural stuff ## maybe have the structure in 1,0
    format and just multiply the by the number of reads that have this excission
    '''
    if ';'.join([str(i) for i in structure]) not in checked_structures:
        for i in range(0, len(structure)):
            if ies.position == i and structure[i] == True:
                ies.add_cntStrc()
                return True

def getCounts_V2(Reads, IESs):
    '''
    goes through the read objects that hold the infos and generates the counts
    '''
    print('. . . generating counts')
    checked_structures={key:{} for key in Reads.keys()}
    totStrc={key:0 for key in Reads.keys()}       ##### maybe make total reads per scaffold; don't know if it will change much but who knows
    totReads={key:0 for key in Reads.keys()}      #####
    for key in Reads.keys():
        print('\tprocessing scaffold: '+str('_'.join(key.split('_')[0:2])))
        if Reads[key] != []:
            for read in Reads[key]:
                print('\t    processing read '+str(Reads[key].index(read)+1)+' of '+str(len(Reads[key])), end='\t\t\r')
                if read.parseCIGAR != []:
                    # if read.structure not in checked_structures:
                    for ies in IESs:
                        if ies.scaff == key:
                            processParseCIGAR_cntNaive_V2(read.parseCIGAR, ies)
                            strucCnt_V2(checked_structures, read.structure, ies)
                    if ';'.join([str(i) for i in read.structure]) in checked_structures[key]:
                        checked_structures[key][';'.join([str(i) for i in read.structure])] += 1
                    else:
                        checked_structures[key][';'.join([str(i) for i in read.structure])] = 1
                else:
                    totReads[key]+=1
                    continue
                totStrc[key]+=1
                totReads[key]+=1
            print('\t    processing read '+str(Reads[key].index(read)+1)+' of '+str(len(Reads[key])))
    # print('\t    processing read '+str(Reads[key].index(read)+1)+' of '+str(len(Reads[key])), end='\t\t\r')
    return totStrc, totReads, checked_structures

def writeOutfile(IESs, total_reads, total_readsStruc, checked_structures, IRSscores):
    '''
    Writes out 2 files; once the calculated scores and a file containing structure
    information and structure count
    '''
    print('. . . writing results to file')
    IRS=open(str(IRSscores)).readlines()
    IRSs={i.split('\t')[0]:'\t'.join(i.split('\t')[2:-4]) for i in IRS}
    with open(str(time.localtime()[0])+'0'+str(time.localtime()[1])+str(time.localtime()[2])+'IES_excision.tsv', 'w') as w:
        w.write('## Output of the IES excision pipeline, tab delimited format\
        \n##IES ID\tscaffold\tposition in PCR fragment\tstart\tstop\ttotal reads\ttotal unique structures\traw naive counts\t\
        raw structural counts\tnaive counts/total reads\tstructural counts/total unique structures\t'+'\t'.join(IRS[0].split('\t')[2:-4])+'\n')
        for ies in IESs:
            if total_reads[ies.scaff] != 0:
                if len(checked_structures[ies.scaff]) != 0:
                    w.write(str(ies.id)+'\t'+str(ies.scaff)+'\t'+str(ies.position)+'\t'+str(ies.start)+\
                    '\t'+str(ies.stop)+'\t'+str(total_reads[ies.scaff])+'\t'+str(total_readsStruc[ies.scaff])+'\t'+str(ies.cntNaive)+\
                    '\t'+str(ies.cntStrc)+'\t'+str(ies.cntNaive/total_reads[ies.scaff])+'\t'+str(ies.cntStrc/len(checked_structures[ies.scaff]))+\
                    '\t'+IRSs[ies.id]+'\n')
                else:
                    w.write(str(ies.id)+'\t'+str(ies.scaff)+'\t'+str(ies.position)+'\t'+str(ies.start)+\
                    '\t'+str(ies.stop)+'\t'+str(total_reads[ies.scaff])+'\t'+str(total_readsStruc[ies.scaff])+'\t'+str(ies.cntNaive)+\
                    '\t'+str(ies.cntStrc)+'\t'+str(ies.cntNaive/total_reads[ies.scaff])+'\tNA'+\
                    '\t'+IRSs[ies.id]+'\n')
            else:
                w.write(str(ies.id)+'\t'+str(ies.scaff)+'\t'+str(ies.position)+'\t'+str(ies.start)+\
                '\t'+str(ies.stop)+'\tNA\tNA\t'+str(ies.cntNaive)+\
                '\t'+str(ies.cntStrc)+'\tNA\tNA\t'+IRSs[ies.id]+'\n') #+str(ies.cntNaive/total_reads[ies.scaff])+'\t'+str(ies.cntStrc/len(checked_structures[ies.scaff]))+\
                #'\t'+IRSs[ies.id]+'\n')
    with open(str(time.localtime()[0])+'0'+str(time.localtime()[1])+str(time.localtime()[2])+'ReadStructures.tsv', 'w') as w:
        w.write('## Readstructures; True = excised; False=retained\n## Scaffold\tStrucutre\tcounts\n')
        for scaf, structures in checked_structures.items():
            w.write('# '+str(scaf)+'\n')
            for structure, count in structures.items():
                w.write(str(scaf)+'\t'+str(structure)+'\t'+str(count)+'\n')

def main(insam, IRS_file, anno):
    '''
    main function that executes the programme
    '''
    IESs, IESsDict = parseIESs(anno)
    getIESPositions(IESsDict, IESs)
    Reads=getSamInfo(str(insam), IESs, IESsDict)
    totalStructure, totalReads, checkedStructures = getCounts_V2(Reads, IESs)
    writeOutfile(IESs, totalReads, totalStructure, checkedStructures, str(IRS_file))

if __name__ == '__main__':
    #parser = ap.ArgumentPareser(description='Program to calculate the IRS for IESs in context to of one single locus / PCR fragment') #### put argparse to have be nice ### NEEDS SAM AND IRS-FILE AS INPUT
    if len(sys.argv[1:]) != 3:
        print('missing arguments!\n>>> USAGE: IESexcissionScript.py insam IRS_file IES_annotationFile')
        sys.exit()

    insam = sys.argv[1]
    IRS_file = sys.argv[2]
    IESanno = sys.argv[3]

    main(insam, IRS_file, IESanno)
