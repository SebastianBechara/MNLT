#! /usr/bin/env python3
import pandas as pd
import os, sys
# FILE='/home/bechara/SSD2/DeepSeqData/SarahIrisPolII_III/1_1735_V1_assembledGCmin40.fa_v_IES.tsv'
#
# df=pd.read_table(FILE, header=None)
# df.columns = ['IES', 'Read', 'alnIden', 'alnLen', 'mismatch', 'gap', 'qStart', 'qEnd', 'qLen', 'sStart', 'sEnd', 'sLen', 'eVal', 'bit']
# dff=df.sort_values(['Read']).iloc[:50,:]
# READ = dff.iloc[0,]['Read']

# def allThere(l, x,y, sort=False):
#     if sort:
#         l.sort()
#     for i in range(x,y):
#         if i in l:
def allThere(l, cnt):
    if l == list(range(cnt)):
        return True
    return False

def getStuff(dff):
    IES_IES_bound =[]
    concatamers=[]
    READ = dff.iloc[0,]['Read']
    l=[[dff.iloc[0,]['IES'] , min(dff.iloc[0,]['sStart'],dff.iloc[0,]['sEnd']), max(dff.iloc[0,]['sStart'],dff.iloc[0,]['sEnd'])]]
    #l=[[row['IES'], min(row['sStart'],row['sEnd']), max(row['sStart'],row['sEnd'])]]
    cnt=0

    for index, row in dff.iterrows():
        # if row['alnIden'] == 100:
        #         pass
        #print(row['Read'])
        if row['Read'] == READ:
            #print(row['alnIden'])
            if row['alnIden'] == 100:# and row['alnLen']==row['qLen']:
                l.append([row['IES'], min(row['sStart'],row['sEnd']), max(row['sStart'],row['sEnd'])])
                #print(row['IES'], min(row['sStart'],row['sEnd']), max(row['sStart'],row['sEnd']), l)
                cnt+=1
                continue
            else: continue
        else:
            if l:
                bound=False
                #print(l)
                l.sort(key=lambda x:x[1])
                conIES=[]
                for pos,item in enumerate(l):
                    #print(pos, len(l)-1)
                    if pos == len(l)-1:#+1 == len(l):
                        #print('YO')
                        if len(l)!=1:
                            #print('YO')
                            #print(f'@pos+1; pos: {l[pos]}; pos-1: {l[pos-1]}; item {item}')
                            if l[pos-1][2] -5 <= item[1] <= l[pos-1][2] +5:
                                bound=True
                                #conIES.append((l[pos-1], pos-1))
                                #conIES.append((item, pos))
                                #print('last')
                                continue #break
                            else: continue
                        else:
                            continue #break
                    if pos == 0:
                        #print('first')
                        if l[pos+1][1]-5 <= item[2] <= l[pos+1][1]+5:
                            conIES.append((item, pos))
                            conIES.append((l[pos+1], pos+1))
                            bound=True
                            continue
                    #print('not last', l)
                    if l[pos+1][1]-5 <= item[2] <= l[pos+1][1]+5:
                        #conIES.append((item, pos))
                        conIES.append((l[pos+1], pos+1))
                        bound=True
                        #print(3)
                        continue
                #print(conIES, [i[1] for i in conIES], cnt)
                if conIES:
                    if allThere([i[1] for i in conIES] ,cnt):
                        if conIES[0][0][0] == conIES[-1][0][0]:
                            concatamers.append(l)
                        pass
                if bound:
                    IES_IES_bound.append(l)
                READ = row['Read']
                l=[[row['IES'], min(row['sStart'],row['sEnd']), max(row['sStart'],row['sEnd'])]]
                cnt=0
                continue
            else:
                READ = row['Read']
                l=[[row['IES'], min(row['sStart'],row['sEnd']), max(row['sStart'],row['sEnd'])]]
                cnt=0
                continue
    return IES_IES_bound, concatamers

def prepDF(FILE):
    df=pd.read_table(FILE, header=None)
    df.columns = ['IES', 'Read', 'alnIden', 'alnLen', 'mismatch', 'gap', 'qStart', 'qEnd', 'qLen', 'sStart', 'sEnd', 'sLen', 'eVal', 'bit']
    dff = dff.sort_values(['Read', 'sStart'])
    for ind, row in dff.iterrows():
        sStart = min(row['sStart'],row['sEnd'])
        sEnd = max(row['sStart'],row['sEnd'])
        row['sStart'] = sStart
        row['sEnd'] = sEnd
    return df

def getReadCount(seqF):
    cnt=0
    for i in open(seqF):
        if i.startswith('>'):
            cnt+=1
    return cnt

def getRes(IES_bound, concats, blast, seqF):
    totReads = getReadCount(seqF)
    print(f'Numbers for: {blast}\nThat many reads contain IES IES junctions: {len(IES_bound)} = {round((len(IES_bound)/totReads)*100,5)}%\nThat may reads seem to be full concatamers: {len(concats)} = {round((len(concats)/totReads)*100,5)}%\n\n{"#"*80}\n\n')
    with open(f'{blast}_Concats.tsv', 'w') as w:
        w.write('\n'.join(['\t'.join([str(i) for i in l]) for l in concats]))
    pass

# def getMeta(FILES, seqFiles):
#     dct={}
#     for blast in FILES:
#         if blast.split('_v_IES')[0] in seqFiles:
#             dct[blast]=
#         else:
#             print('Missing files! Check if you got all neccessarry files...')
#             sys.exit()
#     return dct

def prepDirs(blastPath, seqPath):
    blastFull=[f"{blastPath.rstrip('/')}/{file}" for file in os.listdir(blastPath)]
    blast=[file for file in os.listdir(blastPath)]
    seqFiles=[f"{seqPath.rstrip('/')}/{file}" for file in os.listdir(seqPath) if file.endswith('.fastq.fa')]
    return blast, blastFull, seqFiles

def main(FILES, seqFiles):
    BLASTs, BLASTsFull, SeqFiles = prepDirs(FILES,seqFiles)
    for blast in BLASTsFull:
        IES_bound, concats = getStuff(prepDF(blast))
        for seq in SeqFiles:
            if blast.split('_BlastRef_v_')[0].split('/')[-1] in seq:
                print(f'\n\n...processing {blast.split("_BlastRef_v_")[0].split("/")[-1]}', end='\r')
                getRes(IES_bound, concats, blast.split('/')[-1], seq)
                break

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print('Usage:\n     pol23.py <dir w/ blast reps> <dir w/ corresponding seq files>')
        sys.exit()
    if len(sys.argv[1:]) != 2:
        print('Missing arguments!\n\nUsage:\n     pol23.py <dir w/ blast reps> <dir w/ corresponding seq files>')
        sys.exit()
    blastFiles = sys.argv[1]
    seqFiles = sys.argv[2]
    main(blastFiles, seqFiles)
