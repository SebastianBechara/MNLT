#! /usr/bin/python3

import sys , os , argparse
import matplotlib.pyplot as plt
import numpy as np

def check_folder(title):
    ### maybe make the folders permanent, only let the programm delete the sam/fastq files
    # if os.path.isdir('/home/debian9/host/analysis/' + title) != True:
    #     os.mkdir('/home/debian9/host/analysis/' + title)
    if os.path.isdir('temp') != True:
        os.mkdir('temp')
    if os.path.isdir("temp/mapped") != True:
        os.mkdir("temp/mapped")
    if os.path.isdir("temp/unmapped") != True:
        os.mkdir("temp/unmapped")
    # if os.path.isdir('temp/refs') != True:
    #     os.mkdir('temp/refs')

def seq_files(deepseq_folder, min_size, max_size):
    seq_files=[]
    for file in os.listdir(str(deepseq_folder)):
        if file.endswith('nt.fastq'):
            if int(min_size) <= int(file.split('_')[-1].split('nt')[0]) <= int(max_size):
                print(file)
                seq_files.append(file)
    print(seq_files)
    return seq_files
###########################
###########################

# def runHisat(Refs): ## refs should have the name as string
#     for i in Refs:
#         cmd = ('hisat2-build '+)
#     for i in Refs:
#         cmd = 'hisat2 --un temp/unmapped/'+file+'.'+str(Refs)+'_unmapped.fq --al temp/mapped/'+file+'.'+str(Refs)+'_mapped.fq '

def getNamesList(fastq):
    lst=[i.lstrip('@') for i in open(fastq) if i.startswith('@')]
    with open(fastq+'.Names.list', 'w') as w:
        w.write('\n'.join(lst))
    return fastq+'.Names.list'
############################
############################
def mapSeq(deepseq_folder, min_size, max_size):
    seq_files_list = seq_files(deepseq_folder, min_size, max_size)
    print(seq_files_list)
    cnt=1
    for file in seq_files_list:
        print('Processing file '+str(cnt)+' of '+str(len(seq_files_list))+' files', deepseq_folder+'/'+file)
##############
# bound_cmd w/ samtools
        bound_cmd = 'hisat2 -x /home/bechara/Bioinformatics/References/Hisat2_reference/para_MIC/ies -U '+deepseq_folder+'/'+file+' | samtools view -Sh -F 4 - > temp/'+file+'_boundaries.sam'

##############
        # bound_cmd = 'hisat2 -x /home/bechara/Bioinformatics/References/Hisat2_reference/para_MIC/ies -U '+deepseq_folder+'/'+file+' -S temp/'+file+'_boundaries.sam'
        get_bound_reads = '/home/bechara/Bioinformatics/Scripts/boundarySmallRNAs.py -i temp/'+file+'_boundaries.sam -o temp/mapped/'+file+'.MICboundary_mapped.fastq'

        # filter_reads = '/home/bechara/Programmes/BBTools/bbmap/filterbyname.sh in=' +deepseq_folder+'/'+file+' out=temp/unmapped/'+file+'.MICboundary_unmapped.fastq names='+getNamesList('temp/'+file+'_boundariesMDSIESBoundaryReads.fastq')
##############
        vec_cmd = 'hisat2 --un temp/unmapped/' + file + '.vec_unmapped.fastq \
        --al temp/mapped/' + file + '.vec_mapped.fastq -x /home/bechara/Bioinformatics/References/Hisat2_reference/L4440/vector -U temp/unmapped/'+file+'.MICboundary_unmapped.fastq > temp/shiteVec.sam'

        mac_cmd = 'hisat2 --un temp/unmapped/' + file + '.mac_unmapped.fastq \
        --al temp/mapped/' + file + '.mac_mapped.fastq -x /home/bechara/Bioinformatics/References/Hisat2_reference/para_MAC/mac -U temp/unmapped/' + file +'.vec_unmapped.fastq > temp/shiteMac.sam'

        # mic_all = 'hisat2 -x /home/bechara/Bioinformatics/References/Hisat2_reference/para_MIC/ies -U '+deepseq_folder+'/'+file+' -S temp/'+file+'_MACIES_fullReads.sam'

        ies_cmd = 'hisat2 --un temp/unmapped/' + file + '.ies_unmapped.fastq \
        --al temp/mapped/' + file + '.ies_mapped.fastq -x /home/bechara/Bioinformatics/References/Hisat2_reference/para_MIC/ies -U temp/unmapped/' + file +'.mac_unmapped.fastq > temp/shiteMic.sam'

        oes_cmd = 'hisat2 --un temp/unmapped/' + file + '.oes_unmapped.fastq \
        --al temp/mapped/' + file + '.oes_mapped.fastq -x /home/bechara/Bioinformatics/References/Hisat2_reference/new_para_MIC/oes -U temp/unmapped/' + file +'.ies_unmapped.fastq > temp/shiteOes.sam'

        mito_cmd = 'hisat2 --un temp/unmapped/' + file + '.mito_unmapped.fastq \
        --al temp/mapped/' + file + '.mito_mapped.fastq -x /home/bechara/Bioinformatics/References/Hisat2_reference/mito/mito -U temp/unmapped/' + file +'.oes_unmapped.fastq > temp/shiteMito.sam'

        kleb_cmd = 'hisat2 --un temp/unmapped/' + file + '.kleb_unmapped.fastq \
        --al temp/mapped/' + file + '.kleb_mapped.fastq -x /home/bechara/Bioinformatics/References/Hisat2_reference/klebs/klebs -U temp/unmapped/' + file +'.mito_unmapped.fastq > temp/shiteKlebs.sam'
#         'hisat2 --un temp/unmappendaryFQ(deepseqFolder, file):
# #     os.system('/home/bechara/Bioinformatics/Scripts/boundarySmallRNAs.py -i '+deepseqFolder+'/'+file+' -o temp/mapped/BoundaryReads.fq -o2 OverlappedIESsInfo.tsvd/' + file + '.kleb_unmapped.fastq \
#         --al temp/mapped/' + file + '.kleb_mapped.fastq -x /home/bechara/Bioinformatics/References/Hisat2_reference/klebs/klebs -U temp/unmapped/' + file +'.mito_unmapped.fastq > temp/shiteKlebs.sam'
# ######  filter boudnar reads with script and then continue with cleaned file by seqtk subsample
        # boundcmd = 'hisat2 -x /home/bechara/Bioinformatics/References/Hisat2_reference/para_MIC/ies -U '
        print('bound cmd')
        os.system(bound_cmd)
        print('getbound cmd')
        os.system(get_bound_reads)
        print('filter cmd')
        # os.system(filter_reads)
        os.system('/home/bechara/Programmes/BBTools/bbmap/filterbyname.sh in=' +deepseq_folder+'/'+file+' out=temp/unmapped/'+file+'.MICboundary_unmapped.fastq names='+getNamesList('temp/mapped/'+file+'.MICboundary_mapped.fastq'))
        os.system(vec_cmd)
        os.system(mac_cmd)
        os.system(ies_cmd)
        os.system(oes_cmd)
        os.system(mito_cmd)
        os.system(kleb_cmd)

        cnt+=1
    print('Processing file '+str(cnt)+' of '+str(len(seq_files_list))+' files')

###
# def createBoundaryFQ(deepseqFolder, file):
#     os.system('/home/bechara/Bioinformatics/Scripts/boundarySmallRNAs.py -i '+deepseqFolder+'/'+file+' -o temp/mapped/BoundaryReads.fq -o2 OverlappedIESsInfo.tsv')
###

def count_reads(filename):
    key=filename.split('nt')[0].split('_')[-1]
    cnt=0
    for line in open(filename, 'r'):
        if line.startswith('@'):
            cnt+=1
    return (int(key), int(cnt))

def set_up_dicts(min_size, max_size):

    vec_dict={n:0 for n in range(min_size, max_size+1)}
    mac_dict={n:0 for n in range(min_size, max_size+1)}
###
    mdsiesBound_dict={n:0 for n in range(min_size, max_size+1)}
###
    ies_dict={n:0 for n in range(min_size, max_size+1)}
    oes_dict={n:0 for n in range(min_size, max_size+1)}
    mito_dict={n:0 for n in range(min_size, max_size+1)}
    klebs_dict={n:0 for n in range(min_size, max_size+1)}
    rest_unmapped={n:0 for n in range(min_size, max_size+1)}
    print('Min/Max: ', min_size, max_size)
    print('Vector dict length after initialisation', len(vec_dict))
    for file in os.listdir('temp/mapped/'):
        if file.endswith('fastq'):
            if 'vec_mapped' in file:
                key_cnt_reads=count_reads('temp/mapped/'+file)
                vec_dict[key_cnt_reads[0]] = key_cnt_reads[1]
            if 'mac_mapped' in file:
                key_cnt_reads=count_reads('temp/mapped/'+file)
                mac_dict[key_cnt_reads[0]] = key_cnt_reads[1]
###
            if 'MICboundary_mapped' in file:
                key_cnt_reads=count_reads('temp/mapped/'+file)
                mdsiesBound_dict[key_cnt_reads[0]] = key_cnt_reads[1]
###
            if 'ies_mapped' in file:
                key_cnt_reads=count_reads('temp/mapped/'+file)
                ies_dict[key_cnt_reads[0]] = key_cnt_reads[1]
            if 'oes_mapped' in file:
                key_cnt_reads=count_reads('temp/mapped/'+file)
                oes_dict[key_cnt_reads[0]] = key_cnt_reads[1]
            if 'mito_mapped' in file:
                key_cnt_reads=count_reads('temp/mapped/'+file)
                mito_dict[key_cnt_reads[0]] = key_cnt_reads[1]
            if 'kleb_mapped' in file:
                key_cnt_reads=count_reads('temp/mapped/'+file)
                klebs_dict[key_cnt_reads[0]] = key_cnt_reads[1]

    print('Vector dict length after read counting: ',len(vec_dict))

    for file in os.listdir('temp/unmapped/'):
        if file.endswith('fastq'):
            if 'kleb_unmapped' in file:
                rest_unmapped[count_reads('temp/unmapped/'+file)[0]] = count_reads('temp/unmapped/'+file)[1]
    reads_total = 0
    reads_total += sum(mdsiesBound_dict.values())
    reads_total += sum(vec_dict.values())
    # print('vector: ', reads_total)
    reads_total += sum(mac_dict.values())
    reads_total += sum(ies_dict.values())
    reads_total += sum(oes_dict.values())
    reads_total += sum(mito_dict.values())
    reads_total += sum(klebs_dict.values())
    reads_total += sum(rest_unmapped.values())

    for k,v in vec_dict.items():
        try:
            vec_dict[k] = v / reads_total
        except:
            vec_dict[k] = 0

    for k,v in mac_dict.items():
        try:
            mac_dict[k] = v / reads_total
        except:
            mac_dict[k] = 0
        #mac_dict[k] = v / reads_total
    for k,v in ies_dict.items():
        try:
            ies_dict[k] = v / reads_total
        except:
            ies_dict[k] = 0 #ies_dict[k] = v / reads_total
    #####
    for k,v in mdsiesBound_dict.items():
        try:
            mdsiesBound_dict[k] = v / reads_total
        except:
            mdsiesBound_dict[k] = 0
    #####
    for k,v in oes_dict.items():
        try:
            oes_dict[k] = v / reads_total
        except:
            oes_dict[k] = 0#oes_dict[k] = v / reads_total
    for k,v in klebs_dict.items():
        try:
            klebs_dict[k] = v / reads_total
        except:
            klebs_dict[k] = 0#klebs_dict[k] = v / reads_total
    for k,v in mito_dict.items():
        try:
            mito_dict[k] = v / reads_total
        except:
            mito_dict[k] = 0 #mito_dict[k] = v / reads_total
    for k,v in rest_unmapped.items():
        try:
            rest_unmapped[k] = v / reads_total
        except:
            rest_unmapped[k] = 0

    print('Vector dict length after normalisation: ', len(vec_dict))

    return  vec_dict, mac_dict, ies_dict, mdsiesBound_dict, oes_dict, klebs_dict, mito_dict, rest_unmapped

def make_plot(min_size, max_size,  title,deepseq_folder):

    vec_dict, mac_dict, ies_dict, mdsiesBound_dict, oes_dict, klebs_dict, mito_dict, rest_unmapped = set_up_dicts(min_size, max_size)
    with open(os.getcwd()+'/'+'_'.join(title.split())+'.log', 'w') as w:
        w.write(f'# Used directory: {deepseq_folder}\n#\n#\n\n')
        w.write('sRNA size\tVector\tMDS\tIES\tMDS/IES boundaries\tOES\tKlebsiella\tMitochondrium\tUnmapped\n')
        for i in range(min_size,max_size+1):
            w.write(f'{i}\t{vec_dict[i]}\t{mac_dict[i]}\t{ies_dict[i]}\t{mdsiesBound_dict[i]}\t{oes_dict[i]}\t{klebs_dict[i]}\t{mito_dict[i]}\t{rest_unmapped[i]}\n')
    x_val = np.array(list(range(min_size, max_size+1)))
    print('X values: ', len(x_val))
    vec_val = np.array(list(vec_dict.values()))
    print('Vector keys: ', list(vec_dict.keys()), '\nVector keys list length: ', len(list(vec_dict.keys())))
    print('Vector values: ', list(vec_val), '\nVector key list length: ', len(list(vec_val)))
    mac_val = np.array(list(mac_dict.values()))
    ies_val = np.array(list(ies_dict.values()))
    mdsies_val = np.array(list(mdsiesBound_dict.values()))
    oes_val = np.array(list(oes_dict.values()))
    mito_val = np.array(list(mito_dict.values()))
    klebs_val = np.array(list(klebs_dict.values()))
    rest_val = np.array(list(rest_unmapped.values()))
###########
#    x_val = list(range(min_size, max_size+1))
#    print('X values: ', len(x_val))
#    vec_val = list(vec_dict.values())
#    print('Vector keys: ', list(vec_dict.keys()), '\nVector keys list length: ', len(list(vec_dict.keys())))
#    print('Vector values: ', list(vec_val), '\nVector key list length: ', len(list(vec_val)))
#    mac_val = list(mac_dict.values())
#    ies_val = list(ies_dict.values())
#    oes_val = list(oes_dict.values())
#    mito_val = list(mito_dict.values())
#    klebs_val = list(klebs_dict.values())
#    rest_val = list(rest_unmapped.values())
#############
#### label in plt.bla
    #ind=np.arange(len(x_val))
    plt.figure(figsize=[10,8])
    plt.bar( x_val, vec_val, color='darkviolet')
    plt.bar( x_val, mac_val, color='forestgreen', bottom = vec_val)
    plt.bar( x_val, ies_val, color='red',label='IES', bottom = mac_val+ vec_val)
    #plt.bar(x_val, ies_val, color='red', bottom = mac_val+ vec_val)
    plt.bar( x_val, mdsies_val, color='pink', bottom = ies_val + mac_val + vec_val)
    plt.bar( x_val, oes_val, color='dodgerblue', bottom = mdsies_val+ies_val+ mac_val+vec_val)
    plt.bar( x_val, mito_val, color='orange', bottom = oes_val+mdsies_val+ ies_val+ mac_val+ vec_val)
    plt.bar( x_val, klebs_val, color='black', bottom = mito_val+ oes_val+mdsies_val+ ies_val+ mac_val+ vec_val)
    plt.bar( x_val, rest_val, color='grey', bottom = klebs_val+ mito_val+ oes_val+mdsies_val+ ies_val+ mac_val+ vec_val)

    plt.xlabel('sRNA length [bp]' , fontsize=16)
    plt.ylabel('Reads (normalised to total reads)',fontsize=16)
    lgd=('Vector','MDS','IES','MDS/IES boundaries','OES','Mitochondrial','Klebsiella','Unmapped')
    #plt.legend((p1[0],p2[0],p3[0],p4[0],p5[0],p6[0].p7[0]),('Vector','MDS','IES','OES','Mitochondrial','Klebsiella','Unmapped')) ### dont understand how this works, something with arrayed handles
    plt.legend(lgd, loc='best')#,'MDS')#,'IES','OES','Mitochondrial','klebsiella','Unmapped')
    plt.title(title)
    #plt.show()
    #### saving figures directly dosnt work right know
    plt.savefig(os.getcwd()+'/'+'_'.join(title.split())+'.png', bbox_inches='tight')#'/home/debian9/host/analysis/title+'.png')

def main():
    if len(sys.argv[1:]) != 4:
        print("something didn't work, maybe you messed up the arguments?\n\
        type in following arguments: \nminimum size, maximum size, the folder\
         where your sequencing data is and the title you want to have on your diagramm\
         (also the name of the folder it creates to store the picture)\nin this order!")
        sys.exit()

    min_size = sys.argv[1]
    min_size = int(min_size)
    max_size = sys.argv[2]
    max_size = int(max_size)
    deepseq_folder = sys.argv[3]
    deepseq_folder = str(deepseq_folder)
    title = sys.argv[4]
    title = str(title)
    check_folder(title)    # !!!!!!!!!!!!!!!!!!!!
    # check_hisat_refs()
    mapSeq(deepseq_folder, int(min_size), int(max_size)) # !!!!!!!!!!!!!!!!!!!!!!
    # seq_files(int(min_size), int(max_size), deepseq_folder)
    # mapSeq(deepseq_folder, seq_files)
    make_plot(int(min_size), int(max_size), title, deepseq_folder)
    os.system('rm -rv temp')

main()
