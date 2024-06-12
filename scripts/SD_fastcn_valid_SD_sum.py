import sys,glob,os
dd_contig = {}
dd_ragtag = {}
dd_ragtag_acro = {}
dd_ragtag_autosome = {}
dd_sample = {}
fout = open("SD_summary.autosome.fastcn_valid.using_chm13_acro.txt",'w')
fout.write("sample\tdataset\tsex\thap\ttrio\tpop\tasm_size\tRagTag_size\tRagTag_size_excluding_acro_p\tRagTag_size_autosome\tRagTag_size_autosome_excluding_acro_p\tSD_intra_sum\tSD_inter_sum\n")
sample_info = '../../asm-to-reference-alignment/sample_info.txt'
with open(sample_info,'r') as fp:
    for line in fp:
        line_temp = line.strip().split('\t')
        sample = line_temp[1]
        sex = line_temp[2]
        pop = line_temp[4]
        superpop = line_temp[5]
        trio = line_temp[3]
        if trio == "1":
            trio_idx = "yes"
        else:
            trio_idx = "no"
        dataset = line_temp[0]
        dd_sample[sample] = [dataset,trio_idx,superpop,sex]
'''
Source  Sample  Sex     Trio data available     Population      Superpopulation Family ID       Relation
HPRC    HG01891 Female  1       ACB     AFR     BB05    child
HPRC    HG02109 Female  1       ACB     AFR     BB18    child
'''

flist = glob.glob("../../asm-to-reference-alignment/fasta_gz/*.fai")
for sfile in flist:
    sample_id = sfile.split('/')[-1].split('.')[0]
    hap = sfile.split('/')[-1].split('.')[1]
    dd_contig.setdefault(sample_id,{})
    dd_contig[sample_id].setdefault(hap,{})
    with open(sfile,'r') as fp:
        for line in fp:
            line_temp = line.strip().split("\t")
            dd_contig[sample_id][hap][line_temp[0]] = int(line_temp[1])

flist = glob.glob("../../repeatmask*/fasta/*.RagTag.fasta.fai")
for sfile in flist:
    sample_id = sfile.split('/')[-1].split('.')[0].split('_')[0]
    hap_tmp = sfile.split('/')[-1].split('.')[0].split('_')[1]
    if hap_tmp == "pat":
        hap = "hap1"
    elif hap_tmp == "mat":
        hap = "hap2"
    else:
        hap = hap_tmp
    dd_ragtag.setdefault(sample_id,{})
    dd_ragtag_autosome.setdefault(sample_id,{})
    with open(sfile,'r') as fp:
        for line in fp:
            line_temp = line.strip().split("\t")
            dd_ragtag[sample_id].setdefault(hap,0)
            dd_ragtag[sample_id][hap] += int(line_temp[1])
            if not ("chrX" in line or "chrY" in line or "chrM" in line):
                dd_ragtag_autosome[sample_id].setdefault(hap,0)
                dd_ragtag_autosome[sample_id][hap] += int(line_temp[1])
CHM13 = 0
'''
with open('/net/eichler/vol26/projects/chm13_t2t/nobackups/assemblies/chm13_v1.1_plus38Y_par_masked.fasta.fai','r') as fp:
    for line in fp:
        if ("chrX" in line) or ("chrY" in line) or ("chrM" in line):continue
        line_temp = line.strip().split('\t')
        CHM13 += int(line_temp[1])
'''

flist = glob.glob("acro_p_bed/merged/*.*.acro.ragtag.bed")
for sfile in flist:
    sample_id = sfile.split('/')[-1].split('.')[0]
    hap = sfile.split('/')[-1].split('.')[1].replace("h","hap")
    dd_ragtag_acro.setdefault(sample_id,{})
    with open(sfile,'r') as fp:
        for line in fp:
            line_temp = line.strip().split("\t")
            dd_ragtag_acro[sample_id].setdefault(hap,0)
            dd_ragtag_acro[sample_id][hap] += int(line_temp[2]) - int(line_temp[1])


#flist = glob.glob("./aligned_sum/*.*.contig.bed")
#flist = glob.glob("./fastCN_filtered_ragtag/*.SDs.autosome.intra.bed")#./bed/SDs_merged/*.SDs.autosome.intra.merged.bed")
flist = glob.glob("./acro_p_filtered_ragtag/*.SDs.autosome.intra.merged.bed")#./bed/SDs_merged/*.SDs.autosome.intra.merged.bed")
for sfile in flist:
    SD_intra_region_sum = 0
    SD_inter_region_sum = 0
    RagTag_sum = 0
    RagTag_autosome_sum = 0
    sample_idx = sfile.split('/')[-1].split('.')[0].split('_')
    sample_id = sample_idx[0]
    hap_tmp = sample_idx[1] #"hap"+sample_idx[1]
    if hap_tmp == "pat":
        hap = "hap1"
    elif hap_tmp == "mat":
        hap = "hap2"
    else:
        hap = hap_tmp

    if sample_id not in dd_ragtag_acro:continue
    dataset_info,trio_idx,pop_info,sex = dd_sample[sample_id]
    with open(sfile,'r') as fp:
        for line in fp:
            line_temp = line.strip().split('\t')
            SD_intra_region_sum += int(line_temp[2]) - int(line_temp[1])
    with open(sfile.replace('.intra.','.inter.'),'r') as fp:
        for line in fp:
            line_temp = line.strip().split('\t')
            SD_inter_region_sum += int(line_temp[2]) - int(line_temp[1])
    contig_sum = sum([dd_contig[sample_id][hap][x] for x in dd_contig[sample_id][hap]])
    RagTag_sum = dd_ragtag[sample_id][hap]
    RagTag_sum_excluding_acro = dd_ragtag[sample_id][hap] - dd_ragtag_acro[sample_id][hap]

    RagTag_autosome_sum = dd_ragtag_autosome[sample_id][hap]
    RagTag_autosome_sum_excluding_acro = dd_ragtag_autosome[sample_id][hap] - dd_ragtag_acro[sample_id][hap]
    fout.write(f'{sample_id}\t{dataset_info}\t{sex}\t{hap}\t{trio_idx}\t{pop_info}\t{contig_sum}\t{RagTag_sum}\t{RagTag_sum_excluding_acro}\t{RagTag_autosome_sum}\t{RagTag_autosome_sum_excluding_acro}\t{SD_intra_region_sum}\t{SD_inter_region_sum}\n')
fout.close()


#GM21487.1.contig.bed
