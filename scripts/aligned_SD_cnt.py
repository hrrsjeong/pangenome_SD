import sys,glob

flist = glob.glob("../../SD/sedef/HPRC_chrom/SD_align/bed/SDs_merged_reciprocal90/HG002_pat.SDs.autosome.intra.merged90.bed")

sample_info = '/net/eichler/vol27/projects/hprc/nobackups/analysis/sedef/asm-to-reference-alignment/sample_info.txt'
flist = glob.glob(f"/net/eichler/vol27/projects/hprc/nobackups/analysis/sedef/SD/sedef/*/SD_align/bed/SDs_merged_reciprocal90/*_*.SDs.autosome.in*.merged90.bed")
#os.system(f"mkdir -p aligned/paf/chm13 && mkdir -p aligned/paf/asm")

print ('SD_type\tsample_id\thap\tnum_SD\tnum_alignable_SD\tnum_filtered_SD\tnum_alignable_SD_single_count\tnum_filtered_SD_single_count')
for sfile in flist:
    #print (sfile)
    if "T2T" in sfile:continue
    info_tmp = sfile.split('/')[-1].split('.')
    sample_id,hap_tmp = info_tmp[0].split('_')
    hap = hap_tmp.replace("mat","hap2").replace("pat","hap1")
    SD_bed = sfile
    #print (SD_bed)
    SD_type = sfile.split('/')[-1].split('.')[-3]
    aligned_SD_file = f'./tmp_z2000_aligned/bed/chm13/{sample_id}_{hap}.SD.{SD_type}.chm13.bed'
    filtered_SD_file = f'./chm13_projection/individual/raw/{sample_id}_{hap}.SD.{SD_type}.chm13.bed'
    SD_list = set()
    aligned_SD_list = {}
    filtered_SD_list = {}
    #fout_inter = open(f'SD_count_test/{sample_id}_{hap}.SD.{SD_type}.single_match.txt','w')
    #fout_intra = open(f'SD_count_test/{sample_id}_{hap}.SD.{SD_type}.single_match.txt','w')
    #fout_inter_filtered = open(f'SD_count_test/{sample_id}_{hap}.SD.{SD_type}.single_match.filtered.txt','w')
    #fout_intra_filtered = open(f'SD_count_test/{sample_id}_{hap}.SD.{SD_type}.single_match.filtered.txt','w')
    with open(SD_bed,'r') as fp:
        for line in fp:
            line_temp = line.strip().split('\t')
            SD_list.add(line_temp[3].split('__')[-1])
    with open(aligned_SD_file,'r') as fp:
        for line in fp:
            line_temp = line.strip().split('\t')
            if line_temp[5] not in SD_list:continue
            aligned_SD_list.setdefault(line_temp[5],0)
            aligned_SD_list[line_temp[5]] +=1
    with open(filtered_SD_file,'r') as fp:
        for line in fp:
            line_temp = line.strip().split('\t') #HG002_hap1__97.47945__inter__chr10_RagTag:0-85512
            if line_temp[3].split('__')[-1] not in SD_list:continue
            filtered_SD_list.setdefault(line_temp[3].split('__')[-1],0)
            filtered_SD_list[line_temp[3].split('__')[-1]] +=1
    aligned_SD_list_single = [x for x in aligned_SD_list if aligned_SD_list[x] == 1]
    filtered_SD_list_single = [x for x in filtered_SD_list if filtered_SD_list[x] == 1]
    '''
    with open(aligned_SD_file,'r') as fp:
        for line in fp:
            line_temp = line.strip().split('\t')
            if line_temp[5] in aligned_SD_list_single:
                if SD_type == "intra":
                    fout_intra.write(line)
                if SD_type == "inter":
                    fout_inter.write(line)
            elif line_temp[5] in filtered_SD_list_single:
                if SD_type == "intra":
                    fout_intra_filtered.write(line)
                if SD_type == "inter":
                    fout_inter_filtered.write(line)
    '''
    systemstr = f'{SD_type}\t{sample_id}\t{hap}\t{len(SD_list)}\t{len(aligned_SD_list)}\t{len(filtered_SD_list)}\t{len(aligned_SD_list_single)}\t{len(filtered_SD_list_single)}'
    print (systemstr)

    #fout_inter.close()
    #fout_intra.close()
    #fout_inter_filtered.close()
    #fout_intra_filtered.close()


