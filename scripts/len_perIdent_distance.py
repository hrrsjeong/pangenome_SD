import sys,glob,os

flist = glob.glob("SD_freq_stats/*.stats.bed")
dd = {}
dd_noACRO = {}
for sfile in flist:
    if "above" in sfile:continue
    AF_idx = sfile.split('/')[-1].split('.')[-3]
    SD_type = sfile.split('/')[-1].split('.')[2]
    with open(sfile,'r') as fp:
        for line in fp:
            line_temp = line.strip().split('\t')
            
            sample_info,SD_type_tmp,SD1,SD2,perIdent,tmp = line_temp[3].split('__')
            SD_start,SD_end = [int(x) for x in SD1.split(':')[-1].split('-')]
            SD_start2,SD_end2 = [int(x) for x in SD2.split(':')[-1].split('-')]
            SD_len = str(SD_end - SD_start +1)
            distance_knownSD = line_temp[-1]
            if SD_type == "intra":
                SD_tmp_start,SD_tmp_end = min(SD_start,SD_start2), min(SD_end,SD_end2)
                SD_tmp_start2,SD_tmp_end2 = max(SD_start,SD_start2), max(SD_end,SD_end2)
                SD_flag = SD_tmp_end - SD_tmp_start2
                if SD_flag > 0:
                    if float(SD_flag) / float(SD_len) > 0.2:continue
                
            dd.setdefault(SD_type,{})
            dd[SD_type].setdefault(AF_idx,[])
            dd[SD_type][AF_idx].append((SD_len,perIdent,distance_knownSD))
            if "chr13" in line:continue
            if "chr14" in line:continue
            if "chr15" in line:continue
            if "chr21" in line:continue
            if "chr22" in line:continue
            dd_noACRO.setdefault(SD_type,{})
            dd_noACRO[SD_type].setdefault(AF_idx,[])
            dd_noACRO[SD_type][AF_idx].append((SD_len,perIdent,distance_knownSD))
with open("SD_summary_len_perIdent.txt",'w') as fout:
    fout.write("SD_type\tAF_idx\tSD_length\tSD_perIdent\tdistance_to_knownSD\n")
    for SD_type in dd:
        for AF_idx in dd[SD_type]:
            for i in dd[SD_type][AF_idx]:
                fout.write(SD_type+'\t'+AF_idx+'\t'+'\t'.join(i)+'\n')
    
with open("SD_summary_len_perIdent_noACRO.txt",'w') as fout:
    fout.write("SD_type\tAF_idx\tSD_length\tSD_perIdent\tdistance_to_knownSD\n")
    for SD_type in dd_noACRO:
        for AF_idx in dd_noACRO[SD_type]:
            for i in dd_noACRO[SD_type][AF_idx]:
                fout.write(SD_type+'\t'+AF_idx+'\t'+'\t'.join(i)+'\n')
    



'''head SD_freq_stats/combined.all.inter.hap.cov.AF2-5.stats.bed
chr1	12381475	12383446	HG00512_hap1_nontrio_EAS__inter__chr1_RagTag:12688730-12690773__chr5_RagTag:63172018-63174292__0.9123__chm13_0.93001
chr1	12381475	12383446	HG01952_hap1_trio_AMR__inter__chr1_RagTag:12559564-12561607__chr5_RagTag:67234743-67237017__0.9123__chm13_0.93001
chr1	12381475	12383446	HG03732_hap2_nontrio_SAS__inter__chr1_RagTag:12237067-12239110__chr5_RagTag:67532012-67534286__0.9123__chm13_0.93001
chr1	12381478	12383446	HG02559_hap2_trio_AFR__inter__chr1_RagTag:12282473-12284459__chr5_RagTag:66816215-66818489__0.9093__chm13_0.99395
'''
