import sys,glob,os

dict_query_cnt_chm13 = {}
dict_query_cnt_asm = {}
q_cov_threshold = 0.9
q_cov_threshold_long = 0.8
#q_cov_threshold2 = 0.6
#min_pass_len = 10000
dir_path = "./aligned/bed"
out_path = "./chm13_projection/individual/raw"
#out_path_multimap = "./chm13_projection/individual/multi_mapped"
os.system(f"mkdir -p {out_path}")


def parse_rbed(fname,sname,dd_query_cnt,query_cov_threshold):
    with open(fname,'r') as fp:
        if "intra" in fname:
            SD_idx = "intra"
        elif "inter" in fname:
            SD_idx = "inter"
        elif "all" in fname:
            SD_idx = "all"
        else:
            print ("SD type error",fname)
            sys.exit()

        dd_query_cnt.setdefault(SD_idx,{})
        for line in fp:
            line_temp = line.strip().split('\t')
            if line.startswith("#"):
                title = line_temp
                query_start = title.index("query_start")
                query_end = title.index("query_end")
                query_length = title.index("query_length")
                perID_by_matches = title.index("perID_by_matches")
                query_name = title.index("query_name")
                ref_chr_idx = title.index("#reference_name")
            else:
                if ("chrX" in line) or ("chrY" in line) or ("chrM" in line):continue
                qlen = int(line_temp[query_length])
                q_end = int(line_temp[query_end])
                q_start = int(line_temp[query_start])
                perID = float(line_temp[perID_by_matches])
                ref_chr = line_temp[ref_chr_idx]
                qname = line_temp[query_name]
                q_aligned_len = q_end - q_start + 1
                mapping_chrom = line_temp[5].split('_')[0]
                mapped_chrom  = line_temp[0].split('_')[0]
                if (qlen * q_cov_threshold < q_aligned_len) or (qlen * q_cov_threshold_long < q_aligned_len and q_aligned_len > 10000) or (q_aligned_len > 50000):
                    if perID > 90 and q_aligned_len > 1000:
                        qchrs = [x.split(':')[0].split("_RagTag")[0] for x in qname.split("__")[0:-1]]
                        if ref_chr not in qchrs:continue
                        #if SD_idx == "intra":
                        #    if mapping_chrom != mapped_chrom:
                        #        continue
                        #elif SD_idx == "inter":
                        #    if mapping_chrom == mapped_chrom:
                        #        continue
                        #ref_name = line_temp[0]
                        ref_start,ref_end = line_temp[1],line_temp[2]

                        dd_query_cnt[SD_idx].setdefault(sname,{})
                        dd_query_cnt[SD_idx][sname].setdefault(qname,[])
                        dd_query_cnt[SD_idx][sname][qname].append([ref_chr,ref_start,ref_end,str(round(perID/100,5))])
    return dd_query_cnt
'''
#reference_name	reference_start	reference_end	reference_length	strand	query_name	query_start	query_end	query_length	perID_by_matches	perID_by_events	perID_by_all	matches	mismatches	deletion_events	insertion_events	deletions	insertions
chr10_RagTag	13564	19034	131958902	+	chr10_RagTag:13564-19034	0	5470	5470	100	100	100	5470	0	0	0	0

#reference_name	reference_start	reference_end	reference_length	strand	query_name	query_start	query_end	query_lengthperID_by_matches	perID_by_events	perID_by_all	matches	mismatches	deletion_events	insertion_events	deletions	insertions
chr10	0	81465	134758134	+	chr10_RagTag:0-85329__chr18_RagTag:18-85464__0.9725	4033	85329	85329	99.62536	99.46846	99.05043	80841	304	65	63	320	151
chr18	163478	243123	80542538	+	chr10_RagTag:0-85329__chr18_RagTag:18-85464__0.9725	6083	85329	85329	97.43353	97.04155	91.59808	75017	1976	163	148	2652	2253

'''
sample_info = '/net/eichler/vol27/projects/hprc/nobackups/analysis/sedef/asm-to-reference-alignment/sample_info.txt'
dd_sample = {}

'''
Source	Sample	Sex	Trio data available	Population	Superpopulation	Family ID	Relation
HPRC	HG01891	Female	1	ACB	AFR	BB05	child
'''
with open(sample_info,'r') as fp:
    fp.readline()
    for line in fp:
        line_temp = line.strip().split('\t')
        sample = line_temp[1]
        if line_temp[3] == "1":
            trio_idx = "trio"
        else:
            trio_idx = "nontrio"
        dd_sample[sample] = [line_temp[0],line_temp[2][0],trio_idx,line_temp[4],line_temp[5]]



ref_flist = glob.glob(f"{dir_path}/chm13/*.all.chm13.bed")
for sfile in ref_flist:
    sample_id = sfile.split('/')[-1].split('.')[0]
    dict_query_cnt_chm13 = parse_rbed(sfile,sample_id,dict_query_cnt_chm13,q_cov_threshold)
    #dict_query_cnt_asm = parse_rbed(sfile.replace('chm13','asm'),sample_id,dict_query_cnt_asm,q_cov_threshold)

#os.system(f'mkdir -p {out_path_multimap}/pass')
#os.system(f'mkdir -p {out_path_multimap}/fail')

for SD_type in dict_query_cnt_chm13:
    for sname in dict_query_cnt_chm13[SD_type]:
        sample_id,hap = sname.split('_')
        if sample_id not in dd_sample:continue
        dataset,sex,trio,pop,superpop = dd_sample[sample_id]
        fout = open(f'{out_path}/{sample_id}.{hap}.{dataset}.{sex}.{trio}.{pop}.{superpop}.SD.{SD_type}.chm13.bed','w')
        #fout = open(f'{out_path}/{sname}.SD.{SD_type}.chm13.bed','w')
        #fout_multimap_pass = open(f'{out_path_multimap}/pass/{sname}.SD.{SD_type}.chm13.multimapped_pass.bed','w')
        #fout_multimap_fail = open(f'{out_path_multimap}/fail/{sname}.SD.{SD_type}.chm13.multimapped_fail.bed','w')
        for query in dict_query_cnt_chm13[SD_type][sname]:
            for matched_region_info in dict_query_cnt_chm13[SD_type][sname][query]:
                fout.write('\t'.join(matched_region_info[0:-1])+f'\t{sname}_{trio}_{superpop}__{SD_type}__{query}__chm13_{matched_region_info[-1]}\n')
            #if len(dict_query_cnt_chm13[SD_type][sname][query]) > 1:
            #    for matched_region_info in dict_query_cnt_chm13[SD_type][sname][query]:
            #        fout_multimap_pass.write('\t'.join(matched_region_info[0:-1])+f'\t{sname}__{matched_region_info[-1]}__{SD_type}__'+query+'\n')
            #else:
            #    for matched_region_info in dict_query_cnt_chm13[SD_type][sname][query]:
            #        fout_multimap_fail.write('\t'.join(matched_region_info[0:-1])+f'\t{sname}__{matched_region_info[-1]}__{SD_type}__'+query+'\n')
        fout.close()
        #fout_multimap_fail.close()
        #fout_multimap_pass.close()




