import sys,glob,os

dict_query_cnt = {}
q_cov_threshold = 0.6
perID_threshold = 90
cov_idx = int(q_cov_threshold*100)
dir_path = "./aligned/bed"
dir_path_rm = "./aligned_rm/bed"
#aligned/bed/asm_contig/NA21309_hap1.chm13_gene.protein_coding.asm.contig.bed
multicopy_set = {}
dd_sample = {}

dd_ragtag_SD_completed = {}
flist = glob.glob("/net/eichler/vol28/projects/hprc/nobackups/analysis/sedef/SD/sedef/*/SD_align/fasta/SDs_merged_reciprocal90/*_*.SDs.autosome.intra.merged90.fasta")
for sfile in flist:
    #print (sfile)
    if "T2T" in sfile:continue
    info_tmp = sfile.split('/')[-1].split('.')
    sample_id,hap_tmp = info_tmp[0].split('_')
    hap = hap_tmp.replace("mat","hap2").replace("pat","hap1")
    dd_ragtag_SD_completed.setdefault(sample_id,{})
    dd_ragtag_SD_completed[sample_id][hap] = sfile
    #ragtag_fasta = dd_ragtag_fasta[sample_id][hap]
    #SD_bed = sfile
    #SD_type = sfile.split('/')[-1].split('.')[-3]

sample_info = '../../asm-to-reference-alignment/sample_info.txt'
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

def parse_rbed(fname,fname_rm,sname,dd_query_cnt,query_cov_threshold,query_set,fout_tmp,filtered_bed_file):
    gene_mapping_rm = {}
    gene_mapping_overlap = {}
    filtered_bed_fout = open(filtered_bed_file,'w')
    filtered_bed_fout.write("#reference_name\treference_start\treference_end\treference_length\tstrand\tquery_name\tquery_start\tquery_end\tquery_length\tperID_by_matches\tperID_by_events\tperID_by_all\tmatches\tmismatches\tdeletion_events\tinsertion_events\tdeletions\tinsertions\n")
    with open(fname_rm,'r') as fp:
        for line in fp:
            line_temp = line.strip().split('\t')
            title_temp = "#reference_name\treference_start\treference_end\treference_length\tstrand\tquery_name\tquery_start\tquery_end\tquery_length\tperID_by_matches\tperID_by_events\tperID_by_all\tmatches\tmismatches\tdeletion_events\tinsertion_events\tdeletions\tinsertions"
            title = title_temp.split('\t')
            query_start = title.index("query_start")
            query_end = title.index("query_end")
            query_length = title.index("query_length")
            perID_by_matches = title.index("perID_by_matches")
            query_name = title.index("query_name")
            if line.startswith("#"):
                continue
            else:
                if ("chrX" in line) or ("chrY" in line) or ("chrM" in line):continue
                #if not line_temp[0].startswith("chr"):continue
                qlen = int(line_temp[query_length])
                #if qlen < 1000:continue
                q_end = int(line_temp[query_end])
                q_start = int(line_temp[query_start])
                perID = float(line_temp[perID_by_matches])
                qname = line_temp[query_name]
                gene_mapping_rm.setdefault(qname,{})
                gene_mapping_rm[qname].setdefault(line_temp[0],[])
                gene_mapping_rm[qname][line_temp[0]].append([int(line_temp[1]),int(line_temp[2])])

    with open(fname,'r') as fp:
        for line in fp:
            line_temp = line.strip().split('\t')
            title_temp = "#reference_name\treference_start\treference_end\treference_length\tstrand\tquery_name\tquery_start\tquery_end\tquery_length\tperID_by_matches\tperID_by_events\tperID_by_all\tmatches\tmismatches\tdeletion_events\tinsertion_events\tdeletions\tinsertions"
            title = title_temp.split('\t')
            
            query_start = title.index("query_start")
            query_end = title.index("query_end")
            query_length = title.index("query_length")
            perID_by_matches = title.index("perID_by_matches")
            query_name = title.index("query_name")
            if line.startswith("#"):
                continue
            else:
                if ("chrX" in line) or ("chrY" in line) or ("chrM" in line):continue
                #if not line_temp[0].startswith("chr"):continue
                qlen = int(line_temp[query_length])
                #if qlen < 1000:continue
                q_end = int(line_temp[query_end])
                q_start = int(line_temp[query_start])
                perID = float(line_temp[perID_by_matches])
                qname = line_temp[query_name]
                if qlen * query_cov_threshold < (q_end - q_start) + 1:
                    if perID > perID_threshold:

                        ######check regions are also mapped after repeatmasking
                        
                        flag_rm = 0
                        if qname not in gene_mapping_rm:
                            #print (qname," is not mapped after RM")
                            continue
                        if line_temp[0] not in gene_mapping_rm[qname]:
                            #print (sname)
                            #print (line)
                            #print (gene_mapping_rm[qname])
                            #sys.exit()
                            continue
                        mapped_region_list = gene_mapping_rm[qname][line_temp[0]]
                        for mapped_region in mapped_region_list:
                            if max(int(line_temp[1]),mapped_region[0]) < min(int(line_temp[2]),mapped_region[1]):
                                flag_rm = 1
                                break
                        if flag_rm == 0:continue
                        

                        ####check the region is mapping the same region multiple times####
                        flag_overlap = 0

                        if qname not in gene_mapping_overlap:
                            flag_overlap = 1
                        elif line_temp[0] not in gene_mapping_overlap[qname]:
                            flag_overlap = 1
                        else:
                            tmp_flag = 1
                            for k in gene_mapping_overlap[qname][line_temp[0]]:
                                if min(k[1],int(line_temp[2])) > max(k[0],int(line_temp[1])):
                                    tmp_flag = 0
                                    break
                            if tmp_flag == 1:
                                flag_overlap = 1
                        if flag_overlap == 0:continue

                        gene_mapping_overlap.setdefault(qname,{})
                        gene_mapping_overlap[qname].setdefault(line_temp[0],[])
                        gene_mapping_overlap[qname][line_temp[0]].append([int(line_temp[1]),int(line_temp[2])])
                        
                        dd_query_cnt.setdefault(sname,{})
                        dd_query_cnt[sname].setdefault(qname,0)
                        dd_query_cnt[sname][qname] += 1
                        if sname == "chm13":
                            fout_tmp.write('\t'.join(line_temp[0:3])+'\t'+qname+'\t0\t'+'\t'.join(line_temp[4:])+'\n')
                            continue
                        query_set.setdefault(qname,{})
                        query_set[qname].setdefault(sname,0)
                        query_set[qname][sname] +=1
                        filtered_bed_fout.write(line)
    filtered_bed_fout.close()
    #sys.exit()
    return dd_query_cnt,query_set

'''
#reference_name	reference_start	reference_end	reference_length	strand	query_name	query_start	query_end	query_length	perID_by_matches	perID_by_events	perID_by_all	matches	mismatches	deletion_events	insertion_events	deletions	insertions
chr1	20528	37628	248387328	-	LOC112268260::chr1:20528-37628(-)	0	17100	17100	100	100	100	17100	0	0	0	0	0
chr5	181999855	182017986	182045439	+	LOC112268260::chr1:20528-37628(-)	5	17100	17100	99.37234	99.2959	91.514885	16782	106	11	2	1243	207

'''

fout_chm13_filtered = open(f"tmp_chm13_gene.chm13.filtered.perID_{perID_threshold}.cov{cov_idx}.bed",'w')
#q_flist = glob.glob(f"{dir_path}/asm_contig/*.asm.contig.bed")
q_flist = glob.glob(f"{dir_path}/asm_ragtag_multiexonic_orf200bp/*.chm13_gene.protein_coding.asm.ragtag.bed")
#q_flist = glob.glob(f"{dir_path}/asm_ragtag_rm_repeat/*.asm.ragtag.bed")
#q_flist = glob.glob(f"")
for sfile in q_flist:
    sample_id = sfile.split('/')[-1].split('.')[0]
    sample_idx = sample_id.split('_')[0]
    sfile_rm = f"{dir_path_rm}/asm_ragtag/{sample_id}.chm13_gene.protein_coding.asm.ragtag.bed"
    #HG00673_hap2.chm13_gene.protein_coding.asm.ragtag.bed
    if sample_idx not in dd_sample:
        print ("not in the sample list: ",sfile)
        continue
    if sample_idx not in dd_ragtag_SD_completed:
        print ("not in the SD completed list: ",sfile)
        continue
    if len(dd_ragtag_SD_completed[sample_idx]) != 2:
        print ("only one haplotype has been process..", sfile)
        continue
    filtered_bed = f"{dir_path}/asm_ragtag_filtered/{sample_id}.chm13_gene.protein_coding.asm.filtered.ragtag.bed"
    dict_query_cnt,multicopy_set = parse_rbed(sfile,sfile_rm,sample_id,dict_query_cnt,q_cov_threshold,multicopy_set,fout_chm13_filtered,filtered_bed)

#only one file
ref_flist = glob.glob(f"{dir_path}/chm13_multiexonic_orf200bp/*.chm13_gene.protein_coding.chm13.bed")
#ref_flist = glob.glob(f"{dir_path}/chm13_rm_repeat/*.chm13.bed")
for sfile in ref_flist:
    sample_id = sfile.split('/')[-1].split('.')[0]
    sfile_rm = f"{dir_path_rm}/chm13/{sample_id}.chm13_gene.protein_coding.chm13.bed"
    dict_query_cnt,multicopy_set = parse_rbed(sfile,sfile_rm,sample_id,dict_query_cnt,q_cov_threshold,multicopy_set,fout_chm13_filtered)

fout_chm13_filtered.close()
#sys.exit()



fout = open(f"genecopy_table.chm13.protein_coding.perID_{perID_threshold}.cov{cov_idx}.txt",'w')
fout3 = open(f"genecopy_table.combined_haps.chm13.protein_coding.perID_{perID_threshold}.cov{cov_idx}.txt",'w')

fout_1 = open(f"genecopy_table.chm13.protein_coding.perID_{perID_threshold}.cov{cov_idx}.polymorphic.txt",'w')
fout3_1 = open(f"genecopy_table.combined_haps.chm13.protein_coding.perID_{perID_threshold}.cov{cov_idx}.polymorphic.txt",'w')

fout2 = open(f"sample_specific_multicopy.perID_{perID_threshold}.cov{cov_idx}.txt",'w')


cnnt = {}
candidate_multicopy_set = {}
candidate_multicopy_set_combined_haps = {}
sample_idx_list = sorted(list(set([x.split('_')[0] for x in dict_query_cnt.keys()])))
sample_idx_both_haps = []
for sample_name in sample_idx_list:
    hap1 = sample_name + "_hap1"
    hap2 = sample_name + "_hap2"
    if hap1 in dict_query_cnt.keys():
        if hap2 in dict_query_cnt.keys():
            sample_idx_both_haps.append(sample_name)
sample_idx_both_haps = sorted(sample_idx_both_haps)
dict_query_cnt_combined_hap = {}

for query in multicopy_set:
    for sample_name in sample_idx_both_haps:
        hap1 = sample_name + "_hap1"
        hap2 = sample_name + "_hap2"
        hap1_cnt = 0
        hap2_cnt = 0
        if hap1 in multicopy_set[query]:
            hap1_cnt = multicopy_set[query][hap1]
        if hap2 in multicopy_set[query]:
            hap2_cnt = multicopy_set[query][hap2]
        hap_combined = hap1_cnt + hap2_cnt
        dict_query_cnt_combined_hap.setdefault(sample_name,{})
        dict_query_cnt_combined_hap[sample_name][query] = hap_combined
        if hap_combined > 2:
            candidate_multicopy_set_combined_haps.setdefault(query,{})
            candidate_multicopy_set_combined_haps[query][sample_name] = hap_combined
            
for sample_id in dict_query_cnt:
    for query in dict_query_cnt[sample_id]:
        count = dict_query_cnt[sample_id][query]
        if count > 1:
            candidate_multicopy_set.setdefault(query,{})
            candidate_multicopy_set[query][sample_id] = count
sample_specific_multicopy_set = set()
for query in candidate_multicopy_set:
    if len(candidate_multicopy_set[query].keys()) == 1:
        sample_specific_multicopy_set.add(query)

sample_list2 = sorted(dict_query_cnt_combined_hap.keys())
fout3.write('\t'.join(sample_list2)+'\n')
fout3_1.write('\t'.join(sample_list2)+'\n')
for query in candidate_multicopy_set_combined_haps:
    flag = 0
    if query not in sample_specific_multicopy_set:
        flag = 1
    fout3.write(query)
    if flag:
        fout3_1.write(query)
    for sample_id in sample_list2:
        q_cnt = 0
        if query in dict_query_cnt_combined_hap[sample_id]:
            q_cnt = dict_query_cnt_combined_hap[sample_id][query]
        fout3.write(f'\t{q_cnt}')
        if flag:
            fout3_1.write(f'\t{q_cnt}')
    fout3.write('\n')
    if flag:
        fout3_1.write('\n')

fout3.close()
fout3_1.close()

sample_list = sorted(dict_query_cnt.keys())
fout.write('\t'.join(sample_list)+'\n')
fout_1.write('\t'.join(sample_list)+'\n')
for query in candidate_multicopy_set:
    flag = 0
    if query not in sample_specific_multicopy_set:
        flag = 1
    fout.write(query)
    if flag:
        fout_1.write(query)
    for sample_id in sample_list:
        q_cnt = 0
        if query in dict_query_cnt[sample_id]:
            q_cnt = dict_query_cnt[sample_id][query]
        fout.write(f'\t{q_cnt}')
        if flag:
            fout_1.write(f'\t{q_cnt}')
    fout.write('\n')
    if flag:
        fout_1.write('\n')

    if len(candidate_multicopy_set[query].keys()) == 1:
        for sample_id in candidate_multicopy_set[query]:
            cnnt.setdefault(sample_id,{})
            cnnt[sample_id][query] = 0
fout.close()
fout_1.close()

for sample_id in sample_list:
    cnnnt = 0
    if sample_id in cnnt:
        cnnnt = len(cnnt[sample_id])
    if sample_id == "chm13":continue
    sample_idx, hap = sample_id.split('_')
    info = dd_sample[sample_idx]#[hap]
    fout2.write(sample_id+f'\t{hap}\t'+'\t'.join(info)+f'\t{cnnnt}\n')
fout2.close()

'''
more aligned/bed/asm/NA20509_hap1.chm13_gene.protein_coding.asm.bed
#reference_name	reference_start	reference_end	reference_length	strand	query_name	query_start	query_end	query_length	perID_by_matches	perID_by_eve
nts	perID_by_all	matches	mismatches	deletion_events	insertion_events	deletions	insertions
chr5_RagTag	181950796	181969866	182001041	+	LOC112268260::chr1:20528-37628(-)	5	17100	17100	99.3723	99.301735	87.04741	1678
1	106	10	2	2183	208
chr6_RagTag	170479981	170496747	170527033	+	LOC112268260::chr1:20528-37628(-)	4	17100	17100	99.37663	99.27762	89.46128	15942	100	7	9	724	1054
chr1_RagTag	24983	42878	245507627	-	LOC112268260:
'''
