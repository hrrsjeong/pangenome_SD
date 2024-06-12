import sys,os,glob

perID = 90
cov = 60
dd_cluster = {}
cnt = 0
with open(f"deduplicated_genefamilies.perID_{perID}.cov{cov}.chm13.txt",'r') as fp:
    for line in fp:
        cnt +=1
        line_temp = line.strip().split('\t')
        genes = line_temp[1].split(',')
        gene_fam = f"GeneFamily_{cnt}"
        for gene in genes:
            dd_cluster[gene] = gene_fam

dd_SD_support = {}
with open(f"../SD_mapping_approach/SD_genes/SD_genes.hprc_hgsvc.txt",'r') as fp:
    for line in fp:
        line_temp = line.strip().split('::')[0]
        dd_SD_support[line_temp] = ''

dd_CN_summary = {}
with open(f"summary_CN_variation.txt",'r') as fp:
    summary_title = fp.readline().strip().split('\t')[1:]
    for line in fp:
        line_temp = line.strip().split('\t')
        dd_CN_summary[line_temp[0]] = line_temp[1:]

dd_multiexonic_orf200 = {}
with open("../../repeatmask_t2t/Masked/chm13v2_maskedY_official/multiexonic_protein_coding.orf200.txt",'r') as fp:
    fp.readline()
    for line in fp:
        line_temp = line.strip().split('\t')
        dd_multiexonic_orf200[line_temp[0]] = line_temp[1:]

hap_CN = {}
with open(f"genecopy_table.chm13.protein_coding.perID_{perID}.cov{cov}.txt",'r') as fp:
    hap_CN_title = fp.readline().strip().split('\t')[0:-1]
    for line in fp:
        line_temp = line.strip().split('\t')
        gene = line_temp[0]
        CN = line_temp[1:-1]
        chm13 = line_temp[-1]
        hap_CN[gene] = [chm13,CN]


finput =  f"genecopy_table.combined_haps.chm13.protein_coding.perID_{perID}.cov{cov}.polymorphic.type.txt"
fout = open(f"asm_based_CN_table.perID_{perID}.cov{cov}.txt",'w')
with open(finput,'r') as fp:
    title = fp.readline().strip().split('\t')
    fout.write("Gene_family\tGene\tPos\tGene_length\tNum_Exon\tORF_length_bp\thas_SD\tCN_in_CHM13\tMean_CN_diploid\tMean_CN_AFR\tMean_CN_non-AFR\tstd_CN\tstd_CN_AFR\tstd_CN_non-AFR\tCV_CN\tCV_CN_AFR\tCV_CN_non-AFR\tMannWhitneyU_Stat_AFRvsnon-AFR\tP_value\t"+"\t".join(hap_CN_title)+'\n')
    for line in fp:
        line_temp = line.strip().split('\t')
        gene_name,pos = line_temp[1].split('::')
        if line_temp[1] not in dd_CN_summary:continue
        if line_temp[1] in dd_cluster:
            gene_fam = dd_cluster[line_temp[1]]
        else:
            cnt += 1
            gene_fam = f"GeneFamily_{cnt}"
        start,end = [int(x) for x in pos.split(':')[-1].split('(')[0].split('-')]
        length = str(end-start+1)
        num_exon,orf_len = dd_multiexonic_orf200[gene_name]
        SD_support = "FALSE"
        if gene_name in dd_SD_support:
            SD_support = "TRUE"
        
        fout.write(f"{gene_fam}\t{gene_name}\t{pos}\t{length}\t{num_exon}\t{orf_len}\t{SD_support}\t{hap_CN[line_temp[1]][0]}\t")
        fout.write("\t".join(dd_CN_summary[line_temp[1]])+"\t")
        fout.write("\t".join(hap_CN[line_temp[1]][1])+'\n')
fout.close()
'''
summary_CN_variation.txt
gene	avg_genecopy	avg_genecopy_AFR	avg_genecopy_nonAFR	sd_genecopy	sd_genecopy_AFR	sd_genecopy_nonAFR	cv_genecopy	cv_genecopy_AFR	cv_genecopy_nonAFR	test_stat	p_val
LOC112268260::chr1:20528-37628(-)	10.74	11.73	9.92	1.96	1.95	1.55	0.183	0.167	0.157	467.5	5.312339848884724e-06

more genecopy_table.combined_haps.chm13.protein_coding.perID_90.cov60.poorphic.type.txt
Duplication_Type	Gene	GM19129	GM19320	GM20355	GM21487	HG00096	HG00171	HG002	HG00268	HG00358	HG00438	HG005	HG00512	HG00513	HG0051
8	HG03125	HG03248	HG03371	HG03452	HG03453	HG03456	HG03486	HG03492	HG03516	HG03540	HG03579	HG03683	HG03732	HG03807	HG04036	HG04217	NA1232
9	NA12878	NA18906	NA19238	NA19239	NA19240	NA19317	NA19331	NA19347	NA19650	NA19705	NA19983	NA20129	NA20509	NA20847	NA21309
known	LOC112268260::chr1:20528-37628(-)	11	10	7	16	9	9	10	10	12	8	9	9	88	9	12	11	9	13	10	9	8	10	10	9	9	8	12	11	11	11	98	11	13	11	11	8	12	9	13	11	11	13	9	13	11	11	12	8	12	9	15	9	13	10	14	12	13	11	13	11	12	9	12	9	13	13	17	16	9	12	8	12	12	12	12	11	9	10	9	9	13	10	13	12	11	13	7	10	8	10	12	10	12	8	11
'''
