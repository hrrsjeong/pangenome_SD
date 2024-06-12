import sys,os
import numpy as np
from scipy import stats

perID = 90
cov = 60

dd = {}

with open("fastcn.sample.list",'r') as fp:
    for line in fp:
        line_temp = line.strip().split('/')[-1].split('.')
        if line_temp[0] in dd:continue
        dd[line_temp[0]] = line_temp[2:7]

sample_list = []
AFR_sample_list = []
nonAFR_sample_list = []
fout = open("summary_CN_variation.txt",'w')
fout.write("gene")
fout.write("\tavg_genecopy\tavg_genecopy_AFR\tavg_genecopy_nonAFR")
fout.write("\tsd_genecopy\tsd_genecopy_AFR\tsd_genecopy_nonAFR")
fout.write("\tcv_genecopy\tcv_genecopy_AFR\tcv_genecopy_nonAFR")
fout.write("\ttest_stat\tp_val\n")
##with open(f"genecopy_table.combined_haps.chm13.protein_coding.perID_{perID}.cov{cov}.polymorphic.SD_support.dedup_genefam.txt",'r') as fp:

with open(f"genecopy_table.combined_haps.chm13.protein_coding.perID_{perID}.cov{cov}.polymorphic.txt",'r') as fp:
    title = fp.readline().strip().split('\t')
    for i in title:
        if i not in dd:continue
        sample_list.append(i)
        if dd[i][4] == "AFR":
            AFR_sample_list.append(i)
        else:
            nonAFR_sample_list.append(i)
    sample_idx_list = [title.index(x) for x in sample_list]
    AFR_sample_idx_list = [title.index(x) for x in AFR_sample_list]
    nonAFR_sample_idx_list = [title.index(x) for x in nonAFR_sample_list]

    for line in fp:
        line_temp = line.strip().split('\t')
        gene_name = line_temp[0]
        count_list = np.array([float(x) for x in line_temp[1:]])
        #print (count_list)
        sample_count_list = [count_list[x] for x in sample_idx_list]
        AFR_sample_count_list = [count_list[x] for x in AFR_sample_idx_list]
        nonAFR_sample_count_list = [count_list[x] for x in nonAFR_sample_idx_list]
        mean_all,sd_all = np.mean(sample_count_list),np.std(sample_count_list)
        mean_AFR,sd_AFR = np.mean(AFR_sample_count_list),np.std(AFR_sample_count_list)
        mean_nonAFR,sd_nonAFR = np.mean(nonAFR_sample_count_list),np.std(nonAFR_sample_count_list)
        cv_all = (1+(1/(4*len(sample_idx_list))))*(sd_all/mean_all)
        cv_AFR = (1+(1/(4*len(AFR_sample_idx_list))))*(sd_AFR/mean_AFR)
        cv_nonAFR = (1+(1/(4*len(nonAFR_sample_idx_list))))*(sd_nonAFR/mean_nonAFR)
        if len(list(set(sample_count_list))) == 1:
            MWU_AFRvsnonAFR = [0,"nan"]
        else:
            #ttest_AFRvsnonAFR = stats.ttest_ind(AFR_sample_count_list,nonAFR_sample_count_list)
            MWU_AFRvsnonAFR = stats.mannwhitneyu(AFR_sample_count_list,nonAFR_sample_count_list)
        #print (ttest_AFRvsnonAFR)
        if mean_all < 2:continue
        fout.write(gene_name)
        fout.write('\t'+str(round(mean_all,2)))
        fout.write('\t'+str(round(mean_AFR,2)))
        fout.write('\t'+str(round(mean_nonAFR,2)))
        fout.write('\t'+str(round(sd_all,2)))
        fout.write('\t'+str(round(sd_AFR,2)))
        fout.write('\t'+str(round(sd_nonAFR,2)))
        fout.write('\t'+str(round(cv_all,3)))
        fout.write('\t'+str(round(cv_AFR,3)))
        fout.write('\t'+str(round(cv_nonAFR,3)))
        #fout.write('\t'+str(round(ttest_AFRvsnonAFR[0],3)))
        fout.write('\t'+str(round(MWU_AFRvsnonAFR[0],3)))
        #fout.write('\t'+str(ttest_AFRvsnonAFR[1]))
        fout.write('\t'+str(MWU_AFRvsnonAFR[1]))
        fout.write('\n')
fout.close()
        




'''
../SD_mapping_approach/chm13_projection/individual/fastCN_filtered_2.5/NA20847.hap1.HGSVC.F.nontrio.GIH.SAS.SD.intra.chm13.bed
../SD_mapping_approach/chm13_projection/individual/fastCN_filtered_2.5/NA20847.hap2.HGSVC.F.nontrio.GIH.SAS.SD.intra.chm13.bed


GM19129	GM19320	GM20355	GM21487	HG00096	HG00171	HG002	HG00268	HG00358	HG00438	HG005	HG00512	HG00513	HG00514	HG00621	HG00673	HG00731	HG00732	HG00733	HG00735	HG00741	HG00864HG01071	HG01106	HG01114	HG01123	HG01175	HG01258	HG01352	HG01358	HG01361	HG01457	HG01505	HG01890	HG01891	HG01928	HG01952	HG01978	HG02011	HG02018	HG02055	HG02059	HG02080	HG02106HG02109	HG02145	HG02148	HG02257	HG02486	HG02492	HG02554	HG02559	HG02572	HG02587	HG02622	HG02630	HG02666	HG02717	HG02723	HG02818	HG02886	HG02953	HG03009	HG03065	HG03098	HG03125HG03248	HG03371	HG03452	HG03453	HG03456	HG03486	HG03492	HG03516	HG03540	HG03579	HG03683	HG03732	HG03807	HG04036	HG04217	NA12329	NA12878	NA18906	NA19238	NA19239	NA19240	NA19317NA19331	NA19347	NA19650	NA19705	NA19983	NA20129	NA20509	NA20847	NA21309
LOC112268260::chr1:20528-37628(-)	11	10	7	16	9	9	10	10	12	8	9	9	8	8	9	12	12	10	13	10	9	8	10	10	9	9	8	12	11	11	11	9	8	11	13	11	11	8	12	913	11	11	13	9	13	11	11	12	8	12	9	15	9	13	10	14	12	13	11	13	11	12	9	12	9	13	13	17	16	9	12	8	12	12	12	12	11	9	10	9	9	13	10	13	12	11	13	7	10	8	10	12	10	13	8	11
PERM1::chr1:403906-411933(-)	13	10	7	12	2	14	7	2	3	6	3	2	2	2	2	2	2	9	94	2	4	6	3	12	7	4	6	5	3	3	6	6	2	10	3	3	2	6	2	6	23	2	12	14	7	6	2	5	8	5	6	6	5	9	8	2	7	5	2	8	2	7	11	9	7	10	2	2	2	4	2	2	6	7	2	9	2	2	6	8	2	6	2	2	82	8	2	6	6	2	5	5	6	6
SCNN1D::chr1:709442-723551(+)	2	2	2	2	2	2	2	2	3	2	2	4	3	3	2	2	1	2	13	4	2	2	2	2	2	2	3	3	3	2	3	2	4	2	3	2	2	2	2	2	22	4	2	2	2	2	3	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	23	2	2	2	2	2	2	3	2	2	2	4	2	2	5	2	2	3	2	4	2	2	22	2	2	3	2	2	1	2	3
DVL1::chr1:766862-781044(-)	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	22	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	32	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	22	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	22	2	2	3	2	2	2	2	2
MXRA8::chr1:783234-794097(-)	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	1	2	2	22	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	32	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	22	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	22	2	2	3	2	2	2	2	2
AURKAIP1::chr1:804284-805756(-)	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	1	2	2	22	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	32	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	22	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	22	2	2	3	2	2	2	2	2
CCNL2::chr1:817644-831265(-)	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	1	2	2	22	2	2	2	2	4	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	22	2	2	2	2	2	2	2	2	2	1	2	2	2	2	2	2	2	2	2	2	2	22	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	22	2	2	3	2	2	2	2	2
MRPL20::chr1:833838-839223(-)	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	1	2	2	22	2	2	2	2	4	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	22	2	2	2	2	2	2	2	2	2	1	2	2	2	2	2	2	2	2	2	2	2	21	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	22	2	2	3	2	2	2	2	2
ANKRD65::chr1:850349-853206(-)	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	1	2	2	12	2	2	2	2	4	2	2	2	2	2	2	2	2	2	2	2	1	2	2	2	2	12	1	2	2	2	2	2	1	2	2	1	2	2	2	2	2	2	2	2	2	2	2	21	2	2	2	2	2	2	2	2	2	2	0	2	2	2	2	2	2	2	2	2	2	22	2	2	3	2	2	2	2	2
'''
