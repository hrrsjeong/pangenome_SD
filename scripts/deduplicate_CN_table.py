import sys,os
import numpy as np
perID = 90
cov = 60

dd = {}
dd_info = {}
dd_loc = {}
dd_loc_info = {}
dd_cluster = {}
cluster_flag = set()
with open(f"deduplicated_genefamilies.perID_{perID}.cov{cov}.chm13.txt",'r') as fp:
    for line in fp:
        line_temp = line.strip().split('\t')
        genes = line_temp[1].split(',')
        for gene in genes:
            dd_cluster[gene] = line_temp[0]

finput =  f"genecopy_table.combined_haps.chm13.protein_coding.perID_{perID}.cov{cov}.polymorphic.txt"
fout = open(finput.replace(".txt",".dedup_genefam.txt"),'w')
with open(finput,'r') as fp:
    title = fp.readline().strip()
    fout.write(title+'\n')
    for line in fp:
        line_temp = line.strip().split('\t')
        mean_CN = np.nanmean([float(x) for x in line_temp[1:]])
        if mean_CN < 2:continue
        if line_temp[0].startswith("LOC"):
            continue
        dd.setdefault(line_temp[0],{})
        dd_info[line_temp[0]] = line
        dd[line_temp[0]] = mean_CN

sorted_items = sorted(dd.items(), key=lambda item: item[1],reverse=True)
for i,j in sorted_items:
    if i not in dd_cluster:
        fout.write(dd_info[i])
    else:
        cluster_idx = dd_cluster[i]
        if cluster_idx in cluster_flag:continue
        cluster_flag.add(cluster_idx)
        fout.write(dd_info[i])
    
with open(finput,'r') as fp:
    fp.seek(0)
    title = fp.readline().strip()
    #fout.write(title+'\n')
    for line in fp:
        line_temp = line.strip().split('\t')
        mean_CN = np.nanmean([float(x) for x in line_temp[1:]])
        if mean_CN < 2:continue
        if not line_temp[0].startswith("LOC"):
            continue
        dd_loc.setdefault(line_temp[0],{})
        dd_loc_info[line_temp[0]] = line
        dd_loc[line_temp[0]] = mean_CN

sorted_items = sorted(dd_loc.items(), key=lambda item: item[1],reverse=True)
for i,j in sorted_items:
    if i not in dd_cluster:
        fout.write(dd_loc_info[i])
    else:
        cluster_idx = dd_cluster[i]
        if cluster_idx in cluster_flag:continue
        cluster_flag.add(cluster_idx)
        fout.write(dd_loc_info[i])

fout.close()


dd = {}
dd_info = {}
dd_loc = {}
dd_loc_info = {}
cluster_flag = set()

with open("SD_support.sh",'w') as fout:
    fout.write(f"cat <(head -n 1 genecopy_table.combined_haps.chm13.protein_coding.perID_{perID}.cov{cov}.polymorphic.txt) <(grep -f ../SD_mapping_approach/SD_genes/SD_genes.hprc_hgsvc.txt genecopy_table.combined_haps.chm13.protein_coding.perID_{perID}.cov{cov}.polymorphic.txt) >| genecopy_table.combined_haps.chm13.protein_coding.perID_{perID}.cov{cov}.polymorphic.SD_support.txt\n")
os.system("bash SD_support.sh")

finput =  f"genecopy_table.combined_haps.chm13.protein_coding.perID_{perID}.cov{cov}.polymorphic.SD_support.txt"
fout = open(finput.replace(".txt",".dedup_genefam.txt"),'w')
with open(finput,'r') as fp:
    title = fp.readline().strip()
    fout.write(title+'\n')
    for line in fp:
        line_temp = line.strip().split('\t')
        mean_CN = np.nanmean([float(x) for x in line_temp[1:]])
        if mean_CN < 2:continue
        if line_temp[0].startswith("LOC"):
            continue
        dd.setdefault(line_temp[0],{})
        dd_info[line_temp[0]] = line
        dd[line_temp[0]] = mean_CN

sorted_items = sorted(dd.items(), key=lambda item: item[1],reverse=True)
for i,j in sorted_items:
    if i not in dd_cluster:
        fout.write(dd_info[i])
    else:
        cluster_idx = dd_cluster[i]
        if cluster_idx in cluster_flag:continue
        cluster_flag.add(cluster_idx)
        fout.write(dd_info[i])
    
with open(finput,'r') as fp:
    fp.seek(0)
    title = fp.readline().strip()
    #fout.write(title+'\n')
    for line in fp:
        line_temp = line.strip().split('\t')
        mean_CN = np.nanmean([float(x) for x in line_temp[1:]])
        if mean_CN < 2:continue
        if not line_temp[0].startswith("LOC"):
            continue
        dd_loc.setdefault(line_temp[0],{})
        dd_loc_info[line_temp[0]] = line
        dd_loc[line_temp[0]] = mean_CN

sorted_items = sorted(dd_loc.items(), key=lambda item: item[1],reverse=True)
for i,j in sorted_items:
    if i not in dd_cluster:
        fout.write(dd_loc_info[i])
    else:
        cluster_idx = dd_cluster[i]
        if cluster_idx in cluster_flag:continue
        cluster_flag.add(cluster_idx)
        fout.write(dd_loc_info[i])
    
fout.close()


'''
head genecopy_table.combined_haps.chm13.protein_coding.perID_90.cov60.polymorphic.txt
GM19129	GM19320	GM20355	GM21487	HG00096	HG00171	HG002	HG00268	HG00358	HG00438	HG005	HG00512	HG00513	HG00514	HG00621	HG00673	HG00731	HG00732	HG00733	HG00735	HG00741	HG00864	HG01071	HG01106	HG01114	HG01123	HG01175	HG01258HG01352	HG01358	HG01361	HG01457	HG01505	HG01890	HG01891	HG01928	HG01952	HG01978	HG02011	HG02018	HG02055	HG02059	HG02080	HG02106	HG02109	HG02145	HG02148	HG02257	HG02486	HG02492	HG02554	HG02559	HG02572	HG02587	HG02622	HG02630HG02666	HG02717	HG02723	HG02818	HG02886	HG02953	HG03009	HG03065	HG03098	HG03125	HG03248	HG03371	HG03452	HG03453	HG03456	HG03486	HG03492	HG03516	HG03540	HG03579	HG03683	HG03732	HG03807	HG04036	HG04217	NA12329	NA12878	NA18906NA19238	NA19239	NA19240	NA19317	NA19331	NA19347	NA19650	NA19705	NA19983	NA20129	NA20509	NA20847	NA21309
LOC112268260::chr1:20528-37628(-)	11	10	7	16	9	9	10	10	12	8	9	9	8	8	9	12	12	10	13	10	9	8	10	10	9	9	8	12	11	11	11	9	8	11	13	11	11	8	12	9	13	11	11	13	9	13	11	11	12	8	12	915	9	13	10	14	12	13	11	13	11	12	9	12	9	13	13	17	16	9	12	8	12	12	12	12	11	9	10	99	13	10	13	12	11	13	7	10	8	10	12	10	13	8	11
PERM1::chr1:403906-411933(-)	13	10	7	12	2	14	7	2	3	6	3	2	2	2	2	2	2	9	9	4	2	4	6	3	12	7	4	6	5	3	3	6	6	2	10	3	3	2	6	2	6	2	3	2	12	14	7	6	2	5	8	5	66	5	9	8	2	7	5	2	8	2	7	11	9	7	10	2	2	2	4	2	2	6	7	2	9	2	2	6	82	6	2	2	8	2	8	2	6	6	2	5	5	6	6


head deduplicated_genefamilies.chm13.txt
Cluster_1	ATAD3A::chr1:946665-969189(+),ATAD3B::chr1:906135-943970(+),ATAD3C::chr1:884060-904534(+)
Cluster_2	CDK11A::chr1:1136810-1158725(-),CDK11B::chr1:1069725-1093466(-)
Cluster_3	PRAMEF2::chr1:12401097-12405915(+),PRAMEF14::chr1:12783098-12788341(-),PRAMEF1::chr1:12335684-12340916(+)
Cluster_4	PRAMEF15::chr1:12756817-12763835(+),PRAMEF9::chr1:12566136-12573753(-),PRAMEF25::chr1:12589588-12596703(+),PRAMEF5::chr1:12695430-12704668(+),PRAMEF6::chr1:12482393-12491489(-),PRAMEF4::chr1:12423167-12430160(-),PRAMEF11::chr1:12368615-12375422(-)
'''
