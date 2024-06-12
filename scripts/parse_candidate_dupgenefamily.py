import sys,glob
import numpy as np
perID = 90
cov = 60
dd_region = {}
gene_tmp_set = set()
with open(f"chm13_dupgene_region.perID_{perID}.cov{cov}.collapse.bed",'r') as fp:
    for line in fp:
        line_temp = line.strip().split('\t')
        Region = line_temp[3]
        genes = line_temp[6].split(',')
        gene_names = [x.split('::')[0] for x in genes]
        overlap_frac_list = np.array([float(x) for x in line_temp[-1].split(',')])
        if max(overlap_frac_list) < 0.8:continue
        pass_cov50_indices = list(np.nonzero(overlap_frac_list>0.3)[0])
        pass_cov50_pass_genenames_indices = []
        for x in pass_cov50_indices:
            gene_tmp = gene_names[x].split('-')
            if len(gene_tmp) > 1:
                if len(gene_tmp[-1]) < 3:
                    pass_cov50_pass_genenames_indices.append(x)
            else:
                pass_cov50_pass_genenames_indices.append(x)
        if len(pass_cov50_pass_genenames_indices) < 2:continue
        tmp_gene_list = [genes[x] for x in pass_cov50_pass_genenames_indices]
        dd_region[Region] = tmp_gene_list#[genes[x] for x in pass_cov50_pass_genenames_indices]

Region_list = list(dd_region.keys())
#print (Region_list)
tmp_region_list = list(Region_list)
clusters = {}
cnt = 0
cnnt = 0
fout = open(f"deduplicated_genefamilies.perID_{perID}.cov{cov}.chm13.txt",'w')
for i in range(len(Region_list)):
    if Region_list[i] not in tmp_region_list:
        continue
    if Region_list[i] in tmp_region_list:
        tmp_region_list.remove(Region_list[i])
    group1 = set(dd_region[Region_list[i]])
    tmp_region_list2 = list(tmp_region_list)
    flag = 0

    while flag == 0:
        #print (Region_list[i])
        check_len = 0
        for tmp_region in tmp_region_list2:
            group2 = set(dd_region[tmp_region])
            tmp_group1 = set(group1)
            tmp_intersect = tmp_group1.intersection(group2)
            if len(tmp_intersect) >=1:
                group1 = tmp_group1.union(group2)
                if group1 != tmp_group1:
                    check_len +=1
                if tmp_region in tmp_region_list:
                    tmp_region_list.remove(tmp_region)
        if check_len == 0:
            flag = 1
    #print (len(tmp_region_list))
    cnt +=1
    fout.write(f"Cluster_{cnt}\t"+','.join(list(group1))+'\n')
    cnnt += len(list(group1))
fout.close()
print (cnnt)

        




'''
chr3	161843120	162671965	Region13716	3	+	IQCJ::chr3:161843120-162040655(+),IQCJ-SCHIP1::chr3:161843120-162671965(+),SCHIP1::chr3:162047598-162671965(+)	0.238326,1,0.753298
chr3	187058959	187102316	Region13807	3	+	EEF1AKMT4::chr3:187058959-187068061(+),EEF1AKMT4-ECE2::chr3:187058959-187102316(+),ECE2::chr3:187085300-187102316(+)	0.209931,1,0.392463
chr3	188378368	188387552	Region13824	4	-	LOC124902738::chr7:4699519-4703073(-),LOC107984203::chr10:6824494-6828127(-),LOC124901580::chr7:4708965-4716621(-),LOC105369787::chr12:58306186-5831484
6(-)	0.353332,0.364003,0.801938,0.899281
chr3	198221390	198250329	Region13879	2	-	PPP1R2B::chr5:157369311-157371545(+),PPP1R2::chr3:198221390-198250329(-)	0.048274,1
 hsjeong@ocelot  .../SD_projection_chm13/gene_mapping_approach  cat chm13_dupgene_region.raw.bed | awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$10,$NF/($3-$2)}' | bedtools groupby -i - -g 1,2,3,4,5,6 -c 7,8 -o collapse,collapse > chm13
chm13_dupgene_region.raw.bed   chm13_gene.chm13.filtered.bed
 hsjeong@ocelot  .../SD_projection_chm13/gene_mapping_approach  cat chm13_dupgene_region.raw.bed | awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$10,$NF/($3-$2)}' | bedtools groupby -i - -g 1,2,3,4,5,6 -c 7,8 -o collapse,collapse > chm13_dupgene_region.collapse.bed
 '''
