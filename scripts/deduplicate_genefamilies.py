import sys,os
dd = {}
perID = 90
cov = 60

with open(f"genecopy_table.combined_haps.chm13.protein_coding.perID_{perID}.cov{cov}.polymorphic.txt",'r') as fp:
    for line in fp:
        line_temp = line.strip().split('\t')
        dd[line_temp[0]] = ''

chm13_self_align = f"chm13_gene.chm13.filtered.perID_{perID}.cov{cov}.bed"

with open(chm13_self_align,'w') as fout:
    with open("tmp_"+chm13_self_align,'r') as fp:
        for line in fp:
            line_temp = line.strip().split("\t")
            if line_temp[3] not in dd:continue
            fout.write(line)

with open("tmp_dupgene.sh",'w') as fout:
    #bedtools intersect -a aligned/bed/chm13/chm13.chm13_gene.protein_coding.chm13.bed -b ../../repeatmask_t2t/Masked/chm13v2_maskedY_official/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.merged.bed -f 0.9 -wa | awk -v OFS='\t' '{if ($9 > 1000 && $10 > 90 && $8-$7 > $9*0.6) print $0}' | cut -f 6 | sort | uniq > blacklist_gene.repeat_overlap.txt
    #bedtools intersect -a <(cat chm13_gene.chm13.filtered.perID_90.cov60.bed | sort -k1,1 -k2,2n | bedtools merge -i - -s -c 1,6 -o count,distinct | awk -v OFS='\t' 'START{n=0}{n+=1; if ($4>1) print $1,$2,$3,"Region"n,$4,$5}') -b chm13_gene.chm13.filtered.perID_90.cov60.bed -s -wo >| chm13_dupgene_region.raw.bed
    #cat chm13_dupgene_region.raw.bed | awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$10,$NF/($3-$2)}' | bedtools groupby -i - -g 1,2,3,4,5,6 -c 7,8 -o collapse,collapse > chm13_dupgene_region.collapse.bed
    fout.write("bedtools intersect -a <(cat "+chm13_self_align+" | sort -k1,1 -k2,2n | bedtools merge -i - -s -c 1,6 -o count,distinct | awk -v OFS='\\t' 'START{n=0}{n+=1; if ($4>1) print $1,$2,$3,"+'"Region"n,$4,$5'+"}') -b "+chm13_self_align+f" -s -wo >| chm13_dupgene_region.perID_{perID}.cov{cov}.raw.bed\n")
    fout.write(f"cat chm13_dupgene_region.perID_{perID}.cov{cov}.raw.bed"+" | awk -v OFS='\\t' '{print $1,$2,$3,$4,$5,$6,$10,$NF/($3-$2)}' | bedtools groupby -i - -g 1,2,3,4,5,6 -c 7,8 -o collapse,collapse >|"+f" chm13_dupgene_region.perID_{perID}.cov{cov}.collapse.bed\n")

os.system("bash tmp_dupgene.sh")

