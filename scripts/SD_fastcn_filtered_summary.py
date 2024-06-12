import sys,glob,os
os.system("mkdir -p ./acro_p_filtered_ragtag/tmp")
#flist = glob.glob("chm13_projection/individual/fastCN_filtered_2.5/HG00096.hap1.HGSVC.M.nontrio.GBR.EUR.SD.all.chm13.bed")
flist = glob.glob("chm13_projection/individual/fastCN_filtered_2.5/*.hap*.*.SD.*.chm13.bed")
acro_p_bed = "./acro_p.bed"
for sfile in flist:
    SD_type = sfile.split('/')[-1].split('.')[-3]
    sname,hap = sfile.split('/')[-1].split('.')[0:2]
    acro_p_filtered_tmp = f"./acro_p_filtered_ragtag/tmp/{sname}_{hap}.SDs.autosome.{SD_type}.bed"

    with open(f"tmp_sh/acro_remove.{sname}_{hap}.sh",'w') as fout:
        fout.write(f'cat {sfile} | sort -k1,1 -k2,2n | bedtools subtract -a - -b '+acro_p_bed+f' -A >| {acro_p_filtered_tmp}')
    os.system(f"bash tmp_sh/acro_remove.{sname}_{hap}.sh")
    dd = {}
    acro_p_filtered = acro_p_filtered_tmp.replace("/tmp/","/")
    fout = open(acro_p_filtered,'w')

    with open(acro_p_filtered_tmp,'r') as fp:
        for line in fp:
            line_temp = line.strip().split('\t')[-1].split('__')[2:4]
            for SD_temp in line_temp:
                SD_chr,SD_pos = SD_temp.split(':')
                SD_pos1,SD_pos2 = SD_pos.split('-')
                dd.setdefault(SD_chr,set())
                dd[SD_chr].add((int(SD_pos1),int(SD_pos2)))
        for chrom in dd:
            for pos in dd[chrom]:
                fout.write(f"{chrom}\t{pos[0]}\t{pos[1]}\n")
    fout.close()
    with open(f"tmp_sh/acro_remove.{sname}_{hap}.sh",'w') as fout:
        fout.write(f'cat {acro_p_filtered} | sort -k1,1 -k2,2n | bedtools merge -i - >| {acro_p_filtered.replace(".bed",".merged.bed")}')
    os.system(f"bash tmp_sh/acro_remove.{sname}_{hap}.sh")

            

#chr1	1	140064	HG00096_hap1_nontrio_EUR__all__chr6_RagTag:170577648-170724036__chr1_RagTag:1254724-1400077__0.9909__chm13_0.99094
#chm13_projection/individual/fastCN_filtered_2.5/HG00096.hap1.HGSVC.M.nontrio.GBR.EUR.SD.all.chm13.bed
