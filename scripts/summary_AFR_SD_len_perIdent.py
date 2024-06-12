import os,sys,glob
import time

SD_intra = "tmp_intra_AFR_vs_non-AFR_len_ident.txt"
SD_inter = "tmp_inter_AFR_vs_non-AFR_len_ident.txt"
os.system("echo '' >| "+SD_intra)
os.system("echo '' >| "+SD_inter)
for i in range(1000):
    os.system('cat SD_summary_perIdent_len.intra.txt | grep -v "non-AFR" | shuf -n 100 | bedtools groupby -g 1 -c 4,5 -o mean >> '+SD_intra+' &')
    time.sleep(0.05)
    os.system('cat SD_summary_perIdent_len.intra.txt | grep -e "non-AFR" | shuf -n 100 | bedtools groupby -g 1 -c 4,5 -o mean >> '+SD_intra+' &')
    time.sleep(0.05)

for i in range(1000):
    os.system('cat SD_summary_perIdent_len.inter.txt | grep -v "non-AFR" | shuf -n 100 | bedtools groupby -g 1 -c 4,5 -o mean >> '+SD_inter+' &')
    time.sleep(0.05)
    os.system('cat SD_summary_perIdent_len.inter.txt | grep -e "non-AFR" | shuf -n 100 | bedtools groupby -g 1 -c 4,5 -o mean >> '+SD_inter+' &')
    time.sleep(0.05)
