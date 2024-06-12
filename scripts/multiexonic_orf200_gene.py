import sys,glob,os
import time

'''
head ../../repeatmask_t2t/Masked/chm13v2_maskedY_official/multiexonic_protein_coding.orf200.txt
gene_id	num_exons	orf_bp_len
LOC112268260	5	1602
SAMD11	14	2538
NOC2L	19	2250
KLHL17	12	1929
'''

dd_multiexonic_orf200 = {}
with open("../../repeatmask_t2t/Masked/chm13v2_maskedY_official/multiexonic_protein_coding.orf200.txt",'r') as fp:
    fp.readline()
    for line in fp:
        line_temp = line.strip().split('\t')
        dd_multiexonic_orf200[line_temp[0]] = ''
        
os.system(f"mkdir -p aligned/bed/asm_ragtag_multiexonic_orf200bp/")
os.system(f"mkdir -p aligned/bed/chm13_multiexonic_orf200bp/")

flist = glob.glob("aligned/bed/asm_ragtag/*_hap*.chm13_gene.protein_coding.asm.ragtag.bed")
for sfile in flist:
    print (sfile)
    with open(f"aligned/bed/asm_ragtag_multiexonic_orf200bp/"+sfile.split('/')[-1],'w') as fout:
        with open(sfile,'r') as fp:
            for line in fp:
                if line.startswith("#"):
                    fout.write(line)
                line_temp = line.strip().split('\t')
                gene = line_temp[5].split('::')[0]
                if gene in dd_multiexonic_orf200:
                    fout.write(line)
with open(f"aligned/bed/chm13_multiexonic_orf200bp/chm13.chm13_gene.protein_coding.chm13.bed",'w') as fout:
    with open(f"aligned/bed/chm13/chm13.chm13_gene.protein_coding.chm13.bed",'r') as fp:
        for line in fp:
            if line.startswith("#"):
                fout.write(line)
            line_temp = line.strip().split('\t')
            gene = line_temp[5].split('::')[0]
            if gene in dd_multiexonic_orf200:
                fout.write(line)


        


'''
bedtools intersect -a aligned/bed/chm13/chm13.chm13_gene.protein_coding.chm13.bed -b ../../repeatmask_t2t/Masked/chm13v2_maskedY_official/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.merged.bed -f 0.8 -wa | awk -v OFS='\t' '{if ($9 > 1000 && $10 > 90 && $8-$7 > $9*0.6) print $0}' | cut -f 6 | sort | uniq >| blacklist_gene.repeat_overlap.txt
'''
