import sys,glob,os,random
import argparse

###################
cwd = os.getcwd()
rannum = '{:05}'.format(random.randrange(1, 10**5))
#qsub_list = f"qsub_enrichment_{rannum}.sh"
qsub_list = f"qsub_minimap2.sh"
#conda_env = "liftoff"

systemstr   = 'echo "cd '+cwd
systemstr  += ' && . ~/.bash_profile && . ~/.bashrc_qsub && module load miniconda/4.9.2'
#systemstr  += ' && conda activate '+conda_env
qsub_option = '-q eichler-short.q -pe serial 4 -l h_rt=72:00:00 -l mfree=12G'
#####################

sample_info = '/net/eichler/vol28/projects/hprc/nobackups/analysis/sedef/asm-to-reference-alignment/sample_info.txt'
flist = glob.glob("/net/eichler/vol28/projects/hprc/nobackups/analysis/sedef/SD/sedef/*/SD_align/fasta/SDs_merged_reciprocal90/*_*.SDs.autosome.intra.merged90.fasta")
#../../SD/sedef/HPRC_chrom/SD_align/bed/SDs_merged_reciprocal90/HG00673_mat.SDs.autosome.all.merged90.bed
#os.system(f"mkdir -p aligned/paf/chm13 && mkdir -p aligned/paf/asm")
os.system(f"mkdir -p aligned/bed/chm13 && mkdir -p aligned/bed/asm_contig && mkdir -p aligned/bed/asm_ragtag")
os.system(f"mkdir -p aligned_rm/bed/chm13 && mkdir -p aligned_rm/bed/asm_contig && mkdir -p aligned_rm/bed/asm_ragtag")

chm13 = '/net/eichler/vol28/projects/hprc/nobackups/analysis/sedef/repeatmask_t2t/Masked/chm13v2_maskedY_official/chm13v2.0_maskedY.fa'
chm13_gene_fasta = '/net/eichler/vol28/projects/hprc/nobackups/analysis/sedef/repeatmask_t2t/Masked/chm13v2_maskedY_official/chm13v2.0_maskedY.gene.fasta'
chm13_gene_rm_fasta = '/net/eichler/vol28/projects/hprc/nobackups/analysis/sedef/repeatmask_t2t/Masked/chm13v2_maskedY_official/chm13v2.0_maskedY.gene.rm.fasta'
dd_ragtag_fasta = {}
ragtag_fasta_list = glob.glob("/net/eichler/vol28/projects/hprc/nobackups/analysis/sedef/repeatmask*/fasta/*.RagTag.fasta.fai")
for fasta in ragtag_fasta_list:
    if "RagTag" not in fasta.split('/')[-1]:continue
    info_tmp = fasta.split('/')[-1].split('.RagTag')
    sample_id,hap_tmp = info_tmp[0].split('_')
    hap = hap_tmp.replace("mat","hap2").replace("pat","hap1")
    dd_ragtag_fasta.setdefault(sample_id,{})
    dd_ragtag_fasta[sample_id][hap] = fasta.replace(".fai","")

dd_contig_fasta = {}
contig_fasta_list = glob.glob("/net/eichler/vol28/projects/hprc/nobackups/analysis/sedef/asm-to-reference-alignment/fasta_gz/*.hap*.fasta.gz.fai")
#/net/eichler/vol27/projects/hprc/nobackups/analysis/sedef/repeatmask*/fasta/*.RagTag.fasta.fai
for fasta in contig_fasta_list:
    #if "RagTag" not in fasta.split('/')[-1]:continue
    info_tmp = fasta.split('/')[-1].split('.fasta')
    sample_id,hap = info_tmp[0].split('.')
    #hap = hap_tmp.replace("mat","hap2").replace("pat","hap1")
    dd_contig_fasta.setdefault(sample_id,{})
    dd_contig_fasta[sample_id][hap] = fasta.replace(".fai","")

#/net/eichler/vol27/projects/hprc/nobackups/analysis/sedef/asm-to-reference-alignment/fasta_gz/HG01106.hap1.fasta.gz


with open(qsub_list,'w') as fout:
    command = f' && module load minimap2/2.24 && module load rustybam/0.1.29'
    #command += f' && minimap2 -K 1G -t 4 -cx map-ont -f 5000 -k15 -w10 -p 0.05 -N 200 -m200 -s200 -z10000 --secondary=yes --eqx {chm13} {chm13_gene_fasta} | rb --threads 1 stats --paf >| aligned/bed/chm13/chm13.chm13_gene.protein_coding.chm13.bed'
    command += f' && minimap2 -K 1G -t 4 -cx map-ont -f 5000 -k15 -w10 -p 0.05 -N 200 -m100 -s100 -z20000 --secondary=yes --eqx {chm13} {chm13_gene_rm_fasta} | rb --threads 1 stats --paf >| aligned_rm/bed/chm13/chm13.chm13_gene.protein_coding.chm13.bed'
    ###command += f' && minimap2 -K 1G -t 4 -cx map-ont -f 5000 -k15 -w10 -p 0.05 -N 200 -m200 -s200 -z10000 --secondary=yes --eqx {chm13} {chm13_gene_fasta} >| aligned/paf/chm13/chm13.chm13_gene.protein_coding.chm13.paf'
    ###command += f' && rb --threads 1 stats --paf aligned/paf/chm13/chm13.chm13_gene.protein_coding.chm13.paf >| aligned/bed/chm13/chm13.chm13_gene.protein_coding.chm13.bed'
    run_systemstr  = systemstr + command
    run_systemstr += '" | qsub '+qsub_option
    fout.write(run_systemstr+'\n')
    for sfile in flist:
        if "T2T" in sfile:continue
        print (sfile)
        info_tmp = sfile.split('/')[-1].split('.')
        sample_id,hap_tmp = info_tmp[0].split('_')
        hap = hap_tmp.replace("mat","hap2").replace("pat","hap1")
        ragtag_fasta = dd_ragtag_fasta[sample_id][hap]
        contig_fasta = dd_contig_fasta[sample_id][hap]
        SD_bed = sfile
        SD_type = sfile.split('/')[-1].split('.')[-3]

        command = f' && module load minimap2/2.24 && module load rustybam/0.1.29'
        #command += f' && minimap2 -K 1G -t 4 -cx map-ont -f 5000 -k15 -w10 -p 0.05 -N 200 -m200 -s200 -z10000 --secondary=yes --eqx {ragtag_fasta} {chm13_gene_fasta} | rb --threads 1 stats --paf >| aligned/bed/asm_ragtag/{sample_id}_{hap}.chm13_gene.protein_coding.asm.ragtag.bed'
        #command += f' && minimap2 -K 1G -t 4 -cx map-ont -f 5000 -k15 -w10 -p 0.05 -N 200 -m200 -s200 -z10000 --secondary=yes --eqx {contig_fasta} {chm13_gene_fasta} | rb --threads 1 stats --paf >| aligned/bed/asm_contig/{sample_id}_{hap}.chm13_gene.protein_coding.asm.contig.bed'
        command += f' && minimap2 -K 1G -t 4 -cx map-ont -f 5000 -k15 -w10 -p 0.05 -N 200 -m100 -s100 -z20000 --secondary=yes --eqx {ragtag_fasta} {chm13_gene_rm_fasta} | rb --threads 1 stats --paf >| aligned_rm/bed/asm_ragtag/{sample_id}_{hap}.chm13_gene.protein_coding.asm.ragtag.bed'
        #command += f' && minimap2 -K 1G -t 4 -cx map-ont -f 5000 -k15 -w10 -p 0.05 -N 200 -m50 -s50 -z20000 --secondary=yes --eqx {contig_fasta} {chm13_gene_rm_fasta} | rb --threads 1 stats --paf >| aligned_rm/bed/asm_contig/{sample_id}_{hap}.chm13_gene.protein_coding.asm.contig.bed'
        ###command += f' && minimap2 -K 1G -t 4 -cx map-ont -f 5000 -k15 -w10 -p 0.05 -N 200 -m200 -s200 -z10000 --secondary=yes --eqx {ragtag_fasta} {chm13_gene_fasta} >| aligned/paf/asm/{sample_id}_{hap}.chm13_gene.protein_coding.asm.paf'
        ###command += f' && rb --threads 1 stats --paf aligned/paf/asm/{sample_id}_{hap}.chm13_gene.protein_coding.asm.paf >| aligned/bed/asm/{sample_id}_{hap}.chm13_gene.protein_coding.asm.bed'

        run_systemstr  = systemstr + command
        run_systemstr += '" | qsub '+qsub_option
        fout.write(run_systemstr+'\n')
        #break

with open(qsub_list,'r') as fp:
    for sJob in fp:
        print (sJob)
        #os.system(sJob)

