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
qsub_option = '-q eichler-short.q -pe serial 4 -l h_rt=48:00:00 -l mfree=16G'
#####################

sample_info = '/net/eichler/vol27/projects/hprc/nobackups/analysis/sedef/asm-to-reference-alignment/sample_info.txt'
#flist = glob.glob("/net/eichler/vol27/projects/hprc/nobackups/analysis/sedef/SD/sedef/*/SD_align/fasta/SDs_merged_reciprocal90/*_*.SDs.autosome.*.merged90.fasta")
flist = glob.glob("/net/eichler/vol27/projects/hprc/nobackups/analysis/sedef/SD/sedef/*/SD_align/fasta/SDs/*_*.SDs.autosome.*.fasta")
flist = glob.glob("/net/eichler/vol27/projects/hprc/nobackups/analysis/sedef/SD/sedef/*/SD_align/fasta/SDs/*_*.SDs.autosome.inter.fasta")
#os.system(f"mkdir -p aligned/paf/chm13 && mkdir -p aligned/paf/asm")
os.system(f"mkdir -p aligned/bed/chm13")
#os.system(f"mkdir -p aligned/bed/asm")

chm13 = '/net/eichler/vol27/projects/hprc/nobackups/analysis/sedef/repeatmask_t2t/Masked/chm13v2_maskedY_official/chm13v2.0_maskedY.fa'
dd_ragtag_fasta = {}
ragtag_fasta_list = glob.glob("/net/eichler/vol27/projects/hprc/nobackups/analysis/sedef/repeatmask*/fasta/*.RagTag.fasta.fai")
for fasta in ragtag_fasta_list:
    if "RagTag" not in fasta.split('/')[-1]:continue
    info_tmp = fasta.split('/')[-1].split('.RagTag')
    sample_id,hap_tmp = info_tmp[0].split('_')
    hap = hap_tmp.replace("mat","hap2").replace("pat","hap1")
    dd_ragtag_fasta.setdefault(sample_id,{})
    dd_ragtag_fasta[sample_id][hap] = fasta.replace(".fai","")

with open(qsub_list,'w') as fout:
    for sfile in flist:
        print (sfile)
        #if ".inter." in sfile:continue
        if "T2T" in sfile:continue
        info_tmp = sfile.split('/')[-1].split('.')
        sample_id,hap_tmp = info_tmp[0].split('_')
        hap = hap_tmp.replace("mat","hap2").replace("pat","hap1")
        ragtag_fasta = dd_ragtag_fasta[sample_id][hap]
        SD_bed = sfile
        SD_type = sfile.split('/')[-1].split('.autosome.')[-1].split('.')[0]
        #if os.path.exists(f'aligned/bed/chm13/{sample_id}_{hap}.SD.{SD_type}.chm13.bed'):
        #    if os.path.getsize(f'aligned/bed/chm13/{sample_id}_{hap}.SD.{SD_type}.chm13.bed') > 5000:
        #        continue

        command = f' && module load minimap2/2.24 && module load rustybam/0.1.29'
        command += f' && minimap2 -K 1G -t 4 -cx map-ont -f 5000 -k15 -w10 -p 0.05 -N 200 -m200 -s200 -z10000 --secondary=yes --eqx {chm13} {SD_bed} | rb --threads 1 stats --paf > aligned/bed/chm13/{sample_id}_{hap}.SD.{SD_type}.chm13.bed'
        #command += f' && minimap2 -K 1G -t 4 -cx map-ont -f 5000 -k15 -w10 -p 0.05 -N 200 -m200 -s200 -z10000 --secondary=yes --eqx {ragtag_fasta} {SD_bed} | rb --threads 1 stats --paf > aligned/bed/asm/{sample_id}_{hap}.SD.{SD_type}.asm.bed'

        ##old script#
        #command += f' && minimap2 -K 1G -t 4 -cx map-ont -f 5000 -k15 -w10 -p 0.05 -N 200 -m200 -s200 -z2000 --secondary=yes --eqx {chm13} {SD_bed} > aligned/paf/chm13/{sample_id}_{hap}.SD.{SD_type}.chm13.paf'
        #command += f' && rb --threads 1 stats --paf aligned/paf/chm13/{sample_id}_{hap}.SD.{SD_type}.chm13.paf > aligned/bed/chm13/{sample_id}_{hap}.SD.{SD_type}.chm13.bed'
        #will filter out mismatch of chromosome between chm13 and asm
        #command += f' && minimap2 -K 1G -t 4 -cx map-ont -f 5000 -k15 -w10 -p 0.05 -N 200 -m200 -s200 -z2000 --secondary=yes --eqx {ragtag_fasta} {SD_bed} > aligned/paf/asm/{sample_id}_{hap}.SD.{SD_type}.asm.paf'
        #command += f' && rb --threads 1 stats --paf aligned/paf/asm/{sample_id}_{hap}.SD.{SD_type}.asm.paf > aligned/bed/asm/{sample_id}_{hap}.SD.{SD_type}.asm.bed'

        run_systemstr  = systemstr + command
        run_systemstr += '" | qsub '+qsub_option
        fout.write(run_systemstr+'\n')
        #break

with open(qsub_list,'r') as fp:
    for sJob in fp:
        print (sJob)
        #os.system(sJob)

