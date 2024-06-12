import sys,glob,os,random
import argparse

###################
cwd = os.getcwd()
rannum = '{:05}'.format(random.randrange(1, 10**5))
#qsub_list = f"qsub_enrichment_{rannum}.sh"
qsub_list = f"qsub_SD_combined_chm13_projection.sh"
#conda_env = "liftoff"

systemstr   = 'echo "cd '+cwd
systemstr  += ' && . ~/.bash_profile && . ~/.bashrc_qsub && module load miniconda/4.9.2'
#systemstr  += ' && conda activate '+conda_env
qsub_option = '-q eichler-short.q -pe serial 1 -l h_rt=24:00:00 -l mfree=12G'
#####################

os.system(f"mkdir -p chm13_projection/individual/merged")
os.system(f"mkdir -p chm13_projection/individual/merged_excl_singleton_all")
os.system(f"mkdir -p chm13_projection/individual/merged_excl_singleton_trio")
fout = open(qsub_list,'w')
SD_types = ["intra","inter","all"]
for SD_type in SD_types:
    flist = glob.glob(f'./chm13_projection/individual/merged/*.*.*.*.*.*.*.SD.{SD_type}.chm13.merged.bed')
    flist_trio = glob.glob(f'./chm13_projection/individual/merged/*.*.*.*.trio.*.*.SD.{SD_type}.chm13.merged.bed')
    flist2 = list(flist)
    flist2_trio = list(flist_trio)
    #with open(qsub_list,'w') as fout:
    for sfile in flist:
        print (sfile)
        tmp_flist_pop = [x for x in flist2 if x != sfile]
        out_file1 = f'{sfile.replace("/merged/","/merged_excl_singleton_all/").replace(".chm13.merged.bed",".chm13.excl_singleton.all_samples.bed")}'
        command = f' && bedtools intersect -a {sfile} -b '+ ' '.join(tmp_flist_pop)+f' | sort -k1,1 -k2,2n | bedtools merge -i - >| {out_file1}'
        if ".trio." in sfile:
            tmp_flist_pop2 = [x for x in flist2_trio if x != sfile]
            out_file2 = f'{sfile.replace("/merged/","/merged_excl_singleton_trio/").replace(".chm13.merged.bed",".chm13.excl_singleton.trio.bed")}'
            command += f' && bedtools intersect -a {sfile} -b '+ ' '.join(tmp_flist_pop2)+f' | sort -k1,1 -k2,2n | bedtools merge -i - >| {out_file2}'
        run_systemstr  = systemstr + command
        run_systemstr += '" | qsub '+qsub_option
        fout.write(run_systemstr+'\n')
fout.close()
