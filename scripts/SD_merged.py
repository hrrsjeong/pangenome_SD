import sys,glob,os,random
import argparse

###################
cwd = os.getcwd()
rannum = '{:05}'.format(random.randrange(1, 10**5))
#qsub_list = f"qsub_enrichment_{rannum}.sh"
qsub_list = f"qsub_SD_merge_chm13_projection.sh"
#conda_env = "liftoff"

systemstr   = 'echo "cd '+cwd
systemstr  += ' && . ~/.bash_profile && . ~/.bashrc_qsub && module load miniconda/4.9.2'
#systemstr  += ' && conda activate '+conda_env
qsub_option = '-q eichler-short.q -pe serial 1 -l h_rt=24:00:00 -l mfree=12G'
#####################

#SD_type = 'intra'

sample_info = '/net/eichler/vol28/projects/hprc/nobackups/analysis/sedef/asm-to-reference-alignment/sample_info.txt'
dd_sample = {}

'''
Source	Sample	Sex	Trio data available	Population	Superpopulation	Family ID	Relation
HPRC	HG01891	Female	1	ACB	AFR	BB05	child
'''
with open(sample_info,'r') as fp:
    fp.readline()
    for line in fp:
        line_temp = line.strip().split('\t')
        sample = line_temp[1]
        if line_temp[3] == "1":
            trio_idx = "trio"
        else:
            trio_idx = "nontrio"
        dd_sample[sample] = [line_temp[0],line_temp[2][0],trio_idx,line_temp[4],line_temp[5]]

os.system(f"mkdir -p chm13_projection/individual/raw")
#os.system(f"mkdir -p chm13_projection/individual/fastCN_filtered")
os.system(f"mkdir -p chm13_projection/individual/fastCN_filtered_2.5/tmp")
#os.system(f"mkdir -p chm13_projection/individual/fastCN_filtered_3/tmp")
os.system(f"mkdir -p chm13_projection/individual/merged")
os.system(f"mkdir -p chm13_projection/individual/merged_excl_singleton_all")
os.system(f"mkdir -p chm13_projection/individual/merged_excl_singleton_trio")
'''
os.system(f"mkdir -p chm13_projection/combined")#/incl_singleton")
os.system(f"mkdir -p chm13_projection/accumulated")
os.system(f"mkdir -p chm13_projection/accumulated/tmp_incl_singleton")
os.system(f"mkdir -p chm13_projection/accumulated/tmp_excl_singleton")
'''

#/net/eichler/vol27/projects/human_primate_brain_ips_transcriptomes/nobackups/NPIP_processing_misc/analysis/sd_var/fastcn_dupspace/HG00096.bed
#fastCN_list = glob.glob("/net/eichler/vol27/projects/human_primate_brain_ips_transcriptomes/nobackups/NPIP_processing_misc/analysis/sd_var/fastcn_dupspace/*.bed")
#dd_fastCN = {}
#for fastCN_file in fastCN_list:
#    sample_id = fastCN_file.split('/')[-1].split('.')[0].split("_3_")[0]
#    dd_fastCN.setdefault(sample_id,{})# = fastCN_file
#    if "_3_" in fastCN_file:
#        dd_fastCN[sample_id]["3"] = fastCN_file
#    else:
#        dd_fastCN[sample_id]["2.5"] = fastCN_file

fastCN_list = glob.glob("/net/eichler/vol28/projects/human_primate_brain_ips_transcriptomes/nobackups/NPIP_processing_misc/analysis/sd_var/fastcn_dupspace/20240221/*_2.5*.bed")
dd_fastCN = {}
for fastCN_file in fastCN_list:
    sample_id = fastCN_file.split('/')[-1].split('_2')[0]#.split("_3_")[0]
    dd_fastCN.setdefault(sample_id,{})# = fastCN_file
    dd_fastCN[sample_id]["2.5"] = fastCN_file


#/net/eichler/vol28/projects/human_primate_brain_ips_transcriptomes/nobackups/NPIP_processing_misc/analysis/sd_var/fastcn_dupspace/20240221/*_2.5*


'''
HG00268_hap1.SD.inter.chm13.bed  HG01123_hap2.SD.intra.chm13.bed  HG02145_hap1.SD.all.chm13.bed    HG03065_hap1.SD.inter.chm13.bed  NA12329_hap2.SD.intra.chm13.bed
HG00096.hap1.HGSVC.M.nontrio.GBR.EUR.SD.all.chm13.bed    HG01928.hap2.HPRC.M.trio.PEL.AMR.SD.all.chm13.bed        HG03371.hap1.HGSVC.M.nontrio.ESN.AFR.SD.all.chm13.bed
'''
fout = open(qsub_list,'w')
SD_types = ["intra","inter","all"]
for SD_type in SD_types:
    flist = glob.glob(os.getcwd()+f"/chm13_projection/individual/raw/*.SD.{SD_type}.chm13.bed") #chm13_projection/raw/HG00096_hap1.SD.intra.chm13.bed
    for sfile in flist:
        #info_tmp = sfile.split('/')[-1].split('.')
        sample_id,hap = sfile.split('/')[-1].split('.')[0:2]
        #sample_id,hap = info_tmp[0].split('_')
        if sample_id not in dd_sample:
            print (sample_id,hap,"not in the valid sample list")
            continue
        if sample_id not in dd_fastCN:
            print (sample_id,hap,"not in the fastCN list")
            continue
        if "00514" in sample_id:continue
        if "00733" in sample_id:continue
        if "01890" in sample_id:continue
        if "03452" in sample_id:continue
        if "19240" in sample_id:continue


        if "01596" in sample_id:continue
        if "02769" in sample_id:continue
        if "18534" in sample_id:continue
        if "18989" in sample_id:continue
        if "19036" in sample_id:continue
        if "19384" in sample_id:continue
        dataset,sex,trio,pop,superpop = dd_sample[sample_id]
        out_fastCN_filtered_2_tmp = sfile.replace("/raw/","/fastCN_filtered_2.5/tmp/")
        out_fastCN_filtered_2 = out_fastCN_filtered_2_tmp.replace("/tmp/","/")
        out_filename_merged = sfile.replace("/raw/","/merged/").replace(".chm13.bed",".chm13.merged.bed")
        #out_filename_merged = os.getcwd()+f"/chm13_projection/individual/merged/{sample_id}.{hap}.{dataset}.{sex}.{trio}.{pop}.{superpop}.SD.{SD_type}.chm13.bed"
        if os.path.exists(out_filename_merged):continue
        fastCN_2_file = dd_fastCN[sample_id]["2.5"]
        systemstr2 = f" && cat {sfile} | sort -k1,1 -k2,2n | bedtools intersect -a - -b {fastCN_2_file} -wao | bedtools subtract -a - -b acro_p.bed -A >| {out_fastCN_filtered_2_tmp}" 
        systemstr2 += f" && python fastCN_overlap_filter.py {out_fastCN_filtered_2_tmp}"
        systemstr2 += f" && cat {out_fastCN_filtered_2} | sort -k1,1 -k2,2n | bedtools merge -i - >| {out_filename_merged}" 
        #print (systemstr2)
        #os.system(systemstr2)
        run_systemstr  = systemstr + systemstr2
        run_systemstr += '" | qsub '+qsub_option
        fout.write(run_systemstr+'\n')
        #sys.exit()
fout.close()
'''
    flist = glob.glob(f'./chm13_projection/individual/merged/*.*.*.*.*.*.*.SD.{SD_type}.chm13.bed')
    flist_trio = glob.glob(f'./chm13_projection/individual/merged/*.*.*.*.trio.*.*.SD.{SD_type}.chm13.bed')
    flist2 = list(flist)
    flist2_trio = list(flist_trio)
    #with open(qsub_list,'w') as fout:
    for sfile in flist:
        print (sfile)
        tmp_flist_pop = [x for x in flist2 if x != sfile]
        out_file1 = f'{sfile.replace("/merged/","/merged_excl_singleton_all/").replace(".chm13.bed",".chm13.excl_singleton.all_samples.bed")}'
        command = f' && bedtools intersect -a {sfile} -b '+ ' '.join(tmp_flist_pop)+f' | sort -k1,1 -k2,2n | bedtools merge -i - >| {out_file1}'
        if ".trio." in sfile:
            tmp_flist_pop2 = [x for x in flist2_trio if x != sfile]
            out_file2 = f'{sfile.replace("/merged/","/merged_excl_singleton_trio/").replace(".chm13.bed",".chm13.excl_singleton.trio.bed")}'
            command += f' && bedtools intersect -a {sfile} -b '+ ' '.join(tmp_flist_pop2)+f' | sort -k1,1 -k2,2n | bedtools merge -i - >| {out_file2}'
        run_systemstr  = systemstr + command
        run_systemstr += '" | qsub '+qsub_option
        fout.write(run_systemstr+'\n')
'''
#fout.close()
