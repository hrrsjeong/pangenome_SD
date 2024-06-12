import sys,glob,os,random
import argparse
import time
#SD_type = 'all'

os.system(f"mkdir -p chm13_projection/accumulated_doubleton/tmp")

'''
/net/eichler/vol27/projects/hprc/nobackups/analysis/sedef/SD_projection_chm13/SD_mapping_approach/chm13_projection/individual/merged/GM19129.hap1.HGSVC.F.nontrio.YRI.AFR.SD.intra.chm13.bed
chr1	6430	9484
'''
SD_types = ["intra","all"]
for SD_type in SD_types:
    flist = glob.glob(os.getcwd()+f"/chm13_projection/individual/merged/*.SD.{SD_type}.chm13.merged.bed")
    #HG00438.hap1.HPRC.F.trio.CHS.EAS.SD.all.chm13.merged.bed
    dd_info = {}
    for sfile in flist:
        info_tmp = sfile.split('/')[-1].split('.')
        sample_id,hap,dataset,sex,trio,pop,superpop,SD,SD_type,chm13,xx,bed = info_tmp
        dd_info.setdefault("sample",{})
        dd_info["sample"].setdefault(sample_id,{})
        dd_info["sample"][sample_id][hap] = sfile
        dd_info.setdefault("trio",{})
        dd_info["trio"][sample_id] = trio
        #dd_info["trio"].setdefault(trio,set())
        #dd_info["trio"][trio].add(sample_id)
        dd_info.setdefault("population",{})
        if superpop == "AFR":
            superpop_idx = "AFR"
        else:
            superpop_idx = "non-AFR"
        dd_info["population"].setdefault(superpop,{})
        dd_info["population"][superpop].setdefault(pop,set())
        dd_info["population"][superpop][pop].add(sample_id)
        dd_info.setdefault("population_idx",{})
        dd_info["population_idx"].setdefault(superpop_idx,set())
        dd_info["population_idx"][superpop_idx].add(sample_id)

    #print (sorted(dd_info["population"].keys(),reverse=True))
    #sys.exit()

    #['SAS', 'EUR', 'EAS', 'AMR', 'AFR']
    superpops = sorted(dd_info["population"].keys(),reverse= True)
    #superpops = sorted(dd_info["population_idx"].keys(),reverse= True)
    #superpops = ['SAS', 'AMR', 'EAS', 'EUR', 'AFR']
    print (superpops)
    #trio_idx = sorted(dd_info["trio"].keys(),reverse= True)

    order_idx_trio = []
    order_idx = []
    os.system(f"mkdir -p chm13_projection/accumulated_doubleton/tmp_combined/tmp")

    dir_path = os.getcwd()+"/chm13_projection/accumulated_doubleton/tmp_combined/tmp"
    order_idx = {}
    order_sfile_accumulate = []
    for superpop in superpops:
        for pop in dd_info["population"][superpop]:
            for sample_id in dd_info["population"][superpop][pop]:
                print (sample_id)
                if len(dd_info["sample"][sample_id]) == 1:
                    continue
                trio_idx = dd_info["trio"][sample_id] 
                #if trio_idx != "trio":continue
                idx = (sample_id,pop,superpop,trio_idx)
                idx_join = '.'.join(idx)
                os.system(f'cat {dd_info["sample"][sample_id]["hap1"]} {dd_info["sample"][sample_id]["hap2"]} | sort -k1,1 -k2,2n | bedtools merge -i - >| {dir_path}/../{idx_join}.{SD_type}.bed')
                input_file = f'{dir_path}/../{idx_join}.{SD_type}.bed'#dd_info["sample"][sample_id][hap]
                order_sfile_accumulate.append(input_file)
                order_sfile_idx = len(order_sfile_accumulate)
                out_covfile = f"{dir_path}/{order_sfile_idx}__{idx_join}.{SD_type}.cov.bed"
                join_files = ' '.join(order_sfile_accumulate)
                cov_command = f'cat {join_files} | sort -k1,1 -k2,2n | bedtools genomecov -i - -g ./t2t.autosome.chrom.length -bg >| {out_covfile} &'
                os.system(cov_command)
                time.sleep(0.2)
    time.sleep(80)

    def parse_bed(finp):
        cnt = 0
        region_sum = 0
        with open(finp,'r') as fbed:
            for line in fbed:
                line_temp = line.strip().split('\t')
                region_tmp = int(line_temp[2]) - int(line_temp[1])
                if region_tmp >= 100:cnt+=1
                region_sum += region_tmp
        return cnt,region_sum

    os.system(f'mkdir -p {dir_path}/tmp2')

    for list_idx in reversed(range(1,order_sfile_idx+1)):
        tmp_cov_sh = f'./tmp/{list_idx}_{SD_type}.combined.sh'
        target_cov_file = glob.glob(f"{dir_path}/{list_idx}__*.{SD_type}.cov.bed")[0]
        if len(target_cov_file) < 10:
            print ("glob error",target_cov_file)
            sys.exit()
        target_file_idx = target_cov_file.split('/')[-1].split('.cov.bed')[0]
        if list_idx >= 1:
            #prior_cov_file = glob.glob(f"{dir_path}/{list_idx-1}__*.{SD_type}.cov.bed")[0]
            with open(tmp_cov_sh,'w') as tmp_fout_cov:
                tmp_fout_cov.write("cat "+target_cov_file+" | awk -v OFS='\\t' '{if ($4 > 0) print $0}' >| "+dir_path+"/tmp2/"+target_file_idx+'.singleton.cov.bed\n')
                tmp_fout_cov.write("cat "+target_cov_file+" | awk -v OFS='\\t' '{if ($4 > 1) print $0}' >| "+dir_path+"/tmp2/"+target_file_idx+'.doubleton.cov.bed\n')
                tmp_fout_cov.write("cat "+target_cov_file+" | awk -v OFS='\\t' '{if ($4 > 2) print $0}' >| "+dir_path+"/tmp2/"+target_file_idx+'.polymorphic.cov.bed\n')
                tmp_fout_cov.write("cat "+target_cov_file+" | awk -v OFS='\\t' '{if ($4 >= "+str(list_idx)+"/2) print $0}' >| "+dir_path+"/tmp2/"+target_file_idx+'.maf.cov.bed\n')
                tmp_fout_cov.write("cat "+target_cov_file+" | awk -v OFS='\\t' '{if ($4 == "+str(list_idx)+") print $0}' >| "+dir_path+"/tmp2/"+target_file_idx+'.shared.cov.bed\n')
                tmp_fout_cov.write("cat "+target_cov_file+" | awk -v OFS='\\t' '{if ($4 >= "+str(list_idx)+"*90/100) print $0}' >| "+dir_path+"/tmp2/"+target_file_idx+'.shared90.cov.bed\n')

            os.system("bash "+tmp_cov_sh+' &')
            time.sleep(0.5)

    time.sleep(20)
    fout = open(f"{dir_path}/../../new_SD_accumulated.{SD_type}.combined_haps.txt",'w')
    fout.write("sample\tpop\tsuperpop\ttrio\tSD_type\tnum_SD_singleton\tsum_SD_singleton\tnum_SD_doubleton\tsum_SD_doubleton\tnum_SD_polymorphic\tsum_SD_polymorphic\tnum_SD_shared50\tsum_SD_shared50\tnum_SD_shared90\tsum_SD_shared90\tnum_SD_shared\tsum_SD_shared\n")
    stype_order = ["singleton","doubleton","polymorphic","maf","shared90","shared"]

    for list_idx in range(1,order_sfile_idx+1):
        target_cov_file = glob.glob(f"{dir_path}/{list_idx}__*.{SD_type}.cov.bed")[0]
        target_file_idx = target_cov_file.split('/')[-1].split('.cov.bed')[0]
        info = target_file_idx.split('__')[-1].split('.')
        dd_tmp  = {}
        fout.write('\t'.join(info))
        for tmp_stype in stype_order:
            SD_cnt,SD_sum = parse_bed(f"{dir_path}/tmp2/{target_file_idx}.{tmp_stype}.cov.bed")
            fout.write(f"\t{SD_cnt}\t{SD_sum}")
        fout.write('\n')
    fout.close()


'''
HG00741	PUR	AMR	trio	3469	195828923	3001	169551463	2682	100593922
HG00731	PUR	AMR	nontrio	3494	196951902	3014	170087846	2688	100476622
 ⚙ hsjeong@ocelot  .../SD_projection_chm13/SD_mapping_approach  head -n 2 chm13_projection/accumulated/SD_accumulated.trio.combined_haps.txt
HG04217	ITU	SAS	trio	2777	125669388	2769	124939840	2777	125669388
HG03807	BEB	SAS	trio	2832	131655006	2813	130759305
'''









