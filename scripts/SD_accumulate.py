import sys,glob,os,random
import argparse

#SD_type = 'all'

os.system(f"mkdir -p chm13_projection/combined")
os.system(f"mkdir -p chm13_projection/accumulated")
os.system(f"mkdir -p chm13_projection/accumulated/tmp_incl_singleton")
os.system(f"mkdir -p chm13_projection/accumulated/tmp_excl_singleton")

'''
/net/eichler/vol27/projects/hprc/nobackups/analysis/sedef/SD_projection_chm13/SD_mapping_approach/chm13_projection/individual/merged/GM19129.hap1.HGSVC.F.nontrio.YRI.AFR.SD.intra.chm13.bed
chr1	6430	9484
'''
SD_types = ["intra","inter","all"]
for SD_type in SD_types:
    if SD_type == "intra":
        CHM13_SD = 'chm13v2_maskedY.SDs.autosome.intra.merged.bed'
    elif SD_type == "inter":
        CHM13_SD = 'chm13v2_maskedY.SDs.autosome.inter.merged.bed'
    else:
        CHM13_SD = 'chm13v2_maskedY.SDs.autosome.all.merged.bed'
    print (os.getcwd())
    flist = glob.glob(os.getcwd()+f"/chm13_projection/individual/merged/*.SD.{SD_type}.chm13.merged.bed") 
    #chm13_projection/individual/merged/HG00438.hap2.HPRC.F.trio.CHS.EAS.SD.inter.chm13.merged.bed

    dd_info = {}
    for sfile in flist:
        info_tmp = sfile.split('/')[-1].split('.')
        sample_id,hap,dataset,sex,trio,pop,superpop,SD,SD_type,chm13,mergedx,bed = info_tmp
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

    combined_excl_singleton_all_bed = os.getcwd()+f'/chm13_projection/combined/combined.SD.{SD_type}.all_samples.chm13.excl_singleton.bed'
    combined_excl_singleton_trio_bed = os.getcwd()+f'/chm13_projection/combined/combined.SD.{SD_type}.trio.chm13.excl_singleton.bed'
    combined_incl_singleton_all_bed = os.getcwd()+f'/chm13_projection/combined/combined.SD.{SD_type}.all_samples.chm13.incl_singleton.bed'
    combined_incl_singleton_trio_bed = os.getcwd()+f'/chm13_projection/combined/combined.SD.{SD_type}.trio.chm13.incl_singleton.bed'

    if not (os.path.exists(combined_excl_singleton_all_bed) and os.path.exists(combined_incl_singleton_all_bed)):
        os.system(f"cat "+os.getcwd()+f'/chm13_projection/individual/merged_excl_singleton_all/*.SD.{SD_type}.chm13.excl_singleton.all_samples.bed | sort -k1,1 -k2,2n | bedtools merge -i - >| {combined_excl_singleton_all_bed}')
        os.system(f"cat "+os.getcwd()+f'/chm13_projection/individual/merged_excl_singleton_trio/*.SD.{SD_type}.chm13.excl_singleton.trio.bed | sort -k1,1 -k2,2n | bedtools merge -i - >| {combined_excl_singleton_trio_bed}')
        os.system(f"cat "+os.getcwd()+f'/chm13_projection/individual/merged/*.SD.{SD_type}.chm13.merged.bed | sort -k1,1 -k2,2n | bedtools merge -i - >| {combined_incl_singleton_trio_bed}')
        os.system(f"cat "+os.getcwd()+f'/chm13_projection/individual/merged/*.trio.*.SD.{SD_type}.chm13.merged.bed | sort -k1,1 -k2,2n | bedtools merge -i - >| {combined_incl_singleton_all_bed}')
    #chm13_projection/individual/merged/NA21309.hap2.HPRC.F.trio.MKK.AFR.SD.inter.chm13.bed

    for pop_idx in dd_info["population_idx"]:
        sample_pop_list = dd_info["population_idx"][pop_idx]
        trio_sample_pop_list = [x for x in sample_pop_list if dd_info["trio"][x] == "trio"]
        sample_pop_list_input = [os.getcwd()+f'/chm13_projection/individual/merged_excl_singleton_all/{sample_idx}.hap*.*.SD.{SD_type}.chm13.excl_singleton.all_samples.bed' for sample_idx in sample_pop_list]
        trio_sample_pop_list_input = [os.getcwd()+f'/chm13_projection/individual/merged_excl_singleton_trio/{sample_idx}.hap*.*.SD.{SD_type}.chm13.excl_singleton.trio.bed' for sample_idx in trio_sample_pop_list]
        os.system("cat "+" ".join(sample_pop_list_input)+" | sort -k1,1 -k2,2n | bedtools merge -i - >| "+os.getcwd()+f'/chm13_projection/combined/all_samples.{pop_idx}.SD.{SD_type}.chm13.excl_singleton.bed')
        os.system("cat "+" ".join(trio_sample_pop_list_input)+" | sort -k1,1 -k2,2n | bedtools merge -i - >| "+os.getcwd()+f'/chm13_projection/combined/trio.{pop_idx}.SD.{SD_type}.chm13.excl_singleton.bed')

    #sys.exit()


    order_idx_trio = []
    order_idx = []
    os.system(f"mkdir -p chm13_projection/accumulated/tmp")

    dir_path = os.getcwd()+"/chm13_projection/accumulated/tmp"
    tmp_outfile_excl_singleton_trio_CHM13_SD = dir_path+f"/empty.{SD_type}.trio.excl_singleton_CHM13_SD.bed"
    tmp_outfile_excl_singleton_trio = dir_path+f"/empty.{SD_type}.trio.excl_singleton.bed"
    tmp_outfile_incl_singleton_trio = dir_path+f"/empty.{SD_type}.trio.incl_singleton.bed"
    tmp_outfile_shared_trio = dir_path+f"/empty.{SD_type}.trio.shared.bed"
    tmp_outfile_excl_singleton_all_CHM13_SD = dir_path+f"/empty.{SD_type}.all_samples.excl_singleton_CHM13_SD.bed"
    tmp_outfile_excl_singleton_all = dir_path+f"/empty.{SD_type}.all_samples.excl_singleton.bed"
    tmp_outfile_incl_singleton_all = dir_path+f"/empty.{SD_type}.all_samples.incl_singleton.bed"
    tmp_outfile_shared_all = dir_path+f"/empty.{SD_type}.all_samples.shared.bed"
    os.system(f"echo >| {tmp_outfile_excl_singleton_trio_CHM13_SD}")
    os.system(f"echo >| {tmp_outfile_excl_singleton_trio}")
    os.system(f"echo >| {tmp_outfile_incl_singleton_trio}")
    os.system(f"echo >| {tmp_outfile_shared_trio}")
    os.system(f"echo >| {tmp_outfile_excl_singleton_all_CHM13_SD}")
    os.system(f"echo >| {tmp_outfile_excl_singleton_all}")
    os.system(f"echo >| {tmp_outfile_incl_singleton_all}")
    os.system(f"echo >| {tmp_outfile_shared_all}")
    order_idx = {}
    for superpop in superpops:
        for pop in dd_info["population"][superpop]:
            for sample_id in dd_info["population"][superpop][pop]:
                print (sample_id)
                if len(dd_info["sample"][sample_id]) == 1:
                    continue
                for hap in ["hap1","hap2"]:
                    trio_idx = dd_info["trio"][sample_id] 
                    idx = (sample_id,hap,pop,superpop,trio_idx)
                    idx_join = '.'.join(idx)
                    input_file = dd_info["sample"][sample_id][hap]
                    if trio_idx == "trio":
                        order_idx.setdefault("trio",[])
                        order_idx["trio"].append(idx)
                        if "empty" in tmp_outfile_shared_trio:
                            tmp_outfile_shared_trio = str(input_file)
                        outfile_incl_singleton_trio = f"{dir_path}/{idx_join}.{SD_type}.trio.incl_singleton.accumulated.bed"
                        outfile_excl_singleton_trio = f"{dir_path}/{idx_join}.{SD_type}.trio.excl_singleton.accumulated.bed"
                        outfile_excl_singleton_trio_CHM13_SD = f"{dir_path}/{idx_join}.{SD_type}.trio.excl_singleton_CHM13_SD.accumulated.bed"
                        outfile_shared_trio = f"{dir_path}/{idx_join}.{SD_type}.trio.shared.accumulated.bed"
                        intersect_incl_singleton_trio = f"cat {tmp_outfile_incl_singleton_trio} {input_file} | sort -k1,1 -k2,2n | bedtools merge -i - >| {dir_path}/{idx_join}.{SD_type}.trio.incl_singleton.accumulated.bed"
                        intersect_excl_singleton_trio = f"cat {tmp_outfile_excl_singleton_trio} {input_file} | sort -k1,1 -k2,2n | bedtools merge -i - | bedtools intersect -a {combined_excl_singleton_trio_bed} -b - | sort -k1,1 -k2,2n | bedtools merge -i - >| {dir_path}/{idx_join}.{SD_type}.trio.excl_singleton.accumulated.bed"
                        intersect_excl_singleton_trio_CHM13_SD = f"cat {tmp_outfile_excl_singleton_trio} {input_file} | sort -k1,1 -k2,2n | bedtools merge -i - | bedtools intersect -a {combined_excl_singleton_trio_bed} -b - | sort -k1,1 -k2,2n | bedtools merge -i - | bedtools intersect -a - -b {CHM13_SD} >| {dir_path}/{idx_join}.{SD_type}.trio.excl_singleton_CHM13_SD.accumulated.bed"
                        intersect_shared_trio = f"bedtools intersect -a {tmp_outfile_shared_trio} -b {input_file} | sort -k1,1 -k2,2n | bedtools merge -i - >| {dir_path}/{idx_join}.{SD_type}.trio.shared.accumulated.bed"
                        os.system(intersect_incl_singleton_trio)
                        os.system(intersect_excl_singleton_trio)
                        os.system(intersect_excl_singleton_trio_CHM13_SD)
                        os.system(intersect_shared_trio)
                        tmp_outfile_incl_singleton_trio = str(outfile_incl_singleton_trio)
                        tmp_outfile_excl_singleton_trio = str(outfile_excl_singleton_trio)
                        tmp_outfile_excl_singleton_trio_CHM13_SD = str(outfile_excl_singleton_trio_CHM13_SD)
                        tmp_outfile_shared_trio = str(outfile_shared_trio)

                    order_idx.setdefault("all_samples",[])
                    order_idx["all_samples"].append(idx)
                    if "empty" in tmp_outfile_shared_all:
                        tmp_outfile_shared_all = str(input_file)
                    outfile_incl_singleton_all = f"{dir_path}/{idx_join}.{SD_type}.all_samples.incl_singleton.accumulated.bed"
                    outfile_excl_singleton_all = f"{dir_path}/{idx_join}.{SD_type}.all_samples.excl_singleton.accumulated.bed"
                    outfile_excl_singleton_all_CHM13_SD = f"{dir_path}/{idx_join}.{SD_type}.all_samples.excl_singleton_CHM13_SD.accumulated.bed"
                    outfile_shared_all = f"{dir_path}/{idx_join}.{SD_type}.all_samples.shared.accumulated.bed"
                    intersect_incl_singleton_all = f"cat {tmp_outfile_incl_singleton_all} {input_file} | sort -k1,1 -k2,2n | bedtools merge -i - >| {dir_path}/{idx_join}.{SD_type}.all_samples.incl_singleton.accumulated.bed"
                    intersect_excl_singleton_all = f"cat {tmp_outfile_excl_singleton_all} {input_file} | sort -k1,1 -k2,2n | bedtools merge -i - | bedtools intersect -a {combined_excl_singleton_all_bed} -b - | sort -k1,1 -k2,2n | bedtools merge -i - >| {dir_path}/{idx_join}.{SD_type}.all_samples.excl_singleton.accumulated.bed"
                    intersect_excl_singleton_all_CHM13_SD = f"cat {tmp_outfile_excl_singleton_all} {input_file} | sort -k1,1 -k2,2n | bedtools merge -i - | bedtools intersect -a {combined_excl_singleton_all_bed} -b - | sort -k1,1 -k2,2n | bedtools merge -i - | bedtools intersect -a - -b {CHM13_SD} >| {dir_path}/{idx_join}.{SD_type}.all_samples.excl_singleton_CHM13_SD.accumulated.bed"
                    intersect_shared_all = f"bedtools intersect -a {tmp_outfile_shared_all} -b {input_file} | sort -k1,1 -k2,2n | bedtools merge -i - >| {dir_path}/{idx_join}.{SD_type}.all_samples.shared.accumulated.bed"
                    os.system(intersect_incl_singleton_all)
                    os.system(intersect_excl_singleton_all)
                    os.system(intersect_excl_singleton_all_CHM13_SD)
                    os.system(intersect_shared_all)
                    tmp_outfile_incl_singleton_all = str(outfile_incl_singleton_all)
                    tmp_outfile_excl_singleton_all = str(outfile_excl_singleton_all)
                    tmp_outfile_excl_singleton_all_CHM13_SD = str(outfile_excl_singleton_all_CHM13_SD)
                    tmp_outfile_shared_all = str(outfile_shared_all)

    def parse_bed(finp):
        cnt = 0
        region_sum = 0
        with open(finp,'r') as fbed:
            for line in fbed:
                cnt +=1
                line_temp = line.strip().split('\t')
                region_sum += int(line_temp[2]) - int(line_temp[1])
        return cnt,region_sum

    for data_idx in order_idx:
        fout = open(f"{dir_path.replace('/tmp','')}/SD_accumulated.{SD_type}.{data_idx}.txt",'w')
        fout.write("sample\thap\tpop\tsuperpop\ttrio\tnum_SD\tsum_SD\tnum_SD_polymorphic\tsum_SD_polymorphic\tnum_SD_polymorphic_CHM13\tsum_SD_polymorphic_CHM13\tnum_SD_shared\tsum_SD_shared\n")
        for idx in range(len(order_idx[data_idx])):
            print (idx)
            info = order_idx[data_idx][idx]
            info_join = '.'.join(info)
            fp_incl_singleton = f"{dir_path}/{info_join}.{SD_type}.{data_idx}.incl_singleton.accumulated.bed"
            fp_excl_singleton = f"{dir_path}/{info_join}.{SD_type}.{data_idx}.excl_singleton.accumulated.bed"
            fp_excl_singleton_CHM13_SD = f"{dir_path}/{info_join}.{SD_type}.{data_idx}.excl_singleton_CHM13_SD.accumulated.bed"
            fp_shared = f"{dir_path}/{info_join}.{SD_type}.{data_idx}.shared.accumulated.bed"
            SD_cnt_incl_singleton,SD_sum_incl_singleton = parse_bed(fp_incl_singleton)
            SD_cnt_excl_singleton,SD_sum_excl_singleton = parse_bed(fp_excl_singleton)
            SD_cnt_excl_singleton_CHM13_SD,SD_sum_excl_singleton_CHM13_SD = parse_bed(fp_excl_singleton_CHM13_SD)
            SD_cnt_shared,SD_sum_shared = parse_bed(fp_shared)
            fout.write('\t'.join(info))
            fout.write(f"\t{SD_cnt_incl_singleton}\t{SD_sum_incl_singleton}")
            fout.write(f"\t{SD_cnt_excl_singleton}\t{SD_sum_excl_singleton}")
            fout.write(f"\t{SD_cnt_excl_singleton_CHM13_SD}\t{SD_sum_excl_singleton_CHM13_SD}")
            fout.write(f"\t{SD_cnt_shared}\t{SD_sum_shared}\n")
        fout.close()


'''
HG00741	PUR	AMR	trio	3469	195828923	3001	169551463	2682	100593922
HG00731	PUR	AMR	nontrio	3494	196951902	3014	170087846	2688	100476622
 ⚙ hsjeong@ocelot  .../SD_projection_chm13/SD_mapping_approach  head -n 2 chm13_projection/accumulated/SD_accumulated.trio.combined_haps.txt
HG04217	ITU	SAS	trio	2777	125669388	2769	124939840	2777	125669388
HG03807	BEB	SAS	trio	2832	131655006	2813	130759305
'''









