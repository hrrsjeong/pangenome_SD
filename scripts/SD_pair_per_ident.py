import sys,glob,os
import time

SD_types = ["intra","inter","all"]
for SD_type_name in SD_types:
    flist = glob.glob(f"chm13_projection/individual/merged/*.hap*.*.SD.{SD_type_name}.chm13.merged.bed")
    #chm13_projection/individual/merged/HG00513.hap2.HGSVC.F.nontrio.CHS.EAS.SD.intra.chm13.merged.bed
    dd = {}
    for sfile in flist:
        sample_idx,hap,data_type,sex,trio,subpop,pop,SD,SD_type,x,y,z = sfile.split('/')[-1].split('.')
        fastcn_file = f"fastCN_filtered_ragtag/{sample_idx}_{hap}.SDs.autosome.{SD_type}.bed"
        sedef_file = f"SD_raw/{sample_idx}_{hap}.SDs.autosome.{SD_type}.bed" #fastCN_filtered_ragtag/{sample_idx}_{hap}.SDs.autosome.{SD_type}.bed"
        
        if os.path.exists(fastcn_file):
            dd.setdefault(sample_idx+"__"+hap,{})
            dd[sample_idx+"__"+hap]["fastcn"] = fastcn_file
            dd[sample_idx+"__"+hap]["self"] = sfile.split('/')[-1]
            if os.path.exists(sedef_file):
                dd[sample_idx+"__"+hap]["sedef"] = sedef_file
            else:
                if hap == "hap1":
                    sedef_file = sedef_file.replace("hap1","pat")
                    if os.path.exists(sedef_file):
                        dd[sample_idx+"__"+hap]["sedef"] = sedef_file
                elif hap == "hap2":
                    sedef_file = sedef_file.replace("hap2","mat")
                    if os.path.exists(sedef_file):
                        dd[sample_idx+"__"+hap]["sedef"] = sedef_file
    os.system("mkdir -p SD_pair_fastcn_validated")
    print (len(dd))
    for i in dd:
        fastcn = dd[i]["fastcn"]
        sedef  = dd[i]["sedef"]
        fname = dd[i]["self"]
        systemstr = f"bedtools intersect -a {sedef} -b {fastcn} -f 0.5 -wa >| ./SD_pair_fastcn_validated/{fname.replace('.chm13.merged.bed','.fastcn_valid.bed')} &"
        print (systemstr)
        os.system(systemstr)
        time.sleep(0.3)
    time.sleep(10)
    fout = open("tmp_SD_summary.sh",'w')
    fout.write("echo '' >| SD_summary_perIdent_len."+SD_type_name+".txt\n")
    for i in dd:
        input_name = "SD_pair_fastcn_validated/"+dd[i]["self"].replace(".chm13.merged.bed",".fastcn_valid.bed")
        if ".AFR." in input_name:
            pop = "AFR"
        else:
            pop = "non-AFR"

        fout.write("cat "+input_name+" | awk -v OFS='\\t' '{print"+' "'+pop+'","'+dd[i]["self"].split('.chm13')[0]+'",$1,$15,$24'+"}' >> SD_summary_perIdent_len."+SD_type_name+".txt\n")
    fout.close()
    os.system("bash tmp_SD_summary.sh")


    
        
