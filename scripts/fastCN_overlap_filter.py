import sys
dd = {}
dd_set = set()
fout = open(sys.argv[1].replace("/tmp/","/"),'w')
with open(sys.argv[1],'r') as fp:
    for line in fp:
        line_temp = line.strip().split('\t')
        dd.setdefault(line_temp[3],{"len":int(line_temp[2])-int(line_temp[1]),"overlap":0})
        dd[line_temp[3]]["overlap"] += int(line_temp[-1])
    fp.seek(0)
    for line in fp:
        line_temp = line.strip().split('\t')
        if line_temp[3] in dd_set:continue
        if float(dd[line_temp[3]]["overlap"]) / float(dd[line_temp[3]]["len"]) > 0.5:
            fout.write('\t'.join(line_temp[0:4])+'\n')
fout.close()
        
        


'''head /net/eichler/vol27/projects/hprc/nobackups/analysis/sedef/SD_projection_chm13/SD_mapping_approach/chm13_projection/individual/fastCN_filtered_3/tmp/HG02011_hap1.SD.intra.chm13.bed
chr1    127787  245234  HG02011_hap1__intra__chr1_RagTag:213702109-213816604__chr1_RagTag:86310-203734__0.9737__chm13_0.97906   chr1    0       324723  117447
chr1    127787  245234  HG02011_hap1__intra__chr1_RagTag:86310-203734__chr1_RagTag:213702109-213816604__0.9737__chm13_0.99941   chr1    0       324723  117447
chr1    138513  140860  HG02011_hap1__intra__chr1_RagTag:97020-99367__chr1_RagTag:217738030-217740331__0.9821__chm13_0.99957    chr1    0       324723  2347
chr1    138514  140951  HG02011_hap1__intra__chr1_RagTag:213713058-213715546__chr1_RagTag:217737917-217740339__0.9681__chm13_0.99122    chr1    0       324723  2437
chr1    138531  140841  HG02011_hap1__intra__chr1_RagTag:2177379
'''
