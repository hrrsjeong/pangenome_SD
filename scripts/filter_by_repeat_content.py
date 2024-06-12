import sys
#set(['srpRNA', 'LTR', 'Satellite', 'Retroposon', 'DNA', 'Simple_repeat', 'Unknown', 'scRNA', 'RC', '.', 'snRNA', 'rRNA', 'tRNA', 'LINE', 'SINE', 'Low_complexity'])
finp = sys.argv[1]
foutput = sys.argv[2]
fout = open(foutput,'w')
dd_tmp = {}
#repeat_class = set()
with open(finp,'r') as fp:
    for line in fp:
        line_temp = line.strip().split('\t')
        idx1 = line_temp[5]
        idx2 = '__'.join(line_temp[0:3])
        dd_tmp.setdefault(idx1,{})
        dd_tmp[idx1].setdefault(idx2,{})
        dd_tmp[idx1][idx2]["info"] = line_temp[0:18]
        dd_tmp[idx1][idx2].setdefault(line_temp[24],0)
        dd_tmp[idx1][idx2][line_temp[24]] += int(line_temp[-1])
        #repeat_class.add(line_temp[24])

#print (repeat_class)
for gene in dd_tmp:
    for idx in dd_tmp[gene]:
        idx_chrom,idx_str,idx_end = idx.split('__')
        idx_len = float(idx_end) - float(idx_str)
        if "Satellite" in dd_tmp[gene][idx]:
            if dd_tmp[gene][idx]["Satellite"] / idx_len > 0.8:
                continue
        if "rRNA" in dd_tmp[gene][idx]:
            if dd_tmp[gene][idx]["rRNA"] / idx_len > 0.8:
                continue
        if "LTR" in dd_tmp[gene][idx]:
            if dd_tmp[gene][idx]["LTR"] / idx_len > 0.8:
                continue
        if "Unknown" in dd_tmp[gene][idx]:
            if dd_tmp[gene][idx]["Unknown"] / idx_len > 0.8:
                continue
        fout.write('\t'.join(dd_tmp[gene][idx]["info"])+'\n')


        '''
        for i in dd_tmp[gene][idx]:
            if i == "info":continue
            print (i+":"+str(dd_tmp[gene][idx][i] / idx_len)+" ,"+str(dd_tmp[gene][idx][i]))
        print ('\t'.join(dd_tmp[gene][idx]["info"]))
        '''
            
fout.close()

'''
chr8_RagTag	25558	41565	145587636	-	LOC112268260::chr1:20528-37628(-)	0	17100	17100	99.45366	99.3817	84.33715	15109	83	5	6	815	1908	chr8_RagTag	30485	30554	L1MD3	69	-	LINE	L1	-1.0	3845613 ??
chr8_RagTag	25558	41565	145587636	-	LOC112268260::chr1:20528-37628(-)	0	17100	17100	99.45366	99.3817	84.33715	15109	83	5	6	815	1908	chr8_RagTag	30560	30738	AluYf1	178	-	SINE	Alu	-1.0	3845615 ??
chr8_RagTag	25558	41565	145587636	-	LOC112268260::chr1:20528-37628(-)	0	17100	17100	99.45366	99.3817	84.33715	15109	83	5	6	815	1908	chr8_RagTag	30749	31301	L1MD3	552	-	LINE	L1	-1.0	3845613 ??
chr8_RagTag	25558	41565	145587636	-	LOC112268260::chr1:20528-37628(-)	0	17100	17100	99.45366	99.3817	84.33715	15109	83	5	6	815	1908	chr8_RagTag	31297	31380	L1MD3	83	-	LINE	L1	-1.0	3845613 ??
'''

