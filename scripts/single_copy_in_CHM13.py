import os,sys

with open("tmp_singlecopy_chm13.sh",'w') as fout:
    fout.write("cat genecopy_table.chm13.protein_coding.perID_90.cov60.txt | tail -n +2 | awk '{if ($NF == 1) print $0}' | cut -f 1 >| single_copy_in_chm13.txt\n")

os.system("bash tmp_singlecopy_chm13.sh")
dd = set()
with open("single_copy_in_chm13.txt",'r') as fp:
    for line in fp:
        dd.add(line.strip())

with open("genecopy_table.combined_haps.chm13.protein_coding.perID_90.cov60.polymorphic.type.txt",'w') as fout:
    with open("genecopy_table.combined_haps.chm13.protein_coding.perID_90.cov60.polymorphic.txt",'r') as fp:
        title = "Duplication_Type\tGene\t"+fp.readline()
        fout.write(title)
        for line in fp:
            line_temp = line.strip().split('\t')
            if line_temp[0] in dd:
                fout.write("new\t"+line)
            else:
                fout.write("known\t"+line)

