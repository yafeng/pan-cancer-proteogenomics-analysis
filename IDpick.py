import sys
import operator
import re

input1 = open(sys.argv[1],"r")
output1 = open(sys.argv[2],"w")
output2 = open(sys.argv[3],"w")

header = input1.readline().split("\t")
output1.write("\t".join(header))
loci_col = header.index("Loci")

print ("Loci column number is",loci_col)

loci_dic={}
for line in input1:
    row = line.split("\t")
    row[1] = re.sub(r" ?\([^)]+\)", "", row[1])
    accs = row[1]
    loci = row[loci_col]
    if loci not in loci_dic:
        loci_dic[loci]=[]
    
    if "PGOHUM" in accs:
        acc_list = accs.split(";")
        new_list = []
        for acc in acc_list:
            if "PGOHUM" in acc:
                new_list.append(acc)
                loci_dic[loci].append(acc)
                
        row[1]= ";".join(new_list)
    else:
        acc_list = accs.split(";")
        for acc in acc_list:
            loci_dic[loci].append(acc)
     
    
    output1.write("\t".join(row))
    
input1.close()
output1.close()

id_dic={}
for loci in loci_dic:
    acc_list = loci_dic[loci]
    acc_count = {}
    for acc in set(acc_list):
        acc_count[acc] = acc_list.count(acc)
        
    max_count = max(acc_count.values())
    ids = [key.split("_RF")[0] for key,val in acc_count.items() if val == max_count]
    id_dic[loci] = ";".join(set(ids))

temp = open(sys.argv[2],"r")

header = temp.readline().split("\t")
header.insert(1,"Sequence")
header.insert(2,"Protein")
header.insert(3,"GeneName")
header.insert(4,"Description")

output2.write("\t".join(header))

handle = open("mart_export.txt","r")
enst_dic = {}
annot_dic = {}
for line in handle:
    row=line.strip().split("\t")
    if len(row)>=5:
        enst_dic[row[1]] = row[3]
        annot_dic[row[1]] = row[4]
        
for line in temp:
    row = line.split("\t")
    match_snp = row[6]
    blastp_res = row[7]
    genelist=["NA"]
    annotlist=["NA"]
    if match_snp!="No match in SNP-DB":
        continue;
    if blastp_res=="match to known protein":
        if row[8] not in ["sp|Q9UN81|LORF1_HUMAN","sp|O00370|LORF2_HUMAN"]:
            continue;
    loci = row[loci_col]
    id = id_dic[loci]
    
    try:
        if "PGOHUM" in id:
            genelist=set([enst_dic.get(a.split("_")[1].split(".")[0],id) for a in id.split(";")])
            annotlist=set([annot_dic.get(a.split("_")[1].split(".")[0],id) for a in id.split(";")])
        else:
            genelist=[enst_dic.get(row[8].split("|")[1].split(".")[0],id)]
            annotlist=[annot_dic.get(row[8].split("|")[1].split(".")[0],id)]
    except IndexError:
        print (id)
    seq = re.sub("[\W\d]","",row[0].strip())
    row.insert(1,seq)
    row.insert(2,id)
    row.insert(3,";".join(genelist))
    if row[11] in ["sp|Q9UN81|LORF1_HUMAN","sp|O00370|LORF2_HUMAN"]:
        row.insert(4,row[11])
    else:
        row.insert(4,";".join(annotlist))
    output2.write("\t".join(row))
    
output2.close()