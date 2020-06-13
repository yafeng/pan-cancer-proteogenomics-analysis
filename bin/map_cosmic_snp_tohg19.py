#!/usr/bin/env python3

import sys
import getopt
import os
import re

if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
    print("Warning! wrong command, please read the mannual in Readme.txt.")
    print("Example: python map_cosmic_snp_tohg19.py --input saav_peptable.txt --output saav_pep.coords.txt --cosmic_input CosmicMutantExport.tsv --dbsnp_input snp142CodingDbSnp.txt")
else:
    options, remainder = getopt.getopt(sys.argv[1:],'', ['input=','output=','cosmic_input=','dbsnp_input='])
    for opt, arg in options:
        if opt == '--input': input_file=arg
        elif opt == '--output': output_file=arg
        elif opt == '--cosmic_input': cosmic_file=arg
        elif opt == '--dbsnp_input': dbsnp_file=arg
        else:
            print("Warning! Command-line argument: %s not recognized. Exiting..." % opt); sys.exit()

input=open(input_file,"r")
output=open(output_file,"w")

cosmic_input=open(cosmic_file,"r") # Download at \http://cancer.sanger.ac.uk/cosmic/download
dbsnp_input=open(dbsnp_file,"r") # Download at UCSC table browser, select Variation as group, and snp142CodingDbSnp as table

header = ["#CHROM","POS","ID","INFO"]
output.write("\t".join(header)+"\n")

cosmic_dic={}

print ("reading cosmic database")
header = cosmic_input.readline().strip().split('\t')
geneix = header.index('Gene name')
mrnaix = header.index('Accession Number')
mutdesc = header.index('Mutation Description')
mutcds = header.index('Mutation CDS')
mutaa = header.index('Mutation AA')
mutgenpos = [ix for ix, x in enumerate(header) if re.match('Mutation .*genome position', x)][0]
mutstr = [ix for ix, x in enumerate(header) if re.match('Mutation .*strand', x)][0]

for line in cosmic_input:
    row=line.strip().split("\t")
    if "coding silent" in row[mutdesc]: # correspond to the column "Mutation Description" in the cosmic_input
        continue;
    elif row[mutdesc]=="Unknown":
        continue;
    elif row[mutdesc]=="Whole gene deletion":
        continue;
    else:
        gene=row[geneix]
        mRNA=row[mrnaix]
        dna_mut=row[mutcds] # correspond to the column "Mutation CDS" in the cosmic_input
        aa_mut=row[mutaa] # correspond to the column "Mutation AA" in the cosmic_input
        chr_position=row[mutgenpos] # correspond to the column "Mutation genome position" in the cosmic_input
        chr_strand=row[mutstr] # correspond to the column "Mutation strand" in the cosmic_input
        id="%s:%s:%s" % (mRNA,gene,aa_mut)
        if id not in cosmic_dic:
            cosmic_dic[id]=chr_position

print ("reading cosmic database finished")

print ("reading dbsnp database")
snp_dic={}
cor_dic={}
for line in dbsnp_input:
    row=line.strip().split("\t")
    try:
        rsid=row[4]
    except IndexError:
        print ("no rsid in this line")
        print (line)
    if rsid not in snp_dic:
        snp_dic[rsid]= [row[1],row[3]]
        cor_dic[row[1]+":"+row[3]]=rsid

print ("reading dbsnp database finished")

output_dic={}
output_chr_dic={}

cosmic_inDBsnp_count=0
input.readline()
for line in input:
    row=line.strip().split("\t")
    chr=""
    start=""
    end=""
    try:
        #score=row[12]
        #PSMcount=row[-1]
        pep=re.sub("[\W\d]","",row[0].strip())
        spectrumai_result=row[-1]
        proid = row[1]
    except IndexError:
        print ("INDEX ERROR",line)
        continue;

    acclist=proid.split(";")
    acc=""
    for proid in acclist:
        if "CanProVar_rs" in proid:
            acc=proid
            snpid=acc.split("_")[1]
            if snpid in snp_dic:
                chr=snp_dic[snpid][0]
                start=snp_dic[snpid][1]
                entry="%s\t%s\tID=dbSNP:%s\tSequence=%s;SpectrumAI_result=%s\n" % (chr,start,snpid,pep,spectrumai_result)
                output.write(entry)
            else:
                print (snpid,"not found in dbSNP")
            break;
        elif "COSMIC" in proid:
            acc=proid
            acc_list=acc.split(":")[1:4]
            cosmic_id=":".join(acc_list)
            try:
                chr_position=cosmic_dic[cosmic_id]
            except KeyError:
                print('Entry not found in COSMIC db, skipping {}'.format(cosmic_id))
                continue
            try:
                index1=chr_position.index(":")
                index2=chr_position.index("-")
                chr="chr"+chr_position[:index1].replace("23","X").replace("24","Y") # some of COSMIC entries used chr23 instead of chrX

                chr_start=chr_position[index1+1:index2]
                entry="%s\t%s\tID=COSMIC:%s\tSequence=%s;SpectrumAI_result=%s\n" % (chr,chr_start,cosmic_id,pep,spectrumai_result)
                output.write(entry)
                break;
            except ValueError:
                print ("ValueError",cosmic_id,chr_position)
        else:
            print ("protein id %s are not dnSNP or COSMIC entries, skipped" % proid)

cosmic_input.close()
dbsnp_input.close()
input.close()
output.close()

