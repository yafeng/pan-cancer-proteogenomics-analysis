#!/usr/bin/env python3

import sys
handle = open(sys.argv[1],"r")
output = open(sys.argv[2],"w")

header = handle.readline().split("\t")
output.write("\t".join(header))

pep_col = header.index("Peptide")

for line in handle:
    row = line.split("\t")
    row[pep_col] = row[pep_col][2:-2]
    output.write("\t".join(row))

output.close()
handle.close()

