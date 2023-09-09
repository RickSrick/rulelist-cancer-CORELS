import csv 
import numpy as np
import collections
import math

file_1, file_2 = "data/ov_mutation_frequency_expr_filtered.txt","data/coadread_mutation_frequency_expr_filtered.txt"

rates_1, rates_2 = {}, {}

final_rate = {}

def distance(x,y):
    return math.sqrt(pow(x-y,2))

with open(file_1) as file:
    tsv_file = csv.reader(file, delimiter="\t")

    for line in tsv_file:
        rates_1[line[0]] = (float) (line[1])
    
file.close()

with open(file_2) as file:
    tsv_file = csv.reader(file, delimiter="\t")

    for line in tsv_file:
        rates_2[line[0]] = (float) (line[1])

#print(rates_1, "\n", len(rates_1))
#print(rates_2, "\n", len(rates_2))

for v in rates_1:
    if v in rates_2:
        final_rate[v] = distance(rates_1[v], rates_2[v])
    else:
        final_rate[v] = distance(rates_1[v], 0)

for x in rates_2:
    if x not in final_rate:
        final_rate[x] = distance(rates_2[x], 0)

op = dict(sorted(final_rate.items(), key=lambda item: item[1]))

file = open("tests/output" + ".txt", "a")

print(op)

for j in op:
    file.write(f"{j} -> {op[j]}\n")

#, "PCDH9" NON Si TROVA
my_rule_svn = ["APC", "CHD8", "KRAS", "C11orf30", "TP53"]

for i in my_rule_svn:
    print(f"{i} -> {op[i]}({rates_1[i] if i in rates_1 else 0} - {rates_2[i] if i in rates_2 else 0}) \t (OV: {i in rates_1}, COADRED: {i in rates_2})\n")