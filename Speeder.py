import numpy as np
import pandas as pd
import csv

SNVS_STRICTLY_FILTERED = "data/snvs_strictly_filtered.tsv"
HOTNET2DIR = "data/cancer_type.txt"

cancer_type_one, cancer_type_two = "OV", "COADREAD"

hotnet2_table = pd.read_csv(HOTNET2DIR, sep=" ", names=["Patient", "CancerType"])
hotnet2_table = hotnet2_table[hotnet2_table["CancerType"].isin([cancer_type_one, cancer_type_two])]
hotnet2_table.CancerType[hotnet2_table.CancerType == cancer_type_one] = 1
hotnet2_table.CancerType[hotnet2_table.CancerType == cancer_type_two] = 0

mutations = {}
patients = []
line_num = 0
bound = 3
with open(SNVS_STRICTLY_FILTERED) as file:
    tsv_file = csv.reader(file, delimiter="\t")

    for line in filter(lambda x: x[0] in hotnet2_table['Patient'].values,tsv_file):
        patients.append(line[0])
        for x in line[1:]:
            if x in mutations:
                mutations[x].append(line_num)
            else:
                mutations[x] = [line_num]
        line_num += 1

print(len(patients), len(mutations))

binary_data = np.zeros(shape=(len(hotnet2_table.index), len(mutations)), dtype=np.uint8)
mutation_num = 0

for mut in mutations:
    for ind in mutations[mut]:
        binary_data[ind, mutation_num] = 1
    mutation_num+= 1

bin_mutation = pd.DataFrame(binary_data, columns=mutations).assign(Patient=patients)
bin_mutation = hotnet2_table.merge(bin_mutation, on="Patient",how="left")