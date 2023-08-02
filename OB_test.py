from Test import *

HOTNET2DIR = "data/cancer_type.txt"
SNVS_STRICTLY_FILTERED = "data/snvs_strictly_filtered.tsv"
files = [HOTNET2DIR, SNVS_STRICTLY_FILTERED]

primo = Test(files, "LUAD", "LUSC")
secondo = Test(files, "OV", "COADREAD")

primo.launch_test(0.05, 10000, "curious", 0.01)