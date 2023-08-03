from Test import *
from datetime import * 
import NotificationManager

HOTNET2DIR = "data/cancer_type.txt"
SNVS_STRICTLY_FILTERED = "data/snvs_strictly_filtered.tsv"

files = [HOTNET2DIR, SNVS_STRICTLY_FILTERED]

primo = Test(files, "LUAD", "LUSC")
secondo = Test(files, "OV", "COADREAD")

primo.launch_test(0.05, 100000, "curious", 0.01)
primo.launch_test(0.05, 100000, "bfs", 0.1)

secondo.launch_test(0.05, 1000, "curious", 0.01)
secondo.launch_test(0.01, 10000000, "bfs", 0.1)