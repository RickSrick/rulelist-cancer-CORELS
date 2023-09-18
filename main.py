from Test import *
import NotificationManager

HOTNET2DIR = "data/cancer_type.txt"
SNVS_STRICTLY_FILTERED = "data/snvs_strictly_filtered.tsv"

files = [HOTNET2DIR, SNVS_STRICTLY_FILTERED]

primo = Test(files, "LUAD", "LUSC")
secondo = Test(files, "OV", "COADREAD")

test_type = 2

#10000, 100000, 1000000, 10000000
nodes = [10000000]

#["objective", "lower_bound", "curious", "bfs", "dfs"]
policies = ["bfs"]
#0.01, 0.1, 0.3, 0.4
min_supp = [0.1]

num_test = 3

for node in nodes:
    for policy in policies:
        for min in min_supp:
            for v in range(num_test):
                print(v,":\t",node,"-",policy,"-",min)
                if test_type == 2:
                    secondo.launch_test(0.01, node, policy, min)
                else:
                    primo.launch_test(0.05, node, policy, min)

NotificationManager.notify.send_markdown_text('''*FINE TEST*''')