from Test import *
import NotificationManager
import pstats
import cProfile

HOTNET2DIR = "data/cancer_type.txt"
SNVS_STRICTLY_FILTERED = "data/snvs_strictly_filtered.tsv"

files = [HOTNET2DIR, SNVS_STRICTLY_FILTERED]

primo = Test(files, "LUAD", "LUSC")
secondo = Test(files, "OV", "COADREAD")

NotificationManager.notify.send_markdown_text('''*INIZIO TEST*''')

#with cProfile.Profile() as pr:
secondo.launch_test(0.01, 100000, "objective", 0.01)
secondo.launch_test(0.01, 100000, "lower_bound", 0.01)
secondo.launch_test(0.01, 100000, "curious", 0.01)
secondo.launch_test(0.01, 10000000, "objective", 0.01)
secondo.launch_test(0.01, 10000000, "lower_bound", 0.01)
secondo.launch_test(0.01, 10000000, "curious", 0.01)

primo.launch_test(0.05, 100000, "objective", 0.01)
primo.launch_test(0.05, 100000, "lower_bound", 0.01)
primo.launch_test(0.05, 100000, "curious", 0.01)
primo.launch_test(0.05, 10000000, "objective", 0.01)
primo.launch_test(0.05, 10000000, "lower_bound", 0.01)
primo.launch_test(0.05, 10000000, "curious", 0.01)

secondo.launch_test(0.01, 100000, "objective", 0.03)
secondo.launch_test(0.01, 100000, "lower_bound", 0.03)
secondo.launch_test(0.01, 100000, "curious", 0.03)
secondo.launch_test(0.01, 10000000, "objective", 0.03)
secondo.launch_test(0.01, 10000000, "lower_bound", 0.03)
secondo.launch_test(0.01, 10000000, "curious", 0.03)

primo.launch_test(0.05, 100000, "objective", 0.03)
primo.launch_test(0.05, 100000, "lower_bound", 0.03)
primo.launch_test(0.05, 100000, "curious", 0.03)
primo.launch_test(0.05, 10000000, "objective", 0.03)
primo.launch_test(0.05, 10000000, "lower_bound", 0.03)
primo.launch_test(0.05, 10000000, "curious", 0.03)

NotificationManager.notify.send_markdown_text('''*FINE TEST*''')