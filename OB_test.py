from Test import *
from datetime import * 
import NotificationManager

HOTNET2DIR = "data/cancer_type.txt"
SNVS_STRICTLY_FILTERED = "data/snvs_strictly_filtered.tsv"

files = [HOTNET2DIR, SNVS_STRICTLY_FILTERED]

primo = Test(files, "LUAD", "LUSC")
secondo = Test(files, "OV", "COADREAD")

NotificationManager.notify.send_markdown_text('''*INIZIO TEST*''')

secondo.launch_test(0.01, 10000000, "curious", 0.01)
secondo.launch_test(0.01, 10000000, "lower_bound", 0.01)

NotificationManager.notify.send_markdown_text('''*FINE TEST*''')