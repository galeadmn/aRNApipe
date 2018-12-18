import optparse
import time
import os
import config as config
import refbuild as refbuild

res = refbuild.annotate_gtf("/share/genomes_processed/Macaca_mulatta/genesets.gtf")
print res
