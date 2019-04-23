#!/usr/bin/env python
############################################
## fastqDS.py
## 
#############################################
import sys
import random
import os

if len(sys.argv)<4:
    print("usage: fastqDS.py <path/to/input.fastq.gz> <percentage> <path/to/output.fastq.gz>")
    print("e.g.: fastqDS.py mouse_ChIP_H3K27ac_input.fastq.gz 60 mouse_ChIP_H3K27ac_input_ds.fastq.gz")

srcFastq = sys.argv[1]
percent = int(sys.argv[2])
dstFastq = sys.argv[3]

srcTmp = srcFastq.replace("fastq.gz","fastq")
dstTmp = dstFastq.replace("fastq.gz","fastq")

#os.system("gunzip -c "+srcFastq+" > "+srcTmp)

with open(srcTmp,"r") as src:
    with open(dstTmp, "w") as output:
        for line1 in src:
            line2 = src.readline()
            line3 = src.readline()
            line4 = src.readline()
            if random.randrange(1,101) <= percent:
                output.write(line1)
                output.write(line2)
                output.write(line3)
                output.write(line4)
os.system("gzip "+dstTmp)
os.system("rm "+srcTmp)