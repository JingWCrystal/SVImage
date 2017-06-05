#coding:utf-8 

import pysam
from pysam import VariantFile
import os  
import math
import numpy as np
import gc

zero_file = open("man_made_0.txt", "a")
for sv in open('/mnt/hde/gao/wj/simulate/testSimuResult1/bp.txt'):  
	sv=sv.strip('\n')
	tmp=sv.split(' ')
	bp1=int(tmp[0])+10000
	bp2=int(tmp[1])+10000
	zero_file.write(str(bp1)+" "+str(bp2)+'\n')
zero_file.close()


gc.collect()




