import sys
import os
import pandas as pd
import re
import numpy as np
import scipy as sp
from collections import defaultdict
from itertools import islice
from scipy.stats import ttest_ind
os.environ['KMP_DUPLICATE_LIB_OK']='True'


# In[2]:

class LogFoldScalePlugin:
 def input(self, inputfile):
  self.expfile = inputfile#'input/high-vs-low_metastatic_lines_GSE59857.txt'

 def run(self):
     pass

 def output(self, outputfile):
  exp = pd.read_csv(self.expfile, sep='\t', header=0, index_col=0)
  exp.head()
  #Poorly metastatic: CACO2,COLO201,LS123,SW480,SW1417
  #Highly metastatic: LS174T,COLO320,HCT116,HT29,WIDR,LOVO
  logFC = pd.DataFrame(np.log2(exp.iloc[:,6:12].mean(axis=1) / exp.iloc[:,1:6].mean(axis=1)), columns=['logFC'])
  logFC['pval'] = exp.apply(lambda x: ttest_ind(x[1:6], x[6:12], equal_var=False)[1], axis=1)
  logFC.head()
  logFC['RefSeq'] = exp['RefSeq']
  logFC_r = logFC.groupby(['RefSeq']).agg(np.mean)
  logFC_r.to_csv(outputfile, sep='\t', index=True, index_label='RefSeq')



