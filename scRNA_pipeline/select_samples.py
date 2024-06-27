import os, sys
import random
import scanpy as sc
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

inputdir = sys.argv[1]
adatas = [inputdir + '/' + x for x in os.listdir(inputdir) if x.endswith('filtered.h5')]

refdis = sys.argv[2]
info = pd.read_csv(sys.argv[2], sep='\t')

outputdir = sys.argv[3]
os.makedirs(outputdir, exist_ok=True)

info.set_index('Sample_ID', inplace = True)

adatas = random.sample(adatas, 90)
print(' '.join(adatas))

ages = {'age':[int(info.loc[x.split('/')[-1].split('_')[0], 'Age']) for x in adatas]}
print(len([x for x in ages['age'] if int(x) < 50]))
print(len([x for x in ages['age'] if int(x) > 70]))

ages = pd.DataFrame(ages)
sns.displot(ages, x='age', binwidth=1)
plt.savefig('SCA_age_dist.png')
plt.clf()
