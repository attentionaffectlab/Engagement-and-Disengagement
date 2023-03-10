#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import pandas as pd
from copy import deepcopy
from nltools.utils import get_resource_path
from nltools.file_reader import onsets_to_dm
from nltools.data import Design_Matrix
from nltools.mask import create_sphere
from nltools.mask import expand_mask, collapse_mask
from nltools.data import Adjacency
from scipy.spatial.distance import pdist
from scipy.stats import pearsonr
from scipy import stats
import seaborn as sns
from nltools.stats import correlation_permutation
from scipy.stats import spearmanr
from scipy.spatial import distance
from nltools.stats import matrix_permutation
from nltools.stats import fdr
from scipy.spatial.distance import squareform
from nltools.plotting import plot_brain
from nltools.data import Adjacency
from scipy.stats import rankdata
from mlxtend.evaluate import permutation_test
from dtaidistance import dtw
from dtaidistance import dtw_visualisation as dtwvis
import numpy as np
import csv

# define data repository
path_atlas = '/path/atlas/'
path_fmri = 'path/fmri_files'
output_path = '.../'
path_behavior = '.../'


# In[ ]:


#Templates
mask119 = Brain_Data(path_atlas,'Nsynth119.nii.gz')
allmask119 = expand_mask(mask119)


# In[ ]:


# # for movie watching condition
subs = np.zeros((45, 98, allmask119.shape()[0]))
with open(path_fmri,'fmri_files.txt', 'r') as subject_list:
 
 for s, name in enumerate(subject_list):
     name = name.rstrip("\n")
     dat = Brain_Data(name)
     
     print("Loading subject %s ... " % name)
     print("Extracting time-series")
    
     for i, ma in enumerate(allmask119):
         subs[s,:,i] = dat.extract_roi(ma)
        
     print("Converting to mats...")
     adj = Adjacency()
     for i in range(subs.shape[-1]):
         adj = adj.append(Adjacency(1 - pdist(subs[:,:,i],metric='correlation')))
     adj.write(os.path.join(path,'sad_ISC_119rois_22subj.csv'))


##Transpose output csv file before the next step


# ## Calculate and viusalize brain ISC

# In[ ]:


saddata_119rois = pd.read_csv(os.path.join(path,'sad_ISC_119rois_22subj_transposed.csv'))

sadmean = []

for i in list(range(1,119)):
    r = saddata_200rois['ROI%i' %(i)].mean()
    sadmean.append(r)
    
ISCdata = pd.DataFrame(sadmean)
ISCdata.columns = ['sad average ISC']

print('Average ISC for 21grams movie is:', ISCdata['sad average ISC'].mean())


# In[ ]:


sad_ISC_mean_value = ISCdata['sad average ISC'].values
sad_ISC_mean = mask119.copy()
sad_ISC_mean.data = np.array([x.data*y for (x,y) in zip(allmask119,sad_ISC_mean_value)])
sad_ISC_mean = sad_ISC_mean.sum()

plot_brain(sad_ISC_mean,**{'vmax':0.45}


# # RSA

# ## ISC behavior

# In[ ]:


# compute correlation distance
scale_normalize_raw = pd.read_csv(os.path.join(path,"binary_eng_22subj.csv"))
print(scale_normalize_raw.shape)
scale_normalize = pdist(scale_normalize_raw,metric='hamming')


# correlation: 1-correlation distance 
scale_normalize_r = 1 - scale_normalize
scale_normalize_cor = pd.DataFrame(scale_normalize_r,columns= ["scale_normalize_cor"])


# In[ ]:


#merge behavior and brain
data_119rois = pd.read_csv(os.path.join(path'sad_ISC_119rois_22subj_transposed.csv'))
data_119rois['scale_normalize_cor'] = scale_normalize_cor

scale_cor_square = distance.squareform(scale_normalize_cor['scale_normalize_cor'])

# visualize behavior ISC
mask =np.zeros_like(scale_cor_square)
mask[np.triu_indices_from(mask)] = True
sns.heatmap(scale_cor_square, vmin=.0,vmax=1, mask=mask, cmap = 'Reds')
print(scale_normalize_cor['scale_normalize_cor'])


# ## Mantel test

# In[ ]:


rvals200_1tail = []
pvals200_1tail = []

for i in range(0,119):
    r = matrix_permutation(data_119rois['ROI'+str(i+1)],data_119rois['scale_normalize_cor'],
        n_permute=10000, metric='spearman',tail=1, random_state=None, n_jobs=5)
    print("for ROI %s, the correlation is %s " % (i+1, r))
    rvals200_1tail.append(r['correlation'])
    pvals200_1tail.append(r['p'])


# In[ ]:


rvals200_1tail = pd.DataFrame(rvals200_1tail)
rvals200_1tail.columns=['rvals']
rvals200_1tail.to_csv(os.path.join(sophiepath,'RSA','last','rvals_eng_firstHalf.csv'))

pvals200_1tail = pd.DataFrame(pvals200_1tail)
pvals200_1tail.columns=['pvals']
pvals200_1tail.to_csv(os.path.join(sophiepath,'RSA','last','pvals_eng_firstHalf.csv'))


# In[ ]:


#apply FDR
rvals200_1tail = np.array(rvals200_1tail)
mask = rvals200_1tail < 0
p = np.array(pvals200_1tail)
p[mask] = 1 - p[mask]
fdr(p, q=0.05)

