#!/usr/bin/env python
# coding: utf-8

# # Assignment 1
# ### Richard Wang, Yuanfan Lai, Laurinda He and Kun Liu

# # Part 0: Preliminary Data Analysis

# In[1]:


import pandas as pd
import math


# In[2]:


read_data = pd.read_csv("AB_test_data.csv")


# In[3]:


read_data


# In[4]:


A = read_data[read_data['Variant']=='A']
B = read_data[read_data['Variant']=='B']


# In[5]:


A


# In[6]:


B


# ## Count number of variants

# In[7]:


n_A = A.shape[0]
n_B = B.shape[0]

print("{} {}".format('n_A:',n_A)) 
print("{} {}".format('n_B:',n_B))


# ## Calculate p_A and p_B

# In[9]:


p_A = A[A['purchase_TF'] == True].shape[0] / A.shape[0]
p_B = B[B['purchase_TF'] == True].shape[0] / B.shape[0]

print("{} {}".format('p_A:',p_A)) 
print("{} {}".format('p_B:',p_B))


# # Part 1: A/B test

# ### H0: p_B=p_A
# 
# ### Ha: p_B>p_A

# In[34]:


z = (p_B-p_A)/math.sqrt((p_B*(1-p_B)/n_B))

if z >= 1.64:
    print('z-score is %.4f.>1.64 , Reject null hypothesis. Alternative B improved conversion rates' % (z))
else: 
    print('z-score is %.4f.>1.64 , Fail to reject null hypothesis. Alternative B did not improve conversion rates'% (z))
  


# # Part 2: Calculate optimal sample size

# In[15]:


import scipy.stats as st


# In[17]:


alpha = 0.05
beta = 0.2
pbar = (p_A + p_B) / 2
delta = p_B - p_A

t_0025 = st.norm.ppf(.975)
t_02 = st.norm.ppf(.8)

#use formula given in class to calculate optimal sample size
n_optimal = ((t_0025 * math.sqrt(2*pbar*(1 - pbar)) + t_02*math.sqrt((p_A)*(1-p_A)+p_B*(1-p_B)))**2)/(delta**2)
n_opt = int(round(n_optimal))
n_opt


# <b>We choose the optimal sample size to be 2942.</b>

# # Part 3: Conduct A/B Test 10 times

# ### H0: p_hat=p_A
# 
# ### Ha: p_hat>p_A

# In[19]:


import math


# In[37]:


success = []
samples = []
n=1
#p_A=0.149616

#iterate 10 times to get conduct A/B test 10 times
while n < 11:    
    #randomly select optimal number of samples from variant B
    new_sample = B.sample(n=n_opt)
    
    #append samples of 2942 for each iteration
    samples.append(new_sample.reset_index())
    
    new_success = new_sample['purchase_TF'].sum()
    
    #calculate new p_hat every iteration
    p_hat = new_success/n_opt
    
    #z-score
    z = (p_hat-p_A)/math.sqrt((p_A*(1-p_A)/n_opt))
    
    if z >= 1.64:
        success.append(1)
    else: 
        success.append(0)
    n+=1

print("Out of 10 tests, {} showed significant difference to support that Variant B performs better in improving conversion rate.".format(sum(success)))


# <b>There is enough evidence to support that Variant B (with walkability assessment) is effective in improving conversion rate.</b>

# # Part 4: Sequential Testing on same 10 samples

# In[38]:


import numpy as np


# In[39]:


fA_1 = p_A  #0.149616
fA_0 = 1-p_A
fB_1 = p_B  #0.1766
fB_0 = 1-p_B

#p[timal sample size
n=n_opt

A_bound = np.log(1/alpha)
B_bound = np.log(beta)

test_results = []
iteration_length = []

#use same 10 samples from A/B testing
for sample in samples:
    i = 0
    iter_lambda = 0
    while i < n:
        if sample.purchase_TF[i] == True:
            iter_lambda = iter_lambda + np.log(fB_1 / fA_1)
        else:
            iter_lambda = iter_lambda + np.log(fB_0 / fA_0)
        
        if iter_lambda <= B_bound:
            iteration_length.append(i+1)
            test_results.append("Fail to reject H0, number of trials: {}".format(i+1))
            break
        elif iter_lambda >= A_bound:
            iteration_length.append(i+1)
            test_results.append("Reject H0, number of trials: {}".format(i+1))
            break
        else:
            i = i+1


# In[40]:


print(iteration_length)


# In[45]:


avg_length = sum(iteration_length)/len(iteration_length)
print(avg_length)


# In[48]:


test_results


# In[46]:


print("We were able to stop all 10 tests prior to using the full sample. The average number of iterations required to stop the test is %.0f." % (avg_length))

