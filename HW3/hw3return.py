#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 22 17:03:47 2017

@author: mengranzhang
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 13:17:30 2017

@author: zangz
"""

import pandas as pd
import numpy as np
import sys

if len(sys.argv) < 6:
    print('Usage: HW3.py pricefile K N lam resultname')
    exit(0)
	
df1 = pd.read_csv(sys.argv[1], header = None)
df = df1.iloc[:,0:250]
df.shape
dft = df.T
dfreturn = dft.pct_change()
dfre = dfreturn.iloc[1:249,:]

K = int(sys.argv[2])
N = int(sys.argv[3])
lam = int(sys.argv[4])
res = open(sys.argv[5],"w")
res.write("n %d" % K + "\n")
res.write("\nj_lower_upper_mu\n")
res.write("\n")

L = 0.0
U = max(0.01, 2/float(K))

dnew = dfre.loc[:N-1,:K-1]
dnewp = dnew*100
dfinal = dnewp.T
m = np.mean(dfinal, axis = 1)
cov = np.cov(dfinal)


for i in range(0,K):
	res.write(str(i) + " " + str(L) + " " + str(U) + " " + str('{0:.6f}'.format(m[i])) + "\n")


res.write("\nlambda "+ str(lam) + "\n")

res.write("\ncovariance\n")
res.write("\n")

for i in range(0,K):
	for j in range(0,K):
		res.write(str('{0:.10f}'.format(cov[i,j]))+" ")
	res.write("\n")

res.write("\nEND")

res.close()

print("Done!")


