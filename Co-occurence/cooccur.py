#!/bin/env python3

import sys
import pandas as pd
import numpy as np
import scipy.stats

amr=sys.argv[1]
gen=sys.argv[2]
out=sys.argv[3]

df=pd.read_csv(amr,index_col=0).T
df2=pd.read_csv(gen,index_col=0).T

amr_cols=list(df.keys())
gen_cols=list(df2.keys())

d={}
i=0
with open(out,"w") as coeff:
    coeff.write("var1,var2,coef,pvalue\n")
    for c in amr_cols:
        x=df[c]
        for co in gen_cols:
            y=df2[co]
            cf,p=scipy.stats.spearmanr(x,y)
            if p < 0.01 and cf >0.8:
                coeff.write(f"{c},{co},{cf},{p}\n")
