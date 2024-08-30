import pandas as pd

names=[]
temp=[]
#Column to group data by
col="class"
#Gene Coverage Percentage cutoff
cut=70

main=pd.read_csv("/home/michaelnetherland/Agnes_60_Cattle_Poultry/amrMatrix/amrMatrix.csv")
for s in main["Sample"].unique():
    uni=[]
    name=s

    df=main.loc[main["Sample"]==s]
    
    df2=df.loc[df["% Gene Coverage"]>cut]
    
    #If no hits are above 'cut' -> continue
    if df2.shape[0] == 0:
        names.append(name)
        continue
    
    #If there are duplicate hits -> get the highest coverage value
    ek=list(df2["product_name"].unique())
    
    for u in ek:
        sub=df2.loc[df2["product_name"]==u]
        sub_shape=sub.shape[0]
        if sub_shape > 1:
            tmp=sub.sort_values("% Gene Coverage",ascending=False)
            uni.append(tmp.iloc[[0]])
        else:
            uni.append(sub)
            
    full=pd.concat(uni,ignore_index=True)
    ########################################################
        
    sub_full=full[[col,"Bases Abundance"]]
    Subs=list(sub_full[col].unique())
    sub_dict={"name":name}
    
    if len(Subs)>1:
        for s in Subs:
            tmp_class=sub_full.loc[sub_full[col]==s]
            Sum=tmp_class["Bases Abundance"].sum()
            sub_dict[s]=[Sum]
            
        maked = pd.DataFrame(data=sub_dict)
    
        rep=pd.melt(maked, value_vars=maked.columns, id_vars="name")

    else:
        s=sub_full[col].unique()[0]
        Sum=sub_full["Bases Abundance"].sum()
        sub_dict[s]=[Sum]
            
        maked = pd.DataFrame(data=sub_dict)
        rep=pd.melt(maked, value_vars=s, id_vars="name")
    
    temp.append(rep)
    
df4=pd.concat(temp,ignore_index=True)




zero={}
shape=df4.shape[0]
for c,i in enumerate(range(shape,shape+len(names))):
    
    zero[i]=[names[c],"GLYCOPEPTIDE",0]

Zero=pd.DataFrame.from_dict(zero, orient='index',columns=['name','variable','value'])

df5=pd.concat([df4,Zero],ignore_index=True)

df5.to_csv(f"amr_{col}_melt_{cut}.csv",index=False)

df6=df5.pivot(index='variable', columns='name', values='value').fillna(0)

df6.to_csv(f"amr_{col}_base_abundance_{cut}.csv")
