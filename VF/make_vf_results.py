import glob
import re
import plotly.graph_objects as go
import random
import pandas as pd

main=pd.read_csv("/home/michaelnetherland/Agnes_60_Cattle_Poultry/vfMatrix/vfMatrix.csv")
def rgb_to_hex(temp):
    return '#{:02x}{:02x}{:02x}'.format(temp[0], temp[1], temp[2])


def get_color():
    temp=[]
    for j in range(0,3):
        temp.append(random.randint(0, 255))

    return rgb_to_hex(temp)

base="/home/michaelnetherland/Agnes_60_Cattle_Poultry/test"
states=["Alabama","Tennessee"]
Cut=[40,70]
color2=["#8ce4e1","#309cea","#410e9b","#1a2be3","#f8e331","#000000","#ffcc33","#990000","#000066","#e31a1a"
        ,"#367cb5","#fb7e00","#494e84","#00a65e","#cc58d6","#483396","#ffe326","#14ebf1","#fe66a2","#ffa630"
        ,"#eaffc5","#6f989c","#1e87ff","#840030","#2500d3","#6001df","#9803f0","#f50588","#f49b08","#e34584"
        ,"#034d2f","#c84848","#42bf93","#7da11b","#ceceb7","#906ca4","#2a9d1d"]
Color={}
Sample_Color={}
for cut in Cut:
    check=True
    for s in states:
        total=[]
        above={}
        above2={}
        Set1=set()
        Set2=set()
        Samps=set()
        
        for samp in main["Sample"].unique():
            j=0
            df=main.loc[main["Sample"]==samp]
            #fil=g.split("/")[-1].split(".")[0][1:3]
            #samp=g.split("/")[-1].split(".")[0]
            #Samps.add(samp)
            #trigger=g.split("/")[-1]
            
            fil=samp[1:3]
            Samps.add(samp)
            trigger=samp

            if trigger[0] != s[0]:
                continue


            cov={}
            val={}
            #with open(g,"r") as globe:
            for x in df.index:
                #if re.search("Cov",x):
                    #continue
                line=list(df.loc[x])
                line.remove(samp)
                cover=line[2]
                abund=line[1]
                cat=line[-1]
                #Gene name
                gene=line[0].split()[1].strip("\(\)")
                #Ref/Host name
                name=line[0].split("[")[-1].strip("[]")



                #Parse gene name
                if len(name.split()) > 2:
                    test=name.split()[2]
                    if test == "subsp.":
                        name=" ".join(name.split()[:4])
                    elif re.search(":",test):
                        name=" ".join(name.split()[:3])
                    else:
                        name=" ".join(name.split()[:2])

                if name not in Color.keys():
                    Color[name]=get_color()

                subcat=cat.split(";")[-3]
                if subcat == "":
                    subcat=cat.split(";")[-4]
                if subcat == "":
                    subcat=cat.split(";")[0]

                if gene not in cov.keys():
                    cov[gene]=[cover]
                    val[gene]=[[name,fil,subcat,abund,samp]]
                else:
                    cov[gene].append(cover)
                    val[gene].append([name,fil,subcat,abund,samp])

            with open(f"{base}/vf_gene_base_abundance_{cut}_melt.csv","a") as abb, open(f"{base}/vf_cat_base_abundance_{cut}_melt.csv","a") as catt, open(f"{base}/vf_ref_base_abundance_{cut}_melt.csv","a") as reff:
                #Record per sample
                #filt.write("Gene,Category,Reference,Sample Type,Base Abundance,Coverage\n")

                #Write category and reference matrices for heatmaps
                if check:
                    catt.write("Category,Abundance,Sample\n")
                    reff.write("Reference,Abundance,Sample\n")
                    abb.write("Reference,Abundance,Sample\n")
                    check=False

                #Get max value from duplicate entries
                for key in cov.keys():
                    maxx=0
                    count=0
                    for c,v in enumerate(cov[key]):

                        if float(v) > maxx:
                            maxx=float(v)
                            count=c

                    if maxx > cut:
                        #Write duplication record
                        #if len(cov[key]) > 1:
                            #Sort=sorted(list(cov[key]),reverse=True)
                            #wline=",".join(Sort)
                            #dup.write(f"{samp},{key},{wline}\n")

                        string=""

                        #Entry to use
                        sublist=val[key][count]

                        #Define variables for writing
                        gene=key
                        ref=sublist[0]
                        typ=sublist[2]
                        ab=sublist[3]
                        sa=sublist[4]

                        #Get Sample type
                        for letter in sublist[1]:
                            if letter == "C":
                                string="Cattle"
                            elif letter == "P":
                                string = "Poultry"
                            elif letter == "W":
                                string+="-Water"
                            elif letter == "M":
                                string+="-Manure"
                            elif letter == "S":
                                string+="-Soil"

                        if string not in Sample_Color.keys():
                            Sample_Color[string]=get_color()

                        with open(f"{base}/sankey_record_key_{s}_{cut}.csv","a") as rec:
                            rec.write(f"{gene},{typ},{ref},{string}\n")

                        #filt.write(f"{gene},{typ},{ref},{string},{ab},{maxx}\n")
                        catt.write(f"{typ},{ab},{sa}\n")
                        reff.write(f"{ref},{ab},{sa}\n")
                        abb.write(f"{gene},{ab},{sa}\n")

                        #Create data structures for Sankey plotting
                        change=gene

                        Set1.add(change)
                        Set1.add(ref)
                        Set1.add(string)

                        if gene not in above.keys():
                            above[change]=[ref]
                        else:
                            above[change].append(ref)

                        if ref not in above.keys():
                            above[ref]=[string]
                        else:
                            above[ref].append(string)
                            
        node_color=[]                    

        name_list=list(Set1)
        #print(len(name_list))

        #Write node colors
        for n in name_list:
            if n in Color.keys():
                node_color.append(Color[n])
            elif n in Sample_Color.keys():
                node_color.append(Sample_Color[n])
            else:
                node_color.append("blue")
        
        source=[]
        target=[]
        value=[]
        flow_color=[]

        for key in above.keys():
            temp={}
            for mem in above[key]:
                if mem not in temp.keys():
                    temp[mem]=1
                else:
                    temp[mem]+=1
            #print(key,temp)
            for k in temp.keys():
                sour=name_list.index(key)
                
                tar=name_list.index(k)
                
                source.append(sour)
                target.append(tar)
                value.append(temp[k])
                
                if key in Color.keys():
                    flow_color.append(Color[key])
                elif k in Color.keys():
                    flow_color.append(Color[k])


        fig = go.Figure(data=[go.Sankey(node = dict(pad = 15, thickness = 20, line = dict(color = "black", width = 0.5), label = name_list, color = node_color ), link = dict(source = source, target = target, value = value, color=flow_color))])

        fig.update_layout(title_text=f"Virulence Factors - {s} - {cut}% Coverage", font_size=10,autosize=False, width=1000,height=2500)

        fig.write_image(f"{base}/vf_by_category_{cut}_{s}_sankey_by_gene.pdf")
        fig.write_html(f"{base}/vf_by_category_{cut}_{s}_sankey_by_gene.html")

##### Make Final Dataframes #####
cats=glob.glob(f"{base}/*melt*")
for c in cats:
    fil="_".join(c.split("/")[-1].split("_")[:-1])
    #end=c.split("/")[-1].split("_")[-1]
    if re.search("cat",c):
        df=pd.read_csv(c)
        df2=df.pivot_table(index="Sample",columns="Category",values="Abundance",aggfunc="sum").fillna(0)
    else:
        df=pd.read_csv(c)
        df2=df.pivot_table(index="Sample",columns="Reference",values="Abundance",aggfunc="sum").fillna(0)
        
    cols=list(df2.columns)
    idx2=list(df2.index)
    
    length=len(list(df2.iloc[0]))
    
    df2.to_csv(f"{base}/{fil}.csv")
