library(vegan)

genus=read.csv("genus.normalize.counts.csv",row.names=1)
amr=read.csv("amr_class_full_base_abundance_40.csv",row.names=1)


genus.ord=cmdscale(vegdist(t(genus))
env=envfit(genus.ord, t(amr), perm = 999, scaling = "species")
