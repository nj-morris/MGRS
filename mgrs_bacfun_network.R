library(devtools)
# install_github("zdk123/SpiecEasi") # check gcc, gfortran (install w homebrew)
library(SpiecEasi)
library(phyloseq)
library(igraph)
#devtools::install_github("GraceYoon/SPRING")
# devtools::install_github("stefpeschel/NetCoMi", 
#                         repos = c("https://cloud.r-project.org/",
#                                   BiocManager::repositories()))
library(NetCoMi)
library(microbiome)

ps_bac <- readRDS("~/Library/CloudStorage/OneDrive-UniversityofArizona/Rotation_III/MGRS/physeq_MGRS_16S_unrare_Aug2025.rds")
ps_fun <- readRDS("~/Library/CloudStorage/OneDrive-UniversityofArizona/Rotation_III/MGRS/ps_fungi_unrare_difab.RDS")

# subset to only mgrs fungi/bacteria
ps_bac <- subset_samples(ps_bac, Site == "MG") # 
ps_bac <- subset_samples(ps_bac, !Fecal.ID == "20") # 
ps_fun <- subset_samples(ps_fun, Site == "MG") # 31 
sam_data(ps_fun)

# prune - make sure no zeros, get rid of biologically irrelevant taxa
library(MicEco)
ps_bac <- ps_prune(ps_bac, min.samples = 2, min.reads = 10)
ps_fun <- ps_prune(ps_fun, min.samples = 2, min.reads = 10)
detach("package:MicEco", unload = TRUE) # then detach MicEco - masks clr in spieceasi

# glom to family/genus
ps_bac = tax_glom(ps_bac, taxrank = "Family") #glom by family
# Rename taxonomic table 
ps_bac <- renameTaxa(ps_bac, 
                                  pat = "<name>", 
                                  substPat = "<name>_<subst_name>(<subst_R>)",
                                  numDupli = "Family")
ps_fun = tax_glom(ps_fun, taxrank = "Family") #glom by family
# Rename taxonomic table 
ps_fun <- renameTaxa(ps_fun, 
                     pat = "<name>", 
                     substPat = "<name>_<subst_name>(<subst_R>)",
                     numDupli = "Family")

# filter by abundance
ps_bac<-transform_sample_counts(ps_bac, function(x) x / sum(x) ) #transform to relative abundance
ps_bac<-filter_taxa(ps_bac, function(x) sum(x) > .005, TRUE)

ps_fun<-transform_sample_counts(ps_fun, function(x) x / sum(x) ) #transform to relative abundance
ps_fun<-filter_taxa(ps_fun, function(x) sum(x) > .005, TRUE)


# cross domain interactions
bac_fun <- spiec.easi(list(ps_bac, ps_fun), method='mb', nlambda=40,
                      lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05))
saveRDS(bac_fun, file = "spieceasi_MGRS_family.RDS")

dtype <- c(rep(1,ntaxa(ps_bac)), rep(2,ntaxa(ps_fun)))
plot(adj2igraph(getRefit(bac_fun)), vertex.color=dtype+1, vertex.size=9)


### lauren code ####
# 
getStability(bac_fun) 
getOptMerge(bac_fun)#adjacency matrix with a number between 0-1, corresponding to the fraction of networks, constructed with 80% subsamples, that also contained that edge. 
getOptBeta(bac_fun)#that's the matrix of correlation coeff between taxa
summary(getOptBeta(bac_fun))

ps_combined <- merge_phyloseq(ps_bac, ps_fun)

spiec.graph=adj2igraph(getRefit(bac_fun), vertex.attr=list(name=taxa_names(ps_combined)))

ggnet0=plot_network(spiec.graph, ps_combined, type='taxa', point_size=2, line_color = NULL, color = 
                      "Kingdom", label = "Family")

#count positive edgs + extract reg coefficient matrix
betaMat=as.matrix(symBeta(getOptBeta(bac_fun))) 

#calculate edge lengths
positive=length(betaMat[betaMat>0])/2 
negative=length(betaMat[betaMat<0])/2 
total=length(betaMat[betaMat!=0])/2 

#clustering
clusters=cluster_fast_greedy(spiec.graph)
clusterOneIndices=which(clusters$membership==1)
clusterOneOtus=clusters$names[clusterOneIndices]
clusterTwoIndices=which(clusters$membership==2)
clusterTwoOtus=clusters$names[clusterTwoIndices]
clusterThreeIndices=which(clusters$membership==3)
clusterThreeOtus=clusters$names[clusterThreeIndices]

#summarize clusters
# Extract family names from taxonomy tables
taxa.b <- as.character(tax_table(ps_bac)[, "Family"])
taxa.f <- as.character(tax_table(ps_fun)[, "Family"])
taxa.bf <- c(taxa.b, taxa.f)

#plot with edge colors
otu.ids=colnames(bac_fun$data)
edges=E(spiec.graph)
edge.colors=c()
for(e.index in 1:length(edges)){
  adj.nodes=ends(spiec.graph,edges[e.index])
  xindex=which(otu.ids==adj.nodes[1])
  yindex=which(otu.ids==adj.nodes[2])
  beta=betaMat[xindex,yindex]
  if(beta>0.0000000000000001){
    edge.colors=append(edge.colors,"forestgreen")
  }else if(beta<-0.000000000001){
    edge.colors=append(edge.colors,"red")
  }
}

E(spiec.graph)$color=edge.colors

spiec.graph.b=spiec.graph
nodenames=V(spiec.graph.b)$name
V(spiec.graph.b)$name=getTaxonomy(nodenames, taxa.f, useRownames=TRUE)
E(spiec.graph.b)$arrow.size=5
V(spiec.graph.b)$color="white"
V(spiec.graph.b)$frame.color="black"
tkplot(spiec.graph.b)



# NetCoMi
assoMat1 <- SpiecEasi::symBeta(SpiecEasi::getOptBeta(bac_fun), mode = "ave")
assoMat1 <- as.matrix(assoMat1)

# Get taxa names
taxa.b <- as.character(tax_table(ps_bac)[, "Family"])
taxa.f <- as.character(tax_table(ps_fun)[, "Family"])
taxnames <- c(taxa.b, taxa.f)

colnames(assoMat1) <- rownames(assoMat1) <- taxnames
diag(assoMat1) <- 1

# Network construction (pass association matrices to netConstruct)
# - sparsMethod must be set to "none" because sparsification is already included in SpiecEasi
net_bacfun <- netConstruct(data = assoMat1, 
                                 dataType = "condDependence", 
                                 sparsMethod = "none")

# Network analysis
net1_bacfun <- netAnalyze(net_bacfun, hubPar = "eigenvector")


nodeCols <- c(rep("blue", ntaxa(ps_bac)), rep("lavender", ntaxa(ps_fun)))
names(nodeCols) <- taxnames

plot(net1_bacfun, 
     sameLayout = TRUE, 
     layoutGroup = "union",
     nodeColor = "colorVec", 
     colorVec = nodeCols,
     nodeSize = "eigen", 
     nodeSizeSpread = 1.5,
     labelScale = FALSE,
     cexNodes = 2, 
     cexLabels = 0.5,
     cexHubLabels = 0.7,
     cexTitle = 1,
     )

plotHeat(net_bacfun$assoMat1, textUpp = "none", textLow = "none")

legend(-0.2, 1.2, cex = 3, pt.cex = 4, 
       legend = c("16S", "Fungi"), col = c("lightblue", "orange"), 
       bty = "n", pch = 16) 
####### bacteria only ######
spiec.out=spiec.easi(ps_bac, method="mb",icov.select.params=list(rep.num=50),lambda.min.ratio=0.01) #for timepoint 0, use lambda.min.ratio = 0.1 to get stable network, for 4/6, use 0.01
#rep num = 50 (better to do 999 for pub but network generates the same)

getStability(spiec.out) 
getOptMerge(spiec.out)#adjacency matrix with a number between 0-1, corresponding to the fraction of networks, constructed with 80% subsamples, that also contained that edge. 
getOptBeta(spiec.out)#that's the matrix of correlation coeff between taxa
summary(getOptBeta(spiec.out))

spiec.graph=adj2igraph(getRefit(spiec.out), vertex.attr=list(name=taxa_names(ps_bac)))

ggnet0=plot_network(spiec.graph, ps_bac, type='taxa', label= "Family", point_size=1.5, line_color = NULL)
ggnet0
# NetCoMi
assoMat1 <- SpiecEasi::symBeta(SpiecEasi::getOptBeta(spiec.out), mode = "ave")
assoMat1 <- as.matrix(assoMat1)

# Get taxa names
taxa.b <- as.character(tax_table(ps_bac)[, "Family"])
taxa.f <- as.character(tax_table(ps_fun)[, "Family"])
taxnames <- c(taxa.b, taxa.f)

colnames(assoMat1) <- rownames(assoMat1) <- taxa.b
diag(assoMat1) <- 1

# Network construction (pass association matrices to netConstruct)
# - sparsMethod must be set to "none" because sparsification is already included in SpiecEasi
net_bac <- netConstruct(data = assoMat1, 
                           dataType = "condDependence", 
                           sparsMethod = "none")

# Network analysis
net1_bac <- netAnalyze(net_bac, hubPar = "eigenvector")


nodeCols <- rep("blue", ntaxa(ps_bac))
names(nodeCols) <- taxa.b

plot(net1_bac, 
     sameLayout = TRUE, 
     layoutGroup = "union",
     nodeColor = "colorVec", 
     colorVec = nodeCols,
     nodeSize = "eigen", 
     nodeSizeSpread = 1.5,
     labelScale = FALSE,
     cexNodes = 2, 
     cexLabels = 0.7,
     cexHubLabels = 1,
     cexTitle = 1,
)

plotHeat(net_bac$assoMat1, textUpp = "none", textLow = "none")

##### fungi only ######

spiec.out=spiec.easi(ps_fun, method="mb",icov.select.params=list(rep.num=50),lambda.min.ratio=0.01) #for timepoint 0, use lambda.min.ratio = 0.1 to get stable network, for 4/6, use 0.01
#rep num = 50 (better to do 999 for pub but network generates the same)

getStability(spiec.out) 
getOptMerge(spiec.out)#adjacency matrix with a number between 0-1, corresponding to the fraction of networks, constructed with 80% subsamples, that also contained that edge. 
getOptBeta(spiec.out)#that's the matrix of correlation coeff between taxa
summary(getOptBeta(spiec.out))

spiec.graph=adj2igraph(getRefit(spiec.out), vertex.attr=list(name=taxa_names(ps_fun)))

ggnet0=plot_network(spiec.graph, ps_fun, type='taxa', label= "Family", point_size=1.5, line_color = NULL)
ggnet0
# NetCoMi
assoMat1 <- SpiecEasi::symBeta(SpiecEasi::getOptBeta(spiec.out), mode = "ave")
assoMat1 <- as.matrix(assoMat1)

# Get taxa names
taxa.f <- as.character(tax_table(ps_fun)[, "Family"])

colnames(assoMat1) <- rownames(assoMat1) <- taxa.f
diag(assoMat1) <- 1

# Network construction (pass association matrices to netConstruct)
# - sparsMethod must be set to "none" because sparsification is already included in SpiecEasi
net_fun <- netConstruct(data = assoMat1, 
                        dataType = "condDependence", 
                        sparsMethod = "none")

# Network analysis
net1_fun <- netAnalyze(net_fun, hubPar = "eigenvector")


nodeCols <- rep("lavender", ntaxa(ps_fun))
names(nodeCols) <- taxa.f

plot(net1_fun, 
     sameLayout = TRUE, 
     layoutGroup = "union",
     nodeColor = "colorVec", 
     colorVec = nodeCols,
     nodeSize = "eigen", 
     nodeSizeSpread = 1.5,
     labelScale = FALSE,
     cexNodes = 2, 
     cexLabels = 0.7,
     cexHubLabels = 1,
     cexTitle = 1,
)

plotHeat(net_fun$assoMat1, textUpp = "none", textLow = "none")



# use phyloseq ! 
m <- spiec.easi(ps, method='mb', lambda.min.ratio=1e-2,
                           nlambda=20, pulsar.params=list(rep.num=50))
ig2.mb <- adj2igraph(getRefit(se.mb.amgut2),  vertex.attr=list(name=taxa_names(amgut2.filt.phy)))
plot_network(ig2.mb, amgut2.filt.phy, type='taxa', color="Rank3")



#### tutorial ####
data(amgut1.filt)
depths <- rowSums(amgut1.filt)
amgut1.filt.n  <- t(apply(amgut1.filt, 1, norm_to_total))
amgut1.filt.cs <- round(amgut1.filt.n * min(depths))

d <- ncol(amgut1.filt.cs)
n <- nrow(amgut1.filt.cs)
e <- d

set.seed(10010)
graph <- make_graph('cluster', d, e)
Prec  <- graph2prec(graph)
Cor   <- cov2cor(prec2cov(Prec))

X <- synth_comm_from_counts(amgut1.filt.cs, mar=2, distr='zinegbin', Sigma=Cor, n=n)

se <- spiec.easi(X, method='mb', lambda.min.ratio=1e-2, nlambda=15)
# Applying data transformations...
# Selecting model with pulsar using stars...
# Fitting final estimate with mb...
# done

huge::huge.roc(se$est$path, graph, verbose=FALSE)
stars.pr(getOptMerge(se), graph, verbose=FALSE)
# stars selected final network under: getRefit(se)