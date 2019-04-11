### Multiplot function
#http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
require(Hmisc)
#require(patchwork)
require(ggplot2)
theme_set(theme_bw())
require(reshape2)
require(SNPRelate)
require(gdsfmt)
require(stringr)
require(ggrepel)
require(data.table)
require(ape)
require(directlabels)
require(ggmap)
require(ggtree)
require(phytools)
require(geosphere)
require(gplots)
require(plotrix)
require(ggtree)
require(phangorn)
require(RColorBrewer)
require(vegan)
require(dismo)
require(rasterVis)
require(rgdal)
library(methods)
require(maptools)

tbb1 = c(7027135,7029293)
tbb2 = c(13433925,13436504)

draft='./data/'
setwd(draft)
draft2='./data/outputdir/' ## Directory where figures will be plotted

climate <- getData('worldclim', var='bio', res=2.5)

## Color list copied from by S. Martin twisst plotting script
cols = c(
  "#F0A3FF", #Amethyst
  "#0075DC", #Blue
  "#993F00", #Caramel
  "#4C005C", #Damson
  "#191919", #Ebony
  "#005C31", #Forest
  "#2BCE48", #Green
  "#FFCC99", #Honeydew
  "#808080", #Iron
  "#94FFB5", #Jade
  "#8F7C00", #Khaki
  "#9DCC00", #Lime
  "#C20088", #Mallow
  "#003380", #Navy
  "#FFA405", #Orpiment
  "#FFA8BB", #Pink
  "#426600", #Quagmire
  "#FF0010", #Red
  "#5EF1F2", #Sky
  "#00998F", #Turquoise
  "#E0FF66", #Uranium
  "#740AFF", #Violet
  "#990000", #Wine
  "#FFFF80", #Xanthin
  "#FFFF00", #Yellow
  "#FF5005" #Zinnia
)

###--- LatLon for re-ordering
df = read.table(file='./data/Latlon.txt',header=T)
#colnames(df)[1]="sample"
df$reg=substr(df$sample,1,3)
df$iso=substr(df$sample,5,7)
df$iso[df$reg=='FRG']='FRG'
df$iso[df$lon==-0.4648]='FRA.1' ##- 79
df$iso[df$lon==-0.9449]='FRA.2' ##- 64
df$iso[df$lon==2.8415]='FRA.4' ##-63
df$iso[df$lon==2.1486]='FRA.3' ##- 81
df$iso[df$iso=='WAL']='AUS.1'
df$iso[df$iso=='MAR']='MOR'
df$iso[df$iso=='NAP']='AUS.2' ## Isolate was shipped from Napoli, Italy but comes from Australia (also confirmed by genetic data).
df$iso[df$iso=='KOK']='STA.1'
df$iso[df$iso=='OND']='STA.2'
df$iso[df$iso=='WHR']='STA.3'
df$iso2=substr(df$sample,5,7)
write.table(df,file='recoded_isolates.txt',quote=F,row.names=F)

df$geo=''
df$geo[df$iso2 %in% c("FRA", "ACO", "MAR")]='Mediterranean.Area'
df$geo[df$iso2 %in% c("IND", "WAL","NAP")]='Oceania'
df$geo[df$iso2 %in% c("BRA","FRG")]='South.America'
df$geo[df$iso2 %in% c("WHR","KOK", "OND", "NAM","ZAI")]='Sub.Tropical.Africa'
df$geo[df$iso2 %in% c("CAP","STO","BEN")]='Western.Africa'
df$geo[df$reg=='FRA']='Mediterranean.Area'
df$geo[df$reg=='FRG']='South.America'
df$geo=factor(df$geo)

###-- Coverage Nuclear genome
cov = read.table(file='./data/filt_upd_metadata.txt')
cov$sample = sapply(str_split(cov$V1,"/"),function(x) x[length(x)-1])
cov$mcov = cov$V6
df = merge(df,cov[,c('mcov','sample')],by='sample')

###--- Supplementary table 1
#- mito cov
setwd("./data/MITO/")
coverage=read.table(file='cov_stats.txt',header=F,sep='\t')
colnames(coverage)=c('sample','total','mean',
                     'granular_third_quartile', 'granular_median', 'granular_first_quartile', 
                     'p_bases_above_15')
suptab1 = merge(df,coverage[,c("sample","mean")],by="sample")

#- Lane info
seq = read.table('./data/seq_metadata.out',sep="\t",fill=TRUE)
colnames(seq) = c("sample","Run information")
suptab1 = merge(suptab1,seq,by="sample")

write.csv(suptab1,file=paste0(draft,"supplementary_table1.csv"))

###--- Color definition
clcol = c('#5ab4ac','#ca0020','#d8b365','#252525','#8c510a')
clcol2 = c('#01665e','#ca0020','#d8b365','#252525','#8c510a')

chrom.colors = c(cols[1],cols[12],cols[2],cols[3],cols[20])
chromlist = c('1_Celeg_TT','2_Celeg_TT','3_Celeg_TT','4_Celeg_TT','5_Celeg_TT_axel')

##-- genome size
gensize = 45778363 + 47384193 + 43564237 + 51819793 + 48825595
nsnp = 23868644
nsnp/gensize

####---- MAP

hcon = unique(df[,c('lon','lat')])
h2 = unique(df[,c('lon','lat','geo','iso')])
h2$Cluster = h2$geo
h2$res = "Fenbendazole-resistant and ivermectin-susceptible"
h2$res[h2$iso=='STA.1'|h2$iso=='STA.3'|h2$iso=='AUS.1']="Fenbendazole- and ivermectin-resistant"
h2$res[h2$iso=='FRA.2'|h2$iso=='FRA.4'|h2$iso=='STO'|h2$iso=='ZAI']="Fenbendazole- and ivermectin-susceptible"
h2$res = factor(h2$res)
levels(h2$res) = levels(h2$res)[c(2,1,3)]
h2$Resistance=17
h2$Resistance[h2$res==levels(h2$res)[2]] = 15
h2$Resistance[h2$res==levels(h2$res)[3]] = 19
h2$Resistance = factor(h2$Resistance,labels=levels(h2$res))

h2$iso[h2$iso=='FRG'] = 'GUA'

pdf(file=paste0(draft2,'Figure1A.pdf'),width=14,height=8)
#Using GGPLOT, plot the Base World Map
mp <- NULL
mapWorld <- borders("world", colour="gray50", fill="white") # create a layer of borders
mp <- ggplot() +   mapWorld
#Now Layer the cities on top
mp <- mp+ geom_point(data=h2,aes(x=lon, y=lat ,color=Cluster,shape=Resistance), size=3) +
  xlab("Longitude") + ylab("Latitude") +
  geom_label_repel(data=h2,segment.alpha = 0.8,size=4,
                   aes(x=lon, y=lat, col=Cluster,label=iso),show.legend = FALSE) +
  scale_color_manual(values=clcol) + 
  theme(legend.position=c(.1,.2),text=element_text(size=14))
mp
invisible(dev.off())

####---- Associate climatic information with each isolate
# coords <- data.frame(x=hcon$lon,y=hcon$lat)
# points <- SpatialPoints(coords, proj4string = climate@crs)
# values <- extract(climate,points)
# bioclim <- cbind.data.frame(coordinates(points),values)
# colnames(bioclim)[1]='lon'
# colnames(bioclim)[2]='lat'
# 
# env=merge(df,bioclim,by=c('lat','lon'))
# env$ivm=0
# env$ivm[env$iso2=='WAL']=1
# env$ivm[env$iso2=='KOK'|env$iso2=='WHR']=1
# write.table(env,file='./data/geo_clim_drug_info.txt',quote=F,sep="\t",row.names = FALSE)

env = read.table(file='./data/geo_clim_drug_info.txt',header=T)

#####====== Phylogenetic inference (mtDNA) - Figure 2a ====########
library(ggtree)
require(ape)
require(phytools)
setwd('./data/MITO/')

###-- Rooted Tree
tree0 <- read.tree(file="./pop.muscle.phy_phyml_tree.txt")
plot(density((as.numeric(as.character(tree0$node.label))),na.rm=T))
summary(as.numeric(as.character(tree0$node.label)))
## Rooted tree
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    0.00   14.00   57.00   54.38   98.00  100.00       1 

###-- Unrooted Tree // branch-support is improved
tree <- read.tree(file="./pop_unroot.muscle.phy_phyml_tree.txt")
plot(density((as.numeric(as.character(tree$node.label))),na.rm=T))
summary(as.numeric(as.character(tree$node.label)))
## Unrooted tree
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    0.00   19.75   60.50   57.96  100.00  100.00       1 

##--- Get iso names
tree$tip.label = df$iso[match(tree$tip.label,df$sample)] 
tree$tip.label[tree$tip.label=='FRG']='GUA' ## Modify Guadeloupe encoding
cls = list(c1=c("FRA.1","FRA.2","FRA.3","FRA.4","POR","MOR"),
           c2=c("IND", "AUS.1","AUS.2"),
           c3=c("BRA","GUA"),
           c4=c("STA.1","STA.2","STA.3","NAM","ZAI","ISE"),
           c5=c("BEN", "CAP","STO"))

tree = groupOTU(tree,cls)

nodsize = as.numeric(as.character(tree$node.label))/20 
nodcol=array('#de2d26',length(tree$node.label))
nodcol[which(as.numeric(as.character(tree$node.label))>=70)] = '#41b6c4'

while(0 %in% levels(attributes(tree)$group)){
  attributes(tree)$group[attributes(tree)$group==0] = attributes(tree)$group[which(attributes(tree)$group==0)-1]
  attributes(tree)$group = factor(attributes(tree)$group)
}

pdf(file=paste0(draft2,'Figure2a.pdf'),height=10,width=8) 
print(ggtree(tree,aes(color = group),layout='daylight') +
        #geom_tiplab(size = 2) +
        scale_colour_manual(values = clcol,name = "Location",breaks = cls) +
        geom_treescale(y = -.1,offset = 0.001) +
        geom_nodepoint(color = nodcol, alpha = 0.2, size = nodsize))
dev.off()

rm(tree,tree0)

#####====== Twisst topology plot - Figure 3b ====########
## Twisst plotting function from Martin et al. 2016
setwd("./data/NUC/TW/")

simple.loess.predict <- function(x, y, span, weights = NULL, max = NULL, min = NULL){
  y.loess <- loess(y ~ x, span = span, weights = weights)
  y.predict <- predict(y.loess,x)
  if (is.null(min) == FALSE) {y.predict = ifelse(y.predict > min, y.predict, min)}
  if (is.null(max) == FALSE) {y.predict = ifelse(y.predict < max, y.predict, max)}
  y.predict
}

smooth_df <- function(x, df, span, col.names=NULL, weights=NULL, min=NULL, max=NULL){
  smoothed <- df
  if (is.null(col.names)){col.names=colnames(df)}
  for (col.name in col.names){
    print(paste("smoothing",col.name))
    smoothed[,col.name] <- simple.loess.predict(x,df[,col.name],span = span, max = max, min = min, weights = weights)
  }
  smoothed
}

stack <- function(mat){
  upper <- t(apply(mat, 1, cumsum))
  lower <- upper - mat
  list(upper=upper,lower=lower)
}

interleave <- function(x1,x2){
  output <- vector(length= length(x1) + length(x2))
  output[seq(1,length(output),2)] <- x1
  output[seq(2,length(output),2)] <- x2
  output
}

plot_weights <- function(weights_dataframe,positions=NULL,line_cols=NULL,fill_cols=NULL,xlim=NULL,ylim=c(0,1),stacked=FALSE,cex.lab=2,
                         ylab="Weights", xlab = "Position (Kbp)", main="",xaxt=NULL,yaxt=NULL,bty="n", add=FALSE){
  #get x axis
  x = positions
  #if a two-column matrix is given - plot step-like weights with start and end of each window
  if (is.matrix(x)==TRUE) {
    x = interleave(positions[,1],positions[,2])
    yreps=2
  }
  else {
    if (is.null(x)==FALSE) x = positions
    else x = 1:nrow(weights_dataframe)
    yreps=1
  }

  #set x limits
  if(is.null(xlim)) xlim = c(min(x), max(x))

  #if not adding to an old plot, make a new plot
  if (add==FALSE) plot(0, pch = "", xlim = xlim, ylim=ylim, ylab=ylab, xlab=xlab, main=main,xaxt=xaxt,yaxt=yaxt,bty=bty)

  if (stacked == TRUE){
    y_stacked <- stack(weights_dataframe)
    for (n in 1:ncol(weights_dataframe)){
      y_upper = rep(y_stacked[["upper"]][,n],each=yreps)
      y_lower = rep(y_stacked[["lower"]][,n],each = yreps)
      polygon(c(x,rev(x)),c(y_upper, rev(y_lower)), col = fill_cols[n], border=NA)
    }
  }
  else{
    for (n in 1:ncol(weights_dataframe)){
      y = rep(weights_dataframe[,n],each=yreps)
      polygon(c(x,rev(x)),c(y, rep(0,length(y))), col=fill_cols[n], border=NA)
      lines(x,y, type = "l", col = line_cols[n])
    }
  }
}

options(scipen = 5)

## input data for the beta-tubulin region
#CHROM.1.5800000:8000000.bzr.weights.txt

dir='NAM_FRG_FRA_STO'

##-- Read keep list
kl = read.table(file=paste0('./',dir,'/keep.list'),header=F,sep=' ')
kl$cov = kl$V9

#weights file with a column for each topology
w=readLines(paste0('./',dir,'/CHROM.1.5800000:8000000.bzr.weights.csv.gz'))
w2=w[-c(1,2,3)]
weights = read.csv(textConnection(w2),header=TRUE,sep='\t')

#coordinates file for each window
window_data_file <- paste0('./',dir,"/CHROM.1.5800000:8000000.bzr.raxml_bionj.w50.data.tsv")

#normalise rows so weights sum to 1
weights <- weights / apply(weights, 1, sum)
#retrieve the names of the topologies
topoNames = names(weights)

window_data = read.table(window_data_file, header = T)

s = 7016384-50e3
e = 7039311+50e3

weights = weights[window_data$mid<e & window_data$mid>s,]
window_data = window_data[window_data$mid<e & window_data$mid>s,]

#exclude any rows where data is missing
good_rows = which(is.na(apply(weights,1,sum)) == F)
weights <- weights[good_rows,]
window_data = window_data[good_rows,]

## plot raw data

#plot raw data in "stepped" style, with polygons stacked.
#specify stepped style by providing a matrix of starts and ends for positions
pdf(file = paste0(draft2,"Figure3b.pdf"), width = 10, height = 4)
par(mar = c(4,4,1,1))
plot_weights(weights_dataframe=weights, positions=cbind(window_data$start/1000,window_data$end/1000),
             line_cols=cols, fill_cols=cols, stacked=T,cex.lab=2,xlim =c(s/1000,e/1000))
abline(v=7027384/1000,col='white',lty=2) #- btub
abline(v=7029311/1000,col='white',lty=2) #- btub
dev.off()

## plot topologies
library(ape)
topos = read.tree(file=paste0("./",dir,"/topologies.CHROM.1.5800000:8000000.bzr.trees"))
for (i in 1:length(topos)) topos[[i]] <- root(topos[[i]], "NAM", resolve.root = T)

pdf(file = paste0("topos.rooted.",dir,".pdf"), width = 5, height = 2)
par(mfrow = c(1,3), mar = c(1,1,2,1), xpd=NA)
for (n in 1:length(topos)){
  plot.phylo(topos[[n]], type = "clad", edge.color=cols[n], edge.width=5, label.offset=.1, cex = 1.)
  mtext(side=3,text=paste0("topology",n))
}
dev.off()

######---------------- Same plot for another topology
dir='NAM_FRG_FRA_MOR'
kl = read.table(file=paste0('./',dir,'/keep.list'),header=F,sep=' ')
kl$cov = kl$V9

#weights file with a column for each topology
w=readLines(paste0('./',dir,'/CHROM.1.5800000:8000000.bzr.weights.csv.gz'))
w2=w[-c(1,2,3)]
weights = read.csv(textConnection(w2),header=TRUE,sep='\t')
#weights_file <- "CHROM.1.5800000:8000000.bzr.weights.txt"

#coordinates file for each window
window_data_file <- paste0('./',dir,"/CHROM.1.5800000:8000000.bzr.raxml_bionj.w50.data.tsv")

## read data
#weights = read.table(weights_file, header = T)
#normalise rows so weights sum to 1
weights <- weights / apply(weights, 1, sum)
#retrieve the names of the topologies
topoNames = names(weights)

window_data = read.table(window_data_file, header = T)

weights = weights[window_data$mid<e & window_data$mid>s,]
window_data = window_data[window_data$mid<e & window_data$mid>s,]

#exclude any rows where data is missing
good_rows = which(is.na(apply(weights,1,sum)) == F)
weights <- weights[good_rows,]
window_data = window_data[good_rows,]

## plot raw data

#plot raw data in "stepped" style, with polygons stacked.
#specify stepped style by providing a matrix of starts and ends for positions
pdf(file = paste0(draft2,"supplementary_figure11.pdf"), width = 10, height = 4)
par(mar = c(4,4,1,1))
plot_weights(weights_dataframe=weights, positions=cbind(window_data$start/1000,window_data$end/1000),
             line_cols=cols, fill_cols=cols, stacked=T,cex.lab=2,xlim =c(s/1000,e/1000))
abline(v=7027384/1000,col='white',lty=2) #- btub
abline(v=7029311/1000,col='white',lty=2) #- btub
dev.off()

## plot topologies
library(ape)
topos = read.tree(file=paste0("./",dir,"/topologies.CHROM.1.5800000:8000000.bzr.trees"))
for (i in 1:length(topos)) topos[[i]] <- root(topos[[i]], "NAM", resolve.root = T)

pdf(file = paste0("topos.rooted.",dir,".pdf"), width = 5, height = 2)
par(mfrow = c(1,3), mar = c(1,1,2,1), xpd=NA)
for (n in 1:length(topos)){
  plot.phylo(topos[[n]], type = "clad", edge.color=cols[n], edge.width=5, label.offset=.1, cex = 1.)
  mtext(side=3,text=paste0("topology",n))
}
dev.off()

#####====== Pop CLUSTERING using VQSR SNP - 223 samples ====########
setwd("./data/NUC/PCA/")
chr='WG'

### Creation of a gds file
#ped.fn<-paste("GENO.VQSR.MAF5.",chr,".ped.gz",sep="")
#map.fn<-paste("GENO.VQSR.MAF5.",chr,".fmt.map.gz",sep="")
# snpgdsPED2GDS(ped.fn, map.fn,paste("FLT_VQSR.fmt",chr,".gds",sep=""),
#               family=TRUE, snpfirstdim=FALSE,
#               compress.annotation="ZIP_RA.max", compress.geno="", verbose=TRUE)
### File available on figshare repository: 10.6084/m9.figshare.7923512 
genofile<-snpgdsOpen(paste('FLT_VQSR.fmt',chr,'.gds',sep=''))
mafq = 0.05
callrate = 0.5

#LD based SNP pruning
set.seed(1000)
snpset = snpgdsLDpruning(autosome.only = T,genofile,
                         ld.threshold = .2, maf = mafq, remove.monosnp = T, missing.rate = 1-callrate)
# SNP pruning based on LD:
#   Excluding 0 SNP on non-autosomes
# Excluding 7,708,668 SNPs (monomorphic: TRUE, MAF: 0.05, missing rate: 0.5)
# Working space: 223 samples, 411,574 SNPs
# using 1 (CPU) core
# sliding window: 500,000 basepairs, Inf SNPs
# |LD| threshold: 0.2
# method: composite
# Chromosome 1: 0.17%, 2,596/1,567,226
# Chromosome 2: 0.18%, 2,892/1,569,964
# Chromosome 3: 0.18%, 2,586/1,473,245
# Chromosome 4: 0.17%, 2,976/1,739,874
# Chromosome 5: 0.16%, 2,843/1,769,933
# 13,893 markers are selected in total.
snpid = unlist(snpset)

ccm_pca<-snpgdsPCA(snp.id = snpid, genofile,autosome.only=T,
                   maf = mafq,
                   missing.rate = 1-callrate,
                   num.thread=2)
# Principal Component Analysis (PCA) on genotypes:
#   Excluding 8,017,187 SNPs (non-autosomes or non-selection)
# Excluding 0 SNP (monomorphic: TRUE, MAF: 0.05, missing rate: 0.5)
# Working space: 223 samples, 102,946 SNPs

isolate=df$iso[match(ccm_pca$sample.id,df$sample)]
#isolate[substr(isolate,1,3)=='FRA']=""
Origin=df$geo[match(ccm_pca$sample.id,df$sample)]
npop=length(unique(isolate))
pc.percent <- 100 * ccm_pca$eigenval/sum(ccm_pca$eigenval,na.rm=T)
pc.percent

p1=ggplot(data.frame(ccm_pca$eigenvect),
          aes(x=ccm_pca$eigenvect[,1],
              y=ccm_pca$eigenvect[,2],
              label=ccm_pca$sample.id,
              color=Origin)) +
  geom_point(size=2) +
  geom_text_repel(segment.alpha = 0.2,size=3,aes(x=ccm_pca$eigenvect[,1],y=ccm_pca$eigenvect[,2],col=Origin,label=isolate)) +
  theme(legend.position='none') +
  scale_color_manual(values=clcol) +
  labs(x = paste("PC1 -",round(pc.percent[1],2),' %',sep='')) +
  labs(y = paste("PC2 -",round(pc.percent[2],2),' %',sep='')) +
  theme(text=element_text(size=12))
print(p1)

pdf(file=paste0(draft2,"Supplementary_Figure3.pdf"),width=14,height = 8)
print(p1)
dev.off()
rm(p1)
closefn.gds(genofile)

#####====== Pop CLUSTERING using ANGSD - ALL samples ====########
setwd("./data/PCA/")
chr='GW_ld'

## Data generated using ANGSD (ANGSD_223.sh)
## LD estimation for low coverage data was implemented following Bilton et al. 2018 with ld_lowcov_splitter.py & LD_estim_lowcov.R

# Read input file
covar <- read.table(paste0('./ALL',chr,'.covar'), stringsAsFact=F)
# Read annot file
annot <- read.table('./InfoPop.clst', sep="\t", header=T); # note that plink cluster files are usually tab-separated
# Eigenvalues
eig <- eigen(covar, symm=TRUE);
eig$val <- eig$val/sum(eig$val);
cat(signif(eig$val, digits=3)*100,"\n");
plot(density(RMTstat::rtw(223,beta=1)))

# Tracy Widom test
#test each dimension
m=length(eig$values)
a=eig
res <- matrix(data=NA,nrow=(m-1),ncol=5)
colnames(res) <- c("Dimension","Eigenvalue","nhat","TWstat","P-value")
for (j in 1:(m-1)){
  L1 <- sum(a$values[j:m])
  L2 <- sum(a$values[j:m]^2)
  lambda <- a$values[j]*(m-j)/L1
  nhat <- L1^2/L2
  mu <- (sqrt(nhat-1)+sqrt(m-j))^2/nhat
  sigma <- (sqrt(nhat-1)+sqrt(m-j))/nhat*(1/sqrt(nhat-1)+1/sqrt(m-j))^(1/3)
  twstat <- (lambda-mu)/sigma
  twpvalue <- RMTstat::ptw(twstat,lower.tail=FALSE)
  res[j,1] <- j
  res[j,2] <- lambda
  res[j,3] <- nhat
  res[j,4] <- twstat
  res[j,5] <- twpvalue
}
res[res[,5]<0.05,]
#      Dimension Eigenvalue     nhat   TWstat P-value
# [1,]   1   7.013588 152.7290 22.93302 0
# [2,]   2   5.591693 168.9922 11.68487 0
rm(m,a)
# Parse components to analyze
comp = c(1,2) #<- as.numeric(strsplit(opt$comp, "-", fixed=TRUE)[[1]])

# Plot
PC <- as.data.frame(eig$vectors)
colnames(PC) <- gsub("V", "PC", colnames(PC))
PC$Pop <- factor(annot$CLUSTER)

x_axis = paste0("PC",comp[1])
xtitle=paste0("PC",comp[1]," (",signif(eig$val[comp[1]], digits=3)*100,"%)")
y_axis = paste0("PC",comp[2])
ytitle=paste0("PC",comp[2]," (",signif(eig$val[comp[2]], digits=3)*100,"%)")
PC$Origin=factor(df$geo[match(PC$Pop,df$iso)])
PC$Pop = gsub('FRG','GUA',PC$Pop)

## Figure 1b
pdf(file=paste0(draft2,'Figure1b.pdf'),width=14,height=8)
ggplot(data=PC, aes_string(x=x_axis, y=y_axis, color=PC$Origin)) + 
  theme_bw() + ylab(ytitle) + xlab(xtitle) +
  scale_color_manual(values=clcol) +
  geom_point(size=3) + 
  geom_text_repel(segment.alpha = 0.2,size=3,
                  aes(col=Origin,label=Pop)) +
  #theme(legend.position='none',axis.title = element_text(size=16))
  theme(legend.position=c(.9,0.1),axis.title = element_text(size=16)) +
  labs(colour='Geographical origin')
dev.off()

euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
a=matrix(0,223,223)
for(i in 1:nrow(PC)){
  for(j in 1:nrow(PC)){
    a[i,j]=euc.dist(PC$PC1[i],PC$PC2[j])
  }
}
colnames(a)=df$sample
rownames(a)=df$sample

###------------------ Phylogeography ?
dmx=matrix(0,223,223)

for(i in 1:dim(df)[1]){
  for(j in 1:dim(df)[1]){
    lon1=df$lon[i]
    lat1=df$lat[i]
    lon2=df$lon[j]
    lat2=df$lat[j]
    dmx[i,j] =  distm(c(lon1, lat1), c(lon2, lat2), fun = distHaversine)/1000
  }
}
rownames(dmx)=df$sample
colnames(dmx)=df$sample
mtest1=mantel(dmx,a,method = 'spearman')
mtest1
# Call:
#   mantel(xdis = dmx, ydis = a, method = "spearman") 
# 
# Mantel statistic r: 0.06865 
# Significance: 0.003 
# 
# Upper quantiles of permutations (null model):
#   90%    95%  97.5%    99% 
#   0.0285 0.0355 0.0411 0.0523 
# Permutation: free
# Number of permutations: 999


#####====== ADMIXTURE USING ANGSD - ALL samples ====########
##-- Analysis run with ANGSD (ANGSD_223.sh)
##-- LD estimation for low coverage data was implemented following Bilton et al. 2018 with ld_lowcov_splitter.py & LD_estim_lowcov.R
setwd("./data/NUC/ADMXITURE/")
rep = "./data/NUC/ADMXITURE/"
prefix = "ALL"
k = seq(2,10,1)
k = 3
###--- Admxiture plot for Figure 2c
K = 3
chr='GW_ld'

infile=paste(rep,prefix,".",chr,".",K,".qopt",sep="")

outfile = paste0(draft2,"Figure2c.pdf")
# Admixture
admix <- t(as.matrix(read.table(infile)))
K <- nrow(admix)

pdf(file=outfile,width=14,height=3)

samp<-read.table(paste(rep,"sample.info",sep=""),as.is=T)
pop = df[match(samp$V1,df$sample),c('sample','lon','iso')] 
admix<-admix[,order(pop[,2])]
pop<-pop[order(pop[,2]),]
pop$rn=seq(1:dim(pop)[1])
ord=unique(pop[,3])
bd = pop$rn[which(row.names(pop) %in% row.names(unique(pop[,2:3])))]
lab = array('',dim(pop)[1])
lab[bd]=pop$iso[bd]
colnames(admix)=pop$sample
m=melt(admix)
colnames(m)=c("K","Population","Admixture")
lab = gsub('FRG','GUA',lab)

print(ggplot(m,aes(x=Population,y=Admixture,fill=K,col=K))+
        geom_bar(stat='identity') + 
        scale_fill_manual(values=clcol[c(3,1,4)]) + 
        scale_color_manual(values=clcol[c(3,1,4)]) + 
        geom_vline(xintercept=bd-0.5,colour="white") +
        scale_y_continuous(expand=c(0,0)) +
        theme_classic() + #coord_flip() +
        theme(legend.position='none',
              axis.text.x = element_text(angle = 90, hjust = 1),
              text=element_text(size=18)) + xlab('')+
        scale_x_discrete(breaks=pop$sample,labels=lab))

invisible(dev.off())


##--- Supplementary Figure 5: admixture plots for various K
## This yields 1 pdf file for each K; to be combined into a single plot for Fig. S5.
k=c(2,seq(4,10,1))
chr='GW_ld'
for(K in k){
  infile=paste(rep,prefix,".",chr,".",K,".qopt",sep="")
  outfile=paste(rep,prefix,".",chr,".",K,".ALL.pdf",sep="")
  # Admixture
  admix <- t(as.matrix(read.table(infile)))
  K <- nrow(admix)

  samp<-read.table(paste(rep,"sample.info",sep=""),as.is=T)
  pop = df[match(samp$V1,df$sample),c('sample','lon','iso')]
  admix<-admix[,order(pop[,2])]
  pop<-pop[order(pop[,2]),]
  pop$rn=seq(1:dim(pop)[1])
  ord=unique(pop[,3])
  bd = pop$rn[which(row.names(pop) %in% row.names(unique(pop[,2:3])))]
  lab = array('',dim(pop)[1])
  lab[bd]=pop$iso[bd]
  colnames(admix)=pop$sample
  m = melt(admix)
  colnames(m) = c("K","Population","Admixture")
  lab = gsub('FRG','GUA',lab)
  
  pdf(file=outfile,width=14,height=2)
  print(ggplot(m,aes(x=Population,y=Admixture,fill=K,col=K))+
          geom_bar(stat='identity') + #ggtitle(paste0('K=',K))+
          scale_fill_manual(values=cols) +
          scale_color_manual(values=cols) +
          geom_vline(xintercept=bd-0.5,colour="white") +
          scale_y_continuous(expand=c(0,0)) +
          theme_classic() + #coord_flip() +
          theme(legend.position='none',
                axis.text.x = element_text(angle = 90, hjust = 1),
                text=element_text(size=16)) + xlab('')+
          scale_x_discrete(breaks=pop$sample,labels=lab))
  invisible(dev.off())
}

#####====== ADMIXTURE USING ANGSD - Cross-Validation: K = 3 ====########
setwd("./data/NUC/ADMIXTURE/")
rep="./data/NUC/ADMIXTURE/"
prefix="ALL"
k=seq(2,10,1)
chrl=c('GW_ld.CV1','GW_ld.CV2','GW_ld.CV3','GW_ld.CV4','GW_ld.CV5')
cv1 = array(-1,length(k))
cv2 = array(-1,length(k))
for(K in k){
  for(chr in chrl){
    infile=paste(rep,prefix,".",chr,".",K,".qopt",sep="")
    
    # Admixture
    assign(paste0('admix.',chr),t(as.matrix(read.table(infile))))
  }
  # Med absolute deviation
  cv1[K-1] = mad(c(admix.GW.CV1,admix.GW.CV2,admix.GW.CV3,admix.GW.CV4,admix.GW.CV5))
  # Jackniffing variance
  cv2[K-1] = var(c(admix.GW.CV1,admix.GW.CV2,admix.GW.CV3,admix.GW.CV4,admix.GW.CV5))
  rm(admix.GW.CV1,admix.GW.CV2,admix.GW.CV3,admix.GW.CV4,admix.GW.CV5)
}
pdf(file=paste0(draft2,'Supplementary_Figure6.pdf'))
plot(k,log10(cv1),pch=19,type='b',ylab='Median absolute deviation')
dev.off()

#####====== Diversity and divergence between populations  ====########

###--------- Nucleotide diversity using ANGSD 
rep="./data/NUC/TAJHI"

###---- Summary statistics
theme_set(theme_bw()) 
t.HI=matrix(0,1,7)
t.HI=data.frame(t.HI)
colnames(t.HI)=c('Chr','WinCenter','tP','Tajima','nSites','cov','pop')
t.HI$Chr=factor(t.HI$Chr,levels=c(chromlist))
for(r in c('FRA.1','FRG','NAM')){
  for(c in chromlist){
    temp = read.table(paste0(rep,'/',c,'/',r,'.thetas.idx.pestPG'))
    temp = temp[,c(2,3,5,9,14)]
    colnames(temp) = c('Chr','WinCenter','tP','Tajima','nSites')
    temp$cov = 'HI'
    temp$pop = r
    t.HI = rbind(t.HI,temp)
    rm(temp)
  }
}
t.HI=na.omit(t.HI)
t.HI=t.HI[t.HI$nSites>500,] #- 10% of win width
piMean=aggregate(tP/nSites ~ pop,FUN=mean,data=t.HI)
piMean$std=aggregate(tP/nSites ~ pop,FUN=sd,data=t.HI)[,2]
piMean
#     pop   tP/nSites          std
# 1 FRA.1 0.006519367 0.0005379378
# 2   FRG 0.010595454 0.0007856049
# 3   NAM 0.011322935 0.0007861508

##-- Ne from pi
mu2=2.7e-9
Ne2=piMean$`tP/nSites`/(4*mu2)
Ne2/1e6
#[1] 0.6036451 0.9810606 1.0484199

##-- Overall distribution
wb = 2.47e-4
ovols = 4e-3
ovoln = 1e-3
human = 0.119/100 #Perry et al. Genome Res2012 PMC3317143/
droso = mean(c(0.00531,0.00752,0.01714)) #Langley 2012 Genetics
mpi = max(t.HI$tP/t.HI$nSites)

##-- Supplementary Figure 1
pdf(file=paste0(draft2,'Supplementary_Figure1.pdf'))
t.HI$pop=factor(t.HI$pop)
ggplot(t.HI,aes(x=pop,y=tP/nSites,col=pop)) +
  theme_classic()+theme(text=element_text(size=16)) + 
  xlab('Population') +
  scale_y_continuous(limits=c(0,mpi),breaks=seq(0,mpi,by=0.001)) +
  #geom_hline(yintercept = human,col='blue') + #annotate("Human", x=c(human+0.01),y=200) + 
  geom_hline(yintercept = droso,col='red') + #annotate("D.melanogaster", x=c(droso+0.01),y=200) +
  geom_hline(yintercept = wb,col='forestgreen') + #annotate("W.bancrofti", x=c(wb+0.01),y=200) +
  geom_hline(yintercept = ovols,col='purple') + #annotate("W.bancrofti", x=c(wb+0.01),y=200) +
  geom_hline(yintercept = ovoln,col='pink') + #annotate("W.bancrofti", x=c(wb+0.01),y=200) +
  ylab('Nucleotide diversity') +
  geom_point(size=3) + theme(legend.position='none') +
  scale_color_manual(values=c(clcol[c(1,3,4)]))
dev.off()

rm(t.HI)

###----- Diversity values for ALL populations  
rep = "./data/NUC/TAJ/"
iso5 = unique(df$iso[which(df$iso!='BRA' & df$iso!='FRA.4')]) ## discard pop with less than 5 individuals

###---- Summary statistics
theme_set(theme_bw()) 
t.ALL=matrix(0,1,7)
t.ALL=data.frame(t.ALL)
colnames(t.ALL)=c('Chr','WinCenter','tP','Tajima','nSites','cov','pop')
t.ALL$Chr=factor(t.ALL$Chr,levels=c(chromlist))
for(r in iso5){
  for(c in chromlist){
    temp = read.table(paste0(rep,'/',c,'/',r,'.thetas.idx.pestPG'))
    temp = temp[,c(2,3,5,9,14)]
    colnames(temp) = c('Chr','WinCenter','tP','Tajima','nSites')
    temp$cov = 'HI'
    temp$pop = r
    t.ALL = rbind(t.ALL,temp)
    rm(temp)
  }
}
t.ALL=na.omit(t.ALL)
piMean=aggregate(tP/nSites ~ pop,FUN=mean,data=t.ALL)
piMean$std=aggregate(tP/nSites ~ pop,FUN=sd,data=t.ALL)[,2]
a = aggregate(mcov ~ iso,FUN=mean,data=df)
piMean$mcov = a$mcov[match(piMean$pop,a$iso)]
piMean
#      pop   tP/nSites          std      mcov
# 1    ACO 0.006807788 0.0009633793 1.7920000
# 2  AUS.1 0.006132598 0.0004379338 2.5666667
# 3  AUS.2 0.005093270 0.0003204319 2.7360000
# 4    BEN 0.009389602 0.0023180484 0.3300000
# 5    CAP 0.007210416 0.0010416897 0.7042857
# 6  FRA.1 0.007231135 0.0003690067 3.3248571
# 7  FRA.2 0.005439821 0.0005930883 1.8166667
# 8  FRA.3 0.006230038 0.0008429861 0.4833333
# 9    FRG 0.010777981 0.0003249374 4.3986957
# 10   IND 0.005551998 0.0007544601 4.9500000
# 11   MOR 0.008197893 0.0004382342 3.5000000
# 12   NAM 0.013001924 0.0007489163 5.6153333
# 13 STA.1 0.006692275 0.0018268375 4.0071429
# 14 STA.2 0.004417365 0.0004476453 1.7914286
# 15 STA.3 0.011243171 0.0008985027 2.6376923
# 16   STO 0.009746861 0.0007465504 2.6133333
# 17   ZAI 0.008575769 0.0002907212 2.1080000

summary(piMean$`tP/nSites`)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.004417 0.006133 0.007210 0.007749 0.009390 0.013002 

###-- biased by coverage // strange case of Indonesia
ggplot(piMean,aes(x=`tP/nSites`,y=mcov)) + 
  geom_smooth(method='lm') +
  geom_point(aes(col=pop))

summary.aov(lm(tP/nSites ~ mcov + Chr,data = t.ALL))
#             Df    Sum Sq   Mean Sq F value  Pr(>F)   
# mcov         1 0.0000542 5.418e-05   9.497 0.00283 **
# Chr          4 0.0000180 4.510e-06   0.790 0.53497   
# Residuals   79 0.0004507 5.700e-06   

rcorr(piMean$`tP/nSites`,piMean$mcov)
# x    y
# x 1.00 0.34
# y 0.34 1.00
# 
# n= 17 
# P
# x      y     
# x        0.1783
# y 0.1783  

#####====== BTUB case study ANGSD ALL // coord(tbb1[1], 7029311) ====########

##-- Overall Nucleotide diversity across the genome
rep="./data/NUC/TAJHI/"
theme_set(theme_bw()) 
t.HI = matrix(0,1,8)
t.HI = data.frame(t.HI)
colnames(t.HI) = c('Chr','WinCenter','tP','Tajima','nSites','cov','pop','bin')
temp0 = NULL
temp = NULL
t.HI$Chr = factor(t.HI$Chr,levels=c(chromlist))

###-- First resistant pops
for(r in c('FRA.1','FRG','NAM')){
  for(c in chromlist){
    temp0 = read.table(paste0(rep,'/',c,'/',r,'.thetas.pestPG'))
    temp0 = temp0[,c(2,3,5,9,14)]
    colnames(temp0) = c('Chr','WinCenter','tP','Tajima','nSites')
    temp0$cov = 'HI'
    temp0$pop = r
    temp = rbind(temp,temp0)
    rm(temp0)
  }
  temp$bin = seq(1,nrow(temp),1)
  t.HI = rbind(t.HI,temp)
  temp = NULL
}
t.HI = na.omit(t.HI)
t.HI = t.HI[t.HI$nSites>500,] #- 10% of win width

###-- Add susceptible pops
rep="./data/NUC/TAJ/"
temp0 = NULL
temp = NULL
for(r in c('AUS.2','ZAI')){
  for(c in chromlist){
    temp0 = read.table(paste0(rep,'/',c,'/',r,'.thetas.idx.pestPG'))
    temp0 = temp0[,c(2,3,5,9,14)]
    colnames(temp0) = c('Chr','WinCenter','tP','Tajima','nSites')
    temp0$cov = 'HI'
    temp0$pop = r
    temp = rbind(temp,temp0)
    rm(temp0)
  }
  temp$bin = seq(1,nrow(temp),1)
  t.HI = rbind(t.HI,temp)
  temp = NULL
}
t.HI = na.omit(t.HI)
t.HI = t.HI[t.HI$nSites>500,] #- 10% of win width

t.HI$res='R'
t.HI$res[t.HI$pop=='AUS.2'|t.HI$pop=='ZAI']='S'
t.HI$res = factor(t.HI$res)

t.HI$btub = 'Other'
#t.HI$btub[t.HI$Chr=='1_Celeg_TT' & t.HI$WinCenter>tbb1[1]-250e3 & t.HI$WinCenter<7029311+500e3]='VicinityBtub'
t.HI$btub[t.HI$Chr=='1_Celeg_TT' & t.HI$WinCenter>tbb1[1]-500e3 & t.HI$WinCenter<tbb1[2]+500e3] = 'Vicinity of Hco-tbb-iso-1'
t.HI$group = factor(paste(t.HI$pop,t.HI$res,sep=" - "))

##-- Reduced diversity in R vs. S
t.test(tP/nSites ~ btub,data=t.HI)
# Welch Two Sample t-test
# 
# data:  tP/nSites by btub
# t = 15.075, df = 283.61, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.002994190 0.003893526
# sample estimates:
# mean in group Other mean in group Vicinity of Hco-tbb-iso-1 
# 0.010554549                             0.007110691 
(0.010554549-0.007110691)/0.010554549 #[1] 0.3262913

t.test(Tajima ~ btub,data=t.HI)
# Welch Two Sample t-test
# 
# data:  Tajima by btub
# t = 6.0199, df = 282.5, p-value = 0.000000005407
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.1817867 0.3584244
# sample estimates:
#   mean in group Other mean in group Vicinity of Hco-tbb-iso-1 
# -0.2053464                              -0.4754520 
0.4754520/0.2053464

summary(lm(Tajima ~ btub*pop,data=t.HI))
aggregate(Tajima ~ btub,data=t.HI,FUN=sd)
#                        btub    Tajima
# 1                     Other 0.5967160
# 2 Vicinity of Hco-tbb-iso-1 0.7524641

####------ Plot btub region vs. chromosome I
rep="./data/NUC/TAJHI/"
rm(t.HI)
t.HI = matrix(0,1,8)
t.HI = data.frame(t.HI)
colnames(t.HI) = c('Chr','WinCenter','tP','Tajima','nSites','cov','pop','bin')
temp0 = NULL
temp = NULL
t.HI$Chr = factor(t.HI$Chr,levels=c(chromlist))

###-- First resistant pops
for(r in c('FRA.1','FRG','NAM')){
  for(c in '1_Celeg_TT'){
    temp0 = read.table(paste0(rep,'/',c,'/',r,'.thetas.pestPG'))
    temp0 = temp0[,c(2,3,5,9,14)]
    colnames(temp0) = c('Chr','WinCenter','tP','Tajima','nSites')
    temp0$cov = 'HI'
    temp0$pop = r
    temp = rbind(temp,temp0)
    rm(temp0)
  }
  temp$bin = seq(1,nrow(temp),1)
  t.HI = rbind(t.HI,temp)
  temp = NULL
}
t.HI = na.omit(t.HI)
t.HI = t.HI[t.HI$nSites>500,] #- 10% of win width

###-- Add susceptible pops
rep="./data/NUC/TAJ/"
temp0 = NULL
temp = NULL
for(r in c('AUS.2','ZAI')){
  for(c in chromlist){
    temp0 = read.table(paste0(rep,'/',c,'/',r,'.thetas.pestPG'))
    temp0 = temp0[,c(2,3,5,9,14)]
    colnames(temp0) = c('Chr','WinCenter','tP','Tajima','nSites')
    temp0$cov = 'HI'
    temp0$pop = r
    temp = rbind(temp,temp0)
    rm(temp0)
  }
  temp$bin = seq(1,nrow(temp),1)
  t.HI = rbind(t.HI,temp)
  temp = NULL
}
t.HI = na.omit(t.HI)
t.HI = t.HI[t.HI$nSites>500,] #- 10% of win width

t.HI$res='R'
t.HI$res[t.HI$pop=='AUS.2'|t.HI$pop=='ZAI']='S'
t.HI$res = factor(t.HI$res)

## Replace FRG with GUA
t.HI$pop[t.HI$pop=='FRG']='GUA'

t.HI$btub = 'Other'
t.HI$btub[t.HI$Chr=='1_Celeg_TT' & t.HI$WinCenter>tbb1[1]-500e3 & t.HI$WinCenter<tbb1[2]+500e3] = 'Vicinity of Hco-tbb-iso-1'
t.HI$group = factor(paste(t.HI$res,t.HI$pop,sep=" - "))

p1Divbtub = ggplot(t.HI[t.HI$Chr=='1_Celeg_TT',],aes(x=WinCenter/1e6,y=tP/nSites,col=btub)) +
  theme_classic()+
  theme(legend.position='bottom',text=element_text(size=18))+
  geom_point(size=.5) + scale_color_manual("",values=cols) + 
  facet_wrap(~group,ncol=1)+
  ggtitle('') +
  guides(colour = guide_legend(override.aes = list(size=3)))+
  xlab('Position (Mbp)') + ylab('Nucleotide diversity') 

pdf(file=paste0(draft2,'Supplementary_Figure9.pdf'),width=14,height=8)
print(p1Divbtub)
invisible(dev.off())

##-- Tajima's D with observed and simulated estimates
t.HI$pop[t.HI$pop=='GUA']='FRG'

tHIbtub = t.HI[t.HI$Chr=='1_Celeg_TT' & t.HI$WinCenter>tbb1[1]-5e6 & t.HI$WinCenter<tbb1[2]+5e6,]
tHIbtub$sim = 0
tHIbtub$std = 0

for(pop in c('FRA.1','NAM','FRG')){
  temp=NULL
  c = 1
  temp = rbind(temp,read.table(file=paste0('./data/NUC/MSMC/',pop,'.',c,'.tajimaD.txt')))
  tHIbtub$sim[tHIbtub$pop==pop] = mean(temp$V1)
  tHIbtub$std[tHIbtub$pop==pop] = sd(temp$V1)
  t.HI$sim[t.HI$pop==pop] = mean(temp$V1)
  t.HI$std[t.HI$pop==pop] = sd(temp$V1)
  
  rm(temp)
}

p2Divbtub = ggplot(tHIbtub,aes(x=WinCenter/1e6,y=Tajima,col=btub)) +
  theme_classic()+
  theme(legend.position='bottom',text=element_text(size=18),legend.title = NULL) +
  geom_hline(aes(yintercept = sim),col='gray50',lty=1) + 
  geom_hline(aes(yintercept = sim + 2.58*std),col='gray50',lty=3) + 
  geom_hline(aes(yintercept = sim - 2.58*std),col='gray50',lty=3) + 
  geom_vline(xintercept = (tbb1[1]+(tbb1[2]-tbb1[1])/2)/1e6,col='gray50',lty=2) +
  geom_point(size=.5) + scale_color_manual("",values=cols) + 
  facet_wrap(~group,ncol=1)+
  ggtitle('') + 
  guides(colour = guide_legend(override.aes = list(size=3)))+
  xlab('Position (Mbp)') + ylab("Tajima's D")

pdf(file=paste0(draft2,'Figure3a.pdf'))
print(p2Divbtub)
invisible(dev.off())

## GW plot
pdf(file = paste0(draft2,'Supplementary_Figure8.pdf'),width=14,height=8)
ggplot(t.HI,aes(x=bin,y=Tajima,col=Chr)) +
  theme_classic()+
  theme(legend.position='',text=element_text(size=18),legend.title = NULL) +
  geom_hline(yintercept = 0,col='gray50',lty=3) + 
  geom_vline(xintercept = 700,col='gray50',lty=2) +
  geom_point(size=.5) + scale_color_manual(values=cols) + 
  facet_wrap(~group,ncol=1) +
  ggtitle('') + 
  guides(colour = guide_legend(override.aes = list(size=3)))+
  xlab('Position (Mbp)') + ylab("Tajima's D")
dev.off()

##-- Retrieve BZ haplotypes F167Y E198A F200Y // based on SNP calling
##-- TTC > TAC / GAA > GCA / TTC > TAC
##-- AT-7025426-7029158
##-- Retrieve BZ genotype F167Y E198A F200Y // based on ANGSD Genotype Probability // does not rely on a reference genome
#- 0=A, 1=C, G=2, T=3 
setwd('./data/NUC/')
#pdf(file='BTUB_case_study.pdf')
bzgl = read.table('./data/NUC/bTublike.beagle.gz',header=T)
bzgl = bzgl[bzgl$marker %in% c('1_Celeg_TT_7027536','1_Celeg_TT_7027752','1_Celeg_TT_7027758'),]
bzgl$allele1[bzgl$allele1==0]='A'
bzgl$allele1[bzgl$allele1==1]='C'
bzgl$allele1[bzgl$allele1==2]='G'
bzgl$allele1[bzgl$allele1==3]='T'
bzgl$allele2[bzgl$allele2==0]='A'
bzgl$allele2[bzgl$allele2==1]='C'
bzgl$allele2[bzgl$allele2==2]='G'
bzgl$allele2[bzgl$allele2==3]='T'
genoGL=matrix(0,223,7)
genoGL=data.frame(genoGL)
for(i in seq(1:223)){
  n=0
  for(j in c('1_Celeg_TT_7027536','1_Celeg_TT_7027752','1_Celeg_TT_7027758')){
    n=n+1
    genoGL[i,1] = as.character(df$sample[i])
    gl = c(bzgl[bzgl$marker ==j,3*i+1],bzgl[bzgl$marker ==j,3*i+2],bzgl[bzgl$marker ==j,3*i+3]) ##- GL for maj/maj, maj/min
    gen = c(paste0(bzgl$allele1[n],'/',bzgl$allele1[n]),
            paste0(bzgl$allele1[n],'/',bzgl$allele2[n]),
            paste0(bzgl$allele2[n],'/',bzgl$allele2[n])) ##- mama,mami...
    genoGL[i,2*n] = max(gl)
    genoGL[i,2*n+1] = gen[which.max(gl)]
    rm(gl)
  }
}
colnames(genoGL) = c('sample','GL167','GT167','GL198','GT198','GL200','GT200')
genoGL$pop = df$iso[match(genoGL$sample,df$sample)]
genoGL$mcov = df$mcov[match(genoGL$sample,df$sample)]

table(genoGL$GT167) ## A/A = R 
# A/A T/A T/T 
#   1   1 221 
table(genoGL$GT198) ## C/C = R
# A/A A/C C/C 
# 208   3  12 
table(genoGL$GT200) ## A/A = R
# A/A T/A T/T 
#  52   5 166 

kruskal.test(log(mcov) ~ factor(GT198), data=genoGL)
# data:  log(mcov) by factor(GT198)
# Kruskal-Wallis chi-squared = 5.9078, df = 2, p-value = 0.05214
a=aggregate(log(mcov) ~ GT198,data=genoGL,FUN=mean)
a$std=aggregate(log(mcov) ~ GT198,data=genoGL,FUN=sd)[,2]
a
#   GT198 log(mcov)','std
# 1   A/A  0.626909 0.9976692
# 2   A/C  1.219675 0.4152245
# 3   C/C  1.135981 0.6423307
give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}
p1.1 = ggplot(genoGL,aes(y=mcov,x=GT167,fill=GT167)) + scale_y_log10() +  ylab('Mean coverage') + ggtitle('GL>0') + geom_boxplot()+theme(text=element_text(size=12),legend.position = 'none')+ stat_summary(fun.data = give.n, geom = "text", fun.y = median,position = position_dodge(width = 0.75))
p2.1 = ggplot(genoGL,aes(y=mcov,x=GT198,fill=GT198)) + scale_y_log10() +  ylab('Mean coverage') + ggtitle('GL>0') +geom_boxplot()+theme(text=element_text(size=12),legend.position = 'none')+ stat_summary(fun.data = give.n, geom = "text", fun.y = median,position = position_dodge(width = 0.75))
p3.1 = ggplot(genoGL,aes(y=mcov,x=GT200,fill=GT200)) + scale_y_log10() +  ylab('Mean coverage') + ggtitle('GL>0') +geom_boxplot()+theme(text=element_text(size=12),legend.position = 'none')+ stat_summary(fun.data = give.n, geom = "text", fun.y = median,position = position_dodge(width = 0.75))


summary(lm(log(mcov) ~ GT200, data=genoGL)) ##-- only genotype which we can test
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   1.2962     0.1271  10.198   <2e-16 ***
# GT200T/A     -0.1012     0.4292  -0.236    0.814    
# GT200T/T     -0.8486     0.1457  -5.826    2e-08 ***

# p1=ggplot(genoGL,aes(y=mcov,x=GL167,col=GT167)) + scale_y_log10() + geom_point()+theme(text=element_text(size=16))
# p2=ggplot(genoGL,aes(y=mcov,x=GL198,col=GT198)) + scale_y_log10() + geom_point()+theme(text=element_text(size=16))
# p3=ggplot(genoGL,aes(y=mcov,x=GL200,col=GT200)) + scale_y_log10() + geom_point()+theme(text=element_text(size=16))
# multiplot(p1,p2,p3,cols=1)

rcorr(genoGL$GL167,genoGL$mcov)
rcorr(genoGL$GL198,genoGL$mcov)
rcorr(genoGL$GL200,genoGL$mcov)

##- Correlation of GL with coverage 
genoGL$mGL=rowMeans(genoGL[,c('GL167','GL198','GL200')]) ## GT167 does not contribute
ggplot(genoGL,aes(x=mcov,y=mGL)) + 
  geom_smooth(method='lm') + scale_x_log10() + geom_point()+theme(text=element_text(size=16))

genoGL$status='other'
genoGL$status[genoGL$GT167=='A/A' | genoGL$GT198=='C/C' | genoGL$GT200=='A/A']='BZ-R'
genoGL$status[genoGL$GT167=='T/T' & genoGL$GT198=='A/A' & genoGL$GT200=='T/T']='BZ-S'

##- Filter out samples for which prob are not > 0.6 // GT167 does not contribute
genoGL2=genoGL[genoGL$GL167>=0.6 & genoGL$GL198>=0.6 & genoGL$GL200>=0.6,]
p1.2=ggplot(genoGL2,aes(y=mcov,x=GT167,fill=GT167)) +  scale_y_log10() + ylab('Mean coverage') + ggtitle('GL>0.60') + geom_boxplot()+theme(text=element_text(size=12),legend.position = 'none')+ stat_summary(fun.data = give.n, geom = "text", fun.y = median,position = position_dodge(width = 0.75))
p2.2=ggplot(genoGL2,aes(y=mcov,x=GT198,fill=GT198)) +   scale_y_log10() + ylab('Mean coverage') + ggtitle('GL>0.60') + geom_boxplot()+theme(text=element_text(size=12),legend.position = 'none')+ stat_summary(fun.data = give.n, geom = "text", fun.y = median,position = position_dodge(width = 0.75))
p3.2=ggplot(genoGL2,aes(y=mcov,x=GT200,fill=GT200)) +   scale_y_log10() + ylab('Mean coverage') + ggtitle('GL>0.60') + geom_boxplot()+theme(text=element_text(size=12),legend.position = 'none')+ stat_summary(fun.data = give.n, geom = "text", fun.y = median,position = position_dodge(width = 0.75))

table(genoGL2$GT167) ## A/A = R 
# A/A T/A T/T 
#   1   1  72 
table(genoGL2$GT198) ## C/C = R
# A/A A/C C/C 
#  61   3  10 
table(genoGL2$GT200) ## A/A = R
# A/A T/A T/T 
#  39   4  31 

kruskal.test(log(mcov) ~ factor(GT200), data=genoGL2)
# data:  log(mcov) by factor(GT200)
# Kruskal-Wallis chi-squared = 3.107, df = 2, p-value = 0.2115
a=aggregate(log(mcov) ~ GT200,data=genoGL2,FUN=mean)
a$std=aggregate(log(mcov) ~ GT200,data=genoGL2,FUN=sd)[,2]
a
# GT200 log(mcov) std
# 1   A/A  1.688463 0.6505582
# 2   T/A  1.241839 0.3419151
# 3   T/T  1.735819 0.3616952

kruskal.test(log(mcov) ~ factor(GT198), data=genoGL2)
# data:  log(mcov) by factor(GT198)
# Kruskal-Wallis chi-squared = 3.7472, df = 2, p-value = 0.1536
a=aggregate(log(mcov) ~ GT198,data=genoGL2,FUN=mean)
a$std=aggregate(log(mcov) ~ GT198,data=genoGL2,FUN=sd)[,2]
a
# GT198 log(mcov) std
# 1   A/A  1.368890 0.7291077
# 2   A/C  1.219675 0.4152245
# 3   C/C  1.173507 0.7012118

genoGL3=genoGL[genoGL$GL167>=0.8 & genoGL$GL198>=0.8 & genoGL$GL200>=0.8,]
p1.3=ggplot(genoGL3,aes(y=mcov,x=GT167,fill=GT167)) + scale_y_log10() + ylab('Mean coverage') + ggtitle('GL>0.80') + geom_boxplot()+theme(text=element_text(size=12),legend.position = 'none') + stat_summary(fun.data = give.n, geom = "text", fun.y = median,position = position_dodge(width = 0.75)) 
p2.3=ggplot(genoGL3,aes(y=mcov,x=GT198,fill=GT198)) + scale_y_log10() + ylab('Mean coverage') + ggtitle('GL>0.80') + geom_boxplot()+theme(text=element_text(size=12),legend.position = 'none')+ stat_summary(fun.data = give.n, geom = "text", fun.y = median,position = position_dodge(width = 0.75))
p3.3=ggplot(genoGL3,aes(y=mcov,x=GT200,fill=GT200)) + scale_y_log10() + ylab('Mean coverage') + ggtitle('GL>0.80') + geom_boxplot()+theme(text=element_text(size=12),legend.position = 'none')+ stat_summary(fun.data = give.n, geom = "text", fun.y = median,position = position_dodge(width = 0.75))

pdf(file=paste0(draft2,'Supplementary_Figure12.pdf'))
multiplot(p1.1,p1.2,p1.3,p2.1,p2.2,p2.3,p3.1,p3.2,p3.3,cols=3)
dev.off()

table(genoGL3$GT167) ## A/A = R 
# T/T 
# 31
table(genoGL3$GT198) ## C/C = R
# A/A A/C C/C 
#  25   3  3 
table(genoGL3$GT200) ## A/A = R
# A/A T/A T/T 
#  17   4  10

kruskal.test(log(mcov) ~ factor(GT198), data=genoGL3)
# data:  log(mcov) by factor(GT198)
# Kruskal-Wallis chi-squared = 3.7472, df = 2, p-value = 0.1536

##- Coverage bias on GL // does not impact on the genotype prediction  
ggplot(genoGL2,aes(x=mcov,y=mGL)) +
  geom_smooth(method='lm') + scale_x_log10() + geom_point()+theme(text=element_text(size=16))
rcorr(log(genoGL2$mcov),genoGL2$mGL)

genoGL2$F167Y[genoGL2$GT167=='A/A']='rr'
genoGL2$F167Y[genoGL2$GT167=='T/A']='rS'
genoGL2$F167Y[genoGL2$GT167=='T/T']='SS'

genoGL2$F200Y[genoGL2$GT200=='A/A']='rr'
genoGL2$F200Y[genoGL2$GT200=='T/A']='rS'
genoGL2$F200Y[genoGL2$GT200=='T/T']='SS'

genoGL2$E198A[genoGL2$GT198=='C/C']='rr'
genoGL2$E198A[genoGL2$GT198=='A/C']='rS'
genoGL2$E198A[genoGL2$GT198=='A/A']='SS'

##-- Res level

tmp=reshape2::melt(genoGL2[,c('sample','F167Y','E198A','F200Y')],1)
tmp$iso = df$iso[match(tmp$sample,df$sample)]
tmp$iso[is.na(tmp$iso)]='AUS.2'
ord = order(tmp$sample[1:74])
res = ifelse(tmp$iso[ord] %in% h2$iso[h2$res==levels(h2$res)[2]],'forestgreen','black')

p1 = ggplot(tmp,aes(x=variable,y=sample,fill=value)) +
  scale_fill_brewer(palette='RdPu',direction = -1) + xlab('') + ylab('') +
  geom_tile(col='gray') + theme(legend.position='bottom',text= element_text(size=16),
                                legend.title = element_blank(),
                                axis.text.y = element_text(size=8,colour=res))
tmp2=reshape2::melt(genoGL2[,c('sample','GL167','GL198','GL200')],1)
tmp2$iso = df$iso[match(tmp2$sample,df$sample)]
tmp2$iso[is.na(tmp2$iso)]='AUS.2'
ord = order(tmp2$sample[1:74])
res = ifelse(tmp2$iso[ord] %in% h2$iso[h2$res==levels(h2$res)[2]],'forestgreen','black')

p2=ggplot(tmp2,aes(x=variable,y=sample,fill=value)) + xlab('') + ylab('') +
  geom_tile(color='gray') +theme(legend.position='bottom',text= element_text(size=16),
                                 legend.title = element_blank(),legend.text = element_text(size=8),
                                 axis.text.y = element_text(size=8,colour=res))

pdf(file=paste0(draft2,'Supplementary_Figure13.pdf'),width = 8,height=14)
multiplot(p1,p2,cols=2)
dev.off()

##- EFY samples // 39 samples
efy=genoGL2[genoGL2$GT167=='T/T' & genoGL2$GT198=='A/A' & genoGL2$GT200=='A/A', ]
write.table(efy,file='EFY.txt',quote=F,row.names=F)
table(efy$pop[efy$mcov>=2])

##- BZ-R vs. BZ-S 
genoGL2$status='other'
genoGL2$status[genoGL2$GT167=='A/A' | genoGL2$GT198=='C/C' | genoGL2$GT200=='A/A']='BZ-R'
genoGL2$status[genoGL2$GT167=='T/T' & genoGL2$GT198=='A/A' & genoGL2$GT200=='T/T']='BZ-S'

write.table(genoGL2,file='bz_tot_cov.txt',quote=F,row.names=F)
aggregate(mcov ~ status,data = genoGL2,FUN = min)
# status     mcov
# 1   BZ-R 5.081600
# 2   BZ-S 4.761579
# 3  other 3.800000

table(genoGL2$status[genoGL2$mcov>=5])
# BZ-R  BZ-S other 
# 15     9     1 

table(genoGL2$status[genoGL2$mcov>=2.5],genoGL2$pop[genoGL2$mcov>=2.5])
# BZ-R  BZ-S other 
# 15     9     1 

##-- Haplotype frequency
genoGL2$haplotype=factor(paste0(genoGL2$GT167,'.',genoGL2$GT198,'.',genoGL2$GT200))
table(genoGL2$haplotype)
# A/A.A/A.T/T T/A.A/A.T/T T/T.A/A.A/A T/T.A/A.T/A T/T.A/A.T/T T/T.A/C.T/A T/T.C/C.T/T 
#           1           1          39           1          19           3          10 
tt=data.frame(table(genoGL2$pop,genoGL2$haplotype))
write.table(tt,file='Haplotype_Count_by_pop.txt',quote=F,row.names = F)
p3=ggplot(tt,aes(x=Var1,y=Freq,fill=Var2))+geom_bar(stat='identity')+coord_flip()+xlab('Population')+ylab('Occurrence')+
  theme(legend.title = element_blank(),text=element_text(size=16))
p3

##-- Correlation between haplotype and resistance level for french isolates
fecrt = read.csv(file='./data/btub_fecrt.csv',header=T)
fra = genoGL2[grep('FRA',genoGL2$sample),]
fra$farm = substr(fra$sample,5,7)
fra = merge(fra,fecrt,by='farm')

###----- Study of diversity based on haplotypes (from phased genotypes)
setwd('./data/NUC/')

ghpp = read.csv(file='./data/NUC/haplo_aldiffcount_matrix.tsv',header=T,sep="\t")
ghppm = as.matrix(ghpp[,-1])
ghpp = as.matrix(ghpp[,-1])
colnames(ghppm)=df$iso[match(colnames(ghppm),df$sample)]
rownames(ghppm)=df$iso[match(colnames(ghppm),df$sample)]
# pheatmap::pheatmap(ghppm,fontsize = 8)

pdf(file=paste0(draft2,'Supplementary_Figure10.pdf'),width=14,height=14)
heatmap.2(ghppm,na.rm = T,dendrogram ="col",trace = "none", Colv = "Rowv",scale = "none",
          offsetRow = 0, offsetCol = 0, #colsep = sort(idxc)-1,rowsep=sort(idxr)-1, key=F,
          col = rev(brewer.pal(9,"PuOr")),cexCol = .8,cexRow = .8, key.title = 'Dxy',
          keysize = 1.,
          xlab = "Sample", lhei = c(0.18,0.95),lwid = c(0.18,0.95))
invisible(dev.off())

#####====== Pair-wise FST between populations ANGSD ALL ====########
rep="./data/NUC/FST/"
isolist=c('ACO','AUS.1','AUS.2','BEN','CAP','FRA.1','FRA.2','FRA.3',
          'FRG','IND','MOR','NAM','STA.1','STA.2','STA.3','STO','ZAI')

# i=0
# j=0
# fstpw=matrix(0,1,3)
# fstpw=data.frame(fstpw)
# colnames(fstpw)=c('iso1','iso2','fst')
# for(iso1 in isolist){
#   i=i+1
#   j=0
#   for(iso2 in isolist){
#     j=j+1
#     if(i<j){
#fst = 0
#for(chr in chromlist){
#  tmp = read.table(paste0(rep,chr,"/",iso1,".",iso2,'.fst.txt'),header=T)
#  tmp = tmp[tmp[,3]>=100,4]
#  tmp[tmp<0] = 0
#  fst = c(fst,tmp)
#}
#fst = fst[-1]
#temp = c(iso1,iso2,fst)
#fstpw = rbind(fstpw,temp)
#rm(temp)
#     }
#   }
# }
# fstpw = fstpw[-1,]
# fstpw$fst = as.numeric(as.character(fstpw$fst))
# write.table(fstpw,paste0(rep,'fst_ALL_ANGSD_pw.txt'),quote=F,row.names = F)

fstpw = read.table(paste0(rep,'fst_ALL_ANGSD_pw.txt'),header=T)

fstpw$cov1=df$mcov[match(fstpw$iso1,df$iso)]
fstpw$cov2=df$mcov[match(fstpw$iso2,df$iso)]
fstpw$crosscov=fstpw$cov1*fstpw$cov2
#      iso1  iso2      fst cov1  cov2 crosscov
# 8     ACO   FRG 0.091135 3.70  3.84  14.2080
# 23  AUS.1   FRG 0.108655 2.59  3.84   9.9456
# 37  AUS.2   FRG 0.064266 1.63  3.84   6.2592
# 50    BEN   FRG 0.300531 0.30  3.84   1.1520
# 62    CAP   FRG 0.006139 1.06  3.84   4.0704
# 73  FRA.1   FRG 0.057500 3.83  3.84  14.7072
# 83  FRA.2   FRG 0.056219 0.62  3.84   2.3808
# 92  FRA.3   FRG 0.130263 0.47  3.84   1.8048
# 101   FRG   IND 0.179307 3.84 12.68  48.6912
# 102   FRG   MOR 0.042184 3.84  1.42   5.4528
# 103   FRG   NAM 0.099517 3.84  6.80  26.1120
# 104   FRG STA.1 0.135880 3.84  1.56   5.9904
# 105   FRG STA.2 0.085274 3.84  0.46   1.7664
# 106   FRG STA.3 0.141612 3.84  4.00  15.3600
# 107   FRG   STO 0.053186 3.84  5.01  19.2384
# 108   FRG   ZAI 0.101338 3.84  1.33   5.1072

###--- First plot for Supplementary Figure 19
pSupplFSTANGSD=ggplot(fstpw,aes(x=crosscov,y=fst)) + ggtitle('a') + scale_x_log10() +
  scale_y_continuous(limits=c(0,0.65)) +
  geom_smooth(method='lm') + geom_point() + theme(text=element_text(size=16),legend.position='bottom')+
  xlab('Cross-coverage') + ylab('FST')

rcorr(fstpw$fst,log(fstpw$crosscov),type='pearson') #- r=0.52,p=0


#####====== Pair-wise FST between populations by MAF bin ====########
setwd("./data/NUC/FST/")

complist = unique(read.table(file='./comp_list.txt',header=F))
chrlist=c('1_Celeg_TT','2_Celeg_TT','3_Celeg_TT','4_Celeg_TT','5_Celeg_TT_axel')
maf = seq(1:5)

###------------------------- Aggregate MEAN FST by MAF bin
# tmean = matrix(0,1,3)
# tmean = data.frame(tmean)
# colnames(tmean)=c('MAFbin','MEAN_FST','COMP')
# 
# for(comp in complist$V1){
#   fst = matrix(0,1,4)
#   fst = data.frame(fst)
#   colnames(fst) = c('N_VARIANTS','MEAN_FST','MAFbin','COMP')
#   
#   for(m in maf){
#     for(c in chrlist){
#       tmp = fread(paste('gunzip -c ./FST.',comp,'.',c,'.maf',maf[m],'.windowed.weir.fst.gz',sep=''),header=T)[,c(4,6)]
#       tmp = tmp[which(tmp$N_VARIANTS > 50),] ##-- 50 SNPs min
#       tmp$MEAN_FST[tmp$MEAN_FST<0] = 0
#       tmp$MAFbin = m*0.1
#       tmp$COMP = comp
#       fst=rbind(fst,tmp)
#       rm(tmp)
#     }
#   }
#   tmx = aggregate(MEAN_FST ~ MAFbin,data=fst,FUN=mean)
#   tmx$COMP = comp
#   tmean = rbind(tmean,tmx)
#   rm(tmx)
# }
# Store information
#tmean = tmean[tmean$MAFbin!=0 & tmean$COMP!=0,]
#write.table(tmean,file='FST_pairwise_MEAN_mean.txt',quote=F,row.names = F)

##-- Supplementary Figure 4
tmeanM = read.table(file='FST_pairwise_MEAN_mean.txt',header=T)

pdf(file=paste0(draft2,'Supplementary_Figure4.pdf'),height=28,width=14)
ggplot(tmeanM,aes(x=COMP,y=MEAN_FST,fill=factor(MAFbin))) +
  coord_flip() + guides(fill=guide_legend(title="MAF bin")) +
  geom_bar(stat='identity',position='dodge') +
  theme(legend.position = c(1,1),legend.justification=c(1,1))+
  ylab('FST')+ xlab('Comparison') + scale_y_continuous(limits=c(0,0.6))+
  theme(text = element_text(size=20),axis.text.y = element_text(size=16),
        axis.text.x = element_text(angle = 90, hjust = 1,size=20))
dev.off()

tm = aggregate(MEAN_FST ~ COMP,data=tmeanM,FUN=max)
tm$iso1 = ''
tm$iso2 = ''
il = unique(df$iso[df$iso!='BRA' & df$iso!='FRA.4']) ## at least 5 individuals required
for(i in il){
  v = grep(i,tm$COMP)
  for(k in v){
    if(tm$iso1[k]==''){tm$iso1[k]=i}else{tm$iso2[k]=i}
  }
}

tm$cov1 = df$mcov[match(tm$iso1,df$iso)]
tm$cov2 = df$mcov[match(tm$iso2,df$iso)]
tm$crosscov = tm$cov1*tm$cov2

###-- Supplementary Figure 19
pSupplFSTB = ggplot(tm,aes(x=crosscov,y=MEAN_FST)) + ggtitle('b') + scale_x_log10() +
  xlab('Cross-coverage') + ylab('FST') + scale_y_continuous(limits=c(0,0.65)) +
  geom_smooth(method='lm') + geom_point() + theme(text=element_text(size=16),legend.position='bottom')

pdf(file=paste0(draft2,'Supplementary_Figure19.pdf'),width=14,height=8)
multiplot(pSupplFSTANGSD,pSupplFSTB,cols=2)
dev.off()

rcorr(tm$MEAN_FST,log(tm$crosscov),type='pearson') #- r=-0.05,p=0.5481 n=136

summary.aov(lm(MEAN_FST ~ log(tm$crosscov),data=tm))
#     Df Sum Sq  Mean Sq F value Pr(>F)
# log(tm$crosscov)   1 0.0039 0.003881   0.363  0.548
# Residuals 134 1.4346 0.010706  

##-- Relationship between GUA and others
frg = tm[grep('FRG',tm$COMP),]
frg
#          COMP   MEAN_FST  iso1  iso2 cov1  cov2 crosscov
# 4     ACO.FRG 0.05753186   FRG   ACO 3.84  3.70  14.2080
# 14  AUS.1.FRG 0.11264118 AUS.1   FRG 2.59  3.84   9.9456
# 28  AUS.2.FRG 0.12864988   FRG AUS.2 3.84  1.63   6.2592
# 42    BEN.FRG 0.07552645   BEN   FRG 0.30  3.84   1.1520
# 53    CAP.FRG 0.07237772   CAP   FRG 1.06  3.84   4.0704
# 57  FRA.1.FRG 0.07313125 FRA.1   FRG 3.83  3.84  14.7072
# 62  FRA.2.FRG 0.06560010 FRA.2   FRG 0.62  3.84   2.3808
# 66  FRA.3.FRG 0.10130528 FRA.3   FRG 0.47  3.84   1.8048
# 69    IND.FRG 0.05022449   FRG   IND 3.84 12.68  48.6912
# 75    MOR.FRG 0.05836908   FRG   MOR 3.84  1.42   5.4528
# 82    NAM.FRG 0.18723109   FRG   NAM 3.84  6.80  26.1120
# 92  STA.1.FRG 0.20322438   FRG STA.1 3.84  1.56   5.9904
# 105 STA.2.FRG 0.20639718   FRG STA.2 3.84  0.46   1.7664
# 118 STA.3.FRG 0.12516787   FRG STA.3 3.84  4.00  15.3600
# 125   STO.FRG 0.10382163   FRG   STO 3.84  5.01  19.2384
# 132   ZAI.FRG 0.19798431   FRG   ZAI 3.84  1.33   5.1072

frg$comp='Atlantic'
frg$comp[!(frg$iso1 %in% c('ACO','BEN','CAP','MOR','STO')) & !(frg$iso2 %in% c('ACO','BEN','CAP','MOR','STO'))]='Other'
aggregate(MEAN_FST ~ comp,data=frg,FUN=summary)
ggplot(frg,aes(x=comp,y=MEAN_FST,col=comp)) + geom_point()
wilcox.test(MEAN_FST ~ factor(comp),data=frg)
# data:  MEAN_FST by factor(comp)
# W = 11, p-value = 0.06868
# alternative hypothesis: true location shift is not equal to 0

sto=tm[grep('STO',tm$COMP),]
sto
#','  COMP   MEAN_FST  iso1 iso2  cov1 cov2 crosscov
# 7     ACO.STO 0.35964278   ACO  STO  3.70 5.01  18.5370
# 19  AUS.1.STO 0.24359271 AUS.1  STO  2.59 5.01  12.9759
# 35  AUS.2.STO 0.34107057 AUS.2  STO  1.63 5.01   8.1663
# 47    BEN.STO 0.07431024   BEN  STO  0.30 5.01   1.5030
# 56    CAP.STO 0.39089548   CAP  STO  1.06 5.01   5.3106
# 59  FRA.1.STO 0.15683006 FRA.1  STO  3.83 5.01  19.1883
# 64  FRA.2.STO 0.43596237 FRA.2  STO  0.62 5.01   3.1062
# 68  FRA.3.STO 0.38280964 FRA.3  STO  0.47 5.01   2.3547
# 78    MOR.STO 0.19606765   MOR  STO  1.42 5.01   7.1142
# 84    NAM.STO 0.20108995   NAM  STO  6.80 5.01  34.0680
# 98  STA.1.STO 0.32353887 STA.1  STO  1.56 5.01   7.8156
# 109 STA.2.STO 0.40478384 STA.2  STO  0.46 5.01   2.3046
# 123 STA.3.STO 0.22510935 STA.3  STO  4.00 5.01  20.0400
# 125   STO.FRG 0.10382163   FRG  STO  3.84 5.01  19.2384
# 126   STO.IND 0.19521638   IND  STO 12.68 5.01  63.5268
# 136   ZAI.STO 0.30107363   STO  ZAI  5.01 1.33   6.6633
tm$sto='Other'
tm$sto[grep('STO', tm$COMP)]='STO'
ggplot(tm,aes(x=MEAN_FST,fill=sto)) + geom_density(alpha=.2)
summary.aov(lm(MEAN_FST ~ factor(sto),data=tm))
#Df Sum Sq Mean Sq F value  Pr(>F)    
# factor(sto)   1 0.1124  0.1124   11.36 0.00098 ***
#   Residuals   134 1.3261  0.0099
summary(lm(MEAN_FST ~ factor(sto),data=tm))

###------------------------- subtrop Africa: NAM, ZAI, STA
tm$af='Other'
tm$af[grep('STA', tm$COMP)]='Af'
tm$af[grep('NAM', tm$COMP)]='Af'
tm$af[grep('ZAI', tm$COMP)]='Af'

taf=na.omit(tm[tm$af=='Af',])
taf[grep('FRA',taf$COMP),]
#','    COMP   MEAN_FST  iso1  iso2 cov1 cov2 crosscov    p1    p2   sto af
# 79    NAM.FRA.1 0.23013214 FRA.1   NAM 3.83 6.80  26.0440 FRA.1   NAM Other Af
# 80    NAM.FRA.2 0.34236617 FRA.2   NAM 0.62 6.80   4.2160 FRA.2   NAM Other Af
# 81    NAM.FRA.3 0.19773104 FRA.3   NAM 0.47 6.80   3.1960 FRA.3   NAM Other Af
# 89  STA.1.FRA.1 0.19061544 FRA.1 STA.1 3.83 1.56   5.9748 FRA.1 STA.1    Af Af
# 90  STA.1.FRA.2 0.18991674 FRA.2 STA.1 0.62 1.56   0.9672 FRA.2 STA.1    Af Af
# 91  STA.1.FRA.3 0.19885025 FRA.3 STA.1 0.47 1.56   0.7332 FRA.3 STA.1    Af Af
# 102 STA.2.FRA.1 0.12407060 FRA.1 STA.2 3.83 0.46   1.7618 FRA.1 STA.2    Af Af
# 103 STA.2.FRA.2 0.06482359 FRA.2 STA.2 0.62 0.46   0.2852 FRA.2 STA.2    Af Af
# 104 STA.2.FRA.3 0.13140052 FRA.3 STA.2 0.47 0.46   0.2162 FRA.3 STA.2    Af Af
# 115 STA.3.FRA.1 0.13840276 FRA.1 STA.3 3.83 4.00  15.3200 FRA.1 STA.3    Af Af
# 116 STA.3.FRA.2 0.20435150 FRA.2 STA.3 0.62 4.00   2.4800 FRA.2 STA.3    Af Af
# 117 STA.3.FRA.3 0.19240197 FRA.3 STA.3 0.47 4.00   1.8800 FRA.3 STA.3    Af Af
# 129   ZAI.FRA.1 0.22154924 FRA.1   ZAI 3.83 1.33   5.0939 FRA.1   ZAI Other Af
# 130   ZAI.FRA.2 0.41934695 FRA.2   ZAI 0.62 1.33   0.8246 FRA.2   ZAI Other Af
# 131   ZAI.FRA.3 0.33104480 FRA.3   ZAI 0.47 1.33   0.6251 FRA.3   ZAI Other Af

summary(taf$MEAN_FST[grep('FRA',taf$COMP)])
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.06482 0.16416 0.19773 0.21180 0.22584 0.41935 

taf[grep('AUS',taf$COMP),]
#','     COMP  MEAN_FST  iso1  iso2 cov1 cov2 crosscov    p1    p2   sto af
# 17    AUS.1.NAM 0.3173206 AUS.1   NAM 2.59 6.80  17.6120 AUS.1   NAM Other Af
# 18  AUS.1.STA.2 0.1746202 AUS.1 STA.2 2.59 0.46   1.1914 AUS.1 STA.2    Af Af
# 20    AUS.1.ZAI 0.2545357 AUS.1   ZAI 2.59 1.33   3.4447 AUS.1   ZAI Other Af
# 31    AUS.2.NAM 0.3296348 AUS.2   NAM 1.63 6.80  11.0840 AUS.2   NAM Other Af
# 32  AUS.2.STA.1 0.2405443 AUS.2 STA.1 1.63 1.56   2.5428 AUS.2 STA.1    Af Af
# 33  AUS.2.STA.2 0.1981733 AUS.2 STA.2 1.63 0.46   0.7498 AUS.2 STA.2    Af Af
# 34  AUS.2.STA.3 0.2009236 AUS.2 STA.3 1.63 4.00   6.5200 AUS.2 STA.3    Af Af
# 36    AUS.2.ZAI 0.3098561 AUS.2   ZAI 1.63 1.33   2.1679 AUS.2   ZAI Other Af
# 86  STA.1.AUS.1 0.2121598 AUS.1 STA.1 2.59 1.56   4.0404 AUS.1 STA.1    Af Af
# 112 STA.3.AUS.1 0.1743368 AUS.1 STA.3 2.59 4.00  10.3600 AUS.1 STA.3    Af Af

summary(taf$MEAN_FST[grep('AUS',taf$COMP)])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1743  0.1989  0.2264  0.2412  0.2960  0.3296 

#####====== Pair-wise 1-IBS distance  ====########
setwd("./data/NUC/FST/")
for(chrom in seq(1:5)){
  m=matrix(readBin(paste('hamming.',chrom,'.mdist.bin',sep=''),what='numeric',n=70225),265,265)
  names=read.table(file=paste('hamming.',chrom,'.mdist.id',sep=''))
  colnames(m)=names$V1
  rownames(m)=names$V2
  assign(paste0('d',chrom),m)
}
##- Average 1-Ibs over the whole genome
my_list=list(d1,d2,d3,d4,d5)
d=apply(simplify2array(my_list), 1:2, mean)

##-- Relationship with coverage
md=melt(d)
#md=md[md$value!=0,]
md$c1=df$mcov[match(md$Var1,df$sample)]
md$c2=df$mcov[match(md$Var2,df$sample)]
md$crosscov=md$c1*md$c2

ggplot(md,aes(x=crosscov,y=value))+
  geom_point(size=.1)+xlab('Cross coverage')+ylab('1-IBS') + scale_x_log10() +
  geom_smooth(method='lm')
rcorr(md$crosscov,md$value,type='pearson') ##- r=0.35, p=0
# x    y
# x 1.00 0.31
# y 0.31 1.00
# 
# n
# x     y
# x 49729 49729
# y 49729 70219
# 
# P
# x  y 
# x     0
# y  0   

p2=ggplot(md[md$c1>=2.5 & md$c2>=2.5,],aes(x=crosscov,y=value))+
  geom_point(size=.5)+xlab('Cross coverage')+ylab('1-IBS') + scale_x_log10() +
  geom_smooth(method='lm') + ggtitle('b')

p1=ggplot(md[md$c1<2.5 & md$c2<2.5,],aes(x=crosscov,y=value)) +
  geom_point(size=.5)+xlab('Cross coverage')+ylab('1-IBS') + scale_x_log10() +
  geom_smooth(method='lm') + ggtitle('a')

pdf(file=paste0(draft2,'Supplementary_Figure20.pdf'))
multiplot(p1,p2,cols=1)
dev.off()

rcorr(md$crosscov[md$c1>=2.5 & md$c2>=2.5],md$value[md$c1>=2.5 & md$c2>=2.5],type='pearson')
# x    y
# x 1.00 0.02
# y 0.02 1.00
#
# n= 8281
#
#
# P
# x      y
# x  0.0934
# y 0.0934

#####====== FST vs. Dxy 2.5X ====########
#- Replace lower-triangle of d matrix by average FST 

df3 = df[df$mcov>=2.5 & df$iso!='FRA.4' & df$iso!='BRA',] ## get rid of BRA and FRA.4
tm$p1=tm$iso1
tm$p2=tm$iso2
vec = df3$sample[order(df3$lon)]
o = match(vec,factor(colnames(d)))
#o2 = match(vec,factor(colnames(dmx)))
gd3 = as.matrix(d[o,o])

##-- STO estimates vs. others
md3 = melt(gd3)
md3$sto='Other'
md3$sto[c(grep('STO',md3$Var1),grep('STO',md3$Var2))]='STO'

summary.aov(lm(value ~ sto,data=md3))
#          Df Sum Sq Mean Sq F value Pr(>F)    
#   sto     1  1.568  1.5676   580.9 <2e-16 ***
#   Residuals   7567 20.419  0.0027  

summary(lm(value ~ sto,data=md3))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
#   (Intercept) 0.2703657  0.0006184   437.2   <2e-16 ***
#   stoSTO      0.0572523  0.0023754    24.1   <2e-16 ***

ggplot(md3,aes(x=value,fill=sto))+geom_density()

###--- Phylogeography
# dmx3 = dmx[o2,o2]
# mtest2 = mantel(gd3,dmx3)
# mtest2
# Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel(xdis = gd3, ydis = dmx3) 
# 
# Mantel statistic r: 0.1218 
# Significance: 0.003 
# 
# Upper quantiles of permutations (null model):
# 90%    95%  97.5%    99% 
# 0.0573 0.0769 0.0950 0.1058 
# Permutation: free
# Number of permutations: 999

##-- Plot matrix of FST vs. Dxy
colnames(gd3) = colnames(d)[o]
rownames(gd3) = rownames(d)[o]
for(i in 1:dim(gd3)[1]){ #- rows
  for(j in 1:dim(gd3)[1]){ #- col
    if(j < i){
      isor = df$iso[match(rownames(gd3)[i],df$sample)]
      isoc = df$iso[match(colnames(gd3)[j],df$sample)]
      if(isor!=isoc & isor!='BRA' & isoc!='BRA' & isor!='FRA.4' & isoc!='FRA.4'){ #- lower values
        gd3[i,j] = tm$MEAN_FST[(tm$p2==isor & tm$p1==isoc) | 
                                 (tm$p1==isor & tm$p2==isoc)] }else{
                                   gd3[i,j]=0}
      if(isor=='BRA' | isoc=='BRA' | isor=='FRA.4' | isoc=='FRA.4'){
        gd3[i,j]=NA
      }
    }
  }
}

##--- Rename row and cols
gdr3=df$iso[match(rownames(gd3),df$sample)]
gdc3=df$iso[match(colnames(gd3),df$sample)]
idxr=match(unique(gdr3),gdr3)
idxc=match(unique(gdc3),gdc3)
vec=array('',dim(gd3)[1])
vec[idxc]=gdr3[idxc]
gd2.3=gd3

## Replace FRG with GUA
vec_renamed = gsub('FRG','GUA',vec)
colnames(gd2.3)=vec_renamed
rownames(gd2.3)=vec_renamed

levels(cut(array(gd2.3[lower.tri(gd2.3)]),8)) ##-- FST
# [1] "(-0.00047,0.0588]" "(0.0588,0.118]"    "(0.118,0.176]"     "(0.176,0.235]"     "(0.235,0.294]"     "(0.294,0.353]"    
# [7] "(0.353,0.411]"     "(0.411,0.471]"    

levels(cut(array(gd2.3[upper.tri(gd2.3)]),8)) ##-- 1-IBS
# [1] "(0.107,0.144]" "(0.144,0.18]"  "(0.18,0.216]"  "(0.216,0.252]" "(0.252,0.288]" "(0.288,0.324]" "(0.324,0.36]"  "(0.36,0.397]" 

pdf(file=paste0(draft2,'Figure2b.pdf'),width=14,height=8)
heatmap.2(gd2.3,na.rm = T,dendrogram = "none",Rowv = NULL,Colv = "Rowv",trace = "none",scale = "none",
          offsetRow = 0, offsetCol = 0,colsep = sort(idxc)-1,rowsep=sort(idxr)-1, key=F,
          col = rev(brewer.pal(9,"YlGnBu")),cexRow = 1,cexCol = 1,
          xlab = "Between-group FST", ylab =  "1 - IBS",lhei = c(0.18,0.95),lwid = c(0.18,0.95))
gradient.rect(0.04,0.25,0.06,0.75,nslices = 9,border = F,gradient = "y",col = rev(brewer.pal(9,"YlGnBu")))
text(x = rep(0.035,9), y = seq(0.25,.75,by = 0.06),adj = 1,cex = 0.75,
     labels = c("0.00","0.06","0.12","0.18","0.24","0.29","0.35","0.41","0.47"))
gradient.rect(0.25,0.97,0.75,0.95,nslices = 9,border = F,gradient = "x",col = rev(brewer.pal(9,"YlGnBu")))
text(y = rep(0.99,9), x = seq(0.25,.75,by = 0.06),adj = 1,cex = 0.8,
     labels = c("0.00","0.14","0.18","0.22","0.25","0.29","0.32","0.36","0.40"))
dev.off()

###-- NJ Tree for Supplementary Figure 2
a = nj(d[o,o])
a$tip.label = df$iso[match(a$tip.label,df$sample)]
a$tip.label = gsub('FRG','GUA',a$tip.label)
cls=list(c1=c("FRA.1","FRA.2","FRA.3","FRA.4","POR","MOR"),
         c2=c("IND", "AUS.1","AUS.2"),
         c3=c("BRA","GUA"),
         c4=c("STA.1","STA.2","STA.3","NAM","ZAI"),
         c5=c("BEN", "CAP","STO"))

tree = groupOTU(a,cls)
while(0 %in% levels(attributes(tree)$group)){
  attributes(tree)$group[attributes(tree)$group==0]=attributes(tree)$group[which(attributes(tree)$group==0)+1]
  attributes(tree)$group=factor(attributes(tree)$group)
}
pdf(file=paste0(draft2,'Supplementary_Figure2.pdf'),height=14,width=14)
print(ggtree(tree,aes(color=group),layout='circular') + 
        geom_tiplab2(size=4) +
        scale_colour_manual(values=clcol,name="Location",breaks=cls))
dev.off()

#####====== MSMC Ne estimation ====########

###--- Need multi chrom + 1 per pop + bootstrap
gen <- 40/365 #- Min generation time
mu = 2.7e-9

##-- MSMC analysis
setwd("./data/NUC/MSMC/")
lisestim =list.files(".", pattern="final.txt",recursive = TRUE)

msmc = NULL 
popms = unique(substr(lisestim,1,nchar(lisestim)-16))
for(p0 in popms){
  files = lisestim[grep(p0,lisestim)]
  n=0
  for(f in files){
    n = n+1
    temp = read.table(paste0("./",f,sep=""),header=TRUE)
    temp$rep = n
    temp$POP = p0
    #temp$geo = df$geo[match(p0,df$reg)]
    msmc = rbind(msmc,temp)
  }
}
msmc$Ne = (1/msmc$lambda)/(2*mu)/1e6
msmc$Time = msmc$left_time_boundary/mu*gen
msmc = msmc[msmc$time_index!=0,]

msplot = aggregate(Ne ~ time_index + POP,data=msmc,FUN=mean)
msplot$Time = aggregate(Time ~ time_index + POP,data=msmc,FUN=mean)[,3]
tvec=unique(msmc$time_index)
mi=array(0,length(tvec)*length(popms))
ma=array(0,length(tvec)*length(popms))
i=0
for(p0 in popms){
  for(t in tvec){
    i=i+1
    mi[i]=min(msmc$Ne[msmc$time_index==t & msmc$POP==p0])
    ma[i]=max(msmc$Ne[msmc$time_index==t & msmc$POP==p0])
  }
}
msplot$minNe=mi
msplot$maxNe=ma

##-- Mean Ne estimate
msplot[msplot$time_index==2,]

#### MSMC results from 4 haplotypes 
# time_index   POP    Ne Time  minNe maxNe              Origin
# 2            2 FRA.1 0.112  452 0.1080 0.121  Mediterranean.Area
# 33           2   GUA 0.105  633 0.0859 0.110       South.America
# 64           2   IND 0.233  509 0.2234 0.246             Oceania
# 95           2   MOR 0.252  621 0.2322 0.273  Mediterranean.Area
# 126          2   NAM 0.610  415 0.5604 0.638 Sub.Tropical.Africa
# 157          2 STA.1 0.193  320 0.1660 0.237 Sub.Tropical.Africa

msplot$Origin=df$geo[match(msplot$POP,df$iso)]

##- Get rid of latest estimates that have strange behaviour
msplot = msplot[msplot$time_index>1,]

msplot$POP= factor(gsub('FRG','GUA',msplot$POP))

pdf(file=paste0(draft2,'Figure1c.pdf'))

divcol = c('#2166ac','#006d2c','#b2182b','#8073ac','#542788','#e08214','#2d004b','#fdb863')
ggplot(msplot, aes(x=Time, y=Ne,
                   group=POP)) +
  geom_ribbon(aes(ymin = minNe, ymax = maxNe,fill=POP), alpha=0.2) + 
  geom_line(aes(col=POP)) + 
  scale_y_log10(limits=c(0.05,.8),breaks=c(0.05,0.1,0.2,0.4,0.8)) +
  scale_x_log10(limits=c(200,500000),breaks=c(200,500,1000,2500,5000,10000,25000,50000,125000,250000,500000)) +
  scale_colour_manual(values=divcol) +
  scale_fill_manual(values=divcol) +
  theme(text=element_text(size=16),
        axis.text.x = element_text(angle=45, vjust=0.5)) +
  geom_text(data = subset(msplot, time_index == 2), 
                  aes(label = POP, colour = POP, x = 290, y = Ne), hjust = 1.2) +
  xlab('Years ago') +
  ylab('Ne in millions') +
  theme(legend.position = 'none')
dev.off()

#####====== dadi inferences ====########
setwd("./data/DADI/")
genY = 365/40
mu = 2.7e-9
pairs = list.files(path = "./dadi_sfs/")
##-- Info on n Sites considered: totSites Site Fraction considered by Angsd
nSites = read.table(file = './Summary_nSites.txt',header=F)
nSites$pop = paste0(nSites$V1,'-',nSites$V2)
colnames(nSites)[3] = 'totSites'
colnames(nSites)[4] = 'fracSites'

#pairs = pairs[-length(pairs)] 
dadi = NULL

for(p in pairs){
  temp = read.table(file=paste0('./dadi_sfs/',p,'/Results_Summary_Short.txt'),header=T)
  temp$pop = p
  dadi = rbind(dadi,temp)
}

## Add nsites info
idx = match(dadi$pop,nSites$pop)
dadi$L = gensize*nSites$fracSites[idx]/nSites$totSites[idx] ## L = fraction of Sites considered * gensize
## Discard models with no convergence 
dadi = dadi[dadi$chi.squared!=Inf & dadi$theta > 30 & dadi$log.likelihood < -5,]

###--- Define comparisons of interest
comp = c('AUS.1-AUS.2','AUS.1-FRA.1','AUS.1-STA.1','AUS.1-STA.3',
          'AUS.2-FRA.1','AUS.2-STA.1','AUS.2-STA.3',
          'FRA.1-MOR','FRA.1-NAM','FRA.1-STA.1','FRA.1-STA.3',
          'FRG-MOR','FRG-NAM','FRG-STO',
          'MOR-STO','NAM-STO','STA.1-STO','STA.3-STO')

dadi = dadi[dadi$pop %in% comp,]

###--- First identify pop with simple histories: split & nomig
dadi_firstpass = dadi[dadi$Model %in% c('no_mig','sym_mig','asym_mig'),]

## Select best model
require(dplyr)
retained_firstpass = dadi_firstpass %>% 
  group_by(pop) %>%
  slice(which.min(AIC))

retained_firstpass$Nref = as.numeric(as.character(retained_firstpass$theta))/retained_firstpass$L/4/mu
retained_firstpass$Time1 = NA

v0 = c('asym_mig','sym_mig','no_mig')
retained_firstpass$Time1[retained_firstpass$Model %in% v0] = 2019-2*retained_firstpass$Nref[retained_firstpass$Model %in% v0]/genY*as.numeric(as.character(sapply(str_split(retained_firstpass$optimized_params[retained_firstpass$Model %in% v0],','),function(x) x[length(x)])))

print.data.frame(retained_firstpass[grep('no_mig',retained_firstpass$Model),c('Model','pop','Time1')])

#    Model         pop    Time1
# 1 no_mig AUS.1-AUS.2 1895.089
# 2 no_mig AUS.1-FRA.1 1794.246
# 3 no_mig AUS.2-STA.3 1871.397
# 4 no_mig   FRA.1-NAM 1274.930
# 5 no_mig FRA.1-STA.3 1288.020
# 6 no_mig     FRG-STO 1604.407

simple_pops = retained_firstpass$pop[grep('no_mig',retained_firstpass$Model)]

###--- Scale uncertainty estimates
gim = read.csv(file='./gim_simple_model.txt',header=T,sep='\t')
for(p in gim$pop){
  gim$std.T[gim$pop==p] = 2*retained_firstpass$Nref[retained_firstpass$pop == p]/genY*gim$T[gim$pop==p]
}
gim
#           pop    nu1   nu2       T  theta     std.T
# 1 AUS.1-AUS.2 0.0474  3.74 0.01070  18300  132.5843
# 2 AUS.1-FRA.1 0.9570 15.40 0.33400 144000 7506.7964
# 3 AUS.2-STA.3 0.2680 47.40 0.20200 131000 2981.5775
# 4     FRG-STO 0.0388  6.78 0.00417  28100  167.8498

###--- Pop with migrations
dadi_complex = dadi[!(dadi$pop %in% simple_pops),]

retained_complex = dadi_complex %>% 
  group_by(pop) %>%
  slice(which.min(AIC))

retained_complex$Nref = retained_complex$theta/retained_complex$L/4/mu
retained_complex$Time1 = NA
retained_complex$Time2 = NA
retained_complex$Time3 = NA

v0 = c('sym_mig','asym_mig','no_mig')
v1 = c('sec_contact_sym_mig','sec_contact_asym_mig','anc_asym_mig','anc_sym_mig','founder_sym',
       'sec_contact_sym_mig_size','sec_contact_asym_mig_size','anc_asym_mig_size','anc_sym_mig_size',
       'sym_mig_size','asym_mig_size')
v2 = c('sec_contact_asym_mig_three_epoch','sec_contact_sym_mig_three_epoch',
       'sec_contact_sym_mig_size_three_epoch','sec_contact_asym_mig_size_three_epoch') 

retained_complex$Time1[retained_complex$Model %in% v0] = 2019-2*retained_complex$Nref[retained_complex$Model %in% v0]/genY*as.numeric(as.character(sapply(str_split(retained_complex$optimized_params[retained_complex$Model %in% v0],','),function(x) x[length(x)])))

retained_complex$Time2[retained_complex$Model %in% v1] = 2019-2*retained_complex$Nref[retained_complex$Model %in% v1]/genY*as.numeric(as.character(sapply(str_split(retained_complex$optimized_params[retained_complex$Model %in% v1],','),function(x) x[length(x)])))
retained_complex$Time1[retained_complex$Model %in% v1] = retained_complex$Time2[retained_complex$Model %in% v1]-2*retained_complex$Nref[retained_complex$Model %in% v1]/genY*as.numeric(as.character(sapply(str_split(retained_complex$optimized_params[retained_complex$Model %in% v1],','),function(x) x[length(x)-1])))

retained_complex$Time3[retained_complex$Model %in% v2] = 2019-2*retained_complex$Nref[retained_complex$Model %in% v2]/genY*as.numeric(as.character(sapply(str_split(retained_complex$optimized_params[retained_complex$Model %in% v2],','),function(x) x[length(x)])))
retained_complex$Time2[retained_complex$Model %in% v2] = retained_complex$Time3[retained_complex$Model %in% v2]-2*retained_complex$Nref[retained_complex$Model %in% v2]/genY*as.numeric(as.character(sapply(str_split(retained_complex$optimized_params[retained_complex$Model %in% v2],','),function(x) x[length(x)-1])))
retained_complex$Time1[retained_complex$Model %in% v2] = retained_complex$Time2[retained_complex$Model %in% v2]-2*retained_complex$Nref[retained_complex$Model %in% v2]/genY*as.numeric(as.character(sapply(str_split(retained_complex$optimized_params[retained_complex$Model %in% v2],','),function(x) x[length(x)-2])))

print.data.frame(retained_complex[,c('Model','pop','Time1','Time2','Time3')])
#                                    Model         pop        Time1       Time2    Time3
# 1        sec_contact_sym_mig_three_epoch AUS.1-STA.1   -31486.910   1223.8802 1808.155
# 2                       anc_sym_mig_size AUS.1-STA.3  -424612.737  -3160.7085       NA
# 3                       anc_sym_mig_size AUS.2-FRA.1     1817.574   2017.4245       NA
# 4   sec_contact_sym_mig_size_three_epoch AUS.2-STA.1    -2084.661  -1569.3450 1824.456
# 5  sec_contact_asym_mig_size_three_epoch   FRA.1-MOR   -19458.130 -19123.2591 1843.832
# 6                           sym_mig_size FRA.1-STA.1 -1594604.338 -25365.6822       NA
# 7               sec_contact_sym_mig_size     FRG-MOR     1377.764   1740.2480       NA
# 8                          asym_mig_size     FRG-NAM   -12697.085   1607.5821       NA
# 9        sec_contact_sym_mig_three_epoch     MOR-STO    -1268.039   -833.8575 1876.546
# 10                           anc_sym_mig     NAM-STO   -49483.882   1419.5480       NA
# 11                           anc_sym_mig   STA.1-STO  -530957.568  -3260.1399       NA
# 12                         asym_mig_size   STA.3-STO   -39828.109   1559.8218       NA

###--- Most likely migration models retained for FRA.1-NAM and FRA.1-STA.3 fit better suspected history
retained = dadi[dadi$pop=='FRA.1-NAM'|dadi$pop=='FRA.1-STA.3',] %>% 
  group_by(pop) %>%
  slice(which.min(AIC))

retained$Nref = retained$theta/retained$L/4/mu
retained$Time1 = NA
retained$Time2 = NA
retained$Time3 = NA

v0 = c('sym_mig','asym_mig','no_mig')
v1 = c('sec_contact_sym_mig','sec_contact_asym_mig','anc_asym_mig','anc_sym_mig','founder_sym',
       'sec_contact_sym_mig_size','sec_contact_asym_mig_size','anc_asym_mig_size','anc_sym_mig_size',
       'sym_mig_size','asym_mig_size')
v2 = c('sec_contact_asym_mig_three_epoch','sec_contact_sym_mig_three_epoch',
       'sec_contact_sym_mig_size_three_epoch','sec_contact_asym_mig_size_three_epoch') 

retained$Time1[retained$Model %in% v0] = 2019-2*retained$Nref[retained$Model %in% v0]/genY*as.numeric(as.character(sapply(str_split(retained$optimized_params[retained$Model %in% v0],','),function(x) x[length(x)])))

retained$Time2[retained$Model %in% v1] = 2019-2*retained$Nref[retained$Model %in% v1]/genY*as.numeric(as.character(sapply(str_split(retained$optimized_params[retained$Model %in% v1],','),function(x) x[length(x)])))
retained$Time1[retained$Model %in% v1] = retained$Time2[retained$Model %in% v1]-2*retained$Nref[retained$Model %in% v1]/genY*as.numeric(as.character(sapply(str_split(retained$optimized_params[retained$Model %in% v1],','),function(x) x[length(x)-1])))

retained$Time3[retained$Model %in% v2] = 2019-2*retained$Nref[retained$Model %in% v2]/genY*as.numeric(as.character(sapply(str_split(retained$optimized_params[retained$Model %in% v2],','),function(x) x[length(x)])))
retained$Time2[retained$Model %in% v2] = retained$Time3[retained$Model %in% v2]-2*retained$Nref[retained$Model %in% v2]/genY*as.numeric(as.character(sapply(str_split(retained$optimized_params[retained$Model %in% v2],','),function(x) x[length(x)-1])))
retained$Time1[retained$Model %in% v2] = retained$Time2[retained$Model %in% v2]-2*retained$Nref[retained$Model %in% v2]/genY*as.numeric(as.character(sapply(str_split(retained$optimized_params[retained$Model %in% v2],','),function(x) x[length(x)-2])))

print.data.frame(retained[,c('Model','pop','Time1','Time2','Time3')])
#                                   Model         pop     Time1     Time2      Time3
# 1 sec_contact_asym_mig_size_three_epoch   FRA.1-NAM -122242.6 -121751.9   1893.274
# 2       sec_contact_sym_mig_three_epoch FRA.1-STA.3 -228932.4 -185319.5 -17598.653

##-- output model parameters for Supplementary Table 3
# Add simple models
fout = data.frame(retained_firstpass[grep('no_mig',retained_firstpass$Model),])
# Remove FRA.1-NAM and FRA.1-STA.3 for which more complex models are more plausible
fout = fout[fout$pop!='FRA.1-NAM'& fout$pop!='FRA.1-STA.3',]

# Add more complex models
fout_complex = data.frame(retained_complex)
# Add FRA.1-NAM & FRA.1-STA.3 models
fout_complex = rbind(fout_complex,data.frame(retained))

##-- Output
write.csv(fout,file=paste0(draft2,'supplementary_tableDadi_nomig.csv'),quote=F,row.names=F)
write.csv(fout_complex,file=paste0(draft2,'supplementary_tableDadi_complex.csv'),quote=F,row.names=F)

#####====== MITO tree BEAST  ====########
options(digits = 3)
require(ggtree)
x = read.beast("./data/MITO/HKY_equalrate_nopart_strict_Cele50M_muscle.MCC")

clcol = c('#5ab4ac','#ca0020','#d8b365','#252525','#8c510a')
cols <- ggtree::scale_color(x, by="height")
pop = sapply(str_split(sapply(str_split(fortify(x)$label,"-"), function(x) x[2]),"_"),function(x) x[1])
pop[pop=='ITA']='AUS' ## Rename Australian isolate from Italian lab

cols = clcol[match(df$geo[match(pop,df$reg)],levels(df$geo))]
pop = gsub('FRG','GUA',pop)
pdf(file=paste0(draft2,'Supplementary_Figure7.pdf'),width=8,height=14)
ggtree(x, right=TRUE, mrsd="2011-01-01", color="grey") + theme_tree2() +
  geom_text(aes(x=max(x), label=pop), size=2.5, color=cols, hjust=-.3) +
  scale_x_continuous(breaks=c(-151000,-100000,-50000,0),labels = c(-151000/40,-100000/40,-50000/40,0)) +
  geom_segment(aes(xend=max(x)+200, yend=y), linetype="dotted", size=.1, color=cols) +
  theme(panel.grid.major   = element_line(color="black", size=.2),
        panel.grid.minor   = element_line(color="grey", size=.2),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) 
dev.off()

#####====== XP-CLR for global view on diversifying selection in Africa + Australia ====########
##-- All comparisons
setwd("./data/NUC/XPCLR/")

# xpclr = matrix(0,1,8)
# for(c in seq(1:5)){
#   i=0
#   ##-- AFRICA
#   for(p1 in c("KOK","MAR","NAM","WHR","ZAI")){
#     i=i+1
#     j=0
#     for(p2 in c("KOK","MAR","NAM","WHR","ZAI")){
# j=j+1
# if(j!=i){
#  temp = read.table(file=paste('CHROM',c,'.',p1,'.',p2,'.txt',sep=''),header=F,sep=' ')
#  temp$V1 = c
#  temp = temp[order(temp$V4),]
#  temp$V8 = paste(p1,p2,sep='.')
#  xpclr = rbind(xpclr,temp)
#}
#     }
#   }
# }
# for(c in seq(1:5)){
#   i=0
#   ##-- AUSTRALIA
#   for(p1 in c("AUS.1","AUS.2")){
#     i=i+1
#     j=0
#     for(p2 in c("AUS.1","AUS.2")){
#j=j+1
#if(j!=i){
#  temp = read.table(file=paste('CHROM',c,'.',p1,'.',p2,'.txt',sep=''),header=F,sep=' ')
#  temp$V1 = c
#  temp = temp[order(temp$V4),]
#  temp$V8 = paste(p1,p2,sep='.')
#  xpclr = rbind(xpclr,temp)
#}
#     }
#   }
# }
# xpclr = xpclr[-1,]
# xpclr = xpclr[order(xpclr$V1,xpclr$V4),]
# xpclr$V8= factor(xpclr$V8)
# #xpclr = unique(xpclr)
# chrom_end=array(0,5)
# for(c in seq(1:5)){chrom_end[c+1] = chrom_end[c] + max(xpclr$V4[xpclr$V1==c])}
# for(c in seq(1:5)){
#   xpclr$pos[xpclr$V1==c] = xpclr$V4[xpclr$V1==c] + chrom_end[c]
# }
# write.table(xpclr,file='xpclr_tot.txt',quote=F,row.names=F)

xpclr = read.table(file='xpclr_tot.txt',header=T)
chrom_end = array(0,5)
for(c in seq(1:5)){chrom_end[c+1] = chrom_end[c] + max(xpclr$V4[xpclr$V1==c])}
# 
# RMS <- function(num) cbind(num$V1,sqrt(mean((num$V6)^2)))
# 
#########------- Global hits of diversifying selection across African and Australian populations
# 
# loci = length(unique(xpclr$pos))
# loc = unique(xpclr$pos)
# rms = array(0,loci)
# 
# df_list = split(xpclr[,c('V1','V6')],xpclr$pos)
# rmstALL = lapply(df_list,RMS)
# rmst = lapply(rmstALL,unique)
# 
# loc = as.integer(as.character(names(rmst)))
# chromosome = sapply(rmst,function(x) x[1])
# selection_score = sapply(rmst,function(x) as.numeric(as.character(x[2])))
# score = cbind(loc, chromosome,selection_score)
# score = data.frame(score)
# colnames(score)=c('loc','chromosome','selection_score')
# 
# write.table(score,file='xpclr_selection_scores_tot.txt',row.names=F,quote=F)

score = read.table(file='xpclr_selection_scores_tot.txt',header=T)
candid = NULL
# candid$loc=c(40195769,23169338+chrom_end[3],8721268,
#              16033820+chrom_end[3],3329155+chrom_end[3],40297965+chrom_end[3],
#              9596364,
#              47201208+chrom_end[5])
# candid$chromosome = c(1,3,1,3,3,3,1,5)
# candid$gene = c('glc-4','tax-4','daf-36','snf-9','unc-24','aex-3','B0361.4','pgp-11')

##- Convert to whole genome positions
candid$loc=c(40195769,
             16033820+chrom_end[3],3329155+chrom_end[3],40297965+chrom_end[3],
             9596364,
             47201208+chrom_end[5])
candid$chromosome = c(1,3,3,3,1,5)
candid$gene = c('glc-4','snf-9','unc-24','B0361.4','aex-3','pgp-11')

candid = data.frame(candid)

pXP_ALL = ggplot(score,aes(x=loc/1e6,y=selection_score,col=factor(chromosome))) +
  #geom_vline(xintercept=tbb1[1]/1e6, lty=2,lwd= .1) +
  geom_point(size=.5) +
  geom_hline(yintercept = quantile(score$selection_score,0.999),col='red') +
  scale_colour_manual(values=chrom.colors) +
  geom_label_repel(data=candid,size=1.5,aes(x=loc/1e6,y=30,col=factor(chromosome),label=gene)) +
  geom_segment(data=candid,aes(x=loc/1e6,y=0,xend=loc/1e6,yend=30,col=factor(chromosome)),size=0.2) +
  scale_y_continuous(limits=c(0,max(score$selection_score))) + xlab('Genomic position (Mbp)') +
  scale_x_continuous(limits=c(0,max(score$loc/1e6)),breaks = seq(0,240,10)) +
  theme(legend.position='none',text=element_text(size=8),axis.text=element_text(size=6)) 

tiff(file=paste0(draft2,'Supplementary_Figure14.tif'),width = 120, height = 60,
     units = "mm", res = 600,compression="lzw")
print(pXP_ALL)
dev.off()

###-- Steve's peak
chr5 = score[score$chromosome==5,]
chr5$loc2 = chr5$loc - chrom_end[5]
chr5[chr5$selection_score >quantile(score$selection_score,0.999) & chr5$loc2>38e6 & chr5$loc2<42e6, ]
#               loc chromosome selection_score     loc2
# 1673385 228345418          5        5.906145 39848778
# 1673400 228347418          5        5.906122 39850778
# 1682094 229769418          5        5.777090 41272778
# 1682106 229771418          5        6.730669 41274778
# 1682118 229773418          5        5.633643 41276778
# 1682394 229819418          5        5.476368 41322778
# 1682402 229819803          5        5.532610 41323163
# 1682414 229821803          5        5.532668 41325163
# 1683063 229929805          5        5.852615 41433165
# 1683066 229931418          5        6.644573 41434778
# 1683366 229981418          5        6.281684 41484778
# 1683731 230041743          5        5.980375 41545103

###-- All hits for blast and inference of underlying genes
##- Blast and gene inference is done using the convert_hits_to_gene.sh script 
top = score[score$selection_score>quantile(score$selection_score,0.999),] #- 1740 hits / 840 within genes (534 unique)
c = top$chromosome
up = top$loc-1000-chrom_end[top$chromosome]
dw = top$loc+1000-chrom_end[top$chromosome]
t = cbind(c,up,dw)
t = data.frame(t)
write.table(t,file='./all_hits_allcomp_01p.tsv',row.names=F,sep='\t')

###--- Hcontortus-based GO enrichment analysis for ALL top 0.1% hits 
require(topGO)

# myInterestingGenes = read.table(file='./gene_list.txt',header=F,sep='\t')
# colnames(myInterestingGenes)[1]='V1'
# geneNames = read.table(file='./blastALL01p/HcoGeneList.txt',header=T,sep='\t')
# geneList <- factor(as.integer(geneNames$Gene.stable.ID %in% myInterestingGenes$V1))
# names(geneList) <- geneNames$Gene.stable.ID
# str(geneList)
# geneID2GO = readMappings(file=paste0(draft,'Hco_GO_topGO.txt')) #-- Total GO
# 
# ##- Attach gene score
# hb = read.table(file='./tot_geneid.txt',sep=' ',header=F,fill=T)
# colnames(hb)[1]='loc'
# t$loc[t$c!=5]=paste0('top_',t$c[t$c!=5],'_Celeg_TT:',t$up[t$c!=5],'-',t$dw[t$c!=5],'.fasta.blastn')
# t$loc[t$c==5]=paste0('top_',t$c[t$c==5],'_Celeg_TT_axel:',t$up[t$c==5],'-',t$dw[t$c==5],'.fasta.blastn')
# t$loc=factor(t$loc)
# t = cbind(t,top$selection_score)
# t2 = merge(t,hb,by='loc',all=T)
# colnames(t2)[5]='selection_score'
# 
# a = aggregate(selection_score ~ V3,data=t2,FUN=mean)
# a = a[-1,]
# genescore = a$selection_score[match(a$V3,myInterestingGenes$V1)]
# names(genescore) <- a$V3[match(a$V3,myInterestingGenes$V1)]
# generank = order(genescore,decreasing=T)
# names(generank) = names(genescore)[order(genescore,decreasing=T)]
# 
# ##-- Molecular / Cellular Componenent / Biological Processes //// classicFisher based on geneCounts; KS based on gene scores/ranks
# goterm=c('MF','CC','BP')
# for(i in goterm){
#   GOdata <- new("topGOdata", ontology = i, allGenes = geneList, geneSel = genescore,
#   nodeSize = 10,
#   annot = annFUN.gene2GO, gene2GO = geneID2GO)
#   ##- Fisher
#   resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
#   
#   ##- KS test <> based on gene score + weight01 = account for hierarchy <> stat correction
#   resultKS <- runTest(GOdata, algorithm = "weight01", statistic = "ks")
#   
#   #The elim method was design to be more conservative than the classic method
#   resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
#   
#   #- Table of all results - Classic Fisher
#   topn=0
#   if(sum(score(resultKS) < .05)>1){topn = sum(score(resultKS) < .05)}else{topn=2}
#   allRes1 <- GenTable(GOdata, classicFisher = resultFisher,
#',' classicKS = resultKS, elimKS = resultKS.elim,
#',' orderBy = "classicKS", ranksOf = "classicFisher",
#',' topNodes = topn,numChar=100000)
#   topn=0
#   if(sum(score(resultKS.elim) < .05)>1){topn = sum(score(resultKS.elim) < .05)}else{topn=2}
#   allRes2 <- GenTable(GOdata, classicFisher = resultFisher,
#','  classicKS = resultKS, elimKS = resultKS.elim,
#','  orderBy = "elimKS", ranksOf = "classicFisher",
#','  topNodes = topn,numChar=100000)
#   topn=0
#   if(sum(score(resultFisher) < .05)>1){topn = sum(score(resultFisher) < .05)}else{topn=2}
#   allRes3 <- GenTable(GOdata, classicFisher = resultFisher,
#','  classicKS = resultKS, elimKS = resultKS.elim,
#','  orderBy = "classicFisher", ranksOf = "classicFisher",
#','  topNodes = topn,numChar=100000)
#   
#   write.csv(allRes1,file=paste("./GO_hitsALL_p05_KS",i,".csv",sep="")) ##-- Output x2 significant GO ?
#   write.csv(allRes2,file=paste("./GO_hitsALL_p05_elimKS",i,".csv",sep="")) ##-- Output x2 significant GO ?
#   write.csv(allRes3,file=paste("./GO_hitsALL_p05_Fisher",i,".csv",sep="")) ##-- Output x2 significant GO ?
#   
#   #--- Graph
#   #the subgraph induced by the 5 most signifcant GO
#   #terms as identifed by the elim algorithm. Signifcant nodes are represented as rectangles. The plotted graph
#   #is the upper induced graph generated by these signifcant nodes.
#   printGraph(GOdata, resultKS, firstSigNodes = 10,
#fn.prefix = paste("./tGO_hitsALL_p05_",i,sep=''),
#useInfo = "all", pdfSW = TRUE)
# }

###--- Consider the KS test that partially accounts for multiple testing through parental relationship
allRes1BP = read.csv(file="./GO_hitsALL_p05_KSBP.csv")
allRes1BP$GO.group = 'BP'
allRes1MF = read.csv(file="./GO_hitsALL_p05_KSMF.csv")
allRes1MF$GO.group = 'MF'
allRes1CC = read.csv(file="./GO_hitsALL_p05_KSCC.csv")
allRes1CC$GO.group = 'CC'
allResALL = rbind(allRes1BP,allRes1MF,allRes1CC)
allResALL = allResALL[,-1]
table(allResALL$GO.group)
# BP CC MF 
# 63 19 37

###-- GO terms for which at least 1 gene is annotated
sig = allResALL[allResALL$classicKS < 0.01 ,]
sig = sig[order(sig$classicKS),]

###-- GO terms for which at least 1 gene is annotated
dim(allResALL[allResALL$classicKS < 0.01,])
#[1] 26 10
write.csv2(allResALL[allResALL$classicKS<0.01,],file=paste0(draft,'Table2.csv'),row.names = F,quote=F)
dim(allResALL[allResALL$classicKS<0.01 & allResALL$Significant>0,])
#[1] 14 10

###-- GO terms with p-value < 0.01 + at least 1 gene 
allResALL[allResALL$classicKS<0.01 & allResALL$Significant>0,]
#    GO.ID Term Annotated Significant Expected
# 2   GO:0043412       macromolecule modification      1062    27    28.75
# 4   GO:0042493   response to drug  91     6     2.46
# 7   GO:0043051 regulation of pharyngeal pumping  37     2     1.00
# 8   GO:0006950 response to stress 509    24    13.78
# 9   GO:0010628 positive regulation of gene expression 125     3     3.38
# 11  GO:0046662regulation of oviposition  42     1     1.14
# 12  GO:0001505  regulation of neurotransmitter levels  65     4     1.76
# 13  GO:0006417regulation of translation  55     1     1.49
# 64  GO:0004674       protein serine/threonine kinase activity 240     6     6.99
# 66  GO:0016746 transferase activity, transferring acyl groups 154     4     4.49
# 68  GO:0016627 oxidoreductase activity, acting on the CH-CH group of donors  45     3     1.31
# 70  GO:0016634 oxidoreductase activity, acting on the CH-CH group of donors, oxygen as acceptor  10     1     0.29
# 73  GO:0016757   transferase activity, transferring glycosyl groups 229     6     6.67
# 101 GO:0044428       nuclear part 366     7    10.53


GOdata <- new("topGOdata", ontology = 'BP', allGenes = geneList, geneSel = genescore,
              nodeSize = 10,
              annot = annFUN.gene2GO, gene2GO = geneID2GO)

mm = unlist(genesInTerm(GOdata,"GO:0043412")) # macromolecule modification
myInterestingGenes$V1[myInterestingGenes$V1 %in% mm]
# HCOI00057400	C38C3.4
# HCOI00104200	vhp-1
# HCOI00266800	ufd-2
# HCOI00284300	mes-4
# HCOI00349300	
# HCOI00364200	pmk-2
# HCOI00547800	fzy-1
# HCOI00654300	Y41D4A.5
# HCOI00670400	rib-1
# HCOI00681500	cyl-1
# HCOI00757200	
# HCOI00867600	mlk-1
# HCOI00909200	
# HCOI00922400	prp-19
# HCOI00954100	
# HCOI01113600	
# HCOI01238500	dcaf-1
# HCOI01268000	oga-1
# HCOI01305200	oxi-1
# HCOI01361800	
# HCOI01475600	pigv-1
# HCOI01578200	
# HCOI01710100	rskn-1
# HCOI01710500	csb-1
# HCOI01710500	F53H4.6
# HCOI01957100	ZK524.4
# HCOI02039300	
# HCOI02056000

rd = unlist(genesInTerm(GOdata,"GO:0042493")) ##-  Response to drug 
myInterestingGenes$V1[myInterestingGenes$V1 %in% rd]
#[1] HCOI00391000 HCOI00398300 HCOI00447500 HCOI00489500 HCOI00711600 HCOI00809700
##- gpb-2 prdx-2 T19C4.5 aex-3 fat-3 gar-3

pp = unlist(genesInTerm(GOdata,"GO:0043051")) ##-  pharyngeal pumping
myInterestingGenes$V1[myInterestingGenes$V1 %in% pp]
#[1] HCOI00391000 HCOI00809700
#- gpb-2 gar-3

ro = unlist(genesInTerm(GOdata,"GO:0046662")) ##- regulation of oviposition
myInterestingGenes$V1[myInterestingGenes$V1 %in% ro]
#[1] HCOI00391000
#- gpb-2

sr = unlist(genesInTerm(GOdata,"GO:0006950")) ##-  stress response 
myInterestingGenes$V1[myInterestingGenes$V1 %in% sr]
# HCOI00062600	che-11
# HCOI00104200	vhp-1
# HCOI00127400	ced-1
# HCOI00234900	T24H7.2
# HCOI00289800	exo-3
# HCOI00315700	hid-1
# HCOI00329200	erl-1
# HCOI00364200	pmk-2
# HCOI00398300	prdx-2
# HCOI00417200	dnj-19
# HCOI00525900	
# HCOI00526100	sol-2
# HCOI00661100	tax-4
# HCOI00756200	smrc-1
# HCOI00778000	xpb-1
# HCOI01159100	fan-1
# HCOI01361300	
# HCOI01385300	T06D8.10
# HCOI01426500	rev-1
# HCOI01721500	Y116A8C.13
# HCOI01908800	msh-5
# HCOI01924400	rad-23
# HCOI02091600	dnj-29
# HCOI02140600	gcn-1

##-- Location of genes under enrichment 

##-- GPCR
allResALL[c(57,75),]
# GO.ID','      Term Annotated Significant Expected
# 57 GO:0007186 G-protein coupled receptor signaling pathway','254','    7     6.88
# 75 GO:0004930','   G-protein coupled receptor activity','197','    5     5.74
# Rank.in.classicFisher classicFisher classicKS  elimKS GO.group
# 57     287','  0.829    0.0434 0.01070','BP
# 75     155','  0.845    0.0111 0.04385','MF

gpcr = unlist(genesInTerm(GOdata,"GO:0007186")) ##-  
myInterestingGenes$V1[myInterestingGenes$V1 %in% gpcr]
#[1] HCOI00017000 HCOI00391000 HCOI00403600 HCOI00809700 HCOI01528000 HCOI01571600 HCOI01579400
##- dop-5, gpb-2, H23L24.4, gar-3, npr-3, npr-34, gpa-11
rm(myInterestingGenes,geneID2GO,allRes)

#####------ Check every XP-CLR pair-wise comparison ---####
# i = 0
# ##-- AFRICA
# for(p1 in c("KOK","MAR","NAM","WHR","ZAI")){
#   i = i + 1
#   j = 0
#   rm(temp,top,c,up,dw,t,loci,loc,rmstALL,rmst,score,chromosome,selection_score)
#   for(p2 in c("KOK","MAR","NAM","WHR","ZAI")){
#     j = j + 1
#     if(j>i){
      # temp = xpclr[xpclr$V8==paste0(p1,'.',p2) | xpclr$V8==paste0(p2,'.',p1),]
      # 
      # loci = length(unique(temp$pos))
      # loc = unique(temp$pos)
      # rms = array(0,loci)
      # 
      # df_list = split(temp[,c('V1','V6')],temp$pos)
      # rmstALL = lapply(df_list,RMS)
      # rmst = lapply(rmstALL,unique)
      # 
      # loc = as.integer(as.character(names(rmst)))
      # chromosome = sapply(rmst,function(x) x[1])
      # selection_score = sapply(rmst,function(x) as.numeric(as.character(x[2])))
      # score = cbind(loc, chromosome,selection_score)
      # score = data.frame(score)
      # colnames(score) = c('loc','chromosome','selection_score')
      
      #write.table(score,file=paste0('xpclr_selection_scores_',p1,'.',p2,'.txt'),row.names=F,quote=F)
      
      # pdf(file = paste0('XPCLR.',p1,'_',p2,'.pdf'),width = 14,height = 8)
      # print(ggplot(score,aes(x = loc/1e6,y = selection_score,col = factor(chromosome))) +
      #         geom_vline(xintercept = tbb1[1]/1e6, lty = 2,lwd = .1) +
      #         geom_point(size = .5) +
      #         geom_hline(yintercept = quantile(score$selection_score,0.999),col='red') +
      #         scale_colour_manual(values = chrom.colors) +
      #         scale_y_continuous(limits = c(0,max(score$selection_score))) + xlab('Genomic position (Mbp)') +
      #         scale_x_continuous(limits = c(0,max(score$loc/1e6)+1),breaks = seq(0,max(score$loc/1e6),20)) +
      #         theme(legend.position='none',text=element_text(size=16)) +
      #         ggtitle(paste0('XP-CLR - ',p1,' vs. ',p2)))
      # dev.off()
#       
#       top = score[score$selection_score > quantile(score$selection_score,0.999),]
#       c = top$chromosome
#       up = top$loc - 1000 - chrom_end[top$chromosome]
#       dw = top$loc + 1000 - chrom_end[top$chromosome]
#       sco = top$selection_score
#       t = cbind(c,up,dw,sco)
#       t = data.frame(t)
#       write.table(t,file = paste0('hits_',p1,'.',p2,'_01p.tsv'),row.names=F,quote=F,sep='\t')
#     }
#   }
# }

# ##-- AUSTRALIA
# p1 = 'AUS.1'
# p2 = 'AUS.2'
# temp = xpclr[xpclr$V8 == paste0(p1,'.',p2) | xpclr$V8 == paste0(p2,'.',p1),]
# 
# loci = length(unique(temp$pos))
# loc = unique(temp$pos)
# rms = array(0,loci)
# 
# df_list = split(temp[,c('V1','V6')],temp$pos)
# rmstALL = lapply(df_list,RMS)
# rmst = lapply(rmstALL,unique)
# 
# loc = as.integer(as.character(names(rmst)))
# chromosome = sapply(rmst,function(x) x[1])
# selection_score = sapply(rmst,function(x) as.numeric(as.character(x[2])))
# score = cbind(loc, chromosome,selection_score)
# score = data.frame(score)
# colnames(score) = c('loc','chromosome','selection_score')
# 
# write.table(score,file = paste0('xpclr_selection_scores_',p1,'.',p2,'.txt'),
#       row.names=F,quote=F)
# pdf(file=paste0('XPCLR.',p1,'_',p2,'.pdf'),width=14,height=8)
# ggplot(score,aes(x=loc/1e6,y=selection_score,col=factor(chromosome))) +
#   geom_vline(xintercept=tbb1[1]/1e6, lty=2,lwd= .1) +
#   geom_point(size=.5) +
#   geom_hline(yintercept = quantile(score$selection_score,0.999),col='red') +
#   scale_colour_manual(values=chrom.colors) +
#   scale_y_continuous(limits=c(0,max(score$selection_score))) + xlab('Genomic position (Mbp)') +
#   scale_x_continuous(limits=c(0,max(score$loc/1e6)+1),breaks = seq(0,max(score$loc/1e6),20)) +
#   theme(legend.position='none',text=element_text(size=16)) +
#   ggtitle(paste0('XP-CLR - ',p1,' vs. ',p2))
# dev.off()
# 
###-- Prepare hit files for blast and inference of underlying genes
##- Blast and gene inference is done using the convert_hits_to_gene.sh script 
# top = score[score$selection_score>quantile(score$selection_score,0.999),]
# c = top$chromosome
# up = top$loc-1000-chrom_end[top$chromosome]
# dw = top$loc+1000-chrom_end[top$chromosome]
# sco = top$selection_score
# t = cbind(c,up,dw,sco)
# t = data.frame(t)
# write.table(t,file=paste0('hits_',p1,'.',p2,'_01p.tsv'),row.names=F,quote=F,sep='\t')
# 
# rm(loci,loc,rmstALL,rmst,score,chromosome,selection_score)

###-- GO enrichment for every pair-wise comp
complist=c('AUS.1.AUS.2','KOK.MAR','KOK.NAM','KOK.WHR','KOK.ZAI','MAR.NAM','MAR.WHR','MAR.ZAI','NAM.WHR','NAM.ZAI')
# allResF = matrix(0,1,9)
# allResF = data.frame(allResF)
# colnames(allResF) = c('GO.ID','Term','Annotated','Significant','Expected','classicFisher','classicKS','COMP','GO.group')
allResK = matrix(0,1,10)
allResK = data.frame(allResK)
colnames(allResK) = c('GO.ID','Term','Annotated','Significant','Expected','Rank in classicFisher','classicFisher','classicKS','COMP','GO.group')

for(comp in complist){
  
  myInterestingGenes = read.table(file=paste0('./comp/blast',comp,'/gene_list.txt'),header=F,sep='\t')
  colnames(myInterestingGenes)[1] = 'V1'
  geneNames = read.table(file='./blastALL01p/HcoGeneList.txt',header=T,sep='\t')
  geneList <- factor(as.integer(geneNames$Gene.stable.ID %in% myInterestingGenes$V1))
  names(geneList) <- geneNames$Gene.stable.ID
  str(geneList)
  geneID2GO = readMappings(file=paste0(draft,'Hco_GO_topGO.txt')) #-- Total GO
  
  ##- Attach gene score
  hb = read.table(file=paste0('./comp/blast',comp,'/tot_geneid.txt'),sep=' ',header=F,fill=T)
  colnames(hb)[1] = 'loc'
  t = read.table(file = paste0('hits_',comp,'_01p.tsv'),header=T)
  t$loc[t$c!=5] = paste0('top_',t$c[t$c!=5],'_Celeg_TT:',t$up[t$c!=5],'-',t$dw[t$c!=5],'.fasta.blastn')
  t$loc[t$c==5] = paste0('top_',t$c[t$c==5],'_Celeg_TT_axel:',t$up[t$c==5],'-',t$dw[t$c==5],'.fasta.blastn')
  t$loc = factor(t$loc)
  t2 = merge(t,hb,by='loc',all=T)
  colnames(t2)[5]='selection_score'
  
  a = aggregate(selection_score ~ V3, FUN=mean,data=t2)
  a = a[-1,]
  genescore = a$selection_score[match(a$V3,myInterestingGenes$V1)]
  names(genescore) <- a$V3[match(a$V3,myInterestingGenes$V1)]
  
  #--- Molecular / Cellular Componenent / Biological Processes
  goterm=c('MF','CC','BP')
  topn=0
  for(i in goterm){
    GOdata <- new("topGOdata", ontology = i, allGenes = geneList, geneSel = genescore,
                  nodeSize = 10,
                  annot = annFUN.gene2GO, gene2GO = geneID2GO)
    ##- Fisher
    resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
    
    ##- KS test <> based on gene score + weight01 = account for hierarchy <> stat correction
    resultKS <- runTest(GOdata, algorithm = "weight01", statistic = "ks")
    
    #- Table of all results - Classic Fisher
    if(sum(score(resultKS) < .01)>1){topn = sum(score(resultKS) < .01)}else{topn=2}
    allRes1 <- GenTable(GOdata, classicFisher = resultFisher,
                        classicKS = resultKS,
                        orderBy = "classicKS", ranksOf = "classicFisher",
                        topNodes = topn, numChar=100000)
    
    # if(sum(score(resultFisher) < .01)>1){topn = sum(score(resultFisher) < .01)}else{topn=2}
    # allRes2 <- GenTable(GOdata, classicFisher = resultFisher,
    #classicKS = resultKS,
    #orderBy = "classicFisher", ranksOf = "classicFisher",
    #topNodes = topn,numChar=100000)
    allRes1$COMP = comp
    #allRes2$COMP = comp
    allRes1$GO.group = i
    #allRes2$GO.group = i
    
    #allResF = rbind(allResF,allRes2)
    allResK = rbind(allResK,allRes1)
    
  }
  rm(myInterestingGenes,geneID2GO,allRes1,allRes2,a,t2,hb,t)
}
allResK = allResK[-1,]
write.csv2(allResK,file = "./GO_hits_BYCOMP_p05_KS.csv",quote=F,row.names=F) ##-- Output x2 significant GO ?
write.csv2(allResK[,c(1,2,3,4,5,8,9,10),],file=paste0(draft,'supplementary_table3.csv'),quote=F,row.names=F)

###--- Check congruence between across-comparisons with pair-wise comparisons
allResK = read.csv(file = "./GO_hits_BYCOMP_p05_KS.csv",header=T,row.names=NULL)

common = allResK[(allResK$GO.ID %in% allResALL$GO.ID[allResALL$classicKS<=.01]) & allResK$Significant>0,]
table(common$COMP) ##-- present in all pair-wise comparisons
# AUS.1.AUS.2     KOK.MAR     KOK.NAM     KOK.WHR     KOK.ZAI     MAR.NAM     MAR.WHR     MAR.ZAI     NAM.WHR     NAM.ZAI 
# 3           6           6           4           6           5           5           5           1           4 

##-- Most frequent GO terms
freqgo = data.frame(table(allResK$Term[allResK$Significant>0]))
freqgo[order(freqgo$Freq,decreasing = T),]
#    Var1 Freq
# 5 macromolecule modification    9
# 14 response to stress    9
# 11 protein serine/threonine kinase activity    7
# 6     microtubule cytoskeleton organization    6
# 13   response to drug    3
# 3     AP-type membrane coat adaptor complex    2
# 4  cysteine-type peptidase activity    2
# 8  nuclear part    2
# 10   positive regulation of gene expression    2
# 1   acyl-CoA dehydrogenase activity    1
# 2  anatomical structure homeostasis    1
# 7     neurotransmitter transporter activity    1
# 9  oxidoreductase activity, acting on the CH-CH group of donors    1
# 12    regulation of neurotransmitter levels    1
# 15 transferase activity, transferring acyl groups    1

###-- Anthelmintic-related GO terms found in STA.1-specific comparison (Ivermectin resistant)
goah = c('GO:0042493','GO:0005326','GO:0001505')
common[common$GO.ID %in% goah,]
#         GO.ID                 Term Annotated Significant Expected
# 30  GO:0005326 neurotransmitter transporter activity  42     1     0.19
# 44  GO:0042493  response to drug  91     1     0.32
# 72  GO:0042493  response to drug  91     2     0.34
# 124 GO:0042493  response to drug  91     1     0.26
# 135 GO:0001505 regulation of neurotransmitter levels  65     1     0.18
#     Rank in classicFisher classicFisher classicKS    COMP GO.group
# 30 48       1   0.00065 KOK.MAR MF
# 44 47    0.28   0.00030 KOK.MAR BP
# 72 20   0.081   0.00029 KOK.NAM BP
# 124      79    1.00   5.1e-05 KOK.ZAI BP
# 135      89    1.00   0.00922 KOK.ZAI BP

###-- Retrieved genes associated with anthelmintic-related activity
# complist = c('KOK.MAR','KOK.MAR','KOK.NAM','KOK.NAM','KOK.ZAI')
# ontolist = c('MF','BP','BP','BP','BP')
# golist = c("GO:0005326", "GO:0042493", "GO:0042493", "GO:0042493", "GO:0001505")
# druggene=array(0,6)
# i=0
# 
# for(comp in complist[3]){
#   i = 3 #i+1
#   jj=0
#   myInterestingGenes = read.table(file=paste0('./blast',comp,'/gene_list.txt'),header=F,sep='\t')
#   colnames(myInterestingGenes)[1] = 'V1'
#   geneNames = read.table(file='./blastALL01p/HcoGeneList.txt',header=T,sep='\t')
#   geneList <- factor(as.integer(geneNames$Gene.stable.ID %in% myInterestingGenes$V1))
#   names(geneList) <- geneNames$Gene.stable.ID
#   str(geneList)
#   geneID2GO = readMappings(file=paste0(draft,'Hco_GO_topGO.txt')) #-- Total GO
#   
#   ##- Attach gene score
#   hb = read.table(file=paste0('./blast',comp,'/tot_geneid.txt'),sep=' ',header=F,fill=T)
#   colnames(hb)[1] = 'loc'
#   t = read.table(file = paste0('hits_',comp,'_01p.tsv'),header=T)
#   t$loc[t$c!=5] = paste0('top_',t$c[t$c!=5],'_Celeg_TT:',t$up[t$c!=5],'-',t$dw[t$c!=5],'.fasta.blastn')
#   t$loc[t$c==5] = paste0('top_',t$c[t$c==5],'_Celeg_TT_axel:',t$up[t$c==5],'-',t$dw[t$c==5],'.fasta.blastn')
#   t$loc = factor(t$loc)
#   t2 = merge(t,hb,by='loc',all=T)
#   colnames(t2)[5]='selection_score'
#   
#   a = aggregate(selection_score ~ V3, FUN=mean,data=t2)
#   a = a[-1,]
#   genescore = a$selection_score[match(a$V3,myInterestingGenes$V1)]
#   names(genescore) <- a$V3[match(a$V3,myInterestingGenes$V1)]
#   
#   GOdata <- new("topGOdata", ontology = ontolist[i], allGenes = geneList, geneSel = genescore,
#   nodeSize = 10,
#   annot = annFUN.gene2GO, gene2GO = geneID2GO)
#   ##- KS test <> based on gene score + weight01 = account for hierarchy <> stat correction
#   resultKS <- runTest(GOdata, algorithm = "weight01", statistic = "ks")
#   
#   sr = unlist(genesInTerm(GOdata,golist[i])) ##-  KOK-MAR response to drug
#   tmp = myInterestingGenes$V1[myInterestingGenes$V1 %in% sr]
#   add=F
#   if(length(tmp)==1|is(add)){
#     druggene[i]=as.character(tmp)}else{
# for(jj in 1:length(tmp)){druggene[i+jj-1]=as.character(tmp[jj]);add=T}}
# }

## Gene list with anthelmintic-related function found under selection in comparison involving STA.1
druggene
#[1] "HCOI00389600" "HCOI00032800" "HCOI00032800" "HCOI00243900" "HCOI00032800" "HCOI00489500"
# CRB-SNF-9 (Na symporter) / unc-24 (stomatin-like) / / B0361.4 (aff by sirolimus) / / aex-3 (synaptic vesicle release)

####----- Most likely IVM candidates: aex3, glc4, pgp11 in AUS.1 and AUS.2 
## aex-3                                                                                                                                                                                                                                                                                                                   
#top_1_Celeg_TT:9596364-9598364.fasta.blastn gene_id "HCOI00489500"                                                                                                                                                                                                                                                        

## glc-4                                                                                                                                                                                                                                                                                                                   
#top_1_Celeg_TT:40195769-40197769.fasta.blastn gene_id "HCOI00617300"                                                                                                                                                                                                                                                      
#top_1_Celeg_TT:40199769-40201769.fasta.blastn gene_id "HCOI00617300"                                                                                                                                                                                                                                                      

## pgp-11                                                                                                                                                                                                                                                                                                                  
#top_5_Celeg_TT_axel:47201208-47203208.fasta.blastn gene_id "HCOI00233200"                                                                                                                                                                                                                                                 
#top_5_Celeg_TT_axel:47203208-47205208.fasta.blastn gene_id "HCOI00233200"                                                                                                                                                                                                                                                 

## Nucleotide diversity comparison between Australian isolates
#- limited structure and known status
#- coverage is equivalent 
t.test(mcov ~ iso,data = df[df$iso %in% c('AUS.1','AUS.2'),])
# data:  mcov by iso
# t = -0.29533, df = 13.072, p-value = 0.7724
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -1.407336  1.068669
# sample estimates:
#   mean in group AUS.1 mean in group AUS.2 
# 2.566667            2.736000 

rep = "./data/NUC/TAJ/"
isivm = c('AUS.1','AUS.2') 
t.ALL = NULL
for(r in isivm){
  for(c in c('1_Celeg_TT','5_Celeg_TT_axel')){
    temp = read.table(paste0(rep,c,'/',r,'.thetas.pestPG'))
    temp = temp[,c(2,3,5,9,14)]
    colnames(temp) = c('Chr','WinCenter','tP','Tajima','nSites')
    temp$pop = r
    t.ALL = rbind(t.ALL,temp)
    rm(temp)
  }
}
colnames(t.ALL)=c('Chr','WinCenter','tP','Tajima','nSites','pop')
t.ALL = na.omit(t.ALL)

## Standardize pi values for comparison across isolates
mp = aggregate(tP/nSites ~ pop,FUN = mean,data=t.ALL)
mp$std = aggregate(tP/nSites ~ pop,FUN = sd,data=t.ALL)[,2]
t.ALL$stp = 0
for(p in t.ALL$pop){
  m = mp$`tP/nSites`[mp$pop==p]
  s = mp$std[mp$pop==p]
  t.ALL$stp[t.ALL$pop==p] = ((t.ALL$tP[t.ALL$pop==p]/t.ALL$nSites[t.ALL$pop==p]) - m)/s
}

## Add needed information
t.ALL$res = 'IVM-S'
t.ALL$res[t.ALL$pop %in% c('AUS.1','STA.1')] = 'IVM-R'

t.ALL$geo = 'Africa'
t.ALL$geo[grep ('AUS',t.ALL$pop)] = 'Australia'
t.ALL$res = factor(t.ALL$res)
t.ALL$geo = factor(t.ALL$geo)

## Retain coordinates of interest 
piglc4 = t.ALL[t.ALL$Chr=='1_Celeg_TT'  & t.ALL$WinCenter>40196769-50e3 & t.ALL$WinCenter<40196769+50e3,] 
pipgp11 = t.ALL[t.ALL$Chr=='5_Celeg_TT_axel' & t.ALL$WinCenter>47204208-50e3 & t.ALL$WinCenter<47204208+50e3,] 

colsivm = c('#e31a1c','#41b6c4',
            '#8073ac','#5aae61','#bf812d','#fd8d3c')

p1 = ggplot(piglc4,aes(x=WinCenter/1000,y=tP/nSites,col=pop,shape = res)) +
  scale_color_manual(values = colsivm)+ 
  scale_shape_manual(values = c(17,15)) +
  geom_point(size=3,alpha=.4) + geom_line() +
  scale_x_continuous() +
  #facet_wrap(~ geo) +
  geom_vline(xintercept = 40195.769,lty=3,lwd=.5) +
  geom_vline(xintercept = 40200.769,lty=3,lwd=.5) +
  xlab('Position (Kbp)') + ylab('Nucleotide diversity') +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  #theme(text=element_text(size=16),legend.position ='none',legend.title=element_blank())+
  #theme(legend.position = 'none') +
  theme(legend.title=element_blank())+
  ggtitle('HCOI00617300 (glc-4 ortholog)')

p2 = ggplot(pipgp11,aes(x=WinCenter/1000,y=tP/nSites,col=pop,shape = res)) +
  scale_color_manual(values = colsivm)+ 
  scale_shape_manual(values = c(17,15)) +
  geom_point(size=3,alpha=.4) + geom_line() +
  scale_x_continuous() +
  #facet_wrap(~ geo) +
  geom_vline(xintercept = 47202.208,lty=3,lwd=.5) +
  geom_vline(xintercept = 47204.208,lty=3,lwd=.5) +
  xlab('Position (Kbp)') + ylab('Nucleotide diversity') +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  #theme(text=element_text(size=16),legend.position ='none',legend.title=element_blank())+
  #theme(legend.position = 'none') +
  theme(legend.title=element_blank())+
  ggtitle('HCOI00233200 (pgp-11 ortholog)')

## Test between Australian isolates
summary(lm(tP/nSites ~ geo*res ,data = piglc4))
# Call:
#   lm(formula = tP/nSites ~ geo * res, data = piglc4)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0042288 -0.0015783 -0.0004638  0.0015682  0.0052032 
# 
# Coefficients:
#                       Estimate Std. Error     t value     Pr(>|t|)    
#   (Intercept)            0.0039634  0.0007966   4.975 0.0000065407 ***
#   geoAustralia           0.0043173  0.0011266   3.832     0.000324 ***
#   resIVM-S               0.0048284  0.0009198   5.249 0.0000024369 ***
#   geoAustralia:resIVM-S -0.0097589  0.0014544  -6.710 0.0000000104 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.002519 on 56 degrees of freedom
# Multiple R-squared:  0.4876,	Adjusted R-squared:  0.4601 
# F-statistic: 17.76 on 3 and 56 DF,  p-value: 0.00000003183

summary(lm(tP/nSites ~ geo*res ,data = pipgp11))
# Call:
#   lm(formula = tP/nSites ~ geo * res, data = pipgp11)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0038766 -0.0013659  0.0001058  0.0010976  0.0047591 
# 
# Coefficients:
#                           Estimate   Std. Error t value Pr(>|t|)    
# (Intercept)            0.010488705  0.000674934  15.540  < 2e-16 ***
# geoAustralia          -0.003405393  0.000954501  -3.568 0.000747 ***
# resIVM-S               0.000003261  0.000779346   0.004 0.996676    
# geoAustralia:resIVM-S -0.002484570  0.001232255  -2.016 0.048577 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.002134 on 56 degrees of freedom
# Multiple R-squared:  0.5556,	Adjusted R-squared:  0.5318 
# F-statistic: 23.33 on 3 and 56 DF,  p-value: 6.287e-10

##-- Most likely candidates identified by XP-CLR and STA.1
#"HCOI00389600" snf-9
#blastKOK.MAR/tot_geneid.txt:top_3_Celeg_TT:16033820-16035820.fasta.blastn gene_id "HCOI00389600"
#"HCOI00032800" unc-24
#blastKOK.MAR/tot_geneid.txt:top_3_Celeg_TT:3329155-3331155.fasta.blastn gene_id "HCOI00032800"
#blastKOK.MAR/tot_geneid.txt:top_3_Celeg_TT:3331155-3333155.fasta.blastn gene_id "HCOI00032800"
#blastKOK.NAM/tot_geneid.txt:top_3_Celeg_TT:3329818-3331818.fasta.blastn gene_id "HCOI00032800"
#blastKOK.NAM/tot_geneid.txt:top_3_Celeg_TT:3331818-3333818.fasta.blastn gene_id "HCOI00032800"
#"HCOI00243900" B0361.4
#blastKOK.NAM/tot_geneid.txt:top_3_Celeg_TT:40297965-40299965.fasta.blastn gene_id "HCOI00243900"
#"HCOI00489500" aex-3
#blastKOK.ZAI/tot_geneid.txt:top_1_Celeg_TT:9596364-9598364.fasta.blastn gene_id "HCOI00489500"

rep = "./data/NUC/TAJ/"
isivm = c('STA.1','MOR','NAM','ZAI') 
t.ALL = NULL
for(r in isivm){
  for(c in c('1_Celeg_TT','3_Celeg_TT')){
    temp = read.table(paste0(rep,c,'/',r,'.thetas.pestPG'))
    temp = temp[,c(2,3,5,9,14)]
    colnames(temp) = c('Chr','WinCenter','tP','Tajima','nSites')
    temp$pop = r
    t.ALL = rbind(t.ALL,temp)
    rm(temp)
  }
}
colnames(t.ALL)=c('Chr','WinCenter','tP','Tajima','nSites','pop')
t.ALL = na.omit(t.ALL)

t.ALL$res = 'IVM-S'
t.ALL$res[t.ALL$pop %in% c('AUS.1','STA.1')] = 'IVM-R'

## Retain coordinates of interest 
piaex32 = t.ALL[t.ALL$Chr=='1_Celeg_TT' & t.ALL$WinCenter>9597364-50e3 & t.ALL$WinCenter<9597364+50e3,] 
pisnf9 = t.ALL[t.ALL$Chr=='3_Celeg_TT'  & t.ALL$WinCenter>16034820-50e3 & t.ALL$WinCenter<16034820+50e3,] 
piunc24 = t.ALL[t.ALL$Chr=='3_Celeg_TT'  & t.ALL$WinCenter>3329155-50e3 & t.ALL$WinCenter<3331818+50e3,] 
pib = t.ALL[t.ALL$Chr=='3_Celeg_TT'  & t.ALL$WinCenter>40287965-50e3 & t.ALL$WinCenter<40287965+50e3,] 

###-- Reduced diversity for aex-3
p3=ggplot(piaex32,aes(x=WinCenter/1000,y=tP/nSites,col=pop,shape=res)) +
  scale_color_manual(values = colsivm[3:6])+ 
  scale_shape_manual(values = c(17,15)) +
  geom_point(size=3,alpha=.4) + geom_line() +
  scale_x_continuous() +
  #facet_wrap(~ res) +
  geom_vline(xintercept =9596.364,lty=3,lwd=.5) +
  geom_vline(xintercept = 9598.364,lty=3,lwd=.5) +
  xlab('Position (Kbp)') + ylab('Nucleotide diversity') +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(legend.title=element_blank())+
  ggtitle('HCOI00489500 (aex3 ortholog)')

###-- Reduced diversity for snf-9
p4=ggplot(pisnf9,aes(x=WinCenter/1000,y=tP/nSites,col=pop,shape=res)) +
  scale_color_manual(values = colsivm[3:6])+ 
  scale_shape_manual(values = c(17,15)) +
  geom_point(size=3,alpha=.4) + geom_line() +
  scale_x_continuous() +
  #facet_wrap(~ res) +
  geom_vline(xintercept = 16033.820,lty=3,lwd=.5) +
  geom_vline(xintercept = 16035.820,lty=3,lwd=.5) +
  xlab('Position (Kbp)') + ylab('Nucleotide diversity') +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(legend.title=element_blank())+
  ggtitle('HCOI00389600 (snf-9 ortholog)')

###-- Reduced diversity for unc-24
p5=ggplot(piunc24,aes(x=WinCenter/1000,y=tP/nSites,col=pop,shape=res)) +
  scale_color_manual(values = colsivm[3:6])+ 
  scale_shape_manual(values = c(17,15)) +
  geom_point(size=3,alpha=.4) + geom_line() +
  scale_x_continuous() +
  #facet_wrap(~ res) +
  geom_vline(xintercept = 3329.155,lty=3,lwd=.5) +
  geom_vline(xintercept = 3333.818,lty=3,lwd=.5) +
  xlab('Position (Kbp)') + ylab('Nucleotide diversity') +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(legend.title=element_blank())+
  ggtitle('HCOI00032800 (unc-24 ortholog)')

###-- Reduced diversity for B0361.4
p6=ggplot(pib,aes(x=WinCenter/1000,y=tP/nSites,col=pop,shape=res)) +
  scale_color_manual(values = colsivm[3:6])+ 
  scale_shape_manual(values = c(17,15)) +
  geom_point(size=3,alpha=.4) + geom_line() +
  scale_x_continuous() +
  #facet_wrap(~ res) +
  geom_vline(xintercept = 40297.965,lty=3,lwd=.5) +
  geom_vline(xintercept = 40299.965,lty=3,lwd=.5) +
  xlab('Position (Kbp)') + ylab('Nucleotide diversity') +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(legend.title=element_blank())+
  ggtitle('HCOI00243900 (B0361.4 ortholog)')

pdf(file=paste0(draft2,'Figure3d.pdf'),width=14,height=14)
multiplot(p1,p3,p4,p2,p5,p6,cols=2)
dev.off()

#####====== IVERMECTIN resistance - ANGSD FST ====########

setwd("./data/NUC/FST/")
complist = c('AUS.1.AUS.2','NAM.STA.1','NAM.STA.3','STA.1.ZAI','STA.3.ZAI') #// only comparison for which IVM efficacy data are available
fsttot = matrix(0,1,8)
colnames(fsttot)= c('chr','midPos','Nsites','FST','bin','POS','COMP','cuto')
fsttot = data.frame(fsttot)
topIvmTot = fsttot
cu = array(0,length(complist))
n=0

for(comp in complist){
  n=n+1
  ch = '1_Celeg_TT'
  
  fstivm = read.table(file=paste0('./',ch,'/',comp,'.fst.txt'))
  colnames(fstivm) = c('chr','midPos','Nsites','FST')
  fstivm = fstivm[which(fstivm$Nsites>=500),]
  
  for(ch in chromlist[-1]){
    tmp = read.table(file=paste0('./',ch,'/',comp,'.fst.txt'))
    colnames(tmp) = c('chr','midPos','Nsites','FST')
    tmp = tmp[which(tmp$Nsites>=500),]
    fstivm = rbind(fstivm,tmp)
    rm(tmp)
  }
  fstivm$bin=seq(1:dim(fstivm)[1])
  fstivm = na.omit(fstivm)
  fstivm$FST[fstivm$FST<0]=0
  fstivm$POS=0
  chrom_end=array(0,6)
  for(c in seq(1:5)){chrom_end[c+1] = chrom_end[c] + max(fstivm$midPos[fstivm$chr==chromlist[c]])}
  for(c in seq(1:5)){
    fstivm$POS[fstivm$chr==chromlist[c]] = fstivm$midPos[fstivm$chr==chromlist[c]] + chrom_end[c]
  }
  fstivm = fstivm[(fstivm$POS > chrom_end[1]+0.1e6 & fstivm$POS < chrom_end[2]-0.1e6)| ## Remove noisy signals on chrom edges
                    (fstivm$POS > chrom_end[2]+0.1e6 & fstivm$POS < chrom_end[3]-0.1e6)| 
                    (fstivm$POS > chrom_end[3]+0.1e6 & fstivm$POS < chrom_end[4]-0.1e6)|
                    (fstivm$POS > chrom_end[4]+0.1e6 & fstivm$POS < chrom_end[5]-0.1e6)|
                    (fstivm$POS > chrom_end[5]+0.1e6 & fstivm$POS < chrom_end[6]-0.1e6),]
  fstivm$COMP = comp
  fstivm$cuto = quantile(fstivm$FST,0.995) 
  cu[n] = quantile(fstivm$FST,0.995)
  topIvm = fstivm[fstivm$FST>fstivm$cuto,]
  topIvm$COMP = comp
  fsttot = rbind(fsttot,fstivm)
  topIvmTot = rbind(topIvmTot,topIvm)
  rm(fstivm,topIvm)
}
fsttot = fsttot[-1,]

###--- Figure 3b
afraus=fsttot
afraus$COMP=factor(afraus$COMP)

pdf(file=paste0(draft2,'Figure3b.pdf'),width=14,height=14)
ggplot(afraus,aes(x=POS/1e6,y=FST,col=chr)) +
  xlab('POS (Mbp)') + 
  facet_wrap(~ COMP,ncol=1) +
  scale_color_manual(values=chrom.colors) + geom_point(size=.5) +
  geom_hline(data=afraus,aes(yintercept = cuto), lty=2,lwd=.3) +
  geom_vline(xintercept = tbb1[1]/1e6,lty=3,lwd=.1) +
  geom_vline(xintercept = tbb1[2]/1e6,lty=3,lwd=.1) +
  geom_label_repel(data=candid,size=2.5,aes(x=loc/1e6,y=0.8,col=factor(chromosome),label=gene)) +
  geom_segment(data=candid,aes(x=loc/1e6,y=0.1,xend=loc/1e6,yend=0.8,col=factor(chromosome)),size=0.2) +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.1)) +
  scale_x_continuous(limits=c(0,max(ausivm$POS/1e6)),breaks=seq(0,max(ausivm$POS/1e6),10)) +
  theme(legend.position ='none',text = element_text(size=16),strip.text.x = element_text(size = 18))+ ggtitle('')
dev.off()

##-- Chrom V region
chr5=fsttot[fsttot$chr=='5_Celeg_TT_axel' & fsttot$FST>fsttot$cuto,]
table(chr5$COMP[chr5$midPos>38e6 & chr5$midPos<42e6])
# NAM.STA.1 STA.1.ZAI STA.3.ZAI 
# 11        23        18
table(chr5$COMP[chr5$midPos>45e6 & chr5$midPos<48e6])

##-- Chrom I region
chr1=fsttot[fsttot$chr=='1_Celeg_TT' & fsttot$FST>fsttot$cuto,]
summary(chr1$midPos[chr1$midPos<10e6])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6065000 7072500 7680000 7601974 7977500 9055000 

table(chr1$COMP[chr1$midPos>7e6 & chr1$midPos<7.1e6])
# NAM.STA.1 NAM.STA.3 STA.3.ZAI 
# 2        10         3 

chr1[chr1$midPos>tbb1[1]-5000 & chr1$midPos<tbb1[1]+5000,]
# chr  midPos Nsites      FST bin
# (2020806,2027171)(7000433,7049423)(7000000,7050000) 1_Celeg_TT 7025000   6367 0.426568 691
# (2189660,2199773)(7000598,7049488)(7000000,7050000) 1_Celeg_TT 7025000  10115 0.650556 693
# POS      COMP      cuto
# (2020806,2027171)(7000433,7049423)(7000000,7050000) 7025000 NAM.STA.1 0.4105172
# (2189660,2199773)(7000598,7049488)(7000000,7050000) 7025000 NAM.STA.3 0.4485684


#####====== CLIMATE ADAPTATION  ====########

###--- Pair-Wise FST using ANGSD / Polycomb
setwd("./data/NUC/CLIM/")

clim=read.table(file='./data/koppen_1901-2010.tsv',header=T)
colnames(clim)=c("lon","lat","kop")
clim$kop=as.character(clim$kop)

##- Recode df to get lat and lon to have .25 and .75
dd=unique(df[,c("reg","iso","lat","lon")])
dd$kop=0
v=seq(1,dim(dd)[1])

for(i in v){
  temp_lat = dd$lat[i]
  temp_lon = dd$lon[i]
  temp = clim[which(clim$lon> temp_lon-2 & clim$lon< temp_lon+2 &
                      clim$lat> temp_lat-2 & clim$lat< temp_lat+2),]
  tt=array(0,dim(temp)[1])
  t=seq(1,dim(temp)[1])
  
  for(k in t){
    dlat=abs(temp_lat - temp$lat[k])
    dlon=abs(temp_lon - temp$lon[k])
    d=dlat+dlon
    tt[k] = d
  }
  dd$kop[i]=temp$kop[which.min(tt)]
}
write.table(unique(dd[,c('reg','iso','kop')]),
            file='climate.txt',quote=F,row.names=F)
#     reg   iso      lat      lon kop
# 1   AUS AUS.1 -29.2893 150.9920 Cfa
# 16  BEN   BEN   9.3077   2.3158  Aw
# 21  BRA   BRA -14.2350 -51.9253  Aw
# 23  CAP   CAP  16.5388 -23.0418 BWh
# 30  FRA FRA.1  46.3237  -0.4648 Cfb
# 42  FRA FRA.2  43.1866  -0.9449 Cfb
# 55  FRA FRA.3  43.9251   2.1486 Cfb
# 75  FRA FRA.4  45.6836   2.8415 Cfb
# 87  FRG   FRG  16.2650 -61.5510  Am
# 110 IND   IND  -0.7893 113.9213  Af
# 116 ITA AUS.2 -30.5016 151.6662 Cfb
# 126 MOR   MOR  31.7917  -7.0926 Csa
# 139 NAM   NAM -22.9576  18.4904 BWh
# 154 POR   ACO  37.7412 -25.6756 Cfb
# 159 STA STA.1 -30.5096  29.4063 Cfb
# 173 STA STA.2 -25.6409  28.1708 Cwb
# 187 STA STA.3 -25.3299  31.0163 Cwa
# 200 STO   STO   0.1864   6.6131  Am
# 209 ZAI   ZAI  -4.0383  21.7587  Aw

setwd("./data/NUC/FST/")
###--- TROP vs. ARID: FRG vs. NAM (A vs. B)
###--- TROP vs. TEMP: FRG vs. FRA (A vs. C)
###--- ARID vs. TEMP: NAM vs. FRA (B vs. C)

complist = c('FRA.1.FRG','FRA.1.NAM','FRG.NAM') #// only comparison for which IVM efficacy data are available
fsttot = matrix(0,1,8)
colnames(fsttot)= c('chr','midPos','Nsites','FST','bin','POS','COMP','cuto')
fsttot = data.frame(fsttot)
topIvmTot = fsttot
cu = array(0,length(complist))
n=0

for(comp in complist){
  n=n+1
  ch = '1_Celeg_TT'
  
  fstclim = read.table(file=paste0('./',ch,'/',comp,'.fst.txt'))
  colnames(fstclim) = c('chr','midPos','Nsites','FST')
  fstclim = fstclim[which(fstclim$Nsites>=500),]
  
  for(ch in chromlist[-1]){
    tmp = read.table(file=paste0('./',ch,'/',comp,'.fst.txt'))
    colnames(tmp) = c('chr','midPos','Nsites','FST')
    tmp = tmp[which(tmp$Nsites>=500),]
    fstclim = rbind(fstclim,tmp)
    rm(tmp)
  }
  fstclim$bin=seq(1:dim(fstclim)[1])
  fstclim = na.omit(fstclim)
  fstclim$FST[fstclim$FST<0]=0
  fstclim$POS=0
  chrom_end=array(0,6)
  for(c in seq(1:5)){chrom_end[c+1] = chrom_end[c] + max(fstclim$midPos[fstclim$chr==chromlist[c]])}
  for(c in seq(1:5)){
    fstclim$POS[fstclim$chr==chromlist[c]] = fstclim$midPos[fstclim$chr==chromlist[c]] + chrom_end[c]
  }
  fstclim = fstclim[(fstclim$POS > chrom_end[1]+0.1e6 & fstclim$POS < chrom_end[2]-0.1e6)|
                      (fstclim$POS > chrom_end[2]+0.1e6 & fstclim$POS < chrom_end[3]-0.1e6)|
                      (fstclim$POS > chrom_end[3]+0.1e6 & fstclim$POS < chrom_end[4]-0.1e6)|
                      (fstclim$POS > chrom_end[4]+0.1e6 & fstclim$POS < chrom_end[5]-0.1e6)|
                      (fstclim$POS > chrom_end[5]+0.1e6 & fstclim$POS < chrom_end[6]-0.1e6),]
  fstclim$COMP = comp
  fstclim$cuto = quantile(fstclim$FST,0.995) #mean(fstclim$FST) + 5*sd(fstclim$FST) 
  cu[n] = quantile(fstclim$FST,0.995)
  topIvm = fstclim[fstclim$FST>fstclim$cuto,]
  topIvm$COMP = comp
  fsttot = rbind(fsttot,fstclim)
  topIvmTot = rbind(topIvmTot,topIvm)
  rm(fstclim,topIvm)
}
fsttot = fsttot[-1,]
fsttot = fsttot[fsttot$Nsites>=500 & fsttot$Nsites<30000,]

fsttot$COMP = factor(gsub('FRG','GUA',fsttot$COMP)) ## transform Guadeloupe annotation

####--- Common outliers

setwd("./data/NUC/CLIM/")

fstclimA.climB=fsttot[fsttot$COMP=='FRG.NAM',]
fstclimA.climC=fsttot[fsttot$COMP=='FRA.1.NAM',]
fstclimB.climC=fsttot[fsttot$COMP=='FRA.1.FRG',]

outAB = fstclimA.climB[fstclimA.climB$FST>=fstclimA.climB$cuto,]
outAB$comp = 'AB'
outAC = fstclimA.climC[fstclimA.climC$FST>=fstclimA.climC$cuto,]
outAC$comp = 'AC'
outBC = fstclimB.climC[fstclimB.climC$FST>=fstclimB.climC$cuto,]
outBC$comp = 'BC'
outALL = rbind(outAB,outAC,outBC)

##-- 
tdf = data.frame(table(paste0(outALL$chr,outALL$midPos)))
tdf[tdf$Freq==3,]
#                 Var1 Freq
# 51 1_Celeg_TT6925000    3
# 52 1_Celeg_TT6935000    3
# 53 1_Celeg_TT6945000    3
# 54 1_Celeg_TT6955000    3
# 55 1_Celeg_TT6965000    3
# 56 1_Celeg_TT6975000    3
# 57 1_Celeg_TT6985000    3
# 58 1_Celeg_TT6995000    3
# 75 1_Celeg_TT7165000    3
# 76 1_Celeg_TT7175000    3
# 78 1_Celeg_TT7195000    3
# 79 1_Celeg_TT7205000    3

tdf[tdf$Freq==2,]

o.AB.AC=merge(outAB[,c('chr','midPos')],outAC[,c('chr','midPos')],by=c('chr','midPos'))
o.AB.AC$comp='AB.AC'
dim(o.AB.AC)
#[1] 50 2

o.BC.AC=merge(outBC[,c('chr','midPos')],outAC[,c('chr','midPos')],by=c('chr','midPos'))
o.BC.AC$comp='BC.AC'
dim(o.BC.AC)
#[1] 21 2

o.AB.BC=merge(outAB[,c('chr','midPos')],outBC[,c('chr','midPos')],by=c('chr','midPos'))
o.AB.BC$comp='AB.BC'
dim(o.AB.BC)
#[1] 30 2

tot = rbind(o.AB.AC,o.BC.AC,o.AB.BC)
tot = tot[order(tot$chr,tot$midPos),]
table(tot$chr,tot$comp)
a=merge(o.AB.AC,o.BC.AC,by=c('chr','midPos'))
b=merge(a,o.AB.BC,by=c('chr','midPos'))
b
#           chr  midPos comp.x comp.y  comp
# 1  1_Celeg_TT 6925000  AB.AC  BC.AC AB.BC
# 2  1_Celeg_TT 6935000  AB.AC  BC.AC AB.BC
# 3  1_Celeg_TT 6945000  AB.AC  BC.AC AB.BC
# 4  1_Celeg_TT 6955000  AB.AC  BC.AC AB.BC
# 5  1_Celeg_TT 6965000  AB.AC  BC.AC AB.BC
# 6  1_Celeg_TT 6975000  AB.AC  BC.AC AB.BC
# 7  1_Celeg_TT 6985000  AB.AC  BC.AC AB.BC
# 8  1_Celeg_TT 6995000  AB.AC  BC.AC AB.BC
# 9  1_Celeg_TT 7165000  AB.AC  BC.AC AB.BC
# 10 1_Celeg_TT 7175000  AB.AC  BC.AC AB.BC
# 11 1_Celeg_TT 7195000  AB.AC  BC.AC AB.BC
# 12 1_Celeg_TT 7205000  AB.AC  BC.AC AB.BC

#                 AB.AC AB.BC BC.AC
# 1_Celeg_TT         34    30    14
# 3_Celeg_TT          8     0     4
# 4_Celeg_TT          3     0     3
# 5_Celeg_TT_axel     5     0     0

tot[tot$chr=='4_Celeg_TT',]
# chr   midPos  comp
# 69 4_Celeg_TT 16495000 BC.AC
# 70 4_Celeg_TT 16505000 BC.AC
# 71 4_Celeg_TT 16515000 BC.AC
# 43 4_Celeg_TT 39595000 AB.AC
# 44 4_Celeg_TT 39605000 AB.AC
# 45 4_Celeg_TT 39615000 AB.AC

tot[tot$chr=='5_Celeg_TT_axel',]
#               chr   midPos  comp
# 49 5_Celeg_TT_axel  7345000 AB.AC
# 50 5_Celeg_TT_axel  7355000 AB.AC
# 46 5_Celeg_TT_axel 35855000 AB.AC
# 47 5_Celeg_TT_axel 35865000 AB.AC
# 48 5_Celeg_TT_axel 35875000 AB.AC

#               chr midPos.Min. midPos.1st Qu. midPos.Median midPos.Mean midPos.3rd Qu. midPos.Max.
# 1      1_Celeg_TT     6755000        6955000       7060000     7241795        7412500     8205000
# 2      3_Celeg_TT    10735000       31162500      31265000    27880000       31387500    31415000
# 3      4_Celeg_TT    16495000       16507500      28055000    28055000       39602500    39615000
# 4 5_Celeg_TT_axel     7345000        7355000      35855000    24459000       35865000    35875000

out=tot[,1:2]
out$c = substr(out$chr,1,1)
out$up = out$midPos -5000
out$dw = out$midPos +5000
out$chr = NULL
out$midPos=NULL

##-- Output positions of common hits for blast and gene identification
write.table(out,file='./ANGSD_clim_pw_3pop.txt',quote=F,row.names=F)

##-- Prepare suplementary table 5
##- This file was obtained with convert_hits_to_gene.sh using ANGSD_clim_pw_3pop.txt as an input
hits = read.table(file='./climate_FST_hits.txt',sep="\t",
                header=F,fill=T)
colnames(hits)=c('chr','midPos','GeneID','C.elegans ortholog')
st5 = merge(tot,hits,by=c('chr','midPos'))
st5$comp=as.factor(st5$comp)
levels(st5$comp)=c('Arid vs. Tropical & Arid vs. Temperate',
                   'Arid vs. Tropical & Temperate vs. Tropical',
                   'Temperate vs. Tropical & Arid vs. Temperate')

write.csv(st5,file=paste0(draft,'Supplementary_table5_HitsANGSDclim3pop.csv'),quote=F,row.names=F)


####------ Plot Figure 4b
levels(fsttot$COMP) = c("Temperate vs. Tropical (FRA.1 vs. GUA)",
                        "Temperate vs. Arid (FRA.1 vs. NAM)",
                        "Arid vs. Tropical (NAM vs. GUA)")

# pdf(file=paste0(draft,'Figure4a.pdf'),width=14,height=8)
# ggplot(fsttot,aes(x=POS/1e6,y=FST,col=chr)) +
#   facet_wrap(~ COMP,ncol=1) +
#   geom_hline(data=fsttot,aes(yintercept = cuto), lty=2,lwd=.3) +
#   scale_color_manual(values=chrom.colors) + geom_point(size=.2) +
#   scale_y_continuous(limits=c(0,1),breaks = seq(0,1,.2)) +
#   scale_x_continuous(limits=c(0,240),breaks = seq(0,240,10)) +
#   xlab('Position (Mbp)')+ theme(text=element_text(size=16),legend.position ='none') +
#   ggtitle('')
# dev.off()

fst3 = fsttot[fsttot$chr=='3_Celeg_TT',]
fst3$COMP = factor(fst3$COMP)
levels(fst3$COMP) = c("Temperate vs. Tropical (FRA.1 vs. GUA)",
                      "Temperate vs. Arid (FRA.1 vs. NAM)",
                      "Arid vs. Tropical (NAM vs. GUA)")
ggplot(fst3,aes(x=midPos/1e6,y=FST,col=chr)) +
  facet_wrap(~ COMP,ncol=1) +
  geom_hline(data=fst3,aes(yintercept = cuto), lty=2,lwd=.3) +
  scale_color_manual(values=chrom.colors[3]) + geom_point(size=.4) +
  scale_y_continuous(limits=c(0,1),breaks = seq(0,1,.2)) +
  scale_x_continuous(limits=c(0,45),breaks = seq(0,45,10)) +
  xlab('Position (Mbp)')+ theme(text=element_text(size=20),legend.position ='none') +
  ggtitle('')

####------ Look for pi within each pop

##- get  region boundaries
summary(fst3$midPos[fst3$FST>fst3$cuto & fst3$midPos>30e6 & fst3$midPos<35e6])
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 31135000 31265000 31355000 31331341 31415000 31505000 

piclim = t.HI[t.HI$Chr=='5_Celeg_TT_axel' & t.HI$WinCenter>31155000 & t.HI$WinCenter<31415000 & t.HI$pop %in% c('FRA.1','FRG','NAM'),]

aggregate(tP/nSites ~ pop,data=piclim,FUN=mean)
#     pop tP/nSites
# 1 FRA.1   0.00778
# 2   FRG   0.01336
# 3   NAM   0.01278
piclim$pop = factor(piclim$pop)
m = lm(tP/nSites ~ pop,data=piclim)

summary(m)
# Call:
#   lm(formula = tP/nSites ~ pop, data = piclim)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.007480 -0.004597 -0.000367  0.004733  0.009351 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.00778    0.00104    7.51    2e-10 ***
#   popFRG       0.00558    0.00150    3.73   0.0004 ***
#   popNAM       0.00500    0.00153    3.27   0.0017 ** 

summary(glht(m, mcp(pop="Tukey")))
# Fit: lm(formula = tP/nSites ~ pop, data = piclim)
# 
# Linear Hypotheses:
#                     Estimate Std. Error t value Pr(>|t|)   
#   FRG - FRA.1 == 0  0.005579   0.001496    3.73   0.0012 **
#   NAM - FRA.1 == 0  0.005004   0.001532    3.27   0.0049 **
#   NAM - FRG == 0   -0.000575   0.001562   -0.37   0.9282   

picols = c('#0c2c84','#41ab5d','#cc4c02')
piclim$pop = as.character(piclim$pop)
piclim$pop[piclim$pop=='FRG']='GUA'
piclim$pop = factor(piclim$pop)
pipc = ggplot(piclim,aes(x=WinCenter/1e3,y=tP/nSites,col=pop,fill=pop)) + 
  geom_vline(xintercept = 31179783/1e3,lty=3) +
  geom_vline(xintercept = 31193847/1e3,lty=3) +
  geom_point(size=4,alpha=.4) + geom_line(alpha=.5) +
  scale_color_manual(values=picols) +
  xlab('Position (Kbp)') + ylab('Nucleotide diversity') +
  theme(text=element_text(size=20),legend.position ='top',legend.title=element_blank())


####------ FST over the Pc gene region
setwd("./data/CLIM/")
###--- TROP vs. ARID: FRG vs. NAM (A vs. B)
###--- TROP vs. TEMP: FRG vs. FRA (A vs. C)
###--- ARID vs. TEMP: NAM vs. FRA (B vs. C)

complist = c('FRA.1.FRG','FRA.1.NAM','FRG.NAM') 
fstpc = NULL 

for(c in complist){
  temp = read.table(file=paste0(c,'.fst.txt'))
  temp$COMP = c
  fstpc = rbind(fstpc,temp)
}
colnames(fstpc)= c('chr','midPos','Nsites','FST','COMP')
fstpc$FST[fstpc$FST<0] = 0
fstpc$COMP= factor(gsub('FRG','GUA',fstpc$COMP))
colpc = c('#0c2c84','#41ab5d','#ec7014')
options(digits = 3)
fpc=ggplot(fstpc,aes(x=midPos/1000,y=FST,col=COMP)) +
  scale_color_manual(values = colpc) + 
  geom_point(size=fstpc$FST*5,alpha=fstpc$FST*0.9) +
  scale_x_continuous(breaks = seq(31179783/1000,31193847/1000,2), limits=c(31179783/1000,31193847/1000)) +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
  xlab('Position (Kbp)') + ylab('FST') +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(text=element_text(size=20),legend.position ='top',legend.title=element_blank())

fpc
## Draw gene model
require(Sushi)

## Previous coordinates in V3
# 31179782-31193850
# New coord in latest reference: 31179347-31194248

# load the gtf file for genome V4
gtf.gr <- rtracklayer::import(con ="./data/haemonchus_contortus.PRJEB506.WBPS13.annotations.gtf.gz", format = "gtf" )
gtf.gr <- sortSeqlevels(gtf.gr) # make sure it is sorted
gtf.gr <- sort(gtf.gr)
gtf.gr = data.frame(gtf.gr)
pc.gene = gtf.gr[grep('HCON_00087060',gtf.gr$gene_id),c('seqnames','start','end','gene_name','score','strand','type')]
pc.gene = pc.gene[grep('exon',pc.gene$type),]
pc.gene$strand =1
colnames(pc.gene) = c('chrom','start','stop','gene','score','strand','type')
chromstart =31179000
chromend =31195500
chrom='hcontortus_chr3_Celeg_TT_arrow_pilon'

pdf(file=paste0(draft2,'Figure4b_bottom.pdf'),width=14,height=8)
pg=plotGenes(pc.gene,chrom,chromstart,chromend ,
               types=pc.gene$type,maxrows=1,bheight=0.1,plotgenetype="box",
               bentline=F,labeloffset=.3,fontsize=1.2,arrowlength = 0.01,labeltext=F)
chrom='Chrom. III'
labelgenome( chrom, chromstart,chromend,n=8,scale="Mb",cex=1.4)
dev.off()

pdf(file=paste0(draft2,'Figure4b_top.pdf'),width=14,height=8)
multiplot(pipc,fpc,cols=1)
dev.off()

## GO term enrichment analysis
myInterestingGenes = read.table(file='./gene_list.txt',header=F,sep='\t')
colnames(myInterestingGenes)[1]='V1'
geneNames = read.table(file='../XPCLR/blastALL01p/HcoGeneList.txt',header=T,sep='\t')
geneList <- factor(as.integer(geneNames$Gene.stable.ID %in% myInterestingGenes$V1))
names(geneList) <- geneNames$Gene.stable.ID
str(geneList)
geneID2GO = readMappings(file=paste0(draft,'Hco_GO_topGO.txt')) #-- Total GO

##-- Molecular / Cellular Componenent / Biological Processes //// classicFisher based on geneCounts; KS based on gene scores/ranks
goterm=c('MF','CC','BP')
for(i in goterm){
  GOdata <- new("topGOdata", ontology = i, allGenes = geneList, 
                nodeSize = 10,
                annot = annFUN.gene2GO, gene2GO = geneID2GO)
  ##- Fisher
  resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  
  ##- KS test <> based on gene score + weight01 = account for hierarchy <> stat correction
  resultKS <- runTest(GOdata, algorithm = "weight01", statistic = "ks")
  
  #The elim method was design to be more conservative than the classic method
  resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
  
  #- Table of all results - Classic Fisher
  topn=0
  if(sum(score(resultKS) < .01)>1){topn = sum(score(resultKS) < .01)}else{topn=2}
  allRes1 <- GenTable(GOdata, classicFisher = resultFisher,
                      classicKS = resultKS, elimKS = resultKS.elim,
                      orderBy = "classicKS", ranksOf = "classicFisher",
                      topNodes = topn,numChar=100000)
  topn=0
  if(sum(score(resultKS) < .01)>1){topn = sum(score(resultKS) < .01)}else{topn=2}
  allRes2 <- GenTable(GOdata, classicFisher = resultFisher,
                      classicKS = resultKS, elimKS = resultKS.elim,
                      orderBy = "classicFisher", ranksOf = "classicFisher",
                      topNodes = topn,numChar=100000)
  
  write.csv(allRes1,file=paste("./GO_hitsALL_p05_KS",i,".csv",sep="")) ##-- Output x2 significant GO 
  write.csv(allRes2,file=paste("./GO_hitsALL_p05_Fisher",i,".csv",sep="")) ##-- Output x2 significant GO 
  
  #--- Graph
  # printGraph(GOdata, resultKS, firstSigNodes = 10,
  #            fn.prefix = paste("./tGO_hitsALL_p05_",i,sep=''),
  #            useInfo = "all", pdfSW = TRUE)
}

### BP Fisher
#GO:0003729
GOdata <- new("topGOdata", ontology = 'MF', allGenes = geneList, 
              nodeSize = 10,
              annot = annFUN.gene2GO, gene2GO = geneID2GO)
htr = unlist(genesInTerm(GOdata,"GO:0003729")) ##- transport
myInterestingGenes$V1[which(myInterestingGenes$V1 %in% htr)]
#[1] HCOI00107200 HCOI00456600 // cpb-1

#GO:0006417
GOdata <- new("topGOdata", ontology = 'BP', allGenes = geneList, 
              nodeSize = 10,
              annot = annFUN.gene2GO, gene2GO = geneID2GO)
htr = unlist(genesInTerm(GOdata,"GO:0006417")) ##- transport
myInterestingGenes$V1[which(myInterestingGenes$V1 %in% htr)]
#[1] HCOI00107200 HCOI00456600

#####====== GRADIENT FOREST ANALYSIS = which variable best explain variation ? ====######## 
library(gradientForest)
library(raster)
setwd("./data/NUC/CLIM/")

toad.env.data = data.frame(iso=unique(df$iso[df$iso!='CAP']),lon=unique(df$lon[df$iso!='CAP']),lat=unique(df$lat[df$iso!='CAP']))
files2 <- c(list.files(path="./data/wc2-5", pattern="*?\\.bil$"))
files2 <- paste0("./data/wc2-5/", files2)
r2 <- sapply(files2, function(x) { raster(x) } )
s2 <- stack(r2)
names(s2) <- gsub(".bil", "", names(s2))
names(s2) <- gsub("X.data.wc2.5.", "", names(s2))
toad.env.data <- data.frame(toad.env.data, extract(s2, unique(df[df$iso!='CAP',c('lon','lat')]) ) )
mask = s2
cellIDs = extract(s2,unique(df[df$iso!='CAP',c("lon","lat")]),cellnumbers=T)[,1]
all.env <- data.frame(cell=cellIDs, toad.env.data)

all.env.pred = all.env[all.env$iso %in% c('AUS.2','FRA.1','ZAI','FRG','MOR','NAM','IND','STO'),]

# read in data file with minor allele freqs & env/space variables
# GdF.1-5.in files were prepared from ANGSD mafs using the create_GradientForest_input_ANGSD.py script
gfData <- read.table("GdF.1.in",header=F,sep='\t')
for(i in 2:5){
  tmp = read.table(paste0("GdF.",i,".in"),header=F,sep='\t')
  gfData <- cbind(gfData,tmp[,21:ncol(tmp)])}
dim(gfData)
gfData$V2 = NULL
envGF <- gfData[,2:20] # get climate & MEM variables
for(i in 1:dim(envGF)[2]){colnames(envGF)[i] = paste0('bio',i)}

## build individual SNP datasets
x = seq(21,dim(gfData)[2],1)
SNPs_ref <- gfData[,x] #
candidatesnodata.index <- c() 
##-- Remove SNPs with same obs in less than 8 pop
for (j in (1 : ncol(SNPs_ref)))   {
  if (length(unique(SNPs_ref[ ,j])) < 8){candidatesnodata.index <- append(candidatesnodata.index,j)}
}
SNPs_ref <- SNPs_ref[ , - candidatesnodata.index]

# GRADIENT FOREST MODELING -----------------------------------------------------
maxLevel <- log2(0.368*nrow(envGF)/2) #account for correlations, see ?gradientForest

## Fit gf models for reference SNPs
set.seed(64354) ## reproductibility

## Var reduction
library("corrplot")
require(ade4)
require(factoextra)

cormat = rcorr(as.matrix(envGF))

supplENVA = corrplot(cormat$r, type="upper", order="hclust", 
                     tl.col="black", tl.srt=45)

pc = dudi.pca(envGF, scannf = FALSE, nf = 2)
envGFred =  pc$co

## Compute distance betw. bioclim from coords 
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

a = matrix(0,nrow(pc$co),nrow(pc$co))
for(i in 1:nrow(pc$co)){
  for(j in 1:nrow(pc$co)){
    a[i,j]=euc.dist(pc$co[i,1:2],pc$co[j,1:2])
  }
}
colnames(a) = rownames(pc$co)
rownames(a) = rownames(pc$co)

## Keep bioclim with mean dist > 0.85
colMeans(a)[colMeans(a)>.85]
# bio2      bio4      bio5      bio7     bio10     bio15 
# 1.3229741 1.5042411 1.0588528 1.4446885 0.8677287 1.1947064

## For others, 3 clusters are found; Within each, select variable closest to every other as a summary 
#clus1
which.min(colSums(a)[c('bio6','bio3','bio11','bio9','bio1','bio8')])
#bio9
#clus2
which.min(colSums(a[,c('bio17','bio14','bio19')])[c('bio17','bio14','bio19')])
#bio19
#clus3
which.min(colSums(a[,c('bio12','bio13','bio16','bio18')])[c('bio12','bio13','bio16','bio18')])
#bio13 
# here we pick bio12 instead of bio13 because it is more biologically relevant 
# importance with bio13 gives similar output, i.e. bio13, bio7 most important var.

uncorbio = c(4,7,2,15,9,5,10,17,12)
uncorbion = length(uncorbio)
vecol=rep('black',19)
vecol[uncorbio]='red'

supplENVB = fviz_pca_var(pc, axes = c(1, 2), geom = c("arrow", "text"),
                         label = "all", invisible = "none", labelsize = 5,
                         col.var = vecol, alpha.var = 1, col.quanti.sup = "blue",
                         col.circle = "grey70",
                         select.var = list(name =NULL, cos2 = NULL, contrib = NULL)) + 
  ggtitle('') + theme(legend.position = 'none',text = element_text(size=16))

pdf(file=paste0(draft2,'Supplementary_Figure15a.pdf'),width=14,height=14)
print(supplENVA)
dev.off()

pdf(file=paste0(draft2,'Supplementary_Figure15b.pdf'),width=8,height=8)
supplENVB
dev.off()

##--- importance computation
for(i in 1:length(uncorbio)){uncorbion[i]=paste0("bio",uncorbio[i])}
envGFred2 = envGF[,uncorbion]
n = 500
gfRef1 <- gradientForest(cbind(envGFred2, SNPs_ref), predictor.vars=colnames(envGFred2),
                         response.vars=colnames(SNPs_ref), ntree = n,
                         maxLevel = maxLevel, trace=T, corr.threshold = 0.5)
importance(gfRef1)
# bio12 bio7  bio4  bio2  bio5      bio10      bio15      bio19 bio9 
# 0.02913211 0.02585279 0.02219931 0.02074462 0.01868949 0.01836759 0.01655033 0.01259061 0.01005534 

##- Plot variable importance
gfRef1=c(0.02913211,0.02585279,0.02219931,0.02074462,0.01868949,0.01836759,0.01655033,0.01259061,0.01005534)
names(gfRef1)=c('Annual Precipitation','Temperature Annual Range',
                 'Temperature Seasonality','Mean Diurnal Range','Max Temperature of Warmest Month',
                 'Mean Temperature of Warmest Quarter','Precipitation Seasonality',
                 'Precipitation of Coldest Quarter','Mean Temperature of Driest Quarter')
i = data.frame(gfRef1)
i$var = factor(rownames(i))
i$type = 'Temperature'
i$type[grep('Preci',i$var)] = 'Precipitation'
i$var = factor(i$var,levels=i$var[order(i[,1])])

pdf(file=paste0(draft,'Figure4c.pdf'),width=8,height=8)
ggplot(i,aes(x=var,y=i[,1],fill = type)) +
  coord_flip() + xlab('Climatic variable') + ylab('Proportion of variance explained') +
  scale_y_continuous(limits=c(0,0.03),breaks=seq(0,0.03,0.005)) +
  scale_fill_manual(values=c("#67a9cf","#ef8a62")) +
  theme(legend.position='none',axis.title = element_text(size=16),axis.text= element_text(size=14) ) +
  geom_bar(stat='identity')
dev.off()

##--- importance computation with bio13 instead of bio12
uncorbio2 = c(4,7,2,15,9,5,10,17,13)
for(i in 1:length(uncorbio2)){uncorbion[i]=paste0("bio",uncorbio2[i])}
envGFred2 = envGF[,uncorbion]
n = 500
gfRef1 <- gradientForest(cbind(envGFred2, SNPs_ref), predictor.vars=colnames(envGFred2),
                         response.vars=colnames(SNPs_ref), ntree = n,
                         maxLevel = maxLevel, trace=T, corr.threshold = 0.5)
importance(gfRef1)
# bio13       bio7       bio2       bio5      bio17       bio4      bio10      bio15       bio9 
# 0.02780351 0.02724519 0.02418622 0.02066921 0.02059540 0.02014249 0.01936418 0.01576466 0.01033919 names(importance(gfRef10))

#####====== LFMM ANALYSIS ====########
setwd('./data/NUC/LFMM/')

##-- N unique SNP positions = 200117 // after filtering for missingness and MAF
# complist = list.dirs(path = ".", full.names = F, recursive = TRUE)
# complist = complist[-1]
# complist = data.frame(complist) #unlist(complist))
# colnames(complist)[1]='V1'
# complist$chrom=substr(complist$V1,1,6)
# 
# complist$up=as.integer(as.character(sapply(str_split(substr(complist$V1,8,nchar(as.character(complist$V1))),":"),function(x) x[1])))
# complist$V2=as.integer(as.character(sapply(str_split(substr(complist$V1,8,nchar(as.character(complist$V1))),":"),function(x) x[2])))
# 
# chromlist=unique(substr(complist$chrom,1,6))
# q = 0.05
# 
# for(envar in c(16,21)){ ## envar = bioX+9
#   
#   pvtot = NULL
#   
#   for(i in chromlist){
#     pv = matrix(0,1,2)
#     colnames(pv)=c('pval','pos')
#     
#     tmp=complist[complist$chrom==i,]
#     tmp=tmp[order(tmp$up),]
#     
#     for(k in 1:dim(tmp)[1]){
#       temp = read.table(paste0(i,'.',tmp$up[k],':',tmp$V2[k],'/pval',envar,'.txt'),header=F)
#       colnames(temp)[1]='pval'
#       map = read.table(paste0(i,'.',tmp$up[k],':',tmp$V2[k],'/trimmed.map.file'))
#       temp=cbind(temp,map[,1])
#       colnames(temp)[2]='pos'
#       pv=rbind(pv,temp)
#       rm(temp,map)
#     }
#     pv = pv[-1,] ##- remove init row
#     
#     ##- Aggregate p-values for positions with overlap + remove NA
#     pv = aggregate(pval ~ pos,data=pv,FUN=mean)
#     
#     ## Add chrom info
#     pv$Pos = as.integer(as.character(substr(pv$pos,5,nchar(pv$pos))))
#     pv$Chrom = as.integer(as.character(substr(pv$pos[1],3,3)))
#     
#     ##- Save positions above 5% FDR by CHROM and Bioclim
#     L = length(pv$pval)
#     w = which(sort(pv$pval) < q * (1:L) / L)
#     candidates = order(pv$pval)[w]
#     write.table(pv[candidates,],file=paste0('../',i,'.',colnames(env)[envar],'.txt'),  quote=F,row.names = F)
#  
#     pv$candid = 0
#     pv$candid[candidates] = 1
#     pv$candid = factor(pv$candid)
#     
#     pvtot = rbind(pvtot,pv)
#     rm(pv,candidates)
#   }
#   pvtot$x=seq(1,nrow(pvtot),1)
#   pvtot$Chrom=factor(pvtot$Chrom)
#   rm(pvtot)
# }

###--- Get all significant hits
setwd("./data/NUC/LFMM/")
files = c(list.files(path="./.", pattern="\\.txt$"))
pval = NULL
for(i in files[1:length(files)]){
  if(file.info(i)$size>19 & !(is.na(file.info(i)$size))){
    temp = read.table(i,header=T)
    temp$file = i
    pval = rbind(pval,temp)
  }
}
dim(pval)
#[1] 42  5
length(grep('bio7',pval$file))
#[1] 17
length(grep('bio12',pval$file))
#[1] 25

####------ Plot Figure 4b
pval$chr[pval$Chrom!=5]=paste0(pval$Chrom[pval$Chrom!=5],"_Celeg_TT")
pval$chr[pval$Chrom==5]=paste0(pval$Chrom[pval$Chrom==5],"_Celeg_TT_axel")

chrom_end=c(0,45745000,93100000,136635000,188420000,237215000)
pval$POS = pval$Pos + chrom_end[pval$Chrom]
pval$Variable = 'BIO7'
pval$Variable[grep('bio12',pval$file)] = 'BIO12'

##########------- FIGURE 4a --------------#########
plotLFMM = ggplot(pval,aes(x=POS/1e6,y=-log10(pval),shape=Variable,col=chr))+
  geom_point(size=2) + xlab('') + ylab('-log10(P)') +
  scale_color_manual(values=chrom.colors) +
  scale_y_continuous(limits=c(4,10),breaks = seq(4,10,2)) +
  scale_x_continuous(limits=c(0,240),breaks = seq(0,240,10)) + xlab('Position (Mbp)')+
  theme(text=element_text(size=16),legend.position ='none') #,axis.text.x=element_blank()) 

fsttot$COMP = factor(fsttot$COMP)
levels(fsttot$COMP) = c("Temperate vs. Tropical (FRA.1 vs. GUA)",
                        "Temperate vs. Arid (FRA.1 vs. NAM)",
                        "Arid vs. Tropical (NAM vs. GUA)")

plotFSTangsd = ggplot(fsttot,aes(x=POS/1e6,y=FST,col=chr)) +
  facet_wrap(~ COMP,ncol=1) +
  geom_hline(data=fsttot,aes(yintercept = cuto), lty=2,lwd=.3) +
  scale_color_manual(values=chrom.colors) + geom_point(size=.2) +
  scale_y_continuous(limits=c(0,1),breaks = seq(0,1,.2)) +
  scale_x_continuous(limits=c(0,240),breaks = seq(0,240,10)) +
  xlab('')+ 
  theme(text=element_text(size=16),legend.position ='none') +
  ggtitle('')

library(gridExtra)
gA <- ggplotGrob(plotLFMM)
gB <- ggplotGrob(plotFSTangsd)
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth)
gB$widths[2:5] <- as.list(maxWidth)


pdf(file=paste0(draft2,'Figure4a.pdf'),width=14,height=8)
plot(p4 <- arrangeGrob(
  gB, gA, ncol = 1, heights = c(1.30, 0.30)))
dev.off()
##########------- FIGURE 4a --------------#########

##- Markers found for more than 1 bioclim
nm = names(table(pval$pos)[which(table(pval$pos)>1)])
nmd = pval[pval$pos %in% nm,]
nmd = nmd[order(nmd$pos),]
nmd
# [1] pos   pval  Pos   Chrom file 
# <0 rows> (or 0-length row.names)

##-- Hits by chromosome
table(pval$Chrom)
#  1  2  3  4  5 
# 13  9 11  8  1

table(pval$fil)
# CHROM1.bio12.txt  CHROM1.bio7.txt CHROM2.bio12.txt CHROM3.bio12.txt  CHROM3.bio7.txt CHROM4.bio12.txt  CHROM4.bio7.txt 
# 10  3  9  5  6  1  7 
# CHROM5.bio7.txt 
# 1 
pval$variable=substr(pval$file,8,12)
table(pval$Chrom,pval$variable)
#   bio12 bio7.
# 1    10     3
# 2     9     0
# 3     5     6
# 4     1     7
# 5     0     1

head(pval[order(pval$pval),],n=10)
#             pos              pval      Pos Chrom             file
# 28 M.3.26824240 0.000000002265664 26824240     3  CHROM3.bio7.txt
# 29 M.3.26824197 0.000000138644050 26824197     3  CHROM3.bio7.txt
# 11  M.1.4436114 0.000000602015089  4436114     1  CHROM1.bio7.txt
# 14 M.2.19715876 0.000000848139437 19715876     2 CHROM2.bio12.txt
# 35 M.4.15905191 0.000000886984874 15905191     4  CHROM4.bio7.txt
# 23 M.3.13983070 0.000001179506502 13983070     3 CHROM3.bio12.txt
# 15    M.2.17957 0.000001349731871    17957     2 CHROM2.bio12.txt
# 1   M.1.7178137 0.000001389242038  7178137     1 CHROM1.bio12.txt
# 2   M.1.7178216 0.000002118376537  7178216     1 CHROM1.bio12.txt
# 24 M.3.37491795 0.000002125068323 37491795     3 CHROM3.bio12.txt

## Output locations for BLAST and gene query
c = pval$Chrom
up = pval$Pos - 5000
dw = pval$Pos + 5000
t = cbind(c,up,dw)
write.table(t,file='./tot_hits_lfmm.txt',quote=F,row.names=F) 

# Gene stable ID 	Caenorhabditis elegans (PRJNA13758) gene name 	Gene description
# HCOI00078600 		ISE/inbred ISE genomic scaffold, scaffold_pathogens_Hcontortus_scaffold_1036 [Source:UniProtKB/TrEMBL;Acc:U6NML8]
# HCOI00198200 	anp-1 	Peptidase M1 domain containing protein [Source:UniProtKB/TrEMBL;Acc:U6NQT5]
# HCOI00312500 	ZK892.3 	
# HCOI00312500 	Y57G11C.23 	= Solute Carrier
# HCOI00312500 	C53B4.1 	= SLC + respond to daf-12 (controls dauer) + (GABAergic neuron)
# HCOI00712600 		
# HCOI01130100 		
# HCOI01461000 	Y5F2A.4 	Zinc finger domain containing protein [Source:UniProtKB/TrEMBL;Acc:U6PL64]
# HCOI01509500 	Y37E3.1 	Tetratricopeptide TPR2 domain containing protein [Source:UniProtKB/TrEMBL;Acc:U6PHY9] / TTC5 ortho; stress response DNA dam + heat shock
# HCOI01769300 		
# HCOI02015300 	ugt-45 	UDP-glucuronosyl UDP-glucosyltransferase domain containing protein [Source:UniProtKB/TrEMBL;Acc:U6PY02]
# HCOI02027500 	nlp-40 	= secreted by the intestine and mediates the pacemaker function of the intestine that regulates the rhythmic defecation motor program

###--- Pos and gene id
write.csv(pval,file='LFMM_42pval.csv',quote=F,row.names = F)
head(pval[order(pval$pval),])
# pos              pval      Pos Chrom             file
# 28 M.3.26824240 0.000000002265664 26824240     3  CHROM3.bio7.txt
# 29 M.3.26824197 0.000000138644050 26824197     3  CHROM3.bio7.txt
# 11  M.1.4436114 0.000000602015089  4436114     1  CHROM1.bio7.txt
# 14 M.2.19715876 0.000000848139437 19715876     2 CHROM2.bio12.txt
# 35 M.4.15905191 0.000000886984874 15905191     4  CHROM4.bio7.txt
# 23 M.3.13983070 0.000001179506502 13983070     3 CHROM3.bio12.txt

## 38 unique genes / 106 hits within gene
d = array(-1,nrow(pval))
p = d
for(i in 1:nrow(pval)){
  pos = outALL$midPos[as.integer(as.character(substr(outALL$chr,1,1)))==pval$Chrom[i]]
  subs = outALL[as.integer(as.character(substr(outALL$chr,1,1)))==pval$Chrom[i],]
  j = which.min(abs(pos - pval$Pos[i]))
  d[i] = abs(pos[j] - pval$Pos[i])
  p[i] = paste0(subs$ch[j],'.',subs$midPos[j])
}
d
# [1]     863     784     839     383     904','0   74301      89    1462   74173  128114  109279  597432    2876  630043    2744
# [17] 1443403  273649    2795  274060  188224  273640    9930    3205    4185    9897    2660     240     197      10  797937     172
# [33] 4    7880    4809     485  210808    4859  419651  477259     111     159

p
# [1] "1_Celeg_TT.7179000"','"1_Celeg_TT.7179000"','"1_Celeg_TT.7179000"','"1_Celeg_TT.7184000"      
# [5] "1_Celeg_TT.7179000"','"1_Celeg_TT.7183000"','"1_Celeg_TT.25288000"      "1_Celeg_TT.7183000"      
# [9] "1_Celeg_TT.7179000"','"1_Celeg_TT.25288000"      "1_Celeg_TT.4308000"','"1_Celeg_TT.4308000"      
# [13] "1_Celeg_TT.40128000"      "2_Celeg_TT.19713000"      "2_Celeg_TT.648000"',' "2_Celeg_TT.19713000"     
# [17] "2_Celeg_TT.10400000"      "2_Celeg_TT.33289000"      "2_Celeg_TT.19713000"      "2_Celeg_TT.33289000"     
# [21] "2_Celeg_TT.2196000"','"2_Celeg_TT.33289000"      "3_Celeg_TT.13993000"      "3_Celeg_TT.37495000"     
# [25] "3_Celeg_TT.13993000"      "3_Celeg_TT.13993000"      "3_Celeg_TT.37495000"      "3_Celeg_TT.26824000"     
# [29] "3_Celeg_TT.26824000"      "3_Celeg_TT.26823000"      "3_Celeg_TT.35615000"      "3_Celeg_TT.26824000"     
# [33] "3_Celeg_TT.26823000"      "4_Celeg_TT.21135000"      "4_Celeg_TT.15910000"      "4_Celeg_TT.30390000"     
# [37] "4_Celeg_TT.24052000"      "4_Celeg_TT.15910000"      "4_Celeg_TT.36233000"      "4_Celeg_TT.51288000"     
# [41] "4_Celeg_TT.16545000"      "5_Celeg_TT_axel.41115000"

pval = cbind(pval,p,d)

pval[which(d<10000),]
# pos','  pval      Pos Chrom','      file
# 1   M.1.7178137 1.389242e-06  7178137     1 CHROM1.bio12.txt
# 2   M.1.7178216 2.118377e-06  7178216     1 CHROM1.bio12.txt
# 3   M.1.7178161 3.339127e-06  7178161     1 CHROM1.bio12.txt
# 4   M.1.7183617 4.443068e-06  7183617     1 CHROM1.bio12.txt
# 5   M.1.7178096 5.302831e-06  7178096     1 CHROM1.bio12.txt
# 6   M.1.7183000 1.286546e-05  7183000     1 CHROM1.bio12.txt
# 8   M.1.7183089 3.089454e-05  7183089     1 CHROM1.bio12.txt
# 9   M.1.7177538 4.738724e-05  7177538     1 CHROM1.bio12.txt
# 14 M.2.19715876 8.481394e-07 19715876     2 CHROM2.bio12.txt
# 16 M.2.19715744 4.355767e-06 19715744     2 CHROM2.bio12.txt
# 19 M.2.19715795 2.837455e-05 19715795     2 CHROM2.bio12.txt
# 23 M.3.13983070 1.179507e-06 13983070     3 CHROM3.bio12.txt
# 24 M.3.37491795 2.125068e-06 37491795     3 CHROM3.bio12.txt
# 25 M.3.13988815 6.296489e-06 13988815     3 CHROM3.bio12.txt
# 26 M.3.13983103 1.333593e-05 13983103     3 CHROM3.bio12.txt
# 27 M.3.37492340 2.995079e-05 37492340     3 CHROM3.bio12.txt
# 28 M.3.26824240 2.265664e-09 26824240     3  CHROM3.bio7.txt
# 29 M.3.26824197 1.386441e-07 26824197     3  CHROM3.bio7.txt
# 30 M.3.26823010 2.524743e-06 26823010     3  CHROM3.bio7.txt
# 32 M.3.26824172 7.267905e-06 26824172     3  CHROM3.bio7.txt
# 33 M.3.26823004 1.613099e-05 26823004     3  CHROM3.bio7.txt
# 34 M.4.21127120 3.462805e-06 21127120     4 CHROM4.bio12.txt
# 35 M.4.15905191 8.869849e-07 15905191     4  CHROM4.bio7.txt
# 36 M.4.30389515 3.682638e-06 30389515     4  CHROM4.bio7.txt
# 38 M.4.15905141 2.512507e-05 15905141     4  CHROM4.bio7.txt
# 41 M.4.16544889 3.101342e-05 16544889     4  CHROM4.bio7.txt
# 42 M.5.41114841 3.908132e-06 41114841     5  CHROM5.bio7.txt

cong = pval[which(d<10000),]
write.csv(cong,file='ANGSD_LFMM.csv',row.names = F,quote=F)
cong[order(cong$pval),]
# pos','  pval      Pos Chrom','      file
# 28 M.3.26824240 2.265664e-09 26824240     3  CHROM3.bio7.txt
# 29 M.3.26824197 1.386441e-07 26824197     3  CHROM3.bio7.txt
# 14 M.2.19715876 8.481394e-07 19715876     2 CHROM2.bio12.txt
# 35 M.4.15905191 8.869849e-07 15905191     4  CHROM4.bio7.txt
# 23 M.3.13983070 1.179507e-06 13983070     3 CHROM3.bio12.txt
# 1   M.1.7178137 1.389242e-06  7178137     1 CHROM1.bio12.txt
# 2   M.1.7178216 2.118377e-06  7178216     1 CHROM1.bio12.txt
# 24 M.3.37491795 2.125068e-06 37491795     3 CHROM3.bio12.txt
# 30 M.3.26823010 2.524743e-06 26823010     3  CHROM3.bio7.txt
# 3   M.1.7178161 3.339127e-06  7178161     1 CHROM1.bio12.txt
# 34 M.4.21127120 3.462805e-06 21127120     4 CHROM4.bio12.txt
# 36 M.4.30389515 3.682638e-06 30389515     4  CHROM4.bio7.txt
# 42 M.5.41114841 3.908132e-06 41114841     5  CHROM5.bio7.txt
# 16 M.2.19715744 4.355767e-06 19715744     2 CHROM2.bio12.txt
# 4   M.1.7183617 4.443068e-06  7183617     1 CHROM1.bio12.txt
# 5   M.1.7178096 5.302831e-06  7178096     1 CHROM1.bio12.txt
# 25 M.3.13988815 6.296489e-06 13988815     3 CHROM3.bio12.txt
# 32 M.3.26824172 7.267905e-06 26824172     3  CHROM3.bio7.txt
# 6   M.1.7183000 1.286546e-05  7183000     1 CHROM1.bio12.txt
# 26 M.3.13983103 1.333593e-05 13983103     3 CHROM3.bio12.txt
# 33 M.3.26823004 1.613099e-05 26823004     3  CHROM3.bio7.txt
# 38 M.4.15905141 2.512507e-05 15905141     4  CHROM4.bio7.txt
# 19 M.2.19715795 2.837455e-05 19715795     2 CHROM2.bio12.txt
# 27 M.3.37492340 2.995079e-05 37492340     3 CHROM3.bio12.txt
# 8   M.1.7183089 3.089454e-05  7183089     1 CHROM1.bio12.txt
# 41 M.4.16544889 3.101342e-05 16544889     4  CHROM4.bio7.txt
# 9   M.1.7177538 4.738724e-05  7177538     1 CHROM1.bio12.txt

p[which(d<10000)]
# [1] "1_Celeg_TT.7179000"','"1_Celeg_TT.7179000"','"1_Celeg_TT.7179000"','"1_Celeg_TT.7184000"      
# [5] "1_Celeg_TT.7179000"','"1_Celeg_TT.7183000"','"1_Celeg_TT.7183000"','"1_Celeg_TT.7179000"      
# [9] "2_Celeg_TT.19713000"      "2_Celeg_TT.19713000"      "2_Celeg_TT.19713000"      "3_Celeg_TT.13993000"     
# [13] "3_Celeg_TT.37495000"      "3_Celeg_TT.13993000"      "3_Celeg_TT.13993000"      "3_Celeg_TT.37495000"     
# [17] "3_Celeg_TT.26824000"      "3_Celeg_TT.26824000"      "3_Celeg_TT.26823000"      "3_Celeg_TT.26824000"     
# [21] "3_Celeg_TT.26823000"      "4_Celeg_TT.21135000"      "4_Celeg_TT.15910000"      "4_Celeg_TT.30390000"     
# [25] "4_Celeg_TT.15910000"      "4_Celeg_TT.16545000"      "5_Celeg_TT_axel.41115000"

table(cong$Chrom)
# 1  2  3  4  5 
# 8  3 10  5  1 
length(grep('bio7',cong$file)) #10
length(grep('bio12',cong$file)) #17

##########------- SESSION INFORMATION --------------#########

sessionInfo()
# R version 3.5.0 (2018-04-23)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS  10.14.4
# 
# Matrix products: default
# BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] gridExtra_2.3           Sushi_1.20.0            biomaRt_2.38.0         
# [4] zoo_1.8-4               topGO_2.34.0            SparseM_1.77           
# [7] GO.db_3.7.0             AnnotationDbi_1.44.0    IRanges_2.16.0         
# [10] S4Vectors_0.20.1        Biobase_2.42.0          graph_1.60.0           
# [13] BiocGenerics_0.28.0     dplyr_0.7.8             maptools_0.9-4         
# [16] rgdal_1.3-6             rasterVis_0.45          latticeExtra_0.6-28    
# [19] dismo_1.1-4             vegan_2.5-3             permute_0.9-4          
# [22] RColorBrewer_1.1-2      phangorn_2.4.0          plotrix_3.7-4          
# [25] gplots_3.0.1            geosphere_1.5-7         phytools_0.6-60        
# [28] maps_3.3.0              ggtree_1.14.6           ggmap_2.6.1            
# [31] directlabels_2018.05.22 ape_5.2                 data.table_1.12.0      
# [34] ggrepel_0.8.0           stringr_1.3.1           SNPRelate_1.16.0       
# [37] gdsfmt_1.18.1           reshape2_1.4.3          Hmisc_4.1-1            
# [40] Formula_1.2-3           survival_2.43-3         lattice_0.20-38        
# [43] raster_2.8-4            sp_1.3-1                gradientForest_0.1-17  
# [46] extendedForest_1.6.1    factoextra_1.0.5        ggplot2_3.1.0          
# [49] ade4_1.7-13             corrplot_0.84          
# 
# loaded via a namespace (and not attached):
# [1] backports_1.1.3             fastmatch_1.1-0             plyr_1.8.4                 
# [4] igraph_1.2.2                lazyeval_0.2.1              splines_3.5.0              
# [7] BiocParallel_1.16.5         GenomeInfoDb_1.18.1         digest_0.6.18              
# [10] htmltools_0.3.6             gdata_2.18.0                magrittr_1.5               
# [13] RMTstat_0.3                 checkmate_1.9.1             memoise_1.1.0              
# [16] cluster_2.0.7-1             Biostrings_2.50.2           matrixStats_0.54.0         
# [19] prettyunits_1.0.2           jpeg_0.1-8                  colorspace_1.4-0           
# [22] blob_1.1.1                  xfun_0.4                    crayon_1.3.4               
# [25] RCurl_1.95-4.11             jsonlite_1.6                hexbin_1.27.2              
# [28] bindr_0.1.1                 glue_1.3.0                  gtable_0.2.0               
# [31] zlibbioc_1.28.0             XVector_0.22.0              DelayedArray_0.8.0         
# [34] scales_1.0.0                DBI_1.0.0                   Rcpp_1.0.0                 
# [37] progress_1.2.0              viridisLite_0.3.0           htmlTable_1.13.1           
# [40] tidytree_0.2.1              foreign_0.8-71              bit_1.1-14                 
# [43] mapproj_1.2.6               animation_2.6               httr_1.4.0                 
# [46] htmlwidgets_1.3             acepack_1.4.1               pkgconfig_2.0.2            
# [49] XML_3.98-1.16               nnet_7.3-12                 tidyselect_0.2.5           
# [52] rlang_0.3.1                 munsell_0.5.0               tools_3.5.0                
# [55] RSQLite_2.1.1               yaml_2.2.0                  knitr_1.21                 
# [58] bit64_0.9-7                 caTools_1.17.1.1            purrr_0.2.5                
# [61] RgoogleMaps_1.4.3           bindrcpp_0.2.2              nlme_3.1-137               
# [64] compiler_3.5.0              rstudioapi_0.9.0            png_0.1-7                  
# [67] treeio_1.6.1                clusterGeneration_1.3.4     tibble_2.0.1               
# [70] stringi_1.2.4               Matrix_1.2-15               pillar_1.3.1               
# [73] combinat_0.0-8              bitops_1.0-6                rtracklayer_1.42.1         
# [76] GenomicRanges_1.34.0        R6_2.3.0                    KernSmooth_2.23-15         
# [79] codetools_0.2-16            MASS_7.3-51.1               gtools_3.8.1               
# [82] assertthat_0.2.0            SummarizedExperiment_1.12.0 proto_1.0.0                
# [85] rjson_0.2.20                withr_2.1.2                 GenomicAlignments_1.18.1   
# [88] Rsamtools_1.34.0            mnormt_1.5-5                GenomeInfoDbData_1.2.0     
# [91] hms_0.4.2                   mgcv_1.8-26                 expm_0.999-3               
# [94] quadprog_1.5-5              grid_3.5.0                  rpart_4.1-13               
# [97] tidyr_0.8.2                 coda_0.19-2                 rvcheck_0.1.3              
# [100] numDeriv_2016.8-1           scatterplot3d_0.3-41        base64enc_0.1-3 
