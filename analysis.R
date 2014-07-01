library(vegan)
library(ggplot2)
library(reshape2)
library(plyr)
library(scales)
library(stringr)
library(gplots)
library(RColorBrewer)
library(gridExtra)
library(plotrix)

# in hmp1.v13.hq.otu.counts, remove colname collection and add a new col name to the end (276...)
# the result is .header

occurences = read.csv(file="hmp1.v13.hq.otu.counts.may1", header=T, sep="\t", quote="")
occurences = occurences[1:2799,1:27655]
print(dim(occurences))

# filter columns (species) - remove unclassified and uncounted
lookup = read.csv(file="hmp1.v13.hq.otu.lookup.split", sep="\t")
lookup.classified = subset(lookup, lookup$genus!='unclassified')
occurences <- occurences[,lookup.classified$otu]
print(dim(occurences))

colSum.occurences = colSums(occurences)
print(any(colSum.occurences==0))
occurences <- occurences[,colSum.occurences>0]
colSum.occurences = colSums(occurences)
print(any(colSum.occurences==0))
print(dim(occurences))

rowSum.occurences = rowSums(occurences)
print(any(rowSum.occurences==0))
occurences <- occurences[rowSum.occurences>0,]
rowSum.occurences = rowSums(occurences)
print(any(rowSum.occurences==0))
print(dim(occurences))


# add row (sample) metadata
map = read.csv(file="ppAll_V13_map.no_unavail.txt", sep="\t", row.names = NULL)
map <- map[!duplicated(map$NAP),]
print(dim(occurences))
print(dim(map))
occurences <- merge(map, occurences, by.x="NAP", by.y=0, all.x=F, all.y=F)
print(dim(occurences))
species.cols = seq(from=17,to=dim(occurences)[2])
colSum.occurences = colSums(occurences[,species.cols])
print(any(colSum.occurences==0))
rowSum.occurences = rowSums(occurences[,species.cols])
print(any(rowSum.occurences==0))

# General
sites.df = ddply(occurences, .(HMPBodySite), summarize, n=length(HMPBodySite))
qplot(x=HMPBodySite, y=n, data=sites.df, geom="bar", stat="identity") + 
  theme_bw() + labs(x="Body site", y="# Samples")
ggsave("body site samples distribution.png", width=6, height=4)

ppl.df = ddply(occurences, .(RSID), summarize, n=length(RSID), Sex=unique(Sex)[1], bodysites=length(unique(HMPBodySubsite)), bodysubsites=length(unique(HMPBodySubsite)))
print(paste("Males:",sum(ppl.df$Sex=='male')))
print(paste("Females:",sum(ppl.df$Sex=='female')))
qplot(x=n, data=ppl.df) +
  theme_bw() + labs(x="# Samples per person", y="Count")
ggsave("body site samples distribution per person.png")
qplot(x=bodysites, data=ppl.df) +
  theme_bw() + labs(x="# Body sites sampled per person", y="Count")
ggsave("body site samples distribution per person.png")

# richness
occurences$richness = apply(X=occurences[,species.cols]>0, MARGIN=1, FUN=sum)
qplot(richness, data=occurences, xlab="Species Richness", ylab="Frequnecy", geom="density", fill=Sex, alpha=I(0.5)) + 
  facet_grid(facets=HMPBodySite~.) + theme_bw() + scale_fill_brewer(palette="Set1") + coord_cartesian(xlim=c(0,150))
ggsave("species richness by sex and body site.png")

# rarefaction
## By sex
occur.by.sex = ddply(occurences, .(Sex), function(d) colSums(d[,species.cols]))
row.names(occur.by.sex) = occur.by.sex$Sex
occur.by.sex = subset(occur.by.sex, select=-c(Sex))
step = as.integer(min(rowSums(occur.by.sex)/100))
png("rarefaction by sex.png")
rarecurve(occur.by.sex, step=step, cex=0.8, label=T)
dev.off()

## By site
occur.by.site = ddply(occurences, .(HMPBodySite), function(d) colSums(d[,species.cols]))
row.names(occur.by.site) = occur.by.site$HMPBodySite
occur.by.site = subset(occur.by.site, select=-c(HMPBodySite))
step = as.integer(min(rowSums(occur.by.site)/100))
png("rarefaction by body site.png")
rarecurve(occur.by.site, step=step, cex=0.8)
dev.off()

# diversity
occurences$diversity = exp(diversity(occurences[,species.cols], index="shannon"))
qplot(x=seq_along(diversity), y=diversity, data=occurences, xlab="Sample", ylab="True Diversity", alpha=I(0.5), color=HMPBodySite) + 
  scale_color_brewer(palette="Set1") + theme_bw()
qplot(x=HMPBodySite, y=diversity, data=occurences, color=Sex, position=position_dodge(width=1)) + 
  geom_boxplot(alpha=0, position=position_dodge(width=1)) + scale_color_brewer(palette="Set1") +
  theme_bw() + scale_y_log10() + labs(x="Body Site", y="True Diversity (Shannon)")
ggsave("diversity by sex and body site - boxplot.png")
qplot(diversity, data=occurences, fill=Sex, alpha=I(0.5), xlab="True Diversity (Shannon)", ylab="Frequnecy", geom="density") + 
  facet_grid(facets=HMPBodySite~.) + scale_x_log10() + theme_bw() + scale_fill_brewer(palette="Set1")
ggsave("diversity by sex and body site - distribution.png")

## diversity profile
### the dots are the values, the three lines are min, median, and max
mod = renyi(occur.by.sex, scales=c(0,1,2,3,4,Inf),hill=T)#0, 1/256, 1/128, 1/64, 1/32, 1/16, 1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32, 64, 128, 256, Inf))
png("diversity profile by sex.png")
plot(log(mod), ylab="Log Diversity", xlab="Hill Number", type="p")
dev.off()
mod = renyi(occur.by.site, scales=c(0,1,2,3,4,Inf),hill=T)
png("diversity profile by body site.png")
plot(log(mod), ylab="Log Diversity", xlab="Hill Number", type="p", as.table=F)
dev.off()

# beta diversity
mod = betadiver(occurences[sample(dim(occurences)[1], 100),species.cols])
png("triplot diversity.png")
plot(mod, cex=0.1)
dev.off()

## by site
mod = betadiver(occur.by.site)
png("beta diversity by body sites.png")
plot(mod)
dev.off()

occur.by.subsite = ddply(occurences, .(HMPBodySubsite), function(d) colSums(d[,species.cols]))
row.names(occur.by.subsite) = occur.by.subsite$HMPBodySubsite
occur.by.subsite = subset(occur.by.subsite, select=-c(HMPBodySubsite))
mod = betadiver(occur.by.subsite)
png("beta diversity by body sub sites.png")
plot(mod, cex=0.5)
dev.off() 

## j and sor are similar, with different weights
## z and w  are similar, and are roghly the negative of sor
## z 
betadiver(help=T)

w = betadiver(occurences[, species.cols], "w")
w = betadisper(w, occurences$HMPBodySite)
boxplot(w, main="Whitaker (w) Beta Diversity")

#z = betadiver(occurences[, species.cols], "z")
#boxplot(betadisper(z, occurences$HMPBodySite), main="Z Beta Diversity")

j = betadiver(occurences[, species.cols], "j")
j = betadisper(j, occurences$HMPBodySite)
boxplot(j, main="Jaccard (j) Beta Diversity")

sor = betadiver(occurences[, species.cols], "sor")
sor = betadisper(sor, occurences$HMPBodySite)
boxplot(j, main="Sørensen (sor) Beta Diversity")

r = betadiver(occurences[, species.cols], "r")
r = betadisper(r, occurences$HMPBodySite)
boxplot(r, main="Routledge (r) Beta Diversity")

sim = betadiver(occurences[, species.cols], "sim")
sim = betadisper(sim, occurences$HMPBodySite)
boxplot(r, main="Simpson (sim) Beta Diversity")

gl = betadiver(occurences[, species.cols], "gl")
gl = betadisper(gl, occurences$HMPBodySite)
boxplot(r, main="Lennon (gl) Beta Diversity")


png("beta diversity 3 indx - body site.png", units="in", width=5, height=12, res=300)
def.par <- par(mfrow = c( 3,1 ))
cex.axis = 0.95
#boxplot(w, main="Whitaker", cex.axis=cex.axis)
boxplot(j, main="Jaccard", cex.axis=cex.axis)
#boxplot(sor, main="Sørensen", cex.axis=cex.axis)
#boxplot(r, main="Routledge", cex.axis=cex.axis)
boxplot(sim, main="Simpson", cex.axis=cex.axis)
boxplot(gl, main="Lennon", cex.axis=cex.axis)
par(def.par)
dev.off()


# nMDS
## this takes a long time, so I saved the results
#mds = metaMDS(comm=occurences[,species.cols], k=2, trymax=1)#, plot=T)
#save(mds, file="nMDS.RData.gz", compress="gzip")
load(file="nMDS.RData.gz")
print(paste("Stress",mds$stress))
plot(mds)

## sites
sites.df = cbind(occurences[,1:17], mds$points)
qplot(x=MDS1, y=MDS2, data=sites.df, color=HMPBodySite, size=I(3), shape=Sex) + 
  scale_color_brewer(name="Body Site", palette="Set1") + theme_bw()# + ggtitle("Sites ordination") #facet_grid(.~Sex) +
ggsave("site ordination.png", height=6, width=8)

## species
### all
otus = as.numeric(lapply(names(occurences)[species.cols], function(x) {str_sub(x,2)}))
species.df = cbind(lookup[otus,], data.frame(mds$species))
species.df$frequency = colSums(occurences[,species.cols])/sum(occurences[,species.cols])
q1 = qplot(x=MDS1, y=MDS2, data=species.df, alpha=frequency) + 
  theme_bw() + theme(legend.position="top") + scale_alpha_continuous(name="Frequency") # + ggtitle("Species ordination")
q1
ggsave("species ordination.png", q1)

### frequent phylums
phylum.df = ddply(species.df, .(phylum), summarise, frequency=sum(frequency))
freq.phylums = subset(phylum.df, frequency > 0.0025)$phylum
q2 = qplot(x=MDS1, y=MDS2, data=subset(species.df, phylum %in% freq.phylums), color=phylum) +
  scale_color_brewer(name="Phylum", palette="Set1") + theme_bw() + theme(legend.position="top") #+ ggtitle("Species ordination (frequent phylums)")
q2
ggsave("species ordination - frequent phylums.png", q2)

png("species ordination - 2 panels.png", units="in", width=7, height=14, res=300)
grid.arrange(q1, q2, nrow=2)
dev.off()

qplot(x=phylum, y=-log10(frequency), data=phylum.df, geom="bar", stat="identity") + 
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1))
ggsave("phylum frequencies.png")

# assembly rules
## nestedness
nest = nestedtemp(occurences[1:100,species.cols])
print(nest)
plot(nest)
plot(nest, kind="incid")
png("nestedness.png")
plot(nest)
plot(nest, kind="incid")
dev.off()

## co-occurence

bodysite.order = order(occurences$HMPBodySite, occurences$HMPBodySubsite)
occur.sub = occurences[bodysite.order,]
cor.mat=cor(t(occur.sub[,species.cols]), method="spearman")

png("cooccurence by bodysite.png", units="in", width=12, height=12, res=300)

image(x=1:ncol(cor.mat), y=1:ncol(cor.mat), z=cor.mat, axes=F, xlab="", ylab="")
site.ticks = c(79, 247, 1100, 2197, 2651)
oral.ticks = c(85,256,426,598,760,923,1100,1270,1440)
skin.ticks = c(76,240,405,573)
site.labels = lapply(occur.sub$HMPBodySite[site.ticks], function(x) {stringr::str_sub(string=x, end=1)})
oral.labels =  lapply(occur.sub$HMPBodySubsite[occur.sub$HMPBodySite=="Oral"][oral.ticks], function(x) {stringr::str_sub(string=x, end=3)})
skin.labels =  lapply(occur.sub$HMPBodySubsite[occur.sub$HMPBodySite=="Skin"][skin.ticks], function(x) {stringr::str_sub(string=x, end=3)})
skin.labels = c("LAF", "LRC", "RAF", "RAC")

axis(1, at=site.ticks, labels=site.labels, cex=0.75)
axis(4, at=site.ticks, labels=site.labels, cex=0.75)
axis(3, at=oral.ticks+336, labels=oral.labels, cex=0.75)
axis(3, at=skin.ticks+1865, labels=skin.labels, cex=0.75)

color.legend(xl=-150,xr=-50, yb=ncol(cor.mat)/4, yt=ncol(cor.mat)*3/4, 
             legend=seq(0,1), rect.col=heat.colors(256), cex=1, gradient='y')

dev.off()

