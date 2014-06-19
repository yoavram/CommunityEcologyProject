library(vegan)
library(ggplot2)
library(reshape2)
library(plyr)
library(scales)
library(stringr)

# in hmp1.v13.hq.otu.counts, remove colname collection and add a new col name to the end (276...)
# the result is .header

occurences = read.csv(file="hmp1.v13.hq.otu.counts.may1", header=T, sep="\t", quote="")
occurences = occurences[1:2799,1:27655]
print(dim(occurences))

# filter columns (species) - remove unclassified and uncounted
lookup = read.csv(file="hmp1.v13.hq.otu.lookup.split", sep="\t")
lookup.classified = subset(lookup, lookup$family!='unclassified')
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
ggsave("body site samples distribution.png")

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

## TODO check beta div with less weight on rare species (bray curtis for example)
## j and sor are similar, with different weights
## z and w  are similar, and are roghly the negative of sor
## z 
betadiver(help=T)
w = betadiver(occurences[, species.cols], "w")
w = betadisper(w, occurences$HMPBodySite)
png("beta diversity w - body site.png")
boxplot(w, main="Whitaker Beta Diversity")
dev.off()
#z = betadiver(occurences[, species.cols], "z")
#boxplot(betadisper(z, occurences$HMPBodySite), main="Z Beta Diversity")
j = betadiver(occurences[, species.cols], "j")
j = betadisper(j, occurences$HMPBodySite)
png("beta diversity j - body site.png")
boxplot(j, main="Jaccard Beta Diversity")
dev.off()
#sor = betadiver(occurences[, species.cols], "sor")
#boxplot(betadisper(sor, occurences$HMPBodySite), main="Sorrenson Beta Diversity")
r = betadiver(occurences[, species.cols], "r")
r = betadisper(r, occurences$HMPBodySite)
png("beta diversity r - body site.png")
boxplot(r, main="Routledge Beta Diversity")
dev.off()
#sim = betadiver(occurences[, species.cols], "sim")
#boxplot(betadisper(sim, occurences$HMPBodySite), main="Simpson Beta Diversity")


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
  scale_color_brewer(name="Body Site", palette="Set1") + theme_bw() + ggtitle("Sites ordination") #facet_grid(.~Sex) +
ggsave("site ordination.png")

## species
### all
otus = as.numeric(lapply(names(occurences)[species.cols], function(x) {str_sub(x,2)}))
species.df = cbind(lookup[otus,], data.frame(mds$species))
species.df$frequency = colSums(occurences[,species.cols])/sum(occurences[,species.cols])
qplot(x=MDS1, y=MDS2, data=species.df, alpha=frequency) + 
  theme_bw() + ggtitle("Species ordination") + scale_alpha_continuous(name="Frequency")
ggsave("species ordination.png")

### frequent kingdoms
kingdom.df = ddply(species.df, .(kingdom), summarise, frequency=sum(frequency))
freq.kingdoms = subset(kingdom.df, frequency > 0.0025)$kingdom
qplot(x=MDS1, y=MDS2, data=subset(species.df, kingdom %in% freq.kingdoms), color=kingdom) +
  scale_color_brewer(name="Kingdom", palette="Set1") + theme_bw() + ggtitle("Species ordination (frequent kingdoms)")
ggsave("species ordination - frequent kingdoms.png")

qplot(x=kingdom, y=-log10(frequency), data=kingdom.df, geom="bar", stat="identity") + 
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1))
ggsave("kingdom frequencies.png")

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
cor.mat = cor(t(occurences[bodysite.order, species.cols]), method="spearman")
cor.mat[upper.tri(cor.mat )] <- NA
print(dim(cor.mat))
# http://sebastianraschka.com/Articles/heatmaps_in_r.html
heatmap.2(cor.mat)
