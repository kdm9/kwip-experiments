setwd("/Users/u5264546/sandbox")

library(reshape2)
library(ggplot2)
library(rgl)
library("RColorBrewer")
library(tsne)

##reading the data
#y0<-read.delim("rice_metagenome.dist", header=T)
#kern0<-read.delim("rice_metagenome.kern", header=T)
#stats<-read.delim("rice_metagenome-stats.tab", header=T)
#metadata<-read.delim("rice_metagenome_runinfo.txt", header=T)

## FIELD
#stats<-read.delim("rice_metagenome_field_experiment-stats.tab", header=T)
#y0<-read.delim("rice_metagenome_field_experiment.dist", header=T)
#kern0<-read.delim("rice_metagenome_field_experiment.kern", header=T)
#metadata<-read.delim("rice_metagenome_field_metadata.txt", header=T)

##GRENNHOUSE
#stats<-read.delim("rice_metagenome_greenhouse-stats.tab", header=T)
y0<-read.delim("rice_metagenome_greenhouse.dist", header=T)
kern0<-read.delim("rice_metagenome_greenhouse.kern", header=T)
metadata<-read.delim("rice_metagenome_greenhouse_metadata.txt", header=T)

##TIMECOURSE
#stats<-read.delim("rice_metagenome_TCGHDavis-stats.tab", header=T)
y0<-read.delim("rice_metagenome_TCGHDavis.dist", header=T)
kern0<-read.delim("rice_metagenome_TCGHDavis.kern", header=T)
metadata<-read.delim("rice_metagenome_TC_and_DavisGH_metadata.txt", header=T)

#make matrix
y<-as.matrix(y0[, -1])
kern<-as.matrix(kern0[, -1])

#identity to self (experimental for confidence?)
self_table <- data.frame(colnames(kern),diag(kern))
colnames(self_table) <- c("run_id", "kerndiagonal")



#matching for labels. May look cumbersome, but is supposed to catch if samples are in different order in metadata file and matrix 
m = match(colnames(y), as.character(metadata$SampleID))
#n = paste(metadata$id[m], colnames(y), metadata$name[m], metadata$country[m],  sep="_")
#n = paste(metadata$id[m], metadata$name[m], metadata$country[m],  sep="_")
#biome<-metadata$env_biome_s[m]
#material<-metadata$env_material_s[m]
#location<-metadata$geo_loc_name_s[m]

##FIELD
Site<-metadata$Site[m]
Compartment<-metadata$Compartment[m]
Run<-metadata$Run[m]
Cultivation<-metadata$Cultivation[m]


##GREENHOUSE
Sample<-metadata$SampleID[m]
Site<-metadata$Site[m]
Cultivar<-metadata$Cultivar[m]
Compartment<-metadata$Compartment[m]
Run<-metadata$Run[m]
Description<-metadata$Description[m]

##subsetting the matrix for greenhouse experiment Figure 3 (exclude the soil samples!)
#selection=Compartment=="Rhizosphere"
selection=Cultivar!="Soil"
y<-y[selection, selection]
Sample<-Sample[selection]
Site<-Site[selection]
Cultivar<-Cultivar[selection]
Compartment<-Compartment[selection]
Run<-Run[selection]
Description<-Description[selection]


##TIMECOURSE
Days<-metadata$Days[m]
Cultivar<-metadata$Cultivar[m]
Compartment<-metadata$Compartment[m]
Rep<-metadata$Rep_Run[m]
Description<-metadata$Description[m]




##rename for simplicity (was originally becuase some tools needed the as.dist first)
#data<-as.dist(y)
data<-y

#distance matrix plot
qplot(x=Var1, y=Var2, data=melt(cor(y)), geom="tile", fill=value)

#ordered distance matrix plot
hcl= hclust(as.dist(y))
image(y[hcl$order,hcl$order],col=rainbow(80))

## do PCA
#pca <- prcomp(data, scale=F)
pca <- prcomp(data, scale=T)

#quadratic of sdev should give variance, compare to summary(pca)
#summary(pca)
percent_var=round(100*(pca$sdev)^2 / sum(pca$sdev^2),2) 
barplot(percent_var[1:20], names.arg=c(1:20), xlab="PCA", ylab="% variance explained", main="pca$sdev^2")
plot(cumsum((pca$sdev)^2) / sum(pca$sdev^2), type="b", ylim=c(0,1), main="rice root metagenome, cumulative components")
abline(h=0.8, col="blue")
#abline(v=7, col="blue")
abline(h=0.5, col="dark orange")
#abline(v=3, col="dark orange")



##histograms

########## PLOTS FOR FIELD

#extract the first 9 components into "melted" for easier plotting
#melted <- cbind(group, melt(pca$rotation[,1:9]))
melted <- cbind(Site, Compartment, Cultivation, melt(pca$rotation[,1:9]))
#melted <- cbind(group, melt(pca$x[,1:9]))

Palette1 <- brewer.pal(8, "Paired")
#splot the separation by PC components as histogram
ggplot(data=melted) +
  geom_bar(aes(x=Var1, y=value, fill=factor(Site)), stat="identity") +
  facet_wrap(~Var2) +
  scale_fill_manual(values=Palette1)

##2D PCA plots
scores <- data.frame(Compartment, Cultivation, pca$x[,1:9])

Palette1 <- brewer.pal(8, "Paired")
qplot(x=PC1, y=PC2, data=scores, colour=factor(Site)) +
  scale_colour_manual(values=Palette1) +
  geom_point(size=5) +
  theme(legend.position="right") + ggtitle("rice root metagenome Field")

Palette1 <- c( "blue", "green", "purple")
qplot(x=PC1, y=PC2, data=scores, colour=factor(Compartment)) +
  scale_colour_manual(values=Palette1) +
  geom_point(size=5) +
  theme(legend.position="right") + ggtitle("rice root metagenome Field")

Palette1 <- c( "brown", "orange")
qplot(x=PC3, y=PC2, data=scores, colour=factor(Cultivation)) +
  scale_colour_manual(values=Palette1) +
  geom_point(size=5) +
  theme(legend.position="right") + ggtitle("rice root metagenome Field")


ggplot(aes(x=PC1, y=PC2, colour=Site, shape=Cultivation), data=scores) +
  geom_point(size=5) +
  scale_colour_manual(values=Palette1) +
  theme(legend.position="right") + ggtitle("rice root metagenome Field") +
  theme_classic()
###################################################
###################################################

########## PLOTS FOR GREENHOUSE
#extract the first 9 components into "melted" for easier plotting
#melted <- cbind(group, melt(pca$rotation[,1:9]))
melted <- cbind(Site, Compartment, Cultivar, Run, melt(pca$rotation[,1:9]))

#melted <- cbind(group, melt(pca$x[,1:9]))

#splot the separation by PC components as histogram
Palette1 <- brewer.pal(4, "Paired")
ggplot(data=melted) +
  geom_bar(aes(x=Var1, y=value, fill=factor(Compartment)), stat="identity") +
  facet_wrap(~Var2) +
  scale_fill_manual(values=Palette1)

##2D PCA plots
scores <- data.frame(Sample,Compartment, Cultivar, Site, pca$x[,1:9])

# Palette1 <- brewer.pal(8, "Paired")
# qplot(x=PC1, y=PC2, data=scores, colour=factor(Site)) +
#   scale_colour_manual(values=Palette1) +
#   geom_point(size=5) +
#   theme(legend.position="right") + ggtitle("rice root metagenome Greenhouse")
# 
# Palette1 <- c( "blue", "green", "purple")
# qplot(x=PC1, y=PC2, data=scores, colour=factor(Compartment)) +
#   scale_colour_manual(values=Palette1) +
#   geom_point(size=5) +
#   theme(legend.position="right") + ggtitle("rice root metagenome Greenhouse")
# 
# Palette1 <- c( "brown", "orange")
# qplot(x=PC3, y=PC2, data=scores, colour=factor(Cultivation)) +
#   scale_colour_manual(values=Palette1) +
#   geom_point(size=5) +
#   theme(legend.position="right") + ggtitle("rice root metagenome Greenhouse")

Palette1 <- c("red", "blue", "green", "purple")
Palette1 <- c("#D10014", "#497CAF", "#5FB157", "#8D4EA0")

#PC1 does the compartment
ggplot(aes(x=PC1, y=-PC2, colour=Compartment, shape=Site), data=scores) +
  geom_point(size=5) +
  scale_colour_manual(values=Palette1) +
  theme(legend.position="right") + ggtitle("rice root metagenome Greenhouse") +
  geom_point(size=6) +
  theme_classic()

#1vs3 segregats by site
ggplot(aes(x=PC1, y=PC3, colour=Compartment, shape=Site), data=scores) +
  geom_point(size=5) +
  scale_colour_manual(values=Palette1) +
  theme(legend.position="right") + ggtitle("rice root metagenome Greenhouse") +
  theme_classic()

#for figure 3: separation by Cultivars (without the soil samples!)
#Palette1 <- c("red", "blue", "green", "purple", "darkgreen", "brown")
Palette1 <- rainbow(6)
ggplot(aes(x=PC1, y=PC2, colour=Cultivar, shape=Site), data=scores) +
  geom_point(size=5) +
  scale_colour_manual(values=Palette1) +
  theme(legend.position="right") + ggtitle("rice root metagenome Greenhouse") +
  theme_classic()

###################################################
###################################################

########## PLOTS FOR TIMECOURSE (published timecourse plots: PCoA of the time-series experiment and the greenhouse experiment subsetted to plants growing in Davis soil.
#extract the first 9 components into "melted" for easier plotting
#melted <- cbind(group, melt(pca$rotation[,1:9]))
melted <- cbind(Compartment, Cultivar, Days, Rep, melt(pca$rotation[,1:9]))
#melted <- cbind(group, melt(pca$x[,1:9]))


#splot the separation by PC components as histogram
Palette1 <- brewer.pal(4, "Paired")
ggplot(data=melted) +
  geom_bar(aes(x=Var1, y=value, fill=factor(Compartment)), stat="identity") +
  facet_wrap(~Var2) +
  scale_fill_manual(values=Palette1)

##2D PCA plots
scores <- data.frame(Compartment, Cultivar, Days, Rep, pca$x[,1:9])

# Palette1 <- brewer.pal(8, "Paired")
#  qplot(x=PC1, y=-PC2, data=scores, colour=factor(Days)) +
#    scale_colour_manual(values=Palette1) +
#    geom_point(size=5) +
#    theme(legend.position="right") + ggtitle("Rice Root Metagenome Time-series Analysis by Days")

##compare this to figure 6B
Palette1 <- c("red", "blue", "green", "purple")
qplot(x=-PC1, y=-PC2, data=scores, colour=Compartment) +
  scale_colour_manual(values=Palette1) +
  geom_point(size=6, shape=17) +  
  theme(legend.position="right") + ggtitle("Rice Root Metagenome Time-series Analysis by Compartment")

##compare this to figure 6C
ggplot(scores, aes(x=-PC1, y=-PC2, colour=Days)) +
  scale_colour_gradientn(colours=rainbow(5)) +
  geom_point(size=6, shape=17) +
  theme(legend.position="right") + 
  ggtitle("Rice Root Metagenome Time-series Analysis by Days")

# ggplot(scores, aes(x=-PC1, y=-PC2, shape=Compartment, colour=Days)) + geom_point() + 
#   scale_colour_gradientn(colours=rainbow(4)) +
#   geom_point(size=5) +
#   theme(legend.position="right") + 
#   ggtitle("rice root metagenome Timecourse")
# 
# ggplot(scores, aes(x=PC1, y=PC3, colour=Days)) + geom_point() + 
#   scale_colour_gradientn(colours=rainbow(4)) +
#   geom_point(size=5) +
#   theme(legend.position="right") + 
#   ggtitle("rice root metagenome Timecourse")
# 
# ggplot(scores, aes(x=PC2, y=PC3, colour=Days)) + geom_point() + 
#   scale_colour_gradientn(colours=rainbow(4)) +
#   geom_point(size=5) +
#   theme(legend.position="right") + 
#   ggtitle("rice root metagenome Timecourse")

# 
# Palette1 <- c( "blue", "green", "purple")
# qplot(x=PC1, y=PC2, data=scores, colour=factor(Compartment)) +
#   scale_colour_manual(values=Palette1) +
#   geom_point(size=5) +
#   theme(legend.position="right") + ggtitle("rice root metagenome Greenhouse")
# 
# Palette1 <- c( "brown", "orange")
# qplot(x=PC3, y=PC2, data=scores, colour=factor(Cultivation)) +
#   scale_colour_manual(values=Palette1) +
#   geom_point(size=5) +
#   theme(legend.position="right") + ggtitle("rice root metagenome Greenhouse")

Palette1 <- c("red", "blue", "green", "purple")
Palette1 <- rainbow(8)
#PC1 does the compartment
ggplot(aes(x=PC1, y=PC2, shape=Compartment, colour=as.character(Days)), data=scores) +
  geom_point(size=5) +
  scale_colour_manual(values=Palette1) +
  theme(legend.position="right") + ggtitle("rice root metagenome TIMECOURSE") +
  theme_classic()


###################################################
###################################################




#3D plot for comparisons (with library(rgl))
#setting the colors, not elegant, stole it from somewhere.
color_pallete_function <- colorRampPalette(colors = rainbow(4))
num_colors <- nlevels(Compartment)
subpopulation_colors <- color_pallete_function(num_colors)

#3D plot first 3 dimensions: (needs XQuartz installed)
plot3d(pca$x[,c(1,2,3)], size=10, col=subpopulation_colors[Compartment], type="p", main="rice root metagenome Greenhouse")

#you can plot different dimensions like so:
#plot3d(pca$x[,c(3,4,5)], size=10, col=subpopulation_colors[metadata$group][m], main="rice root metagenome")

#you may need to resize the window and then plot the legend (again)
legend3d("topleft", legend = levels(material), col = subpopulation_colors, pch = 19, cex = 1.5 )

## stick labels on.
#text3d(pca$x[,c(1,2,3)], texts=colnames(y),  adj=c(-0.1,0.1), cex=.7)
text3d(pca$x[,c(1,2,3)], texts=Run, adj=c(-1,1), cex=.7)
text3d(pca$x[,c(1,2,3)], texts=Cultivar, adj=c(-1,1), cex=.7)



##################### THIS IS FULLY EXPERIMENTAL ##################

#Mclust
pcat<-as.data.frame(pca$x[,1:7])
#fit <- Mclust(pcat[2:4], G=4, modelNames = "EEV")
fit <- Mclust(pcat[2:4], G=4)
plot3d(pcat[2:4], col = fit$classification, size=10) 
plot(pcat[2:3], cex = fit$uncertainty*10)


#t-sne
tsne_test<-tsne(data, k=4)
open3d()
plot3d(tsne_test[,c(1,2,3)], size=10, col=subpopulation_colors[material], main="rice root metagenome")
text3d(tsne_test[,c(1,2,3)], texts=colnames(y), adj=c(-0.2,0.2), cex=.7)

#################
#Contrained PCOA code from Edwards et al
# Compartment
library(vegan)
wuf.cap.comp <- capscale(as.dist(y) ~ Compartment + Condition(Site + Cultivar + Run), add = T)
wuf.cap.comp <- capscale(as.dist(y) ~ Compartment + Site + Cultivar + Run, add = T)

#wuf.cap.comp <- capscale(as.dist(y) ~ Compartment + Site + Condition(Cultivar + Run), add = T)
anova(wuf.cap.comp)
wuf.cap.comp.axes <- data.frame(scores(wuf.cap.comp)$sites)
#wuf.cap.comp.axes$Compartment <- factor(wuf.cap.comp.axes$Compartment, levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
percent_explained <- wuf.cap.comp$CCA$eig / sum(wuf.cap.comp$CCA$eig) * 100
comp.col <- c("#984EA3", "#4DAF4A", "#377EB8")
ggplot(wuf.cap.comp.axes, aes(x = CAP1, y = CAP2, color = Compartment, shape=Site)) +
  geom_vline(x = 0, alpha = 0.9) +
  geom_hline(y = 0, alpha = 0.9) +
  geom_point(size = 9, alpha = 0.75) +
  theme_classic() +
  labs(x = "Constrained PCo1", y = "Constrained PCo2") +
  scale_color_manual(values = comp.col) +
  theme(text = element_text(size = 30))

#Genotype
wuf.cap.geno <- capscale(as.dist(y) ~ Cultivar + Condition(Site + Compartment + Run), add = T)
anova(wuf.cap.geno)
wuf.cap.geno.axes <- data.frame(scores(wuf.cap.geno)$sites)
#wuf.cap.geno.axes$Cultivar <- factor(wuf.cap.geno.axes$Cultivar, levels = c("Glab_B", "Glab_E", "93-11", "IR50", "M104", "Nipponbare"))
percent_explained <- wuf.cap.geno$CCA$eig / sum(wuf.cap.geno$CCA$eig) * 100
cult.cols <- c("#006600", "#B2D1B2", "#FF9933", "#FFD6AD", "#660066", "#D1B2D1")
ggplot(wuf.cap.geno.axes, aes(x = CAP1, y = -CAP2, color = Cultivar)) +
  geom_vline(x = 0, alpha = 0.9) +
  geom_hline(y = 0, alpha = 0.9) +
  geom_point(size = 9, alpha = 0.9) +
  theme_classic() +
  labs(x = "Constrained PCo1 (38.7%)", y = "Constrained -PCo2 (26.6%)") +
  scale_color_manual(values = cult.cols) +
  theme(text = element_text(size = 30))



