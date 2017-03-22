
#weighted unifraq
library(tidyr)
library(vegan)
library(ggplot2)
wun <- read.table("weighted.gh.unifrac", header = T, row.names = 1)
gh <- read.table("greenhouse_metadata.txt", header = T, sep = "\t")
wun2 <- wun[match(gh$SampleID, row.names(wun)), match(gh$SampleID, colnames(wun))]
#w.pc <- capscale(as.dist(wun2) ~ 1)
w.pc <- cmdscale(as.dist(wun2))
w.axes <- cbind(gh, scores(w.pc)$sites)
w.axes %>% 
  ggplot(aes(MDS1, MDS2, color = Compartment, shape = Site)) +
  geom_point()

#unweighted unifraq
library(tidyr)
library(vegan)
library(ggplot2)
wun <- read.table("unweighted.gh.unifrac", header = T, row.names = 1)
gh <- read.table("greenhouse_metadata.txt", header = T, sep = "\t")
wun2 <- wun[match(gh$SampleID, row.names(wun)), match(gh$SampleID, colnames(wun))]
#w.pc <- capscale(as.dist(wun2) ~ 1)
w.pc <- cmdscale(as.dist(wun2))
w.axes <- cbind(gh, scores(w.pc)$sites)
w.axes %>% 
  ggplot(aes(-MDS1, MDS2, color = Compartment, shape = Site)) +
  geom_point(size=3)
