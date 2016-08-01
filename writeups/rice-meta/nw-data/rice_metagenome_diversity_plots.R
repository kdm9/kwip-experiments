
library(reshape2)
library(ggplot2)
library(rgl)
library("RColorBrewer")
library(tsne)

##reading the data
diversity<-read.delim("table_for_Diversity_plots.txt", header=T)
names(diversity)
plot(diversity$num_reads, diversity$num_kmers)

greenhouse<-subset(diversity, Dataset=="greenhouse")
field<-subset(diversity, Dataset=="field")
time_series<-subset(diversity, Dataset=="time_series")

#ggplot(greenhouse, aes(x=Site, y=Entropy, fill=Compartment)) + 
#  geom_boxplot(position=position_dodge(.8)) +
#  theme_classic()

#ggplot(greenhouse, aes(x=Site, y=Entropy*num_kmers/num_reads, fill=Compartment)) + 
#  geom_boxplot(position=position_dodge(.8)) +
#  theme_classic()


ggplot(greenhouse, aes(x=Compartment,y=Entropy*num_kmers/num_reads, group=Compartment)) + 
  geom_boxplot(aes(fill=Compartment), position=position_dodge(.8)) +
  facet_grid(. ~Site) +
  theme(axis.text.x=element_blank())

+
  title("Rice Root Metagenome, alpha-Diversity")

# 
# ggplot(greenhouse, aes(x=Site, y=num_reads, fill=Compartment)) + 
#   geom_boxplot(position=position_dodge(.8)) +
#   theme_classic()
# 
# ggplot(greenhouse, aes(x=Site, y=num_kmers, fill=Compartment)) + 
#   geom_boxplot(position=position_dodge(.8)) +
#   theme_classic()
# 
# ggplot(greenhouse, aes(x=Site, y=num_kmers/num_reads, fill=Compartment)) + 
#   geom_boxplot(position=position_dodge(.8)) +
#   theme_classic()
# 
# ggplot(greenhouse, aes(x=Site, y=WeightedSampleSum, fill=Compartment)) + 
#   geom_boxplot(position=position_dodge(.8)) +
#   theme_classic()
# 
# 
# 
# 
