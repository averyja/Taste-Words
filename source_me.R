library(igraph)
library(gtools)
library(aricode)
library(rlist)
library(reshape2)
library(ggplot2)
library(ggsignif)
library(plyr)
library(dplyr)
library(stringr)
library(RColorBrewer)

## ========================================================================================== ##
## VARIABLES
## ========================================================================================== ##
groups = c("emotion","taste","color","emotion_taste","emotion_color","taste_color")
group_labels = c("Emotion","Taste","Color","Taste-Emotion","Color-Emotion","Taste-Color")

taste_concepts = read.csv(sprintf("concept_lists/%s_concepts.txt","taste"),row.names = 1, header=TRUE,sep="\t"); taste_concepts=as.character(taste_concepts$Concept);
color_concepts = read.csv(sprintf("concept_lists/%s_concepts.txt","color"),row.names = 1, header=TRUE,sep="\t"); color_concepts=as.character(color_concepts$Concept); 
emotion_concepts = read.csv(sprintf("concept_lists/%s_concepts.txt","emotion"),row.names = 1, header=TRUE,sep="\t"); emotion_concepts=as.character(emotion_concepts$Concept);

order = c(emotion_concepts,taste_concepts,color_concepts)


r_colors = c(brewer.pal(9,"Reds")[6],"gold",brewer.pal(9,"Blues")[5],brewer.pal(9,"Oranges")[5],brewer.pal(9,"Purples")[5],brewer.pal(9,"Greens")[6])
r_adj_colors = c(brewer.pal(9,"Reds")[8],"goldenrod4",brewer.pal(9,"Blues")[8],brewer.pal(9,"Oranges")[7],brewer.pal(9,"Purples")[8],brewer.pal(9,"Greens")[8])

color_community<-c("red","yellow","blue","orange","purple","green",'#5d8aa8',
                   '#f0f8ff','#e32636','#efdecd','#e52b50','#ffbf00','#ff033e','#9966cc',
                   '#a4c639','#f2f3f4','#cd9575','#915c83','#faebd7','#008000','#8db600',
                   '#fbceb1','#00ffff','#7fffd4','#4b5320','#e9d66b','#b2beb5','#87a96b',
                   '#ff9966','#a52a2a','#a67b5b','#fad6a5')


## ========================================================================================== ##
## FUNCTIONS
## ========================================================================================== ##
# Annotate p-values with significance labels
sig_label = function(p){
  ifelse(p>0.05,"NS",ifelse(p<0.001,"***",ifelse(p<0.01,"**","*")))
}
## ========================================================================================== ##
# remove certain characters in concept names
clean_names <- function(x) {
  new_names <- str_remove(x,pattern="\\s\\(.*")
  return(new_names)
}
## ========================================================================================== ##
## Function to change the range of an input vector or matrix to between 0..1
range <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
## ========================================================================================== ##
## Transform edge strengths to normal distribution, between 0 and 1, set 0's to NA
trans_range <- function(x){range(ifelse(x > 0,log(x),NA), na.rm = TRUE)}
## ========================================================================================== ##
## Function for permutation tests for linear models
perm.lm.test <- function(uni_glm,iter) {
  names(uni_glm) = c("edge_vec","valence","arousal")
  true_index = (summary(lm(edge_vec ~ valence + arousal, data = uni_glm)))$r.squared
  perm_array <- replicate(iter, sample(uni_glm$edge_vec))
  perm.dist <- apply(perm_array, 2, function(x) (summary(lm(x ~ uni_glm$valence + uni_glm$arousal)))$r.squared)
  pval = mean(perm.dist > true_index)
  return(pval)
}
## ========================================================================================== ##
# Aggregate data to get summary stats table
summary_stats = function(dframe,factor,index){
  dframe = dframe[,c(factor,index)]; names(dframe) = c(factor,"idx")
  stats_table = ddply(dframe, c(factor), summarise, N = length(idx[!is.na(idx)]), mean = mean(idx, na.rm = T),
          sd = sd(idx, na.rm = T), se = sd/sqrt(N), max_height = mean + se, T = mean/se, 
          p.val = round(pt(-abs(T),df=N-1), digits = 3), label = sig_label(p.val))
  return(stats_table)
}
## ========================================================================================== ##
# Function for generating a dataframe full of corrected paired t-test results 
paired_ttests = function(factor,concept_data,index){
  all_factors = as.character(unique(concept_data[,factor]))
  df=t(combn(all_factors,2)) # list out all combinations of paired samples
  dset = list() #empty list for data
  for(i in 1:nrow(df)){
    a = df[i,1]; b = df[i,2] 
    aset = concept_data[concept_data[,factor] == a,]; rownames(aset) = aset$pair
    bset = concept_data[concept_data[,factor] == b,]; rownames(bset) = bset$pair
    set_pairs = intersect(rownames(aset),rownames(bset)) 
    test_set = data.frame(cbind(aset[set_pairs,c(index)],bset[set_pairs,c(index)]))
    names(test_set)=c(a,b)
    test=t.test(test_set[,a],test_set[,b],paired = T, na.action("na.omit"))
    stats=c(a,b,mean(test_set[,a],na.rm=T),mean(test_set[,b],na.rm=T),test$estimate,test$parameter,test$stderr,test$statistic,test$p.value)
    dset[[i]]=stats
  }
  dset = data.frame(do.call("rbind", dset)) # bind to dataframe
  names(dset) = c("A","B","A.mean","B.mean","Diff.AB","DoF","STDERR","T","p.val") # add names
  for(i in 3:9){dset[,i]=as.numeric(as.character(dset[,i]))} # convert back to numeric
  dset$p.adj=p.adjust(dset$p.val,method="fdr",n=length(dset$p.val))
  dset$label=sig_label(dset$p.adj)
  return(dset)
}
## ========================================================================================== ##
# Function for generating a dataframe full of corrected unpaired t-test results 
unpaired_ttests = function(factor,dframe,index){
  all_factors = as.character(unique(dframe[,factor]))
  df = t(combn(all_factors,2)) # list out all combinations of paired samples
  dframe = dframe[,c(factor,index)]; names(dframe) = c("group_label","edge")
  dset = list() #empty list for data
  for(i in 1:nrow(df)){
    a = df[i,1]; b = df[i,2] 
    test = with(dframe, t.test(edge[group_label == a],edge[group_label == b]))
    stats = c(a,b,mean(dframe[dframe$group_label == a,"edge"],na.rm=T),
              mean(dframe[dframe$group_label == b,"edge"],na.rm=T),0,test$parameter,test$stderr,test$statistic,test$p.value)
    dset[[i]]=stats
  }
  dset = data.frame(do.call("rbind", dset)) # bind to dataframe
  names(dset) = c("A","B","A.mean","B.mean","Diff.AB","DoF","STDERR","T","p.val") # add names
  for(i in 3:9){dset[,i]=as.numeric(as.character(dset[,i]))} # convert back to numeric
  dset$Diff.AB = dset$A.mean - dset$B.mean
  dset$p.adj=p.adjust(dset$p.val,method="fdr",n=length(dset$p.val))
  dset$label=sig_label(dset$p.adj)
  return(dset)
}
## ========================================================================================== ##
