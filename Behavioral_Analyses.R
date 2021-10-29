# Behavioral Analyses for Avery et al. Taste-Words study

#source all the good stuff here
source("source_me.R")

working_dir = "behavioral_data"
## ========================================================================================== ##
## Read in mturk data csv files and output behavioral similarity matrix

for (k in c(1,2)){
  infile = sprintf("%s/mturk_day%s_data.csv",working_dir,k) #set the input file
  outfile = sprintf("%s/mturk_day%s_pval_mat.csv",working_dir,k) #set the output file
  
  data = data.frame(read.csv(infile));
  word1 = data$word1; word2 = data$word2; word3 = data$word3; answer = data$answer
  
  # Assemble list of all concept words in dataset
  words = unique(c(levels(word1),levels(word2),levels(word3)))
  words = words[!words %in% c("PLUS","MINUS","EQUAL")]
  
  # Now generate p-value matrix using behavioral triplet data
  pval_mat=array(0,c(length(words),length(words)))
  row.names(pval_mat)=words; colnames(pval_mat)=words;
  Y = pval_mat;
  
  for (i in words){
    for (j in words){
      if (i == j){pval_mat[i,j] = 1; Y[i,j] = NA; next}
      # return indices of triplets containing both words
      triplets = which((word1 == i | word2 == i | word3 == i) & (word1 == j | word2 == j | word3 == j))
      answers = answer[triplets] # subset answer column by those indices
      nval = length(answers) # how many triples are in this set
      num_targets = sum(answers == j | answers == i) #find the number of times either i or j is in the list of answers
      pval_mat[i,j] = 1 - (num_targets)/nval # proportion we're looking for
      Y[i,j] = nval # store this value in a separate matrix
    }
  }

  print(pval_mat)
  # write resulting p-value matrix to a csv file
  write.csv(pval_mat, file = outfile, quote = TRUE, eol = "\n", na = "NA", row.names = TRUE)
}
## ========================================================================================== ##
## Combine matrices from both days
mturk_data1 = read.csv(sprintf("%s/mturk_day%s_data.csv",working_dir,1), row.names = 1)
mturk_data2 = read.csv(sprintf("%s/mturk_day%s_data.csv",working_dir,2), row.names = 1)
cor.test(mturk_data1[upper.tri(mturk_data1)],mturk_data2[upper.tri(mturk_data2)])
beh_data = (mturk_data1 + mturk_data2)/2

# sort the columns and rows of the matrix into the correct order
base_concepts = clean_names(c(emotion_concepts,taste_concepts,color_concepts))
beh_data = beh_data[base_concepts,base_concepts]
row.names(beh_data) = c(emotion_concepts,taste_concepts,color_concepts); colnames(beh_data) = row.names(beh_data)

## ========================================================================================== ##
## ========================================================================================== ##
## Now import valence and arousal ratings from database to generate V & A similarity matrices
full_war=data.frame(read.table("BRM-emot-submit.csv", header = TRUE, sep = ",", row.names = 2))
# Also import concreteness ratings for another analysis
BRM_ratings = read.csv("Concreteness_ratings_Brysbaert_et_al_BRM.csv",header = T, row.names = 1)

all_concepts = toupper(row.names(full_war))
full_war=data.frame(scale(full_war[,c(2,5)])) # subset to Valence, Arousal, and Dominance values only

## ========================================================================================== ##
# Check if the words we're using are in the full_war list
color_words = clean_names(color_concepts); color_words=gsub("GREY","GRAY",color_words,ignore.case = TRUE)
emotion_words = clean_names(emotion_concepts)
taste_words = clean_names(taste_concepts[1:12])

taste_words %in% all_concepts #All there except for 'cooked' and 'unripe'
color_words %in% all_concepts #All there, 'GREY' is 'GRAY'
emotion_words %in% all_concepts #All there

## ========================================================================================== ##
## Assemble a data frame for some later analyses of concept database stats 
con_words = c(emotion_words,taste_words,color_words)
con_group = c(rep("Emotion",24),rep("Taste",12),rep("Color",13))
val_order = c(emotion_concepts,taste_concepts[1:12],color_concepts)
 
war_stats=full_war[tolower(con_words),] #grab stats for all our words
row.names(war_stats) = val_order #rename rows to conform to overall scheme
BRM_stats = BRM_ratings[tolower(con_words),"Conc.M"] 
row.names(BRM_stats) = val_order #rename rows to conform to overall scheme

con_frame = data.frame(con_group,war_stats,BRM_stats)
names(con_frame) = c("group","valence","arousal","concreteness")
con_frame$group = factor(con_frame$group, levels=c("Emotion","Taste","Color"))
# Scale concept stats to range of 0..1, for ease of plotting
con_frame$valence = range(con_frame$valence);con_frame$arousal = range(con_frame$arousal);con_frame$concreteness = range(con_frame$concreteness)

write.csv(con_frame, file = sprintf("%s/concept_database_stats.csv",working_dir), quote = TRUE, eol = "\n", na = "NA",  row.names = TRUE)

## ========================================================================================== ##
# create distance matrices for correlation analyses, scaled to a range of 0..1
val_dist=data.matrix(range(dist(war_stats[1], method = "euclidean", diag = TRUE, upper = TRUE), na.rm =T))
ar_dist=data.matrix(range(dist(war_stats[2], method = "euclidean", diag = TRUE, upper = TRUE), na.rm =T))

# subtract 1 to transform to similarity matrices
val_sim = 1 - val_dist; ar_sim = 1 - ar_dist;
val_vec = val_sim[upper.tri(val_sim)]; ar_vec = ar_sim[upper.tri(ar_sim)]
## ========================================================================================== ##
## Create a dataframe in long-format for analyses and plotting
#rename and subset rows and columns
beh_mat = beh_data[val_order,val_order]
diag(beh_mat) = NA; beh_mat[lower.tri(beh_mat)] = NA
df = melt(as.matrix(beh_mat), na.rm = TRUE)
names(df) = c("row","col","edge")

df$index = 0

df[which(df$row %in% emotion_concepts&df$col %in% emotion_concepts),"index"] = 1
df[which(df$row %in% taste_concepts&df$col %in% taste_concepts),"index"] = 2
df[which(df$row %in% color_concepts&df$col %in% color_concepts),"index"] = 3
df[which(df$row %in% emotion_concepts&df$col %in% taste_concepts),"index"] = 4
df[which(df$row %in% emotion_concepts&df$col %in% color_concepts),"index"] = 5
df[which(df$row %in% taste_concepts&df$col %in% color_concepts),"index"] = 6

df$group = sapply(df$index,function(x) groups[x])
df$group_label = sapply(df$index,function(x) group_labels[x])
df$valence = val_vec; df$arousal = ar_vec;
df$type = "within"; df[df$group %in% groups[4:6],"type"] = "between"

## Write behavioral data to a file
write.csv(df, file = sprintf("%s/behavioral_data_table.csv",working_dir), quote = TRUE, eol = "\n", na = "NA",  row.names = TRUE)
df = read.csv(sprintf("%s/behavioral_data_table.csv",working_dir), row.names = 1) #read in if needed
df$group = factor(df$group, levels=groups)
df$group_label = factor(df$group_label, levels=group_labels)
## ========================================================================================== ##
# Aggregate to get summary stats table - see source_me.R for function details
edge_stats = summary_stats(df,'group_label',"edge")

## ========================================================================================== ##
## ========================================================================================== ##
# Some statistical tests
summary(aov(edge ~ group+type,data = df)) #overall effect of concept type

# Run all pairwise comparisons
beh_stats = unpaired_ttests("group_label",df,"edge")
write.csv(beh_stats, file = sprintf("%/beh_t-tests.csv",working_dir), quote = TRUE, eol = "\n", na = "NA",  row.names = TRUE)

## Test relationship with valence and arousal
summary(aov(edge ~ group*valence*arousal, data = df))
summary(lm(edge ~ group*valence, data = df[df$group==groups[4]|df$group==groups[5],]))
summary(aov(edge ~ group*valence*arousal, data = df[df$group==groups[4]|df$group==groups[5],]))
