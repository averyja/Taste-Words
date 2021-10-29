# CLICS Colexification Analyses for Avery et al. Taste-Words study
#source all the good stuff here
source("source_me.R")

## ========================================================================================== ##
## ========================================================================================== ##
working_dir = "graph_files/ETC_concepts"

family = "UniversalNet" # read in universal net to get concept names
full_graph <- read.graph(sprintf("%s/%s.gml",working_dir,family),format="gml")
concepts = vertex_attr(full_graph)$emonames # read in full list of concepts
adj_matrix = matrix(as_adj(full_graph,type="both",sparse = F),nrow=length(concepts),ncol=length(concepts), dimnames = list(concepts,concepts))
adj_matrix[adj_matrix > 0] = 1 #just in case, turn all non-zero values into 1
adj_matrix = adj_matrix[order,order]
diag(adj_matrix) = NA; adj_matrix[lower.tri(adj_matrix)] = NA

## ========================================================================================== ##
# create a df for edges and grab pair names
df = melt(as.matrix(adj_matrix), na.rm = TRUE)
names(df) = c("row","col","edge")
concept_pairs = paste(df$row,df$col,sep = ":")
row.names(df) = concept_pairs

df$index = 0

df[which(df$row %in% emotion_concepts&df$col %in% emotion_concepts),"index"] = 1
df[which(df$row %in% taste_concepts&df$col %in% taste_concepts),"index"] = 2
df[which(df$row %in% color_concepts&df$col %in% color_concepts),"index"] = 3
df[which(df$row %in% emotion_concepts&df$col %in% taste_concepts),"index"] = 4
df[which(df$row %in% emotion_concepts&df$col %in% color_concepts),"index"] = 5
df[which(df$row %in% taste_concepts&df$col %in% color_concepts),"index"] = 6

df$group = sapply(df$index,function(x) groups[x])
df$group = factor(df$group, levels=groups)
df$group_label = sapply(df$index,function(x) group_labels[x])
df$group_label = factor(df$group_label, levels=group_labels)

df$type = "within"; df[df$group %in% groups[4:6],"type"] = "between"

## ========================================================================================== ##
## ========================================================================================== ##
## ========================================================================================== ##
# list out family files in directory
gml_files = list.files(working_dir, pattern = "[a-z].gml") 
families = sub(".gml","",gml_files)
families = families[!families %in% c("FamilySimNet","UniversalNet")]

# write out edges from graphs
edge_file = sprintf("%s/ETC_family_edges.csv",working_dir)

partitions=array(0,c(nrow(df),length(families)),dimnames = list(c(1:nrow(df)),families))
row.names(partitions) = concept_pairs
# Grab all the edges from the adjacency matrices for each language family 
for (j in 1:length(families)){
  family = families[j]
  full_graph <- read.graph(sprintf("%s/%s.gml",working_dir,family),format="gml")
  edges <- get_all_edges(full_graph,order)
  partitions[,family] = edges # put family concept edges in array
}
drops = names(which(colSums(partitions) == 0)) # remove any families with zero edges b/t all concepts
partitions = partitions[,families[which(!families %in% drops)]]

all_edges = cbind(df,partitions)
write.csv(all_edges, file = edge_file, quote = TRUE, eol = "\n", na = "NA",  row.names = TRUE)

## ========================================================================================== ##
## Now run loop to create matrix subgraphs, graph partitions, and then compute pairwise JI across families.
index = "JI"; test = jaccard
# Create Big DF with all concept data for all groups
concept_data = list()

for (i in c(1:6)) {
  # Set names of concept pairs and read in concepts
  group = groups[[i]]; set = strsplit(group,"_")[[1]]
  print(sprintf("#===== Begin Set: %s ========#",group))

  index_file = sprintf("%s/%s_edge_%s_mat.csv",working_dir,group,index)
  dframe_file = sprintf("%s/%s_edge_%s_adjusted_dframe.csv",working_dir,group,index)

  set_parts = partitions[which(df$group == group),]
  drops = names(which(colSums(set_parts) == 0)) # remove any families with zero edges b/t all concepts
  set_parts = set_parts[,colnames(set_parts)[!colnames(set_parts) %in% drops]]
  
  index_mat = par_calc(set_parts,test) # now calculate similarity of family edges
  pval_mat = perm_test(set_parts,test,1000) # adjust calculations for chance
  
  dframe = make_dframe(index_mat,pval_mat,group,index) # write out to a dataframe for comparison
  
  write.csv(index_mat, file = index_file, quote = TRUE, eol = "\n", na = "NA", row.names = TRUE)
  write.csv(dframe, file = dframe_file, quote = TRUE, eol = "\n", na = "NA", row.names = TRUE)
  
  concept_data[[group]] = dframe
  
}
## ========================================================================================== ##
# Create Big DF with all concept data for all groups
# Now combine concept list into big DF
concept_data = do.call("rbind", concept_data) 
write.csv(concept_data, file = sprintf("%s/Colex_ETC_%s_data_full.csv",working_dir,index), row.names = TRUE)

# Use below if not starting from scratch
# concept_data = read.csv(sprintf("%s/Colex_ETC_%s_data_full.csv",working_dir,"JI"),row.names = 1)
concept_data$group = factor(concept_data$group, levels=groups)
concept_data$pair = as.factor(concept_data$pair)

concept_data$group_label = "emotion"
for(i in 1:6){
  concept_data[concept_data$group==groups[i],"group_label"] = group_labels[i]
}
concept_data$group_label = factor(concept_data$group_label, levels=group_labels)

concept_data$type = "within"; concept_data[concept_data$group %in% groups[4:6],"type"] = "between"
concept_data$type = as.factor(concept_data$type)

## ========================================================================================== ##
# Aggregate to get summary stats table - see source_me.R for function details
D = summary_stats(concept_data,'group_label',"JI")
# Now do same for adjusted ji
E = summary_stats(concept_data,'group_label',"JI.adj")
## Run paired t-tests for concept data - see source_me.R for function details
JI_tests = paired_ttests("group_label",concept_data,"JI")
JI.adj_tests = paired_ttests("group_label",concept_data,"JI.adj")

## ========================================================================================== ##
write.csv(JI_tests, file = sprintf("%s/Colex_ETC_%s_t-tests.csv",working_dir,index), row.names = TRUE)
write.csv(JI.adj_tests, file = sprintf("%s/Colex_ETC_%s.adj_t-tests.csv",working_dir,index), row.names = TRUE)

write.csv(D, file = sprintf("%s/Colex_ETC_%s_summary_stats.csv",working_dir,index), row.names = TRUE)
write.csv(E, file = sprintf("%s/Colex_ETC_%s.adj_summary_stats.csv",working_dir,index), row.names = TRUE)

## ========================================================================================== ##
## Run some statistical tests
ji_model = aov(JI ~ group, data = concept_data)
ji_adj_model = aov(JI.adj ~ group, data = concept_data)

summary(aov(JI.adj ~ type+group + Error(pair), data = concept_data))

#########################################
## Test for relationship with geographic distance

ji_model = aov(JI ~ group*distance, data = concept_data)
ji_adj_model = aov(JI.adj ~ group*distance, data = concept_data)

summary(ji_model)
summary(ji_adj_model)

## Now, regress out the effect of distance, create new JI.resid values from the residuals
## And test groups as before

ji_dist_model = lm(JI ~ distance, data = concept_data)

concept_data$JI.resid = range(ji_dist_model$residuals)

F = summary_stats(concept_data,'group_label',"JI.resid")
JI_resid_tests = paired_ttests("group_label",concept_data,"JI.resid")

