#source all the good stuff here
source("source_me.R")

## ========================================================================================== ##
## Behavioral Data Figures ====
working_dir = "behavioral_data"
## ========================================================================================== ##
## First import database stats to create Concreteness, Valence, and Arousal Comparisons
## Figure S1
con_frame = read.csv(sprintf("%s/concept_database_stats.csv",working_dir), row.names = 1)

values = c("valence","arousal","concreteness")
for(i in values){
  bar_plot = ggplot(con_frame, aes_string("group", i, fill = "group")) + scale_fill_manual(values=r_colors[1:3]) +
    stat_summary(geom = "bar", fun = mean, position = "dodge", colour="black") + 
    stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width = 0.2)+
    geom_hline(yintercept=0) +
    theme(strip.text.x = element_text(color="black", size=14, face="bold")) +
    geom_signif(comparisons=list(c("Taste", "Emotion"),c("Color", "Emotion"),c("Color", "Taste")), 
                test = "t.test", map_signif_level = TRUE, y_position = c(.8, .9, 1.0), tip_length = 0, vjust=0.2) +
    ggtitle(sprintf("%s by Concept Category",str_to_sentence(i))) + theme(plot.title = element_text(hjust = 0.5, color="black", size=20, face="bold.italic")) +
    theme(axis.text.x = element_text(color="black",size = 12,hjust=NULL, face="bold")) +
    xlab("Concept Set") + theme(axis.title.x = element_blank()) +
    ylab(sprintf("Average %s - scaled",str_to_sentence(i))) + theme(axis.title.y = element_text(hjust = 0.5, color="black", size=14, face="bold")) +
    theme(legend.position = "none")
  plot(bar_plot)
  ggsave(sprintf("%s/%s_comparisons.pdf",working_dir,str_to_sentence(i)),dpi=300)
}

## ========================================================================================== ##
## Combine matrices from both days
mturk_data1 = read.csv(sprintf("%s/mturk_day%s_data.csv",working_dir,1), row.names = 1)
mturk_data2 = read.csv(sprintf("%s/mturk_day%s_data.csv",working_dir,2), row.names = 1)
beh_data = (mturk_data1 + mturk_data2)/2

# Read in behavioral p-values
df = read.csv(sprintf("%s/behavioral_data_table.csv",working_dir), row.names = 1) #read in if needed
## ========================================================================================== ##
## Plot beh data, For Figure 1d
start_height=0.8; bar_size = 0; text_size = 4; ypos = 1

beh_plot = ggplot(df, aes(group_label, edge, fill = group_label)) + scale_fill_manual(values=r_colors[1:6]) +
  stat_summary(geom = "bar", fun = mean, position = "dodge", colour="black") + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width = 0.2)+
  geom_hline(yintercept=0) +
  theme(strip.text.x = element_text(color="black", size=14, face="bold")) +
  geom_signif(comparisons=list(c("Taste-Emotion", "Color-Emotion"),c("Color", "Emotion")), test = "t.test", map_signif_level = TRUE, y_position = c(0.30, 0.75), tip_length = 0, vjust=0.2) +
  ggtitle("Semantic Similarity by Concept Category") + theme(plot.title = element_text(hjust = 0.5, color="black", size=20, face="bold.italic")) +
  theme(axis.text.x = element_text(size = 8,hjust=NULL, face="bold")) +
  xlab("Concept Set") + theme(axis.title.x = element_text(hjust = 0.5, color="black", size=14, face="bold")) +
  ylab("Average P-value") + theme(axis.title.y = element_text(hjust = 0.5, color="black", size=14, face="bold")) +
  labs(fill = "Concept Set")
plot(beh_plot)
ggsave(sprintf("%s/Behavioral_Data.pdf",working_dir),dpi=300)

test_frame = filter(df,group %in% groups[c(4,5)])

beh_plot = ggplot(test_frame, aes(group_label, edge, fill = group_label)) + scale_fill_manual(values=r_colors[4:5]) +
  stat_summary(geom = "bar", fun = mean, position = "dodge", colour="black") + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width = 0.2)+
  geom_hline(yintercept=0) +
  theme(strip.text.x = element_text(color="black", size=14, face="bold")) +
  geom_signif(comparisons=list(c("Taste-Emotion", "Color-Emotion")), test = "t.test", map_signif_level = TRUE, y_position = c(0.30), tip_length = 0, vjust=0.2) +
  ggtitle("Semantic Similarity by Concept Category") + theme(plot.title = element_text(hjust = 0.5, color="black", size=20, face="bold.italic")) +
  theme(axis.text.x = element_text(color="black",size = 12,hjust=NULL, face="bold")) +
  xlab("Concept Set") + theme(axis.title.x = element_blank()) +
  ylab("Average P-value") + theme(axis.title.y = element_text(hjust = 0.5, color="black", size=14, face="bold")) +
  labs(fill = "Concept Set")
plot(beh_plot)
ggsave(sprintf("%s/Beh_TE_vs_CE.pdf",working_dir),dpi=300)

## ========================================================================================== ##
## Create Scatterplots for each concept group
## Used for Figure 2b,c
for (i in 1:6){
  group = groups[i]; set = strsplit(group,"_")[[1]]
  group_label = group_labels[i]
  
  group_frame = df[df$group == group,]
  
  uni_glm = data.frame(group_frame$edge,group_frame$val,group_frame$ar); names(uni_glm) = c("edge_vec","valence","arousal")
  uni_model = lm(edge_vec ~ valence + arousal, data = uni_glm)
  uni_sum = summary(uni_model)
  uni_glm$model_fit = range(uni_model$fitted.values)
  R = cor.test(uni_glm$edge_vec,uni_glm$model_fit)
  #Run permutation test on uni_glm to get permutation-tested p-value for model r-squared
  model_pval = perm.lm.test(uni_glm,1000)
  
  val_stats = round(uni_sum$coefficients[2,c(1,4)], digits = 3) # grab coefficients
  ar_stats = round(uni_sum$coefficients[3,c(1,4)], digits = 3) # grab coefficients
  val_stats[3] = sig_label(val_stats[2]); ar_stats[3] = sig_label(ar_stats[2])
  
  # Plot a scatterplot of glm fit of valence + arousal vs. scaled edge strength
  scatterplot = ggplot(uni_glm,aes(x=model_fit)) + 
    geom_point(aes(y=edge_vec,size=edge_vec),color=r_colors[i],alpha=.7) + geom_point(aes(y=edge_vec,size=edge_vec),shape = 1,colour = "black") +
    geom_smooth(aes(y=edge_vec),method="lm",se=T,size=.5,fill="black",color='black') +
    scale_x_continuous(expand=c(0.02,0)) +
    scale_y_continuous(expand=c(0.02,0)) +
    scale_shape_manual(values=c(18, 20)) +
    theme(legend.position="top")+
    coord_cartesian(xlim = c(0,1),ylim = c(-.20,1))+
    xlab("Scaled Model Fit") + theme(axis.title.x = element_text(hjust = 0.5, color="black", size=18, face="bold")) +
    ylab("Behavioral p-values") + theme(axis.title.y = element_text(hjust = 0.5, color="black", size=18, face="bold")) +
    annotate("text", label = sprintf("%s",group_label), size = 8, x = 0.20, y = 0.95) +
    annotate("text", label = sprintf("Model: R^2 = %s, p = %s; %s",round(uni_sum$r.squared, digits = 2),round(model_pval, digits = 3), sig_label(model_pval)), size = 6, x = 0.5, y = -0.02) +
    annotate("text", label = sprintf("Valence: Beta = %s, p = %s; %s",val_stats[1],val_stats[2], val_stats[3]), size = 6, x = 0.5, y = -0.09) +
    annotate("text", label = sprintf("Arousal: Beta = %s, p = %s; %s",ar_stats[1],ar_stats[2], ar_stats[3]), size = 6, x = 0.5, y = -0.17) +
    theme_classic() + theme(axis.title.x=element_blank(),axis.title.y=element_blank()) +
    theme(legend.position="none",plot.margin=unit(c(0,0,0,0),"points"))
  
  plot(scatterplot)
  ggsave(sprintf("%s/%s_concepts_beh_data_scatterplot.pdf",working_dir,group),dpi=300)
}

## ========================================================================================== ##
## Colex Data Figures ====
## Make Colex graphs and plots
working_dir = "graph_files/ETC_concepts"

## ========================================================================================== ##
## create circle network graph plots for all families
## In Figure 3a,b
# list out family files in directory
gml_files = list.files(working_dir, pattern = "[a-z].gml") 
families = sub(".gml","",gml_files)
families = families[!families %in% c("FamilySimNet")]

picdir = sprintf("%s/Bin_Graphs/",working_dir)
dir.create(picdir)

# Plot all the graphs to a picture directory
for (j in 1:length(families)){
  family = families[j]
  full_graph <- read.graph(sprintf("%s/%s.gml",working_dir,family),format="gml")
  net = graph_prep(full_graph,family,order)
  pdf(sprintf("%s/%s.pdf",picdir,family),10,10)
  graph_plot(net,family)
  dev.off()
}

for (i in 1:6){
  group = groups[i]
  picdir = sprintf("%s/Bin_Graphs/%s_graphs",working_dir,group)
  dir.create(picdir)
  # Plot all the graphs to a picture directory
  for (j in 1:length(families)){
    family = families[j]
    full_graph <- read.graph(sprintf("%s/%s.gml",working_dir,family),format="gml")
    net = graph_prep(full_graph,family,order)
    net = delete_edges(net, which(E(net)$index != i))
    pdf(sprintf("%s/%s.pdf",picdir,family),10,10)
    graph_plot(net,family)
    dev.off()
  }
}
## ========================================================================================== ##
## ========================================================================================== ##
## Read in Concept Data csv, run stats, and generate bar plots for Colex data
## In Figure 3d,e and supplemental figures
index = "JI"; test = jaccard
concept_data = read.csv(sprintf("%s/Colex_ETC_%s_data_full.csv",working_dir,"JI"),row.names = 1)
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
## ========================================================================================== ##
## Create Bar plots, use annotation from paired t-tests for significance labels
start_height=max(E$max_height); bar_size = 0; text_size = 4; ypos = 1

bar_plot = ggplot(concept_data, aes(group_label, JI, fill = group_label)) + scale_fill_manual(values=r_colors[1:6]) +
  stat_summary(geom = "bar", fun = mean, position = "dodge", colour="black") + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width = 0.2)+
  geom_hline(yintercept=0) +
  theme(strip.text.x = element_text(color="black", size=14, face="bold")) +
  geom_signif(comparisons=list(c("Taste-Emotion", "Emotion"),c("Color-Emotion", "Emotion"),c("Color-Emotion", "Taste-Emotion")),
              test = "t.test", annotations = c(JI_tests[3,"label"],JI_tests[4,"label"],JI_tests[13,"label"]), 
              y_position = c(start_height+0.02, start_height+.04, 0.25), tip_length = 0, vjust=0.2) +
  ggtitle(sprintf("%s comparisons across concept sets",index)) + 
  theme(plot.title = element_text(hjust = 0.5, color="black", size=16, face="bold.italic")) +
  theme(axis.text.x = element_text(color="black",size = 12,hjust=NULL, face="bold")) +
  xlab("Concept Set") + theme(axis.title.x = element_blank()) +
  ylab(sprintf("Average %s",index)) + theme(axis.title.y = element_text(hjust = 0.5, color="black", size=14, face="bold")) +
  theme(legend.position = "none")

# Print it
plot(bar_plot)
ggsave(sprintf("%s/Average_%s_values_per_concept_group.pdf",working_dir,index), dpi = 300)
## ========================================================================================== ##
## Create Bar plots of perm-adjusted values, use annotation from paired t-tests for significance labels
start_height=max(E$max_height); bar_size = 0; text_size = 4; ypos = 1

bar_plot = ggplot(concept_data, aes(group_label, JI.adj, fill = group_label)) + scale_fill_manual(values=r_adj_colors[1:6]) +
  stat_summary(geom = "bar", fun = mean, position = "dodge", colour="black") + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width = 0.2)+
  geom_hline(yintercept=0) +
  theme(strip.text.x = element_text(color="black", size=14, face="bold")) +
  geom_signif(comparisons=list(c("Taste-Emotion", "Emotion"),c("Color-Emotion", "Emotion"),c("Color-Emotion", "Taste-Emotion")),test = "t.test", annotations = c(JI.adj_tests[3,"label"],JI.adj_tests[4,"label"],JI.adj_tests[13,"label"]), y_position = c(start_height+0.02, start_height+.04, 0.25), tip_length = 0, vjust=0.2) +
  ggtitle(sprintf("Adjusted %s comparisons",index)) + theme(plot.title = element_text(hjust = 0.5, color="black", size=16, face="bold.italic")) +
  theme(axis.text.x = element_text(color="black",size = 12,hjust=NULL, face="bold")) +
  xlab("Concept Set") + theme(axis.title.x = element_blank()) +
  ylab(sprintf("Average %s - Perm Adjusted",index)) + theme(axis.title.y = element_text(hjust = 0.5, color="black", size=14, face="bold")) +
  theme(legend.position = "none")

# Print it
plot(bar_plot)
ggsave(sprintf("%s/Average_%s_values_per_perm_adjusted.pdf",working_dir,index), dpi = 300)

## ========================================================================================== ##
## Plot panel for within concept data

within_data = filter(concept_data,group %in% groups[1:3])

bar_plot = ggplot(within_data, aes(group_label, JI, fill = group_label)) + scale_fill_manual(values=r_colors[1:6]) +
  stat_summary(geom = "bar", fun = mean, position = "dodge", colour="black") + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width = 0.2)+
  geom_hline(yintercept=0) +
  theme(strip.text.x = element_text(color="black", size=14, face="bold")) +
  geom_signif(comparisons=list(c("Color", "Emotion"),c("Taste", "Emotion")),test = "t.test", annotations = c(JI_tests[2,"label"],JI_tests[1,"label"]), y_position = c(0.25, 0.2), tip_length = 0, vjust=0.2) +
  ggtitle(sprintf("%s comparisons across concept sets",index)) + theme(plot.title = element_text(hjust = 0.5, color="black", size=16, face="bold.italic")) +
  theme(axis.text.x = element_text(color="black",size = 12,hjust=NULL, face="bold")) +
  xlab("Concept Set") + theme(axis.title.x = element_blank()) +
  ylab(sprintf("Average %s",index)) + theme(axis.title.y = element_text(hjust = 0.5, color="black", size=14, face="bold")) +
  labs(fill = "Concept Set") + theme(legend.text= element_text(color="black", size=12, face="bold"))
plot(bar_plot)
ggsave(sprintf("%s/Average_%s_values_Within_Set_Concepts.pdf",working_dir,index),dpi=300)

## Plot panel for Emo-related concept data
emo_data = filter(concept_data,group %in% groups[c(1,4,5)])

bar_plot = ggplot(emo_data, aes(group_label, JI, fill = group_label)) + scale_fill_manual(values=r_colors[c(1,4,5)]) +
  stat_summary(geom = "bar", fun = mean, position = "dodge", colour="black") + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width = 0.2)+
  geom_hline(yintercept=0) +
  theme(strip.text.x = element_text(color="black", size=14, face="bold")) +
  geom_signif(comparisons=list(c("Taste-Emotion", "Emotion"),c("Color-Emotion", "Emotion"),c("Color-Emotion", "Taste-Emotion")),
              test = "t.test", annotations = c(JI_tests[3,"label"],JI_tests[4,"label"],JI_tests[13,"label"]), 
              y_position = c(0.15, 0.20, 0.25), tip_length = 0, vjust=0.2) +
  ggtitle(sprintf("%s comparisons across concept sets",index)) + 
  theme(plot.title = element_text(hjust = 0.5, color="black", size=16, face="bold.italic")) +
  theme(axis.text.x = element_text(color="black",size = 12,hjust=NULL, face="bold")) +
  xlab("Concept Set") + theme(axis.title.x = element_blank()) +
  ylab(sprintf("Average %s",index)) + theme(axis.title.y = element_text(hjust = 0.5, color="black", size=14, face="bold")) +
  labs(fill = "Concept Set") + theme(legend.text= element_text(color="black", size=12, face="bold"))
plot(bar_plot)
ggsave(sprintf("%s/Average_%s_values_Emo_Set_Concepts.pdf",working_dir,index),dpi=300)

## ========================================================================================== ##
### Make Taste-Emotion JI heatmap figure
### Figure 3c
i = 4
group = groups[[i]];
temo_data = concept_data[concept_data$group_label == "Taste-Emotion",]

# Heatmap 
ggplot(data = temo_data, aes(Fam1, Fam2, fill = JI))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "yellow", mid = "red", 
                       limit = c(0,1), space = "Lab", 
                       name="Jaccard\nIndex") +
  theme_minimal()+ 
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size = 10,hjust=NULL, face="bold")) +
  coord_fixed()
ggsave(sprintf("%s/TE_JI_matrix_heatmap.pdf",working_dir), dpi = 300)

## ========================================================================================== ##
## Plot figures for relationship b/t JI values and distance
## If desired

for (i in 1:6){
  group = groups[i]
  group_label = group_labels[i]
  #dframe = concept_data[concept_data$group == group,]
  dframe = filter(concept_data,group == groups[i])
  R=cor.test(dframe$JI,dframe$distance)
  R.adj=cor.test(dframe$JI.adj,dframe$distance)
  
  #Scatter
  scatterplot = ggplot(dframe,aes(x=distance)) + 
    geom_point(aes(y=JI,size=JI),color=r_colors[i],alpha=.7) + geom_point(aes(y=JI,size=JI),shape = 1,colour = "black") +
    geom_smooth(aes(y=JI),method="lm",se=T,size=.5,fill=r_colors[i],color='black') +
    scale_x_continuous(expand=c(0.02,0)) +
    scale_y_continuous(expand=c(0.02,0)) +
    scale_shape_manual(values=c(18, 20)) +
    theme(legend.position="top")+
    coord_cartesian(xlim = c(1000,19000),ylim = c(-.20,1.05))+
    annotate("text", label = sprintf("%s: JI vs. Geo Distance: r = %s, p = %s; %s",group_label,round(R$estimate, digits = 2),round(R$p.value, digits = 2), sig_label(R$p.value)), size = 6, x = 10000, y = -0.05) +
    theme_classic() +
    theme(legend.position="none",plot.margin=unit(c(0,0,0,0),"points"))
  
  plot(scatterplot)
  ggsave(sprintf("%s/%s_%s_vs_Distance_scatterplot.pdf",working_dir,group,index), dpi = 300)
}

