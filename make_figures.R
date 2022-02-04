#source all the good stuff here
source("source_me.R")

## ========================================================================================== ##
## Behavioral Data Figures ====
working_dir = "behavioral_data"
## ========================================================================================== ##
## First import database stats to create Concreteness, Valence, and Arousal Comparisons
## Figure S1 =====
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
## Figure 1d =====
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
## Figure 2b,c =====
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
