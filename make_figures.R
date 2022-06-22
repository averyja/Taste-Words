#source all the good stuff here
setwd("/Users/averyja/Documents/GitHub/Taste-Words")
source("source_me.R")

## ========================================================================================== ##
## Behavioral Data Figures ====
working_dir = "behavioral_data"
## ========================================================================================== ##
## Heatmaps for Figure 1
## Combine matrices from both days
mturk_data1 = read.csv(sprintf("%s/mturk_day%s_pval_mat.csv",working_dir,1), row.names = 1)
mturk_data2 = read.csv(sprintf("%s/mturk_day%s_pval_mat.csv",working_dir,2), row.names = 1)
beh_data = (mturk_data1 + mturk_data2)/2

display.brewer.pal(9,"YlOrRd")
mypalette=brewer.pal(9,"YlOrRd")
## ========================================================================================== ##
## Heatmap using ggplot

df=melt(as.matrix(beh_data[clean_names(order),clean_names(order)]), na.rm = T)

concept_frame = data.frame(clean_names(order),c(rep("Emotion",24),rep("Taste",14),rep("Color",13)))
names(concept_frame) = c("concept","group_label")
concept_frame$color = c(rep(r_adj_colors[1],24),rep(r_adj_colors[2],14),rep(r_adj_colors[3],13))

ggplot(df,aes(x=Var1,y=Var2,fill=value)) +
  geom_tile(color="black") + scale_fill_gradientn(colors = brewer.pal(9,"YlOrRd"))+
  theme(axis.title.y = element_blank()) + 
  theme(axis.text.y = element_text(size=10, face="bold", color = concept_frame$color)) + 
  theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())+
  labs(title="Semantic Similarity Matrix") + theme(plot.title = element_text(size=32, face="bold", hjust = 0.5)) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2), panel.spacing.x = unit(0,"line")) +
  theme(legend.text = element_text(face="bold")) + 
  theme(legend.title = element_blank()) + 
  labs(fill = "Connection\nStrength") + theme(legend.position = "none")
ggsave(sprintf("%s/Semantic Similarity Matrix ggplot.png",working_dir),width=8,height=8,dpi=300)

## ========================================================================================== ##
## ggplot version of 1C heatmaps
temo_heatmap=heatmap(as.matrix(beh_data[clean_names(taste_concepts),clean_names(emotion_concepts)]), col=mypalette)
temo_rows = clean_names(taste_concepts)[temo_heatmap$rowInd]
temo_cols = clean_names(emotion_concepts)[temo_heatmap$colInd]
temo_frame=melt(as.matrix(beh_data[temo_rows,temo_cols]), na.rm = T)

ggplot(temo_frame,aes(x=Var1,y=Var2,fill=value)) +
  geom_tile(color="black") + scale_fill_gradientn(colors = brewer.pal(9,"YlOrRd"))+
  theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank()) + 
  theme(axis.text.y = element_text(size=10, face="bold", color = r_adj_colors[1])) + 
  theme(axis.text.x = element_text(size=10, face="bold", color = r_adj_colors[2]))+
  labs(title="Taste-Emotion") + theme(plot.title = element_text(size=32, face="bold", hjust = 0.5)) + 
  # theme(panel.border = element_rect(colour = "black", fill=NA, size=2), panel.spacing.x = unit(0,"line")) +
  theme(legend.text = element_text(face="bold")) + 
  theme(legend.title = element_blank()) + 
  labs(fill = "Connection\nStrength") + theme(legend.position = "none")
ggsave(sprintf("%s/taste-emotion heatmap ggplot.png",working_dir),width=10,height=5,dpi=300)

## ========================================================================================== ##
cemo_heatmap=heatmap(as.matrix(beh_data[clean_names(color_concepts),clean_names(emotion_concepts)]), col=mypalette)
cemo_rows = clean_names(color_concepts)[cemo_heatmap$rowInd]
cemo_cols = clean_names(emotion_concepts)[cemo_heatmap$colInd]
cemo_frame=melt(as.matrix(beh_data[cemo_rows,cemo_cols]), na.rm = T)

ggplot(cemo_frame,aes(x=Var1,y=Var2,fill=value)) +
  geom_tile(color="black") + scale_fill_gradientn(colors = brewer.pal(9,"YlOrRd"))+
  theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank()) + 
  theme(axis.text.y = element_text(size=10, face="bold", color = r_adj_colors[1])) + 
  theme(axis.text.x = element_text(size=10, face="bold", color = r_adj_colors[3]))+
  labs(title="Color-Emotion") + theme(plot.title = element_text(size=32, face="bold", hjust = 0.5)) + 
  # theme(panel.border = element_rect(colour = "black", fill=NA, size=2), panel.spacing.x = unit(0,"line")) +
  theme(legend.text = element_text(face="bold")) + 
  theme(legend.title = element_blank()) + 
  labs(fill = "Connection\nStrength") + theme(legend.position = "none")
ggsave(sprintf("%s/color-emotion heatmap ggplot.png",working_dir),width=10,height=5,dpi=300)

## ========================================================================================== ##

# Read in behavioral p-values
df = read.csv(sprintf("%s/behavioral_data_table.csv",working_dir), row.names = 1) #read in if needed
df$group = factor(df$group, levels=groups)
df$group_label = factor(df$group_label, levels=group_labels)

test_frame = filter(df,group %in% groups[c(4,5)])

## ========================================================================================== ##
## Create Box-Violin plot for Figure 2c
start_height=0.8; bar_size = 0; text_size = 4; ypos = 1

ggplot(test_frame, aes(group_label, edge, fill = group_label)) + 
  scale_fill_manual(values=r_colors[4:5]) + 
  geom_violin() + geom_boxplot(fill = "white", alpha = 0.75)+
  # stat_summary(geom = "bar", fun = mean, position = "dodge", colour="black") + 
  # stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width = 0.2)+
  geom_hline(yintercept=0) +
  # geom_jitter()+
  theme(strip.text.x = element_text(color="black", size=14, face="bold")) +
  geom_signif(comparisons=list(c("Taste-Emotion", "Color-Emotion")), test = "t.test", map_signif_level = TRUE, 
              y = 0.80, tip_length = 0, vjust=0.2) +
  ggtitle("Semantic Similarity by Concept Category") + theme(plot.title = element_text(hjust = 0.5, color="black", size=12, face="bold.italic")) +
  theme(axis.text.x = element_text(color="black",size = 12,hjust=NULL, face="bold")) +
  xlab("Concept Set") + theme(axis.title.x = element_blank()) +
  ylab("Average Similarity Value") + theme(axis.title.y = element_text(hjust = 0.5, color="black", size=12, face="bold")) +
  labs(fill = "Concept Set")+ theme(legend.position="none")
ggsave(sprintf("%s/Beh_TE_vs_CE_violin_box.png",working_dir),width=4, height = 4.5,dpi=300)

## ========================================================================================== ##
## Behavioral Data Figures ====
working_dir = "behavioral_data"
val_order = c(emotion_concepts,taste_concepts[1:12],color_concepts)

data_table = read.csv(sprintf("behavioral_data/behavioral_data_table.csv", row.names = 1))
data_table$row = factor(data_table$row, levels = val_order)
data_table$col = factor(data_table$col, levels = val_order)
# create half-mat for new Figure 3

ggplot(data_table,aes(x=row,y=col,fill=edge)) +
  geom_tile(color="black") + scale_fill_gradientn(colors = brewer.pal(9,"YlOrRd"))+
  theme(axis.title.y = element_blank()) + 
  theme(axis.text.y = element_blank()) + 
  theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.spacing.x = unit(0,"line")) +
  theme(legend.text = element_text(face="bold")) + 
  theme(legend.title = element_blank()) + 
  labs(fill = "Connection\nStrength") + theme(legend.position = "none")
ggsave(sprintf("%s/half mat ggplot.png",working_dir),width=6,height=6,dpi=300)

# create valence half-mat for new Figure 3
ggplot(data_table,aes(x=row,y=col,fill=valence)) +
  geom_tile(color="black") + scale_fill_gradientn(colors = brewer.pal(9,"YlOrRd"))+
  theme(axis.title.y = element_blank()) + 
  theme(axis.text.y = element_blank()) + 
  theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.spacing.x = unit(0,"line")) +
  theme(legend.text = element_text(face="bold")) + 
  theme(legend.title = element_blank()) + 
  labs(fill = "Connection\nStrength") + theme(legend.position = "none")
ggsave(sprintf("%s/valence half mat ggplot.png",working_dir),width=6,height=6,dpi=300)

# create arousal half-mat for new Figure 3
ggplot(data_table,aes(x=row,y=col,fill=arousal)) +
  geom_tile(color="black") + scale_fill_gradientn(colors = brewer.pal(9,"YlOrRd"))+
  theme(axis.title.y = element_blank()) + 
  theme(axis.text.y = element_blank()) + 
  theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.spacing.x = unit(0,"line")) +
  theme(legend.text = element_text(face="bold")) + 
  theme(legend.title = element_blank()) + 
  labs(fill = "Connection\nStrength") + theme(legend.position = "none")
ggsave(sprintf("%s/arousal half mat ggplot.png",working_dir),width=6,height=6,dpi=300)

## ========================================================================================== ##
## Create Scatterplots for each concept group
## Figure 3b,c =====
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
  ggsave(sprintf("%s/%s_concepts_beh_data_scatterplot.pdf",working_dir,group),width=6,height=6,dpi=300)
}
