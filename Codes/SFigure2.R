##### Supplemental figure 2
### use these variables from Figure1.R

# mutation.data_PANCAN


### Supplemental Figure 2 (A) - Boxplot
####
# distinctive colors for the box plots of COCA subtypes
scale_fillscolor_iclust <- function(...){
  ggplot2:::manual_scale(
    c('fill','color') ,
    values = setNames(c(iCluster$iclust.colors), iCluster$group), 
    ...
  )
}

#### numbers (x.axis)

iclust.num <- data.frame(iclust=as.factor(names(table(na.omit(reorder(mutation.data_PANCAN$iCluster,mutation.data_PANCAN$CCS_ct,median))))),
                  sum=as.numeric(table(sort(reorder(mutation.data_PANCAN$iCluster,mutation.data_PANCAN$CCS_ct,median)))))
iclust.num$iclust <- paste('C',iclust.num$iclust,sep = '')

SFig2A.boxplot <-  ggplot(na.omit(mutation.data_PANCAN), aes(x = reorder(iCluster,CCS_ct,median), y = rescale(CCS_ct,newrange = c(0,1)), fill = reorder(iCluster,CCS_ct,median),order=iCluster)) +
  geom_boxplot(outlier.fill = NULL,outlier.size = 3,show.legend = F,alpha=0.9,outlier.colour = NULL,outlier.alpha = 0,outlier.color = NULL) + 
  geom_jitter(show.legend = F,size=0.8,aes(color=reorder(iCluster,CCS_ct,median),fill=reorder(iCluster,CCS_ct,median))) +
  stat_boxplot(geom='errorbar',width=0.5) +
  geom_boxplot(outlier.fill = NULL,outlier.size = 3,show.legend = F,alpha=0.9,outlier.colour = NULL,outlier.alpha = 0,outlier.color = NULL)  + 
  scale_fillscolor_iclust() +
  scale_x_discrete(labels=paste(iclust.num$iclust," (",iclust.num$sum,")",sep="")) +
  labs(x='iCluster groups') +
  scale_y_continuous(breaks = c(0.25,0.5,0.75,1),sec.axis = sec_axis(trans = (~.),name="",CCS_ticks,labels=c('Low','Inter','High')),name = expression(atop(bold("Cell-cycle score"))))+
  
  ggtitle(label = 'CCS score across iCluster subgroups') +
  geom_hline(aes(yintercept=quantile(rescale(newrange = c(0,1),mutation.data_PANCAN$CCS_ct),probs = seq(0,1,length.out = 100))[c(34)]), color="forestgreen", linetype="dashed",size=1)+
  geom_hline(aes(yintercept=quantile(rescale(newrange = c(0,1),mutation.data_PANCAN$CCS_ct),probs = seq(0,1,length.out = 100))[c(67)]), color="tomato", linetype="dashed",size=1)+
  coord_cartesian(ylim=c(0,0.85)) +
  theme(legend.key = element_rect(fill='white'),
        strip.text.x =element_text(size = 24),
        axis.title.x = element_text(size = 24,face='bold'),
        axis.title.y = element_text(size = 24,face='bold'),
        plot.title = element_text(size = 29, face = "bold",hjust = c(0.5,0.5)),
        plot.subtitle = element_text(size = 26,hjust = c(0.5,0.5)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.line.y.right =  element_blank(), 
        axis.text.y.right=element_text(size=24,angle = 90,face = 'bold',colour = c('black','darkgray','gold4'),hjust = c(0,0.4,0.5)),
        axis.text.y.left=element_text(size=20,face = 'bold',colour = 'black'),
        axis.title.y.left=element_text(size=24,face = 'italic',colour = 'black'),
        axis.text.x=element_text(size=20,angle = 70,hjust = c(1,1)),axis.ticks.y = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.7,0.5),'cm'),
        legend.key.size = unit(1.5, "cm"),legend.title = element_text(size = 15,face='bold'))


######## Supplemental Figure 2 (B) - heatmap

### USE PANCAN_CCS , colors (generated in Figure1.R)

# function to color the iCluster groups
color.iclust <- function(iclust) { as.character(factor(iclust, levels = as.character(levels(PANCAN_ccs$iClust_code)),
                                                       labels = as.character(iCluster$iclust.colors)))}
### icluster group colors
color.iclust_ccs <- sapply(PANCAN_ccs$iClust_code,color.iclust)

# colors for cancer types, CCS and iCluster groups
Col.matrix_ccs <- cbind(color.CCS_ccs,color.cancer.type_ccs,color.iclust_ccs)
colnames(Col.matrix_ccs) <- c("Cell cycle score","Cancer types","iCluster")

# Heatmap
pdf(file= "~/Projects/PANCAN/graphs/Heatmaps/CCS_heatmap.PANCAN.sfigx.pdf",height = 8.72, width = 11.69, onefile = T, paper = "special")
SFigure2B.heatmap <- heatmap.3((log2(as.matrix(Biobase::exprs(PANCAN_ccs) + 1))),
          key = TRUE,na.rm= TRUE,notecex=0.8,keysize = 1,KeyValueName="Expression level",
          side.height.fraction = 1.5,hclustfun = "ward.D2",
          distfun = function(x) as.dist((1-cor(t(x)))/2),
          labCol = FALSE, labRow = FALSE, Rowv = NA, Colv = NA,
          scale = 'row', dendrogram = 'none',  breaks = col_breaks,col=hm_color,
          density.info = 'none', trace = 'none',
          RowSideColors = t(as.matrix(color.Cellcycle_ccs)),
          ColSideColors = Col.matrix_ccs.coca)

