####### 
# Figure 1
#######

### load(PANCAN.RData)

######## Boxplots including Figure 1A, Figure 1B

#Subsetting CCS score, cancer type, iCluster, COCA subtypes, Mutational loads

mutation.data_PANCAN <- (subset(pData(PANCAN), select = c(type,CCS_ct,CCS,CCS_intra,iCluster,iCluster.annot,COCA_k32_annot,COCA_K32.colors,COCA_k32,770:902)))
mutation.data_PANCAN$COCA_shortCODE <- sub(":.*","",x = mutation.data_PANCAN$COCA_k32_annot);mutation.data_PANCAN$CCS_ct <- log2(mutation.data_PANCAN$CCS_ct+1)


### Figure 1 (A) - Boxplot
# distinctive colors for the box plots of cancer types

scale_fillscolor_cancer <- function(...){
  ggplot2:::manual_scale(
    c('fill','color') ,
    values = setNames(c(cancer_names$cancer.colors), cancer_names$abbs), 
    ...
  )
}
### dashlines  on the boxplots CCS continues - catagorical (Low-inter-high)
CCS_ticks <- c(as.numeric(quantile(rescale(mutation.data_PANCAN$CCS_ct,newrange = c(0,1)),probs = seq(0,1,length.out = 100))[c(2,50,97)]))
### numbers (x.axis)
cancer.num <- data.frame(type=as.factor(names(table(na.omit(reorder(mutation.data_PANCAN$type,mutation.data_PANCAN$CCS_ct,median))))),
                  sum=as.numeric(table(sort(reorder(mutation.data_PANCAN$type,mutation.data_PANCAN$CCS_ct,median)))))

Fig1A.boxplot <-  ggplot((mutation.data_PANCAN), aes(x = reorder(type,CCS_ct,median), y = rescale(CCS_ct,newrange = c(0,1)), fill = reorder(type,CCS_ct,median),order=type)) +
  geom_boxplot(outlier.fill = NULL,outlier.size = NULL,show.legend = F,alpha=0.9,outlier.colour = NULL,outlier.alpha = 0,outlier.color = NULL,outlier.shape = NULL) + 
  geom_jitter(show.legend = F,size=0.8,aes(color=reorder(type,CCS_ct,median),fill=reorder(type,CCS_ct,median))) +
  stat_boxplot(geom='errorbar',width=0.5)+ 
  geom_boxplot(outlier.fill = NULL,outlier.size = 3,show.legend = F,alpha=0.9,outlier.colour = NULL,outlier.alpha = 0,outlier.color = NULL,outlier.shape = NULL)  + 
  scale_fillscolor_cancer() +
  scale_y_continuous(breaks = c(0.25,0.5,0.75,1),sec.axis = sec_axis(trans = (~.),name="",CCS_ticks,labels=c('Low','Inter','High')),name = expression(atop(bold("Cell-cycle score"))))+
  labs(x='Cancer types') +
  ggtitle(label = 'CCS score across cancer types') +
  geom_hline(aes(yintercept=quantile(rescale(newrange = c(0,1),mutation.data_PANCAN$CCS_ct),probs = seq(0,1,length.out = 100))[c(34)]), color="forestgreen", linetype="dashed",size=1)+
  geom_hline(aes(yintercept=quantile(rescale(newrange = c(0,1),mutation.data_PANCAN$CCS_ct),probs = seq(0,1,length.out = 100))[c(67)]), color="tomato", linetype="dashed",size=1)+
  scale_x_discrete(labels=paste(cancer.num$type," (",cancer.num$sum,")",sep="")) +
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

### Figure 1 (B) - Boxplot
# distinctive colors for the box plots of COCA subtypes

COCA.colors <- na.omit(unique(mutation.data_PANCAN[,8:9]))
scale_fillscolor_coca <- function(...){
  ggplot2:::manual_scale(
    c('fill','color') ,
    values = setNames(c(COCA.colors$COCA_K32.colors), as.factor(COCA.colors$COCA_k32)), 
    ...
  )
}
#### numbers (x.axis)
coca.num <- data.frame(coca=as.factor(names(table(na.omit(reorder(mutation.data_PANCAN$COCA_shortCODE,mutation.data_PANCAN$CCS_ct,median))))),
                  sum=as.numeric(table(sort(reorder(mutation.data_PANCAN$COCA_shortCODE,mutation.data_PANCAN$CCS_ct,median)))))


Fig1B.boxplot <-  ggplot(na.omit(mutation.data_PANCAN), aes(x = reorder(COCA_k32,CCS_ct,median), y = rescale(CCS_ct,newrange = c(0,1)), fill = reorder(COCA_k32,CCS_ct,median),order=COCA_k32)) +
  geom_boxplot(outlier.fill = NULL,outlier.size = 3,show.legend = F,alpha=0.9,outlier.colour = NULL,outlier.alpha = 0,outlier.color = NULL) + 
  geom_jitter(show.legend = F,size=0.8,aes(color=reorder(COCA_k32,CCS_ct,median),fill=reorder(COCA_k32,CCS_ct,median))) +
  stat_boxplot(geom='errorbar',width=0.5) +
  geom_boxplot(outlier.fill = NULL,outlier.size = 3,show.legend = F,alpha=0.9,outlier.colour = NULL,outlier.alpha = 0,outlier.color = NULL)  + 
  scale_fillscolor_coca()+
  scale_x_discrete(labels=paste(coca.num$coca," (",coca.num$sum,")",sep="")) +
  labs(x='COCA subtypes') +
  scale_y_continuous(breaks = c(0.25,0.5,0.75,1),sec.axis = sec_axis(trans = (~.),name="",CCS_ticks,labels=c('Low','Inter','High')),name = expression(atop(bold("Cell-cycle score"))))+
  ggtitle(label = 'CCS score across COCA subtypes') +
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



##### Figure 1C - heatmap


## subseting a PANCAN data for Cell-cycle genes

PANCAN_ccs <- PANCAN[,order(PANCAN$CCS)]
fData(PANCAN_ccs)$phase <- factor(fData(PANCAN_ccs)$phase,levels = c('G1','G1/S','S','G2','G2/M','M'))
PANCAN_ccs <- PANCAN_ccs[order(fData(PANCAN_ccs)$phase),]
PANCAN_ccs <- PANCAN_ccs[,!is.na(PANCAN$type)]
levels(PANCAN_ccs$COCA_k32_annot) <- as.factor(unique(CoCA_subtypes$COCA_k32_annot))

#### Colors 
hm_color = colorRampPalette(c("#0083FF", "black", "yellow"))(n = 299)
col_breaks <- c(seq(-6.5,-0.2, length = 100), seq(-0.1, 0.1, length = 100), seq(0.2, 6.5,length = 100)) 

## color for cell-cycle phases 
color.CC <-function(CC){if(CC =="G1") "blue"
  else if (CC=="G1/S") "purple"
  else if (CC=="S") "chartreuse3"
  else if (CC=="G2") "orange"
  else if (CC=="G2/M") "red"
  else if (CC=="M") "pink"}

color.tum <- function(tum) { as.character(factor(tum, levels = as.character(levels(PANCAN_ccs$type)),
                                                 labels = as.character(cancer_names$cancer.colors)))}
color.coca <- function(coca) { as.character(factor(coca, levels = as.character(levels(PANCAN_ccs$COCA_k32_annot)),
                                                   labels = as.character(unique(CoCA_subtypes$COCA_K32.colors))))}


### CCS colors
color.CCS_ccs <- (ifelse((PANCAN_ccs$CCS == '1'),'black',
                         ifelse((PANCAN_ccs$CCS == '2'),'gray',
                                ifelse((PANCAN_ccs$CCS == '3'),'gold','black'))))

### cell cycle phases colors
color.Cellcycle_ccs <- unlist(lapply(fData(PANCAN_ccs)$phase, color.CC))
### cancer type colors
color.cancer.type_ccs <- sapply(PANCAN_ccs$type,color.tum)
### coca group colors
color.coca_ccs <- sapply(PANCAN_ccs$COCA_k32_annot,color.coca)
# colors for cancer types, CCS and COCA groups
Col.matrix_ccs.coca <- cbind(color.CCS_ccs,color.cancer.type_ccs,color.coca_ccs)
colnames(Col.matrix_ccs.coca) <- c("Cell cycle score","Cancer types","COCA subtypes")

# Heatmap
Figure1C.heatmap <- heatmap.3((log2(as.matrix(Biobase::exprs(PANCAN_ccs) + 1))),
          key = TRUE,na.rm= TRUE,notecex=0.8,keysize = 1,
          side.height.fraction = 1.5,hclustfun = "ward.D2",
          distfun = function(x) as.dist((1-cor(t(x)))/2),
          labCol = FALSE, labRow = FALSE, Rowv = NA, Colv = NA,
          scale = 'row', dendrogram = 'none',  breaks = col_breaks,col=hm_color,
          density.info = 'none', trace = 'none',
          RowSideColors = t(as.matrix(color.Cellcycle_ccs)),
          ColSideColors = Col.matrix_ccs)
