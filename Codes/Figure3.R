### use these variables from Figure2.R

# CCS*.barplot;CCS*.gain.barplot and CCS*.loss.barplot, Gains

### Boxplots (Figure 3A:D)
### Figure 3A - top 50 genes

CCS.mut.boxplot <- cbind(c(t(sort(colSums(CCS1.barplot[,-c(1,ncol(CCS1.barplot))]),decreasing = T)[1:50])),
                          c(t(sort(colSums(CCS2.barplot[,-c(1,ncol(CCS2.barplot))]),decreasing = T)[1:50])),
                          c(t(sort(colSums(CCS3.barplot[,-c(1,ncol(CCS3.barplot))]),decreasing = T)[1:50])))
colnames(CCS.mut.boxplot) <- c("Low","Inter","High")
CCS.mut.boxplot_melted <-  melt(CCS.mut.boxplot)
colnames(CCS.mut.boxplot_melted) <- c('Genes','CCS','value')
CCS.mut.boxplot_melted$CCS <- factor(CCS.mut.boxplot_melted$CCS,levels = c("Low","Inter","High"))
CCS.mut.boxplot_melted <- cbind(CCS.mut.boxplot_melted,score='Inter-Cancer')

Fig3A.boxplot <- ggplot(melt((CCS.mut.boxplot_melted)), aes(x=CCS, y=log(value))) + 
  geom_boxplot(aes(col=CCS),outlier.fill = NULL,outlier.size = 3,show.legend = F,outlier.alpha = 1) +
  geom_jitter(show.legend = F,alpha= 0.9,size=3,aes(col=CCS,fill=CCS)) +
  stat_boxplot(geom='errorbar',width=0.3) +
  geom_boxplot(aes(fill=CCS),outlier.shape = NA,show.legend = F,outlier.colour = NULL,fatten=1) +
  scale_fill_manual(values=c("black","ivory3","gold3"),name="Cell Cycle Score",labels=c('Low','Inter','High')) + 
  scale_color_manual(values=c("black","ivory3","gold3"),name="Cell Cycle Score",labels=c('Low','Inter','High')) +
  scale_x_discrete(name="Cell Cycle Score") +
  scale_y_continuous(name="log (Number of mutations)",breaks = c(4,6,8),limits = c(3,10)) +
  stat_summary(geom = "crossbar", width=0.74, fatten=3, color="white", 
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))})+
  ggtitle(label = 'Oncogene mutation frequency in CCS subgroups') +
  theme(legend.key = element_rect(fill='white'),
        axis.title.x = element_text(size = 29, face = "bold"),
        axis.title.y = element_text(size = 29, face = "italic"),
        plot.title = element_text(size = 35, face = "bold",hjust = c(0.5,0.5)),
        plot.subtitle = element_text(size = 12.5,hjust = c(0.5,0.5)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x =element_text(size=34,face = 'bold',angle = 45,hjust = c(1,1)),
        axis.text.y =element_text(size=34,face = 'bold',hjust = c(1,1)),
        plot.margin = unit(c(0.5,0.5,0,0.5),'cm'),
        strip.background = element_rect(fill='gray100'),strip.text = element_text(face = 'italic',size = 34))

### TUKEYHSD test

Fig3A.boxplot.stats <- TukeyHSD(aov(CCS.mut.boxplot_melted$value~CCS.mut.boxplot_melted$CCS))

### Figure 3B - Chromosomal gains
CCS.gain.boxplot <- t(rbind.fill(as.data.frame(t(colSums(CCS1.gain.barplot[,-c(1,ncol(CCS1.gain.barplot))]))),
                                  as.data.frame(t(colSums(CCS2.gain.barplot[,-c(1,ncol(CCS2.gain.barplot))]))),
                                  as.data.frame(t(colSums(CCS3.gain.barplot[,-c(1,ncol(CCS1.gain.barplot))]))))); colnames(CCS.gain.boxplot) <- c("Low","Inter","High")

CCS.gain.boxplot_melted <-  melt(CCS.gain.boxplot);colnames(CCS.gain.boxplot_melted) <- c('Genes','CCS','value');CCS.gain.boxplot_melted$CCS <- factor(CCS.gain.boxplot_melted$CCS,levels = c("Low","Inter","High"))
CCS.gain.boxplot_melted <- cbind(CCS.gain.boxplot_melted,score='Inter-Cancer')

Fig3B.boxplot <- ggplot(melt((CCS.gain.boxplot_melted)), aes(x=CCS, y=log(value))) + 
  geom_boxplot(aes(col=CCS),outlier.fill = NULL,outlier.size = 3,show.legend = F,outlier.alpha = 1) +
  geom_jitter(show.legend = F,alpha= 0.9,size=3,aes(col=CCS,fill=CCS)) +
  stat_boxplot(geom='errorbar',width=0.2) +
  geom_boxplot(aes(fill=CCS),outlier.shape = NA,show.legend = F,outlier.colour = NULL,fatten=1) +
  scale_fill_manual(values=c("black","ivory3","gold3"),name="Cell Cycle Score",labels=c('Low','Inter','High')) + 
  scale_color_manual(values=c("black","ivory3","gold3"),name="Cell Cycle Score",labels=c('Low','Inter','High')) +
  scale_x_discrete(name="Cell Cycle Score") +
  scale_y_continuous(name="log (Number of arm-level gains)",breaks = c(4,6,8),limits = c(3,10)) +
  stat_summary(geom = "crossbar", width=0.74, fatten=3, color="white", 
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))})+
  ggtitle(label = 'Arm-level gains frequency in CCS subgroups') +
  theme(legend.key = element_rect(fill='white'),
        axis.title.x = element_text(size = 29, face = "bold"),
        axis.title.y = element_text(size = 29, face = "italic"),
        plot.title = element_text(size = 35, face = "bold",hjust = c(0.5,0.5)),
        plot.subtitle = element_text(size = 12.5,hjust = c(0.5,0.5)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x =element_text(size=34,face = 'bold',angle = 45,hjust = c(1,1)),
        axis.text.y =element_text(size=34,face = 'bold',hjust = c(1,1)),        
        plot.margin = unit(c(0.5,0.5,0,0.5),'cm'),
        strip.background = element_rect(fill='gray100'),strip.text = element_text(face = 'italic',size = 34))

Fig3B.boxplot.stats <- TukeyHSD(aov(CCS.gain.boxplot_melted$value~CCS.gain.boxplot_melted$CCS))


### Figure 3C - Chromosomal losss
CCS.loss.boxplot <- t(rbind.fill(as.data.frame(t(colSums(CCS1.loss.barplot[,-c(1,ncol(CCS1.loss.barplot))]))),
                                 as.data.frame(t(colSums(CCS2.loss.barplot[,-c(1,ncol(CCS2.loss.barplot))]))),
                                 as.data.frame(t(colSums(CCS3.loss.barplot[,-c(1,ncol(CCS1.loss.barplot))]))))); colnames(CCS.loss.boxplot) <- c("Low","Inter","High")

CCS.loss.boxplot_melted <-  melt(CCS.loss.boxplot);colnames(CCS.loss.boxplot_melted) <- c('Genes','CCS','value');CCS.loss.boxplot_melted$CCS <- factor(CCS.loss.boxplot_melted$CCS,levels = c("Low","Inter","High"))
CCS.loss.boxplot_melted <- cbind(CCS.loss.boxplot_melted,score='Inter-Cancer')

Fig3C.boxplot <- ggplot(melt((CCS.loss.boxplot_melted)), aes(x=CCS, y=log(value))) + 
  geom_boxplot(aes(col=CCS),outlier.fill = NULL,outlier.size = 3,show.legend = F,outlier.alpha = 1) +
  geom_jitter(show.legend = F,alpha= 0.9,size=3,aes(col=CCS,fill=CCS)) +
  stat_boxplot(geom='errorbar',width=0.2) +
  geom_boxplot(aes(fill=CCS),outlier.shape = NA,show.legend = F,outlier.colour = NULL,fatten=1) +
  scale_fill_manual(values=c("black","ivory3","gold3"),name="Cell Cycle Score",labels=c('Low','Inter','High')) + 
  scale_color_manual(values=c("black","ivory3","gold3"),name="Cell Cycle Score",labels=c('Low','Inter','High')) +
  scale_x_discrete(name="Cell Cycle Score") +
  scale_y_continuous(name="log (Number of arm-level losses)",breaks = c(4,6,8),limits = c(3,10)) +
  stat_summary(geom = "crossbar", width=0.74, fatten=3, color="white", 
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))})+
  ggtitle(label = 'Arm-level losses frequency in CCS subgroups') +
  theme(legend.key = element_rect(fill='white'),
        axis.title.x = element_text(size = 29, face = "bold"),
        axis.title.y = element_text(size = 29, face = "italic"),
        plot.title = element_text(size = 35, face = "bold",hjust = c(0.5,0.5)),
        plot.subtitle = element_text(size = 12.5,hjust = c(0.5,0.5)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x =element_text(size=34,face = 'bold',angle = 45,hjust = c(1,1)),
        axis.text.y =element_text(size=34,face = 'bold',hjust = c(1,1)),        
        plot.margin = unit(c(0.5,0.5,0,0.5),'cm'),
        strip.background = element_rect(fill='gray100'),strip.text = element_text(face = 'italic',size = 34))

Fig3C.boxplot.stats <- TukeyHSD(aov(CCS.loss.boxplot_melted$value~CCS.loss.boxplot_melted$CCS))

### Figure 3D - Aneuploidy score

CCS.Aneu.boxplot <- melt(Gains[,c(2,3)],id.vars = c("Aneuploidy Score"))

Fig3D.boxplot <-  ggplot(CCS.Aneu.boxplot, aes(x = value, y = `Aneuploidy Score`)) + 
  geom_boxplot(aes(col=value),outlier.fill = NULL,outlier.size = 0.0001,show.legend = F,outlier.alpha = 0) + 
  geom_jitter(show.legend = F,alpha= 0.7,size=0.5,aes(col=value,fill=value)) +
  stat_boxplot(geom='errorbar',width=0.3) +
  geom_boxplot(aes(fill=value),outlier.shape = NA,show.legend = F,outlier.colour = NULL,fatten=1) +
  scale_fill_manual(values=c("black","ivory3","gold3"),name="Cell Cycle Score",labels=c('Low','Inter','High')) + 
  scale_color_manual(values=c("black","ivory3","gold3"),name="Cell Cycle Score",labels=c('Low','Inter','High')) +
  scale_x_discrete(name="Cell Cycle Score",labels=c('Low','Inter','High')) +
  stat_summary(geom = "crossbar", width=0.74, fatten=3, color="white", 
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))})+
  ggtitle(label = 'Aneuploidy score in CCS subgroups') +
  ylab('Aneuploidy Score') +
  coord_cartesian(ylim=c(0,58)) +
  theme(legend.key = element_rect(fill='white'),
        axis.title.x = element_text(size = 29, face = "bold"),
        axis.title.y = element_text(size = 29, face = "italic"),
        plot.title = element_text(size = 35, face = "bold",hjust = c(0.5,0.5)),
        plot.subtitle = element_text(size = 12.5,hjust = c(0.5,0.5)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x =element_text(size=34,face = 'bold',angle = 45,hjust = c(1,1)),
        axis.text.y =element_text(size=34,face = 'bold',hjust = c(1,1)),
        plot.margin = unit(c(0.5,0.5,0,0.5),'cm'),
        strip.background = element_rect(fill='gray100'),strip.text = element_text(face = 'italic',size = 34))

Fig3D.boxplot.stats <- TukeyHSD(aov(CCS.Aneu.boxplot$`Aneuploidy Score`~as.factor(CCS.Aneu.boxplot$value)))

######
# Figure 3E - KM curve

Fig3E.KM <- survplot(Surv(PANCAN$Surv_15_cens,PANCAN$PFI_15_cens) ~ (PANCAN$CCS),hr.pos=NA,show.nrisk = T, data= PANCAN, 
         main = bquote(atop(bold('All patients\n (n = 9515)'),italic('\nPan-cancer'))),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","ivory3","gold3"),snames=c("Low","Inter","High"),legend.pos = 'bottomleft',stitle=NA,mark=20,cex.lab=1.1, cex.main=1.0)

CCS.15.pval.PANCAN = coxph(Surv(PANCAN$Surv_15_cens,PANCAN$PFI_15_cens) ~ PANCAN$CCS_ct, data = PANCAN)

pval=ifelse(formatC(format='f',summary(CCS.15.pval.PANCAN)$sctest[3], digits=3) < '0.001', '< 0.001',formatC(format='f',summary(CCS.15.pval.PANCAN)$sctest[3], digits=3))
legend("topright", bty="n", legend = paste("p ", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)




