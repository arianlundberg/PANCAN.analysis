#### Figure 2
### Used MAFTOOLS to download mutation information (gene level)

### FIGURE 2A (mutations)

### CCS low

CCS1.barplot <- pData(PANCAN) %>% group_by(type) %>% filter(CCS==1) %>% dplyr::select('TP53':ncol(.)) %>% summarise_all(sum)  
CCS1.barplot$col<- as.character(cancer_names[match(CCS1.barplot$type,cancer_names$abbs),]$cancer.colors)
CCS1.barplot_melted <- reshape::melt(transform(CCS1.barplot[c('type','col',names(sort(colSums(CCS1.barplot[,-c(1,ncol(CCS1.barplot))]),decreasing = T))[1:15])],type=type),id.vars = c('type','col'))

### CCS Intermediate

CCS2.barplot <- pData(PANCAN) %>% group_by(type) %>% filter(CCS==2) %>% dplyr::select('TP53':ncol(.)) %>% summarise_all(sum)  
CCS2.barplot$col<- as.character(cancer_names[match(CCS2.barplot$type,cancer_names$abbs),]$cancer.colors)
CCS2.barplot_melted <- reshape::melt(transform(CCS2.barplot[c('type','col',names(sort(colSums(CCS2.barplot[,-c(1,ncol(CCS2.barplot))]),decreasing = T))[1:15])],type=type),id.vars = c('type','col'))

### CCS High

CCS3.barplot <- pData(PANCAN) %>% group_by(type) %>% filter(CCS==3) %>% dplyr::select('TP53':ncol(.)) %>% summarise_all(sum)  
CCS3.barplot$col<- as.character(cancer_names[match(CCS3.barplot$type,cancer_names$abbs),]$cancer.colors)
CCS3.barplot_melted <- reshape::melt(transform(CCS3.barplot[c('type','col',names(sort(colSums(CCS3.barplot[,-c(1,ncol(CCS3.barplot))]),decreasing = T))[1:15])],type=type),id.vars = c('type','col'))


CCS1.bp <- ggplot(CCS1.barplot_melted, aes(x=reorder(variable,-value), y=value, fill=type)) + 
  geom_bar(stat="identity") + 
  geom_hline(yintercept=c(500), linetype="dashed",color='black',cex=1.1)+
  scale_fill_manual(values=CCS1.barplot$col,name="Cancer types")+
  scale_y_continuous(name=expression(bold(paste("Number of mutations")))) + scale_x_discrete(name="Genes") +
  ggtitle(label = 'CCS Low\n(n = 3145)') +
  theme(legend.key = element_rect(fill='white'),legend.position = 'None',
        strip.text.x =element_text(size = 25),
        axis.title.x = element_text(size = 25,face='bold'),
        axis.title.y = element_text(size = 25,face='bold'),
        plot.title = element_text(size = 35, face = "bold",hjust = c(0.5,0.5),vjust=(-10)),
        plot.subtitle = element_text(size = 27,hjust = c(0.5,0.5)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.text.x=element_text(size = 24,angle = 45,hjust = c(1,1)),
        axis.text.y=element_text(size = 30,face='bold'),
        plot.margin = unit(c(0.5,0.5,0.7,0.5),'cm'),
        axis.line.x = element_line(color="gray", size = 0.5),
        axis.line.y = element_line(color="gray", size = 0.5),
        legend.key.size = unit(1.3, "cm"),legend.title = element_text(size = 24,face='bold'),legend.text = element_text(size = 22))+ coord_cartesian(ylim = c(0,1800))

CCS2.bp <- ggplot(CCS2.barplot_melted, aes(x=reorder(variable,-value), y=value, fill=type)) + 
  geom_bar(stat="identity") +    
  geom_hline(yintercept=c(500), linetype="dashed",color='black',cex=1.1)+ 
  scale_fill_manual(values=CCS2.barplot$col,name="Cancer types")+
  scale_y_continuous(name=expression(bold(paste("Number of mutations")))) + scale_x_discrete(name="Genes") +
  ggtitle(label = 'CCS Intermediate\n(n = 3184)') +
  theme(legend.key = element_rect(fill='white'),legend.position = 'None',
        strip.text.x =element_text(size = 25),
        axis.title.x = element_text(size = 25,face='bold'),
        axis.title.y = element_text(size = 25,face='bold'),
        plot.title = element_text(size = 35, face = "bold",hjust = c(0.5,0.5),vjust=(-10)),
        plot.subtitle = element_text(size = 27,hjust = c(0.5,0.5)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.text.x=element_text(size = 24,angle = 45,hjust = c(1,1)),
        axis.text.y=element_text(size = 30,face='bold'),
        plot.margin = unit(c(0.5,0.5,0.7,0.5),'cm'),
        axis.line.x = element_line(color="gray", size = 0.5),
        axis.line.y = element_line(color="gray", size = 0.5),
        legend.key.size = unit(1.3, "cm"),legend.title = element_text(size = 24,face='bold'),legend.text = element_text(size = 22))+ coord_cartesian(ylim = c(0,1800))

CCS3.bp <- ggplot(CCS3.barplot_melted, aes(x=reorder(variable,-value), y=value, fill=type)) + 
  geom_bar(stat="identity") + 
  geom_hline(yintercept=c(500), linetype="dashed",color='black',cex=1.1) +
  scale_fill_manual(values=CCS3.barplot$col,name="Cancer types")+
  scale_y_continuous(name=expression(bold(paste("Number of mutations")))) + scale_x_discrete(name="Genes") +
  ggtitle(label = 'CCS High\n(n = 3186)') +
  theme(legend.key = element_rect(fill='white'),legend.position = 'None',
        strip.text.x =element_text(size = 25),
        axis.title.x = element_text(size = 25,face='bold'),
        axis.title.y = element_text(size = 25,face='bold'),
        plot.title = element_text(size = 35, face = "bold",hjust = c(0.5,0.5),vjust=(-10)),
        plot.subtitle = element_text(size = 27,hjust = c(0.5,0.5)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.text.x=element_text(size = 24,angle = 45,hjust = c(1,1)),
        axis.text.y=element_text(size = 30,face='bold'),
        plot.margin = unit(c(0.5,0.5,0.7,0.5),'cm'),
        axis.line.x = element_line(color="gray", size = 0.5),
        axis.line.y = element_line(color="gray", size = 0.5),
        legend.key.size = unit(1.3, "cm"),legend.title = element_text(size = 24,face='bold'),legend.text = element_text(size = 22))+ coord_cartesian(ylim = c(0,1800))

Fig2A.barplots <- arrangeGrob(CCS1.bp,CCS2.bp,CCS3.bp, nrow = 1, ncol=3)

### FIGURE 2B (Gains)

Gains <- pData(PANCAN) %>% group_by(type) %>% dplyr::select('Aneuploidy Score','CCS','1p':'22 (22q)') %>% as.data.frame();Gains[,4:ncol(Gains)] <- ifelse(Gains[,4:ncol(Gains)]>0,1,0)

# CCS1
CCS1.gain.barplot <- Gains %>% group_by(type) %>% filter(CCS==1) %>% dplyr::select(4:ncol(.)) %>% summarise_all(sum,na.rm=TRUE)
CCS1.gain.barplot$col<- as.character(cancer_names[match(CCS1.gain.barplot$type,cancer_names$abbs),]$cancer.colors)
CCS1.gain.barplot_melted <- melt(transform(CCS1.gain.barplot[c('type','col',names(sort(colSums(CCS1.gain.barplot[,-c(1,ncol(CCS1.gain.barplot))]),decreasing = T))[1:15])],type=type),id.vars = c('type','col'))
CCS1.gain.barplot_melted$variable <- as.factor(gsub(pattern = 'X',replacement = '',x = CCS1.gain.barplot_melted$variable))

# CCS2
CCS2.gain.barplot <- Gains %>% group_by(type) %>% filter(CCS==2) %>% dplyr::select(4:ncol(.)) %>% summarise_all(sum,na.rm=TRUE)
CCS2.gain.barplot$col<- as.character(cancer_names[match(CCS2.gain.barplot$type,cancer_names$abbs),]$cancer.colors)
CCS2.gain.barplot_melted <- melt(transform(CCS2.gain.barplot[c('type','col',names(sort(colSums(CCS2.gain.barplot[,-c(1,ncol(CCS2.gain.barplot))]),decreasing = T))[1:15])],type=type),id.vars = c('type','col'))
CCS2.gain.barplot_melted$variable <- as.factor(gsub(pattern = 'X',replacement = '',x = CCS2.gain.barplot_melted$variable))
CCS2.gain.barplot_melted$variable <- as.factor(gsub(pattern = '13..13q.',replacement = '13q',x = CCS2.gain.barplot_melted$variable))

# CCS3
CCS3.gain.barplot <- Gains %>% group_by(type) %>% filter(CCS==3) %>% dplyr::select(4:ncol(.)) %>% summarise_all(sum,na.rm=TRUE)
CCS3.gain.barplot$col<- as.character(cancer_names[match(CCS3.gain.barplot$type,cancer_names$abbs),]$cancer.colors)
CCS3.gain.barplot_melted <- melt(transform(CCS3.gain.barplot[c('type','col',names(sort(colSums(CCS3.gain.barplot[,-c(1,ncol(CCS3.gain.barplot))]),decreasing = T))[1:15])],type=type),id.vars = c('type','col'))
CCS3.gain.barplot_melted$variable <- as.factor(gsub(pattern = 'X',replacement = '',x = CCS3.gain.barplot_melted$variable))


CCS1.gain.bp <- ggplot(CCS1.gain.barplot_melted, aes(x=reorder(variable,-value), y=value, fill=type)) + 
  geom_bar(stat="identity") + 
  geom_hline(yintercept=c(500), linetype="dashed",color='black',cex=1) +
  scale_fill_manual(values=CCS1.gain.barplot$col,name="Cancer types")+
  scale_y_continuous(name=expression(bold(paste("Number of arm-level gains")))) + scale_x_discrete(name="Chromosome location") +
  ggtitle(label = 'CCS Low\n(n = 3145)') +
  theme(legend.key = element_rect(fill='white'),legend.position = 'None',
        strip.text.x =element_text(size = 25),
        axis.title.x = element_text(size = 25,face='bold'),
        axis.title.y = element_text(size = 25,face='bold'),
        plot.title = element_text(size = 35, face = "bold",hjust = c(0.5,0.5),vjust=(-10)),
        plot.subtitle = element_text(size = 27,hjust = c(0.5,0.5)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.text.x=element_text(size = 24,angle = 50,hjust = c(1,1),face='bold'),
        axis.text.y=element_text(size = 30,face='bold'),
        plot.margin = unit(c(0.5,0.5,0.7,0.5),'cm'),
        axis.line.x = element_line(color="gray", size = 0.5),
        axis.line.y = element_line(color="gray", size = 0.5),
        legend.key.size = unit(1.3, "cm"),legend.title = element_text(size = 24,face='bold'),legend.text = element_text(size = 22))+ coord_cartesian(ylim = c(0,1500),xlim=c(1,14.77))

CCS2.gain.bp <- ggplot(CCS2.gain.barplot_melted, aes(x=reorder(variable,-value), y=value, fill=type)) + 
  geom_bar(stat="identity") + 
  geom_hline(yintercept=c(500), linetype="dashed",color='black',cex=1) +
  scale_fill_manual(values=CCS2.gain.barplot$col,name="Cancer types")+
  scale_y_continuous(name=expression(bold(paste("Number of arm-level gains")))) + scale_x_discrete(name="Chromosome location") +
  ggtitle(label = 'CCS Intermediate\n(n = 3184)') +
  theme(legend.key = element_rect(fill='white'),legend.position = 'None',
        strip.text.x =element_text(size = 25),
        axis.title.x = element_text(size = 25,face='bold'),
        axis.title.y = element_text(size = 25,face='bold'),
        plot.title = element_text(size = 35, face = "bold",hjust = c(0.5,0.5),vjust=(-10)),
        plot.subtitle = element_text(size = 27,hjust = c(0.5,0.5)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.text.x=element_text(size = 24,angle = 50,hjust = c(1,1),face='bold'),
        axis.text.y=element_text(size = 30,face='bold'),
        plot.margin = unit(c(0.5,0.5,0.7,0.5),'cm'),
        axis.line.x = element_line(color="gray", size = 0.5),
        axis.line.y = element_line(color="gray", size = 0.5),
        legend.key.size = unit(1.3, "cm"),legend.title = element_text(size = 24,face='bold'),legend.text = element_text(size = 22))+ coord_cartesian(ylim = c(0,1500),xlim=c(1,14.77))


CCS3.gain.bp <- ggplot(CCS3.gain.barplot_melted, aes(x=reorder(variable,-value), y=value, fill=type)) + 
  geom_bar(stat="identity") + 
  geom_hline(yintercept=c(500), linetype="dashed",color='black',cex=1) +
  scale_fill_manual(values=CCS3.gain.barplot$col,name="Cancer types")+
  scale_y_continuous(name=expression(bold(paste("Number of arm-level gains")))) + scale_x_discrete(name="Chromosome location") +
  ggtitle(label = 'CCS High\n(n = 3186)') +
  theme(legend.key = element_rect(fill='white'),legend.position = 'None',
        strip.text.x =element_text(size = 25),
        axis.title.x = element_text(size = 25,face='bold'),
        axis.title.y = element_text(size = 25,face='bold'),
        plot.title = element_text(size = 35, face = "bold",hjust = c(0.5,0.5),vjust=(-10)),
        plot.subtitle = element_text(size = 27,hjust = c(0.5,0.5)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.text.x=element_text(size = 24,angle = 50,hjust = c(1,1),face='bold'),
        axis.text.y=element_text(size = 30,face='bold'),
        plot.margin = unit(c(0.5,0.5,0.7,0.5),'cm'),
        axis.line.x = element_line(color="gray", size = 0.5),
        axis.line.y = element_line(color="gray", size = 0.5),
        legend.key.size = unit(1.3, "cm"),legend.title = element_text(size = 24,face='bold'),legend.text = element_text(size = 22))+ coord_cartesian(ylim = c(0,1500),xlim=c(1,14.77))

Fig2B.barplots <- arrangeGrob(CCS1.gain.bp,CCS2.gain.bp,CCS3.gain.bp, nrow = 1, ncol=3)



### FIGURE 2B (Deletion)

Losses <- pData(PANCAN) %>% group_by(type) %>% dplyr::select('CCS','1p':'22 (22q)') %>% as.data.frame();Losses[,3:ncol(Losses)] <- ifelse(Losses[,3:ncol(Losses)]<0,1,0)

# CCS1
CCS1.loss.barplot <- Losses %>% group_by(type) %>% filter(CCS==1) %>% dplyr::select(3:ncol(.)) %>% summarise_all(sum,na.rm=TRUE)
CCS1.loss.barplot$col<- as.character(cancer_names[match(CCS1.loss.barplot$type,cancer_names$abbs),]$cancer.colors)
CCS1.loss.barplot_melted <- melt(transform(CCS1.loss.barplot[c('type','col',names(sort(colSums(CCS1.loss.barplot[,-c(1,ncol(CCS1.loss.barplot))]),decreasing = T))[1:15])],type=type),id.vars = c('type','col'))
CCS1.loss.barplot_melted$variable <- as.factor(gsub(pattern = 'X',replacement = '',x = CCS1.loss.barplot_melted$variable))
CCS1.loss.barplot_melted$variable <- as.factor(gsub(pattern = '13..13q.',replacement = '13q',x = CCS1.loss.barplot_melted$variable))
CCS1.loss.barplot_melted$variable <- as.factor(gsub(pattern = '14..14q.',replacement = '14q',x = CCS1.loss.barplot_melted$variable))
CCS1.loss.barplot_melted$variable <- as.factor(gsub(pattern = '21..21q.',replacement = '21q',x = CCS1.loss.barplot_melted$variable))
CCS1.loss.barplot_melted$variable <- as.factor(gsub(pattern = '22..22q.',replacement = '22q',x = CCS1.loss.barplot_melted$variable))

# CCS2
CCS2.loss.barplot <- Losses %>% group_by(type) %>% filter(CCS==2) %>% dplyr::select(3:ncol(.)) %>% summarise_all(sum,na.rm=TRUE)
CCS2.loss.barplot$col<- as.character(cancer_names[match(CCS2.loss.barplot$type,cancer_names$abbs),]$cancer.colors)
CCS2.loss.barplot_melted <- melt(transform(CCS2.loss.barplot[c('type','col',names(sort(colSums(CCS2.loss.barplot[,-c(1,ncol(CCS2.loss.barplot))]),decreasing = T))[1:15])],type=type),id.vars = c('type','col'))
CCS2.loss.barplot_melted$variable <- as.factor(gsub(pattern = 'X',replacement = '',x = CCS2.loss.barplot_melted$variable))
CCS2.loss.barplot_melted$variable <- as.factor(gsub(pattern = '13..13q.',replacement = '13q',x = CCS2.loss.barplot_melted$variable))
CCS2.loss.barplot_melted$variable <- as.factor(gsub(pattern = '13..13q.',replacement = '13q',x = CCS2.loss.barplot_melted$variable))
CCS2.loss.barplot_melted$variable <- as.factor(gsub(pattern = '14..14q.',replacement = '14q',x = CCS2.loss.barplot_melted$variable))
CCS2.loss.barplot_melted$variable <- as.factor(gsub(pattern = '21..21q.',replacement = '21q',x = CCS2.loss.barplot_melted$variable))
CCS2.loss.barplot_melted$variable <- as.factor(gsub(pattern = '22..22q.',replacement = '22q',x = CCS2.loss.barplot_melted$variable))

# CCS3
CCS3.loss.barplot <- Losses %>% group_by(type) %>% filter(CCS==3) %>% dplyr::select(3:ncol(.)) %>% summarise_all(sum,na.rm=TRUE)
CCS3.loss.barplot$col<- as.character(cancer_names[match(CCS3.loss.barplot$type,cancer_names$abbs),]$cancer.colors)
CCS3.loss.barplot_melted <- melt(transform(CCS3.loss.barplot[c('type','col',names(sort(colSums(CCS3.loss.barplot[,-c(1,ncol(CCS3.loss.barplot))]),decreasing = T))[1:15])],type=type),id.vars = c('type','col'))
CCS3.loss.barplot_melted$variable <- as.factor(gsub(pattern = 'X',replacement = '',x = CCS3.loss.barplot_melted$variable))
CCS3.loss.barplot_melted$variable <- as.factor(gsub(pattern = '13..13q.',replacement = '13q',x = CCS3.loss.barplot_melted$variable))
CCS3.loss.barplot_melted$variable <- as.factor(gsub(pattern = '14..14q.',replacement = '14q',x = CCS3.loss.barplot_melted$variable))
CCS3.loss.barplot_melted$variable <- as.factor(gsub(pattern = '21..21q.',replacement = '21q',x = CCS3.loss.barplot_melted$variable))
CCS3.loss.barplot_melted$variable <- as.factor(gsub(pattern = '22..22q.',replacement = '22q',x = CCS3.loss.barplot_melted$variable))


CCS1.loss.bp <- ggplot(CCS1.loss.barplot_melted, aes(x=reorder(variable,-value), y=value, fill=type)) + 
  geom_bar(stat="identity") + 
  geom_hline(yintercept=c(500), linetype="dashed",color='black',cex=1) +
  scale_fill_manual(values=CCS1.loss.barplot$col,name="Cancer types")+
  scale_y_continuous(name=expression(bold(paste("Number of arm-level losses")))) + scale_x_discrete(name="Chromosome location") +
  ggtitle(label = 'CCS Low\n(n = 3145)') +
  theme(legend.key = element_rect(fill='white'),legend.position = 'None',
        strip.text.x =element_text(size = 25),
        axis.title.x = element_text(size = 25,face='bold'),
        axis.title.y = element_text(size = 25,face='bold'),
        plot.title = element_text(size = 35, face = "bold",hjust = c(0.5,0.5),vjust=(-10)),
        plot.subtitle = element_text(size = 27,hjust = c(0.5,0.5)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.text.x=element_text(size = 24,angle = 50,hjust = c(1,1),face='bold'),
        axis.text.y=element_text(size = 30,face='bold'),
        plot.margin = unit(c(0.5,0.5,0.7,0.5),'cm'),
        axis.line.x = element_line(color="gray", size = 0.5),
        axis.line.y = element_line(color="gray", size = 0.5),
        legend.key.size = unit(1.3, "cm"),legend.title = element_text(size = 24,face='bold'),legend.text = element_text(size = 22))+ coord_cartesian(ylim = c(0,1500),xlim=c(1,14.77))

CCS2.loss.bp <- ggplot(CCS2.loss.barplot_melted, aes(x=reorder(variable,-value), y=value, fill=type)) + 
  geom_bar(stat="identity") + 
  geom_hline(yintercept=c(500), linetype="dashed",color='black',cex=1) +
  scale_fill_manual(values=CCS2.loss.barplot$col,name="Cancer types")+
  scale_y_continuous(name=expression(bold(paste("Number of arm-level losses")))) + scale_x_discrete(name="Chromosome location") +
  ggtitle(label = 'CCS Intermediate\n(n = 3184)') +
  theme(legend.key = element_rect(fill='white'),legend.position = 'None',
        strip.text.x =element_text(size = 25),
        axis.title.x = element_text(size = 25,face='bold'),
        axis.title.y = element_text(size = 25,face='bold'),
        plot.title = element_text(size = 35, face = "bold",hjust = c(0.5,0.5),vjust=(-10)),
        plot.subtitle = element_text(size = 27,hjust = c(0.5,0.5)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.text.x=element_text(size = 24,angle = 50,hjust = c(1,1),face='bold'),
        axis.text.y=element_text(size = 30,face='bold'),
        plot.margin = unit(c(0.5,0.5,0.7,0.5),'cm'),
        axis.line.x = element_line(color="gray", size = 0.5),
        axis.line.y = element_line(color="gray", size = 0.5),
        legend.key.size = unit(1.3, "cm"),legend.title = element_text(size = 24,face='bold'),legend.text = element_text(size = 22))+ coord_cartesian(ylim = c(0,1500),xlim=c(1,14.77))


CCS3.loss.bp <- ggplot(CCS3.loss.barplot_melted, aes(x=reorder(variable,-value), y=value, fill=type)) + 
  geom_bar(stat="identity") + 
  geom_hline(yintercept=c(500), linetype="dashed",color='black',cex=1) +
  scale_fill_manual(values=CCS3.loss.barplot$col,name="Cancer types")+
  scale_y_continuous(name=expression(bold(paste("Number of arm-level losses")))) + scale_x_discrete(name="Chromosome location") +
  ggtitle(label = 'CCS High\n(n = 3186)') +
  theme(legend.key = element_rect(fill='white'),legend.position = 'None',
        strip.text.x =element_text(size = 25),
        axis.title.x = element_text(size = 25,face='bold'),
        axis.title.y = element_text(size = 25,face='bold'),
        plot.title = element_text(size = 35, face = "bold",hjust = c(0.5,0.5),vjust=(-10)),
        plot.subtitle = element_text(size = 27,hjust = c(0.5,0.5)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.text.x=element_text(size = 24,angle = 50,hjust = c(1,1),face='bold'),
        axis.text.y=element_text(size = 30,face='bold'),
        plot.margin = unit(c(0.5,0.5,0.7,0.5),'cm'),
        axis.line.x = element_line(color="gray", size = 0.5),
        axis.line.y = element_line(color="gray", size = 0.5),
        legend.key.size = unit(1.3, "cm"),legend.title = element_text(size = 24,face='bold'),legend.text = element_text(size = 22))+ coord_cartesian(ylim = c(0,1500),xlim=c(1,14.77))

Fig2C.barplots <- arrangeGrob(CCS1.loss.bp,CCS2.loss.bp,CCS3.loss.bp, nrow = 1, ncol=3)


