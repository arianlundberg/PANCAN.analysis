##### Supplemental figure 3 and Supplemental figure 4

registerDoParallel(cores = detectCores(all.tests = FALSE, logical = TRUE)-2)
`%notin%` <- Negate(`%in%`)
PFI_surv_15_PANCAN= Surv((PANCAN$Surv_15_cens),PANCAN$PFI_15_cens)


# Multivariate analysis For CCS (low,intermediate,high):
### Pathological stage - releveled/ Stage I as a reference group
PT_fun <- function(PT){if(PT  %in% c("Stage I","Stage IA","Stage IB")) "I"
  else if (PT %in% c("Stage II","Stage IIA","Stage IIB","Stage IIC")) "II"
  else if (PT %in% c("Stage III","Stage IIIA","Stage IIIB","Stage IIIC")) "III"
  else if (PT %in% c("Stage IV","Stage IVA","Stage IVB","Stage IVC")) "IV"
  else if (PT %in% c("","[Discrepancy]","[Unknown]","I or II NOS","Stage 0","Stage IS","Stage X",NA)) NA}
PANCAN$CCS_inter.relevel <- as.numeric(PANCAN$CCS)
PANCAN$CCS_intra.relevel <- as.numeric(PANCAN$CCS_intra)
PANCAN$pathologic_stage_relevel <- unlist(lapply(PANCAN$pathologic_stage, PT_fun))
PANCAN$pathologic_stage_relevel <- factor(PANCAN$pathologic_stage_relevel,c('I','II','III','IV'))

PANCAN$age_group2 <- as.numeric(cut(PANCAN$age_at_initial_pathologic_diagnosis,breaks=c('-inf',30,40,50,60,70,'inf'), include.lowest=T))

PANCAN$radiation_therapy_relevel <- factor(ifelse(PANCAN$radiation_therapy %in% c("","[Discrepancy]","[Unknown]",NA),NA,PANCAN$radiation_therapy),levels=c('YES','NO'))

###### P.VALUES (UNI AND MULTIVARIATE) ADJUSTED WITH FDR CORRECTION
#### EXCLUDE CANCERS THAT MAY NOT ADJUSTED FOR EVERYTHING
adjusted.all.cancers <- levels(PANCAN$type)[(levels(PANCAN$type)) %notin% c('CESC','DLBC','GBM','KIRC','LGG','LUAD','OV','PCPG','PRAD','SARC','THYM','UCEC','UCS')]
adjusted.not.path.cancers <- c('CESC','DLBC','GBM','LGG','LUAD','OV','PCPG','PRAD','SARC','THYM','UCEC','UCS')

### CCS-INTER

uni.p.inter <- data.frame(foreach(i = levels(PANCAN$type), .packages = 'survival', .combine = 'rbind') %dopar% {
  uni.p  <- summary(coxph(Surv((PANCAN$Surv_15_cens),PANCAN$PFI_15_cens) ~
                            (PANCAN$CCS_inter.relevel), subset = (PANCAN$type %in% i), data = PANCAN))$sctest[3]},row.names=levels(PANCAN$type))


### Adjusted for all
multi.p.all <- data.frame(foreach(i = adjusted.all.cancers, .packages = 'survival', .combine = 'rbind') %dopar% {
  p.all <- summary(coxph(Surv((PANCAN$Surv_15_cens),PANCAN$PFI_15_cens) ~ 
                           age_group2 + as.factor(gender) +
                           radiation_therapy_relevel + pathologic_stage_relevel + 
                           CCS_inter.relevel, data=PANCAN,subset=(PANCAN$type %in% i)))$coeff[,'Pr(>|z|)']['CCS_inter.relevel']},row.names = adjusted.all.cancers)

### Adjusted for all except Pathological stage
multi.not.path <- data.frame(foreach(i = adjusted.not.path.cancers, .packages = 'survival', .combine = 'rbind') %dopar% {
  not.path <- summary(coxph(Surv((PANCAN$Surv_15_cens),PANCAN$PFI_15_cens) ~ 
                              age_group2 + as.factor(gender) +
                              radiation_therapy_relevel + 
                              # pathologic_stage_relevel + 
                              CCS_inter.relevel, data=PANCAN,subset=(PANCAN$type %in% i)))$coeff[,'Pr(>|z|)']['CCS_inter.relevel']},row.names= adjusted.not.path.cancers)


### Adjusted for all except radiotherapy
multi.not.radio <- data.frame(CCS_inter.relevel=summary(coxph(Surv((PANCAN$Surv_15_cens),PANCAN$PFI_15_cens) ~ 
                                                                age_group2 + as.factor(gender) +
                                                                #radiation_therapy_relevel + 
                                                                pathologic_stage_relevel + 
                                                                CCS_inter.relevel, data=PANCAN,subset=(PANCAN$type %in% 'KIRC')))$coeff[,'Pr(>|z|)']['CCS_inter.relevel'],row.names = 'KIRC')


####### GATHER all p. values (uni and multi) to adjust for multiple testing using BH method

multi.p <- data.frame(multi.p=rbind(multi.p.all,multi.not.path,multi.not.radio))
p.tables <- merge(uni.p.inter,multi.p,by='row.names',sort=F);rownames(p.tables) <- p.tables$Row.names;p.tables <- p.tables[,-1];colnames(p.tables) <- c('uni.p','multi.p')

p.tables <- p.tables %>% tibble::add_column(uni.BH = p.adjust(.$uni.p, method = 'BH'),.after = 'multi.p')
p.tables.inter <- p.tables %>% tibble::add_column(multi.BH = p.adjust(.$multi.p, method = 'BH'),.after = 'uni.BH')


#### CCS-intra


uni.p.intra <- data.frame(foreach(i = levels(PANCAN$type), .packages = 'survival', .combine = 'rbind') %dopar% {
  uni.p  <- summary(coxph(Surv((PANCAN$Surv_15_cens),PANCAN$PFI_15_cens) ~
                            (PANCAN$CCS_intra.relevel), subset = (PANCAN$type %in% i), data = PANCAN))$sctest[3]},row.names=levels(PANCAN$type))
### Adjusted for all
multi.p.all <- data.frame(foreach(i = adjusted.all.cancers, .packages = 'survival', .combine = 'rbind') %dopar% {
  p.all <- summary(coxph(Surv((PANCAN$Surv_15_cens),PANCAN$PFI_15_cens) ~ 
                           age_group2 + as.factor(gender) +
                           radiation_therapy_relevel + pathologic_stage_relevel + 
                           CCS_intra.relevel, data=PANCAN,subset=(PANCAN$type %in% i)))$coeff[,'Pr(>|z|)']['CCS_intra.relevel']},row.names = adjusted.all.cancers)

### Adjusted for all except Pathological stage
multi.not.path <- data.frame(foreach(i = adjusted.not.path.cancers, .packages = 'survival', .combine = 'rbind') %dopar% {
  not.path <- summary(coxph(Surv((PANCAN$Surv_15_cens),PANCAN$PFI_15_cens) ~ 
                              age_group2 + as.factor(gender) +
                              radiation_therapy_relevel + 
                              # pathologic_stage_relevel + 
                              CCS_intra.relevel, data=PANCAN,subset=(PANCAN$type %in% i)))$coeff[,'Pr(>|z|)']['CCS_intra.relevel']},row.names= adjusted.not.path.cancers)


### Adjusted for all except radiotherapy
multi.not.radio <- data.frame(CCS_intra.relevel=summary(coxph(Surv((PANCAN$Surv_15_cens),PANCAN$PFI_15_cens) ~ 
                                                                age_group2 + as.factor(gender) +
                                                                #radiation_therapy_relevel + 
                                                                pathologic_stage_relevel + 
                                                                CCS_intra.relevel, data=PANCAN,subset=(PANCAN$type %in% 'KIRC')))$coeff[,'Pr(>|z|)']['CCS_intra.relevel'],row.names = 'KIRC')


####### GATHER all p. values (uni and multi) to adjust for multiple testing using BH method

multi.p <- data.frame(multi.p=rbind(multi.p.all,multi.not.path,multi.not.radio))
p.tables <- merge(uni.p.intra,multi.p,by='row.names',sort=F);rownames(p.tables) <- p.tables$Row.names;p.tables <- p.tables[,-1];colnames(p.tables) <- c('uni.p','multi.p')

p.tables <- p.tables %>% tibble::add_column(uni.BH = p.adjust(.$uni.p, method = 'BH'),.after = 'multi.p')
p.tables.intra <- p.tables %>% tibble::add_column(multi.BH = p.adjust(.$multi.p, method = 'BH'),.after = 'uni.BH')



##### SUPPLEMENTAL FIGURE 3

#### KM - CCS - ALL CANCERS


pdf(file= "destination folder",width = 10.72, height = 18.19,
    onefile = T, paper = "special")
par(mfrow=c(5,3),xpd=NA,
    omi=c(0.85, 0.85, 0.85, 0.85), mai=c(0.7,0.55,0.7,0.55))

#### CCS-inter


## ACC

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='ACC'),
         main = bquote(atop(bold('ACC'),'\n(n = 76)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[1], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[1], digits=3))
legend("topright", bty="n", legend = paste("p =", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)
mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[1], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[1], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p =", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+1.2,y=legend("topright", bty="n", legend = paste("p =", pval),
                                                                                      lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p* =", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

## BLCA

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='BLCA'),
         main = bquote(atop(bold('BLCA'),'\n(n = 398)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[2], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[2], digits=3))

legend("topright", bty="n", legend = paste("p =", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[2], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[2], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p =", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+1.3,y=legend("topright", bty="n", legend = paste("p =", pval),
                                                                                      lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p* =", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)



legend(y =1.9,x=-3,legend = c("Low","Intermediate","High"),x.intersp = 0.2,
       text.width = c(2.5,2.5,5),fill=c("black","gray","gold3"),
       horiz = T, bty="n",border=FALSE, cex=1.5)
mtext(outer = T,line=3,expression(bold("Cell cycle score")~"("~italic("Pan-cancer")~")"),cex=1)


## BRCA

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='BRCA'),
         main = bquote(atop(bold('BRCA'),'\n(n = 1038)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[3], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[3], digits=3))

legend("topright", bty="n", legend = paste("p =", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[3], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[3], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p =", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+1.4,y=legend("topright", bty="n", legend = paste("p =", pval),
                                                                                      lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p* =", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)



## CESC

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='CESC'),
         main = bquote(atop(bold('CESC'),'\n(n = 291)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("gray","gold3"),snames=c("Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[4], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[4], digits=3))

legend("topright", bty="n", legend = paste("p =", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[4], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[4], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p =", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+1.2,y=legend("topright", bty="n", legend = paste("p =", pval),
                                                                                      lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p** =", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)



## CHOL

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='CHOL'),
         main = bquote(atop(bold('CHOL'),'\n(n = 36)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[5], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[5], digits=3))

legend("topright", bty="n", legend = paste("p =", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[5], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[5], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p =", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+0.5,y=legend("topright", bty="n", legend = paste("p =", pval),
                                                                                      lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p* =", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)



## COAD

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='COAD'),
         main = bquote(atop(bold('COAD'),'\n(n = 428)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[6], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[6], digits=3))

legend("topright", bty="n", legend = paste("p =", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[6], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[6], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p =", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+1.2,y=legend("topright", bty="n", legend = paste("p =", pval),
                                                                                      lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p* =", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)


## DLBC

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='DLBC'),
         main = bquote(atop(bold('DLBC'),'\n(n = 47)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("gray","gold3"),snames=c("Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[7], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[7], digits=3))

legend("topright", bty="n", legend = paste("p =", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[7], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[7], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p =", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+1.3,y=legend("topright", bty="n", legend = paste("p =", pval),
                                                                                      lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p** =", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)


## ESCA

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='ESCA'),
         main = bquote(atop(bold('ESCA'),'\n(n = 161)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[8], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[8], digits=3))

legend("topright", bty="n", legend = paste("p =", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[8], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[8], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p =", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+1,y=legend("topright", bty="n", legend = paste("p =", pval),
                                                                                    lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p* =", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)



## GBM

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='GBM'),
         main = bquote(atop(bold('GBM'),'\n(n = 151)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[9], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[9], digits=3))

legend("topright", bty="n", legend = paste("p =", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[9], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[9], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p =", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+0.35,y=legend("topright", bty="n", legend = paste("p =", pval),
                                                                                       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p** =", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

## HNSC

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='HNSC'),
         main = bquote(atop(bold('HNSC'),'\n(n = 503)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[10], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[10], digits=3))

legend("topright", bty="n", legend = paste("p =", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[10], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[10], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p =", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+1.3,y=legend("topright", bty="n", legend = paste("p =", pval),
                                                                                      lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p* =", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)



## KICH

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='KICH'),
         main = bquote(atop(bold('KICH'),'\n(n = 65)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray"),snames=c("Low","Inter"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[11], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[11], digits=3))


mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[11], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[11], digits=3))
legend(x=legend("bottomright", bty="n", legend = paste("p* =", mpval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+1.7,y=legend("bottomright", bty="n", legend = paste("p* =", mpval),
                                                                                      lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y+0.11, bty="n", legend = paste("p ", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

## KIRC

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='KIRC'),
         main = bquote(atop(bold('KIRC'),'\n(n = 480)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[12], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[12], digits=3))


mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[12], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[12], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p ", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+1.2,y=legend("topright", bty="n", legend = paste("p ", pval),
                                                                                      lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p^ =", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)



## KIRP

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='KIRP'),
         main = bquote(atop(bold('KIRP'),'\n(n = 280)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[13], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[13], digits=3))


mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[13], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[13], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p ", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+1.55,y=legend("topright", bty="n", legend = paste("p ", pval),
                                                                                       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p* =", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)



## LGG

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='LGG'),
         main = bquote(atop(bold('LGG'),'\n(n = 507)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[14], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[14], digits=3))

mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[14], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[14], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p ", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+1.24,y=legend("topright", bty="n", legend = paste("p ", pval),
                                                                                       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p** ", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)



## LIHC

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='LIHC'),
         main = bquote(atop(bold('LIHC'),'\n(n = 355)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[15], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[15], digits=3))

legend("topright", bty="n", legend = paste("p =", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[15], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[15], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p =", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+0.9,y=legend("topright", bty="n", legend = paste("p =", pval),
                                                                                      lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p* =", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)



## LUAD

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='LUAD'),
         main = bquote(atop(bold('LUAD'),'\n(n = 498)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[16], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[16], digits=3))

legend("topright", bty="n", legend = paste("p =", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[16], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[16], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p =", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+1.2,y=legend("topright", bty="n", legend = paste("p =", pval),
                                                                                      lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p** =", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)


## LUSC

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='LUSC'),
         main = bquote(atop(bold('LUSC'),'\n(n = 479)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[17], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[17], digits=3))

legend("topright", bty="n", legend = paste("p =", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[17], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[17], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p =", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+1.32,y=legend("topright", bty="n", legend = paste("p =", pval),
                                                                                       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p* =", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)


## MESO

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='MESO'),
         main = bquote(atop(bold('MESO'),'\n(n = 81)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[18], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[18], digits=3))

legend("topright", bty="n", legend = paste("p =", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[18], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[18], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p =", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+0.6,y=legend("topright", bty="n", legend = paste("p =", pval),
                                                                                      lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p* =", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)



## OV

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='OV'),
         main = bquote(atop(bold('OV'),'\n(n = 289)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[19], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[19], digits=3))

legend("topright", bty="n", legend = paste("p =", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[19], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[19], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p =", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+1.26,y=legend("topright", bty="n", legend = paste("p =", pval),
                                                                                       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p** =", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)


## PAAD

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='PAAD'),
         main = bquote(atop(bold('PAAD'),'\n(n = 158)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[20], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[20], digits=3))

legend("topright", bty="n", legend = paste("p =", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[20], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[20], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p =", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+0.7,y=legend("topright", bty="n", legend = paste("p =", pval),
                                                                                      lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p* =", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)


## PCPG

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='PCPG'),
         main = bquote(atop(bold('PCPG'),'\n(n = 160)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray"),snames=c("Low","Inter"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[21], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[21], digits=3))
mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[21], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[21], digits=3))
legend(x=legend("bottomright", bty="n", legend = paste("p** =", mpval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+2.7,y=legend("bottomright", bty="n", legend = paste("p** =", mpval),
                                                                                      lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y+0.11, bty="n", legend = paste("p =", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)


## PRAD

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='PRAD'),
         main = bquote(atop(bold('PRAD'),'\n(n = 471)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray"),snames=c("Low","Inter"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[22], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[22], digits=3))
mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[22], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[22], digits=3))
legend(x=legend("bottomright", bty="n", legend = paste("p** =", mpval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+2.5,y=legend("bottomright", bty="n", legend = paste("p** =", mpval),
                                                                                      lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y+0.11, bty="n", legend = paste("p =", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)




## READ

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='READ'),
         main = bquote(atop(bold('READ'),'\n(n = 154)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[23], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[23], digits=3))

legend("topright", bty="n", legend = paste("p =", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[23], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[23], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p =", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+1,y=legend("topright", bty="n", legend = paste("p =", pval),
                                                                                    lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p* =", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)



## SARC

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='SARC'),
         main = bquote(atop(bold('SARC'),'\n(n = 242)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[24], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[24], digits=3))

mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[24], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[24], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p ", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+1.36,y=legend("topright", bty="n", legend = paste("p ", pval),
                                                                                       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p** =", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)


## SKCM

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='SKCM'),
         main = bquote(atop(bold('SKCM'),'\n(n = 458)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[25], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[25], digits=3))

legend("topright", bty="n", legend = paste("p =", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[25], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[25], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p =", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+1.32,y=legend("topright", bty="n", legend = paste("p =", pval),
                                                                                       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p* =", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)


## STAD

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='STAD'),
         main = bquote(atop(bold('STAD'),'\n(n = 404)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[26], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[26], digits=3))

legend("topright", bty="n", legend = paste("p =", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[26], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[26], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p =", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+1,y=legend("topright", bty="n", legend = paste("p =", pval),
                                                                                    lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p* =", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)


## TGCT

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='TGCT'),
         main = bquote(atop(bold('TGCT'),'\n(n = 133)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[27], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[27], digits=3))

legend("topright", bty="n", legend = paste("p =", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[27], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[27], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p =", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+1.32,y=legend("topright", bty="n", legend = paste("p =", pval),
                                                                                       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p* =", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)


## THCA

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='THCA'),
         main = bquote(atop(bold('THCA'),'\n(n = 462)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

THCA.15.pval.PANCAN = coxph(Surv((PANCAN$Surv_15_cens),PANCAN$PFI_15_cens) ~
                              PANCAN$CCS_inter.relevel, subset = (PANCAN$type=='THCA'), data = PANCAN)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[28], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[28], digits=3))
mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[28], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[28], digits=3))
legend(x=legend("bottomright", bty="n", legend = paste("p* =", mpval),col = NULL,
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+1.9,y=legend("bottomright", bty="n", legend = paste("p* =", mpval),col=NULL,
                                                                                      lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y+0.11, bty="n", legend = paste("p =", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)


## THYM

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='THYM'),
         main = bquote(atop(bold('THYM'),'\n(n = 103)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[29], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[29], digits=3))
mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[29], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[29], digits=3))
legend(x=legend("bottomright", bty="n", legend = paste("p** =", mpval),col = NULL,
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+1.9,y=legend("bottomright", bty="n", legend = paste("p** =", mpval),col=NULL,
                                                                                      lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y+0.11, bty="n", legend = paste("p =", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

## UCEC

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='UCEC'),
         main = bquote(atop(bold('UCEC'),'\n(n = 514)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[30], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[30], digits=3))

mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[30], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[30], digits=3))
legend(x=legend("bottomright", bty="n", legend = paste("p** =", mpval),col = NULL,
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+2.2,y=legend("bottomright", bty="n", legend = paste("p** =", mpval),col=NULL,
                                                                                      lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y+0.11, bty="n", legend = paste("p =", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

## UCS

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='UCS'),
         main = bquote(atop(bold('UCS'),'\n(n = 56)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("gray","gold3"),snames=c("Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[31], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[31], digits=3))

legend("topright", bty="n", legend = paste("p =", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[31], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[31], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p =", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+0.99,y=legend("topright", bty="n", legend = paste("p =", pval),
                                                                                       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p** =", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)



## UVM

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='UVM'),
         main = bquote(atop(bold('UVM'),'\n(n = 80)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray"),snames=c("Low","Inter"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

UVM.15.pval.PANCAN = coxph(Surv((PANCAN$Surv_15_cens),PANCAN$PFI_15_cens) ~
                             PANCAN$CCS_inter.relevel, subset = (PANCAN$type=='UVM'), data = PANCAN)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[32], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[32], digits=3))

mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[32], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[32], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p ", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+0.69,y=legend("topright", bty="n", legend = paste("p ", pval),
                                                                                       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p* =", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

dev.off()

############ 

#### SUPPLEMENTAL FIGURE 4

#### KM - CCS-inter vs Intra / CANCER type

pdf(file= "folder",width = 10.72, height = 18.19,title = "Inter vs Intra-CCS",
    #height = 8.72, width = 11.69,
    onefile = T, paper = "special")
#par(mfrow=c(2,3),mar=c(5,5,5,5))
par(mfrow=c(4,2),xpd=NA,
    omi=c(0.85, 0.85, 0.85, 0.85), mai=c(0.7,0.55,0.7,0.55))


## KIRC - Inter

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='KIRC'),
         main = bquote(atop(bold('KIRC'),'\n(n = 480)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

KIRC.15.pval.PANCAN = coxph(Surv((PANCAN$Surv_15_cens),PANCAN$PFI_15_cens) ~
                              PANCAN$CCS_inter.relevel, subset = (PANCAN$type=='KIRC'), data = PANCAN)


pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[12], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[12], digits=3))


KIRC_CCS.inter.uni.LR.PANCAN <- data.frame("n"= KIRC.15.pval.PANCAN$n,"events"=KIRC.15.pval.PANCAN$nevent,"cIndex"= summary(KIRC.15.pval.PANCAN)$concordance[1],"low.CI95%"=summary(KIRC.15.pval.PANCAN)$conf.int[3],"upp.CI95%"=summary(KIRC.15.pval.PANCAN)$conf.int[4],"LR.X2"=signif(anova(KIRC.15.pval.PANCAN)['PANCAN$CCS_inter.relevel',2], digits =3),"p.value"=signif(anova(KIRC.15.pval.PANCAN)['PANCAN$CCS_inter.relevel',4], digits =3), row.names = "KIRC_CCS.inter")

KIRC.15.mpval.PANCAN <- coxph(Surv((PANCAN$Surv_15_cens),PANCAN$PFI_15_cens) ~
                                age_group2 + as.factor(gender) +
                                #radiation_therapy_relevel + 
                                pathologic_stage_relevel +
                                CCS_inter.relevel, data=PANCAN,subset=(PANCAN$type=='KIRC'))

mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[12], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[12], digits=3))

legend(x=legend("topright", bty="n", legend = paste("p ", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+0.7,y=legend("topright", bty="n", legend = paste("p ", pval),
                                                                                      lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p^ =", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)


KIRC_CCS.inter.multi.LR.PANCAN <- data.frame("n"= KIRC.15.mpval.PANCAN$n,"events"=KIRC.15.mpval.PANCAN$nevent,"LR.X2"=signif(anova(KIRC.15.mpval.PANCAN)['CCS_inter.relevel',2], digits =3),
                                             "p.value"=signif(anova(KIRC.15.mpval.PANCAN)['CCS_inter.relevel',4], digits =3), row.names = "KIRC_CCS.inter")


## KIRC - Intra

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_intra.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='KIRC'),
         main = bquote(atop(bold('KIRC'),'\n(n = 480)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

KIRC.15.pval.PANCAN = coxph(Surv((PANCAN$Surv_15_cens),PANCAN$PFI_15_cens) ~
                              PANCAN$CCS_intra.relevel, subset = (PANCAN$type=='KIRC'), data = PANCAN)


pval=ifelse(formatC(format='f',p.tables.intra$uni.BH[12], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.intra$uni.BH[12], digits=3))


KIRC_CCS.intra.uni.LR.PANCAN <- data.frame("n"= KIRC.15.pval.PANCAN$n,"events"=KIRC.15.pval.PANCAN$nevent,"cIndex"= summary(KIRC.15.pval.PANCAN)$concordance[1],"low.CI95%"=summary(KIRC.15.pval.PANCAN)$conf.int[3],"upp.CI95%"=summary(KIRC.15.pval.PANCAN)$conf.int[4],"LR.X2"=signif(anova(KIRC.15.pval.PANCAN)['PANCAN$CCS_intra.relevel',2], digits =3),"p.value"=signif(anova(KIRC.15.pval.PANCAN)['PANCAN$CCS_intra.relevel',4], digits =3), row.names = "KIRC_CCS.intra")

KIRC.15.mpval.PANCAN <- coxph(Surv((PANCAN$Surv_15_cens),PANCAN$PFI_15_cens) ~
                                age_group2 + as.factor(gender) +
                                #radiation_therapy_relevel + 
                                pathologic_stage_relevel +
                                CCS_intra.relevel, data=PANCAN,subset=(PANCAN$type=='KIRC'))

mpval=ifelse(formatC(format='f',p.tables.intra$multi.BH[12], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.intra$multi.BH[12], digits=3))

legend(x=legend("topright", bty="n", legend = paste("p =", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+0.7,y=legend("topright", bty="n", legend = paste("p =", pval),
                                                                                      lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p^ =", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)


KIRC_CCS.intra.multi.LR.PANCAN <- data.frame("n"= KIRC.15.mpval.PANCAN$n,"events"=KIRC.15.mpval.PANCAN$nevent,"LR.X2"=signif(anova(KIRC.15.mpval.PANCAN)['CCS_intra.relevel',2], digits =3),
                                             "p.value"=signif(anova(KIRC.15.mpval.PANCAN)['CCS_intra.relevel',4], digits =3), row.names = "KIRC_CCS.intra")




legend(y =1.6,x=-10,legend = c("Low","Intermediate","High"),x.intersp = 0.2,
       text.width = c(2.5,2.5,4),fill=c("black","gray","gold3"),
       horiz = T, bty="n",border=FALSE, cex=1.5)
mtext(outer = T,line=3,expression(bold("Cell cycle score")),cex=1.5)
mtext(outer = T,line=-0.5,expression(italic("Pan-cancer")),cex=1.2,adj = 0.23)
mtext(outer = T,line=-0.5,expression(italic("Intra-cancer")),cex=1.2,adj = 0.8)



## KIRP -- Inter


KIRP.15.pval.PANCAN = coxph(Surv((PANCAN$Surv_15_cens),PANCAN$PFI_15_cens) ~
                              PANCAN$CCS_inter.relevel, subset = (PANCAN$type=='KIRP'), data = PANCAN)



KIRP_CCS.inter.uni.LR.PANCAN <- data.frame("n"= KIRP.15.pval.PANCAN$n,"events"=KIRP.15.pval.PANCAN$nevent,"cIndex"= summary(KIRP.15.pval.PANCAN)$concordance[1],"low.CI95%"=summary(KIRP.15.pval.PANCAN)$conf.int[3],"upp.CI95%"=summary(KIRP.15.pval.PANCAN)$conf.int[4],"LR.X2"=signif(anova(KIRP.15.pval.PANCAN)['PANCAN$CCS_inter.relevel',2], digits =3),"p.value"=signif(anova(KIRP.15.pval.PANCAN)['PANCAN$CCS_inter.relevel',4], digits =3), row.names = "KIRP_CCS.inter")


KIRP.15.mpval.PANCAN <- coxph(Surv((PANCAN$Surv_15_cens),PANCAN$PFI_15_cens) ~ 
                                age_group2 + as.factor(gender) +
                                radiation_therapy_relevel + pathologic_stage_relevel + 
                                CCS_inter.relevel, data=PANCAN,subset=(PANCAN$type=='KIRP'))


KIRP_CCS.inter.multi.LR.PANCAN <- data.frame("n"= KIRP.15.mpval.PANCAN$n,"events"=KIRP.15.mpval.PANCAN$nevent,"LR.X2"=signif(anova(KIRP.15.mpval.PANCAN)['CCS_inter.relevel',2], digits =3),
                                             "p.value"=signif(anova(KIRP.15.mpval.PANCAN)['CCS_inter.relevel',4], digits =3), row.names = "KIRP_CCS.inter")



## KIRP -- Intra


KIRP.15.pval.PANCAN = coxph(Surv((PANCAN$Surv_15_cens),PANCAN$PFI_15_cens) ~
                              PANCAN$CCS_intra.relevel, subset = (PANCAN$type=='KIRP'), data = PANCAN)


KIRP_CCS.intra.uni.LR.PANCAN <- data.frame("n"= KIRP.15.pval.PANCAN$n,"events"=KIRP.15.pval.PANCAN$nevent,"cIndex"= summary(KIRP.15.pval.PANCAN)$concordance[1],"low.CI95%"=summary(KIRP.15.pval.PANCAN)$conf.int[3],"upp.CI95%"=summary(KIRP.15.pval.PANCAN)$conf.int[4],"LR.X2"=signif(anova(KIRP.15.pval.PANCAN)['PANCAN$CCS_intra.relevel',2], digits =3),"p.value"=signif(anova(KIRP.15.pval.PANCAN)['PANCAN$CCS_intra.relevel',4], digits =3), row.names = "KIRP_CCS.intra")


KIRP.15.mpval.PANCAN <- coxph(Surv((PANCAN$Surv_15_cens),PANCAN$PFI_15_cens) ~ 
                                age_group2 + as.factor(gender) +
                                radiation_therapy_relevel + pathologic_stage_relevel + 
                                CCS_intra.relevel, data=PANCAN,subset=(PANCAN$type=='KIRP'))

mpval=ifelse(formatC(format='f',p.tables.intra$multi.BH[13], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.intra$multi.BH[13], digits=3))

KIRP_CCS.intra.multi.LR.PANCAN <- data.frame("n"= KIRP.15.mpval.PANCAN$n,"events"=KIRP.15.mpval.PANCAN$nevent,"LR.X2"=signif(anova(KIRP.15.mpval.PANCAN)['CCS_intra.relevel',2], digits =3),
                                             "p.value"=signif(anova(KIRP.15.mpval.PANCAN)['CCS_intra.relevel',4], digits =3), row.names = "KIRP_CCS.intra")


## LGG - Inter

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='LGG'),
         main = bquote(atop(bold('LGG'),'\n(n = 507)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

LGG.15.pval.PANCAN = coxph(Surv((PANCAN$Surv_15_cens),PANCAN$PFI_15_cens) ~
                             PANCAN$CCS_inter.relevel, subset = (PANCAN$type=='LGG'), data = PANCAN)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[14], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[14], digits=3))

legend("topright", bty="n", legend = paste("p ", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)


LGG_CCS.inter.uni.LR.PANCAN <- data.frame("n"= LGG.15.pval.PANCAN$n,"events"=LGG.15.pval.PANCAN$nevent,"cIndex"= summary(LGG.15.pval.PANCAN)$concordance[1],"low.CI95%"=summary(LGG.15.pval.PANCAN)$conf.int[3],"upp.CI95%"=summary(LGG.15.pval.PANCAN)$conf.int[4],"LR.X2"=signif(anova(LGG.15.pval.PANCAN)['PANCAN$CCS_inter.relevel',2], digits =3),"p.value"=signif(anova(LGG.15.pval.PANCAN)['PANCAN$CCS_inter.relevel',4], digits =3), row.names = "LGG_CCS.inter")


LGG.15.mpval.PANCAN <- coxph(Surv((PANCAN$Surv_15_cens),PANCAN$PFI_15_cens) ~ 
                               age_group2 + as.factor(gender) +
                               radiation_therapy_relevel +
                               #pathologic_stage_relevel + 
                               CCS_inter.relevel, data=PANCAN,subset=(PANCAN$type=='LGG'))

mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[14], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[14], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p ", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+0.7,y=legend("topright", bty="n", legend = paste("p ", pval),
                                                                                      lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p** ", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

LGG_CCS.inter.multi.LR.PANCAN <- data.frame("n"= LGG.15.mpval.PANCAN$n,"events"=LGG.15.mpval.PANCAN$nevent,"LR.X2"=signif(anova(LGG.15.mpval.PANCAN)['CCS_inter.relevel',2], digits =3),
                                            "p.value"=signif(anova(LGG.15.mpval.PANCAN)['CCS_inter.relevel',4], digits =3), row.names = "LGG_CCS.inter")

## LGG - Intra


survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_intra.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='LGG'),
         main = bquote(atop(bold('LGG'),'\n(n = 507)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

LGG.15.pval.PANCAN = coxph(Surv((PANCAN$Surv_15_cens),PANCAN$PFI_15_cens) ~
                             PANCAN$CCS_intra.relevel, subset = (PANCAN$type=='LGG'), data = PANCAN)

pval=ifelse(formatC(format='f',p.tables.intra$uni.BH[14], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.intra$uni.BH[14], digits=3))

legend("topright", bty="n", legend = paste("p ", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)


LGG_CCS.intra.uni.LR.PANCAN <- data.frame("n"= LGG.15.pval.PANCAN$n,"events"=LGG.15.pval.PANCAN$nevent,"cIndex"= summary(LGG.15.pval.PANCAN)$concordance[1],"low.CI95%"=summary(LGG.15.pval.PANCAN)$conf.int[3],"upp.CI95%"=summary(LGG.15.pval.PANCAN)$conf.int[4],"LR.X2"=signif(anova(LGG.15.pval.PANCAN)['PANCAN$CCS_intra.relevel',2], digits =3),"p.value"=signif(anova(LGG.15.pval.PANCAN)['PANCAN$CCS_intra.relevel',4], digits =3), row.names = "LGG_CCS.intra")


LGG.15.mpval.PANCAN <- coxph(Surv((PANCAN$Surv_15_cens),PANCAN$PFI_15_cens) ~ 
                               age_group2 + as.factor(gender) +
                               radiation_therapy_relevel +
                               #pathologic_stage_relevel + 
                               CCS_intra.relevel, data=PANCAN,subset=(PANCAN$type=='LGG'))

mpval=ifelse(formatC(format='f',p.tables.intra$multi.BH[14], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.intra$multi.BH[14], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p ", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+0.7,y=legend("topright", bty="n", legend = paste("p ", pval),
                                                                                      lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p** ", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

LGG_CCS.intra.multi.LR.PANCAN <- data.frame("n"= LGG.15.mpval.PANCAN$n,"events"=LGG.15.mpval.PANCAN$nevent,"LR.X2"=signif(anova(LGG.15.mpval.PANCAN)['CCS_intra.relevel',2], digits =3),
                                            "p.value"=signif(anova(LGG.15.mpval.PANCAN)['CCS_intra.relevel',4], digits =3), row.names = "LGG_CCS.intra")



## SARC -- Inter

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='SARC'),
         main = bquote(atop(bold('SARC'),'\n(n = 242)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

SARC.15.pval.PANCAN = coxph(Surv((PANCAN$Surv_15_cens),PANCAN$PFI_15_cens) ~
                              PANCAN$CCS_inter.relevel, subset = (PANCAN$type=='SARC'), data = PANCAN)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[24], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[24], digits=3))

legend("topright", bty="n", legend = paste("p ", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)


SARC_CCS.inter.uni.LR.PANCAN <- data.frame("n"= SARC.15.pval.PANCAN$n,"events"=SARC.15.pval.PANCAN$nevent,"cIndex"= summary(SARC.15.pval.PANCAN)$concordance[1],"low.CI95%"=summary(SARC.15.pval.PANCAN)$conf.int[3],"upp.CI95%"=summary(SARC.15.pval.PANCAN)$conf.int[4],"LR.X2"=signif(anova(SARC.15.pval.PANCAN)['PANCAN$CCS_inter.relevel',2], digits =3),"p.value"=signif(anova(SARC.15.pval.PANCAN)['PANCAN$CCS_inter.relevel',4], digits =3), row.names = "SARC_CCS.inter")


SARC.15.mpval.PANCAN <- coxph(Surv((PANCAN$Surv_15_cens),PANCAN$PFI_15_cens) ~ 
                                age_group2 + as.factor(gender) +
                                radiation_therapy_relevel +
                                #pathologic_stage_relevel + 
                                CCS_inter.relevel, data=PANCAN,subset=(PANCAN$type=='SARC'))

mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[24], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[24], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p ", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+0.8,y=legend("topright", bty="n", legend = paste("p ", pval),
                                                                                      lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p** =", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

SARC_CCS.inter.multi.LR.PANCAN <- data.frame("n"= SARC.15.mpval.PANCAN$n,"events"=SARC.15.mpval.PANCAN$nevent,"LR.X2"=signif(anova(SARC.15.mpval.PANCAN)['CCS_inter.relevel',2], digits =3),
                                             "p.value"=signif(anova(SARC.15.mpval.PANCAN)['CCS_inter.relevel',4], digits =3), row.names = "SARC_CCS.inter")



## SARC -- Intra

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_intra.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='SARC'),
         main = bquote(atop(bold('SARC'),'\n(n = 242)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

SARC.15.pval.PANCAN = coxph(Surv((PANCAN$Surv_15_cens),PANCAN$PFI_15_cens) ~
                              PANCAN$CCS_intra.relevel, subset = (PANCAN$type=='SARC'), data = PANCAN)

pval=ifelse(formatC(format='f',p.tables.intra$uni.BH[24], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.intra$uni.BH[24], digits=3))

legend("topright", bty="n", legend = paste("p =", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)


SARC_CCS.intra.uni.LR.PANCAN <- data.frame("n"= SARC.15.pval.PANCAN$n,"events"=SARC.15.pval.PANCAN$nevent,"cIndex"= summary(SARC.15.pval.PANCAN)$concordance[1],"low.CI95%"=summary(SARC.15.pval.PANCAN)$conf.int[3],"upp.CI95%"=summary(SARC.15.pval.PANCAN)$conf.int[4],"LR.X2"=signif(anova(SARC.15.pval.PANCAN)['PANCAN$CCS_intra.relevel',2], digits =3),"p.value"=signif(anova(SARC.15.pval.PANCAN)['PANCAN$CCS_intra.relevel',4], digits =3), row.names = "SARC_CCS.intra")


SARC.15.mpval.PANCAN <- coxph(Surv((PANCAN$Surv_15_cens),PANCAN$PFI_15_cens) ~ 
                                age_group2 + as.factor(gender) +
                                radiation_therapy_relevel +
                                #pathologic_stage_relevel + 
                                CCS_intra.relevel, data=PANCAN,subset=(PANCAN$type=='SARC'))

mpval=ifelse(formatC(format='f',p.tables.intra$multi.BH[24], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.intra$multi.BH[24], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p =", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+0.8,y=legend("topright", bty="n", legend = paste("p =", pval),
                                                                                      lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p** =", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

SARC_CCS.intra.multi.LR.PANCAN <- data.frame("n"= SARC.15.mpval.PANCAN$n,"events"=SARC.15.mpval.PANCAN$nevent,"LR.X2"=signif(anova(SARC.15.mpval.PANCAN)['CCS_intra.relevel',2], digits =3),
                                             "p.value"=signif(anova(SARC.15.mpval.PANCAN)['CCS_intra.relevel',4], digits =3), row.names = "SARC_CCS.intra")

## UVM -- Inter

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_inter.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='UVM'),
         main = bquote(atop(bold('UVM'),'\n(n = 80)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray"),snames=c("Low","Inter"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

UVM.15.pval.PANCAN = coxph(Surv((PANCAN$Surv_15_cens),PANCAN$PFI_15_cens) ~
                             PANCAN$CCS_inter.relevel, subset = (PANCAN$type=='UVM'), data = PANCAN)

pval=ifelse(formatC(format='f',p.tables.inter$uni.BH[32], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$uni.BH[32], digits=3))

legend("topright", bty="n", legend = paste("p ", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

UVM_CCS.inter.uni.LR.PANCAN <- data.frame("n"= UVM.15.pval.PANCAN$n,"events"=UVM.15.pval.PANCAN$nevent,"cIndex"= summary(UVM.15.pval.PANCAN)$concordance[1],"low.CI95%"=summary(UVM.15.pval.PANCAN)$conf.int[3],"upp.CI95%"=summary(UVM.15.pval.PANCAN)$conf.int[4],"LR.X2"=signif(anova(UVM.15.pval.PANCAN)['PANCAN$CCS_inter.relevel',2], digits =3),"p.value"=signif(anova(UVM.15.pval.PANCAN)['PANCAN$CCS_inter.relevel',4], digits =3), row.names = "UVM_CCS.inter")


UVM.15.mpval.PANCAN <- coxph(Surv((PANCAN$Surv_15_cens),PANCAN$PFI_15_cens) ~ 
                               age_group2 + as.factor(gender) +
                               radiation_therapy_relevel + 
                               pathologic_stage_relevel + 
                               CCS_inter.relevel, data=PANCAN,subset=(PANCAN$type=='UVM'))

mpval=ifelse(formatC(format='f',p.tables.inter$multi.BH[32], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.inter$multi.BH[32], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p ", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+0.5,y=legend("topright", bty="n", legend = paste("p ", pval),
                                                                                      lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p* =", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

UVM_CCS.inter.multi.LR.PANCAN <- data.frame("n"= UVM.15.mpval.PANCAN$n,"events"=UVM.15.mpval.PANCAN$nevent,"LR.X2"=signif(anova(UVM.15.mpval.PANCAN)['CCS_inter.relevel',2], digits =3),
                                            "p.value"=signif(anova(UVM.15.mpval.PANCAN)['CCS_inter.relevel',4], digits =3), row.names = "UVM_CCS.inter")


## UVM -- Intra

survplot(PFI_surv_15_PANCAN ~ PANCAN$CCS_intra.relevel,hr.pos=NA,show.nrisk = T, data= PANCAN, subset = (PANCAN$type=='UVM'),
         main = bquote(atop(bold('UVM'),'\n(n = 80)')),xlab=expression(bold('Time (years)')),ylab=expression(bold("PFI"))
         ,lwd=3,col=c("black","gray","gold3"),snames=c("Low","Inter","High"),legend.pos = NA,stitle=NA,mark=20,cex.lab=1, cex.main=1.2)

UVM.15.pval.PANCAN = coxph(Surv((PANCAN$Surv_15_cens),PANCAN$PFI_15_cens) ~
                             PANCAN$CCS_intra.relevel, subset = (PANCAN$type=='UVM'), data = PANCAN)

pval=ifelse(formatC(format='f',p.tables.intra$uni.BH[32], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.intra$uni.BH[32], digits=3))

legend("topright", bty="n", legend = paste("p =", pval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)


UVM_CCS.intra.uni.LR.PANCAN <- data.frame("n"= UVM.15.pval.PANCAN$n,"events"=UVM.15.pval.PANCAN$nevent,"cIndex"= summary(UVM.15.pval.PANCAN)$concordance[1],"low.CI95%"=summary(UVM.15.pval.PANCAN)$conf.int[3],"upp.CI95%"=summary(UVM.15.pval.PANCAN)$conf.int[4],"LR.X2"=signif(anova(UVM.15.pval.PANCAN)['PANCAN$CCS_intra.relevel',2], digits =3),"p.value"=signif(anova(UVM.15.pval.PANCAN)['PANCAN$CCS_intra.relevel',4], digits =3), row.names = "UVM_CCS.intra")


UVM.15.mpval.PANCAN <- coxph(Surv((PANCAN$Surv_15_cens),PANCAN$PFI_15_cens) ~ 
                               age_group2 + as.factor(gender) +
                               radiation_therapy_relevel + 
                               #    pathologic_stage_relevel + 
                               CCS_intra.relevel, data=PANCAN,subset=(PANCAN$type=='UVM'))

mpval=ifelse(formatC(format='f',p.tables.intra$multi.BH[32], digits=3) < '0.001', '< 0.001',formatC(format='f',p.tables.intra$multi.BH[32], digits=3))
legend(x=legend("topright", bty="n", legend = paste("p =", pval),
                lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$x+0.4,y=legend("topright", bty="n", legend = paste("p =", pval),
                                                                                      lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)$text$y-0.11, bty="n", legend = paste("p* =", mpval),
       lty = FALSE, xjust = 0.5, yjust = 0.5, cex = 1.2)

UVM_CCS.intra.multi.LR.PANCAN <- data.frame("n"= UVM.15.mpval.PANCAN$n,"events"=UVM.15.mpval.PANCAN$nevent,"LR.X2"=signif(anova(UVM.15.mpval.PANCAN)['CCS_intra.relevel',2], digits =3),
                                            "p.value"=signif(anova(UVM.15.mpval.PANCAN)['CCS_intra.relevel',4], digits =3), row.names = "UVM_CCS.intra")


dev.off()


