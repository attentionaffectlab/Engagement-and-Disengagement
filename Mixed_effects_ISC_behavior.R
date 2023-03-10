packages <- c('ggplot2','readr','tidyr','Hmisc','plyr','RColorBrewer','reshape2','lme4','effects','emmeans','sjPlot','LMERConvenienceFunctions','ggsignif','standardize','scales')
lapply(packages, require, character.only = TRUE)



setwd('/media/melanni/DATADRIVE1/Echoes/scripts/Rcodes')
subj_isc=read.csv('csv_files/isc_subj.csv', sep=',');
dur_episode= read.csv('csv_files/duration_episodes.csv', sep=','); #duration of episodes
isc_cont= read.csv('csv_files/isc_continuous.csv', sep=',');


#------------------------------ISC Sujective emotional arousal annotations--------------------------------------------

#Mixed effects model (MEM)
isc_cont$isc_scaled= scale(isc_cont$isc)
modelC<-lmer(isc_scaled~movie+(1|Subject),data=isc_cont) #model

dfEC=as.data.frame(effect('movie', modelC)) # save MEM effect to data frame
dfEC$isc_scaled=dfEC$fit

tab_model(modelC,show.intercept = F,show.std = T,show.est = F) #MEM results in table


#--------------------------------ISC of eng-diseng --------------------------------------------------------------
# scale, center isc scores
subj_isc$isc_scaled_engDiseng= scale(subj_isc$isc)

#Mixed effects model
modelED<-lmer(isc_scaled_engDiseng~condition*phase+number+(1|Subject), data= subj_isc)
modelNum<-lmer(isc_scaled_engDiseng~phase*number+(1|Subject), data= subj_isc)
tab_model(modelNum,show.intercept = F,show.std = T,show.est = F)

dfE=as.data.frame(effect('condition*phase', modelED)) #interaction model
dfE$isc_scaled_engDiseng=dfE$fit
dfE$phase<-factor(dfE$phase)

plot(effect('phase*number', modelED),layout = c(1,1))

# ## post-hoc contrasts
posthoci=emmeans(modelED,pairwise ~ phase|condition | number,var='condition',  infer = TRUE)
dfC=as.data.frame(posthoci$contrasts)
posthocN=emmeans(modelNum,pairwise ~ phase | number,var='phase',  infer = TRUE)
dfCN=as.data.frame(posthocN$contrasts)

# ## format table to 3 digits after comma
is.num <- sapply(dfC, is.numeric)
dfC[is.num] <- lapply(dfC[is.num], round, 3)
dfC
summary(modelED)

# MEM model
tab_model(modelED,show.intercept = F,show.std = T,show.est = F)

##################

#change order in ggplot:
#Turn 'treatment' column into a character vector
#dfN$condition <- as.character(dfN$condition)
#dfN$condition <- factor(dfN$condition, levels=unique(dfN$condition))

#---models------
modelED_augmented <- augment(modelED) #this is for the newer version of emmeans works ~.~
modelC_augmented <- augment(modelC) 
############################## PLOT DATA #########################################################

library(raincloudplots)
#source("/media/melanni/fc8c7ab1-5327-4e16-ade6-bf26d19b2ebf/Echoes/scripts/Rcodes/geom_flat_violin.R")
GroupCol = c('darkgreen','#a6611a')

# Engagement-disengagement rainclouds including estimated MEM effects and raw data
p1<-ggplot(data=modelED_augmented,aes(x = condition, y = isc_scaled_engDiseng, color = phase, fill = phase))+ 
  
  geom_flat_violin(aes(fill = phase),position = position_nudge(x = .2, y = 0),
                   adjust = 1.5, trim = FALSE, alpha = .3, colour = NA)+
  
  geom_point(position = position_jitter(width = .15), size = 2,alpha=.1,aes(color=phase,fill=phase))+
  theme_classic()+
  geom_boxplot(aes(fill=phase,color=phase), color = 'black', outlier.shape=NA, width = 0.32,alpha = 0.3)+
  geom_errorbar(data=dfE,aes(ymin=lower, ymax=upper, color=phase), 
                width=.2,position = position_nudge(x = .35, y = 0)) +
  stat_summary(fun = "mean", geom = "point", shape = 23, size = 2, fill = NA,
                position = position_dodge(width = .33)) +
  scale_color_manual(values=GroupCol)+
  guides(col = guide_legend(nrow = 1),size=.2) + 
  scale_fill_manual(values=GroupCol, labels = c("Dis-engagement", "Engagement"))+
  ylab('ISC (subjective arousal annotations)')+
  theme(text = element_text(size=10,  family="sans"),strip.background=element_blank(),
        legend.text = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        axis.title.x = element_text(size=10),
        axis.text.x =  element_text(size=10),
        axis.text.y =  element_text(size=10),
        strip.text.x = element_text(size=10),
        axis.ticks.x = element_blank())


p1

ggsave(filename='/output_path/rainplots_subj.pdf',p1,
       width = 1.8, height = 2, units ='in')

#ISC continuous 

GroupCol=c('gray17','gray17')
p2<-ggplot(data=modelC_augmented,aes(x = movie, y = isc_scaled, fill = movie))+ 
  
  geom_flat_violin(aes(fill = movie),position = position_nudge(x = .2, y = 0),
                   adjust = 1.5, trim = FALSE, alpha = 0.4, colour = NA)+
  
  geom_point(position = position_jitter(width = .12), size = 3,alpha=.3,aes(color=movie,fill=movie))+
  theme_classic()+
  geom_boxplot(aes(fill=movie,color=movie), color = 'black', outlier.shape=NA, width = 0.28,alpha = 0.4)+
  geom_errorbar(data=dfEC,aes(ymin=lower, ymax=upper, color=movie,fill=movie), 
                width=.2,position = position_nudge(x = .35, y = 0)) +
  stat_summary(fun = "mean", geom = "point", shape = 23, size = 2, fill = NA,
               width=.2, position = position_dodge(width = .33)) +
  
  scale_color_manual(values=GroupCol)+
  
  guides(col = guide_legend(nrow = 1),size=.2) + 
  scale_fill_manual(values=GroupCol, labels = c("Dis-engagement", "Engagement"))+
  #ylab('ISC (subjective arousal annotations)')+
  theme(text = element_text(size=10,  family="sans"),strip.background=element_blank(),
        legend.text = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        axis.title.x = element_text(size=10),
        axis.text.x =  element_text(size=10),
        axis.text.y =  element_text(size=10),
        strip.text.x = element_text(size=10),
        axis.ticks.x = element_blank())


p2

