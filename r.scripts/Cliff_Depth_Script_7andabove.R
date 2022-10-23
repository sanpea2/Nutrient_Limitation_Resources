############################ Libraries
library(MuMIn)
library(lme4)
library(DescTools)
library(PMCMRplus)
library(psych)
library(MASS)
library(rcompanion)
library(ggplot2)
library(car)
library(sjPlot)
library(bbmle)
library(emmeans)
library(nlme)
library(fitdistrplus)
library(ggpubr)
library(dplyr)
############################
setwd("C:/users/aldoa/Desktop/NDS_MEGA/Cliff_August//")
nds.data = read.csv("Cliff_August_5m_DONE.csv", header = T)

tail(nds.data)
############################
#:::::: boxplot for 6 treatments  ::::::::
all.data.c = nds.data[nds.data$sample == "Control",]
all.data.p = nds.data[nds.data$sample == "PO4",]
all.data.no3 = nds.data[nds.data$sample == "NO3",]
all.data.no3p = nds.data[nds.data$sample == "NO3+P",]
all.data.nh4 = nds.data[nds.data$sample == "NH4",]
all.data.nh4p = nds.data[nds.data$sample == "NH4+P",]

boxplot((all.data.c$chla_corrected), (all.data.p$chla_corrected), (all.data.no3$chla_corrected), (all.data.no3p$chla_corrected), (all.data.nh4$chla_corrected), (all.data.nh4p$chla_corrected),
frame = F, axes = FALSE, col = c("white", "gray", "gray", "black", "gray", "black"), ylim = c(0,2.2))
axis(side = 2, at = seq(from = 0, to = 2.5, by = .5))
mtext(c("Control", "PO4", "Nitrate", "Nitrate_PO4", "NH4","NH4_PO4" ), side = 1, at = 1:6, cex = 1.5)
mtext(c("5 Meter NDS Sep7-26" ), side = 3, at = 3.5, cex = 1.5)
mtext(expression(paste("chlorphyll-a (ug cm"^"-2"*")")), side = 2, line = 2, cex = 1.5)
#mtext(expression(paste("chlorphyll-a (ug cm"^"-2"*"day"^"-1"*")")), side = 2, line = 2, cex = 1.5)
############################
plotNormalHistogram(nds.data$chla_corrected)
plotNormalHistogram(log10(nds.data$chla_corrected))
nds.data_tuk = transformTukey(nds.data$chla_corrected,plotit=TRUE) #Tukey transformation
plotNormalHistogram(nds.data_tuk)
############################
qqPlot(nds.data$chla_corrected)
qqPlot(notransformmodel)
qqPlot(log10(nds.data$chla_corrected))
qqPlot(nds.data_tuk)

############################
notransformmodel = lm(nds.data$chla_corrected ~ nds.data$sample, 
                      data=nds.data)
residuals.notransform = residuals(notransformmodel)
###
logtransformmodel = lm(log10(nds.data$chla_corrected) ~ nds.data$sample,
                       data=nds.data)
residuals.logtransform = residuals(logtransformmodel)
###
all.data_tukmodel = lm(nds.data_tuk ~ nds.data$sample,
           data=nds.data)
residuals.all.data_tuk = residuals(all.data_tukmodel)
#############################
shapiro.test(nds.data$chla_corrected)
shapiro.test(residuals.notransform)
shapiro.test(log10(nds.data$chla_corrected))
shapiro.test(residuals.logtransform)
shapiro.test(nds.data_tuk)
shapiro.test(residuals.all.data_tuk)

###############################################################
###############################################################
###############################################################
#Anova(M2Var, type="3")
# use Type I only when there is a serious theoretical reason for it, use Type II when there is no interaction, use Type III when there is interaction.

nds.data$tuktransform <- nds.data_tuk
nds.data$logtransform <- log10(nds.data$chla_corrected)

M2<-gls(chla_corrected~sample, data=nds.data)
M2Var<-update (M2,weights = varIdent(form=~ 1 | sample))
qqPlot(residuals(M2))
shapiro.test(resid(M2))
Anova(M2)
Anova(M2Var)

M2Tuk<-gls(tuktransform~sample, data=nds.data)
M2TukVar<-update (M2Tuk,weights = varIdent(form=~ 1 | sample))
qqPlot(residuals(M2Tuk))
shapiro.test(resid(M2Tuk))
Anova(M2Tuk)
Anova(M2TukVar)


M2Log<-gls(logtransform~sample, data=nds.data)
M2LogVar<-update (M2Log,weights = varIdent(form=~ 1 | sample))
qqPlot(residuals(M2Log))
shapiro.test(resid(M2Log))
Anova(M2Log)
Anova(M2LogVar)

AICtab(M2Log,M2LogVar,M2TukVar, M2Tuk, M2, M2Var, base=T, sort=T, weights=T)
#choose lowest score for best fit model

##
plotNormalHistogram(residuals(M2Var))
plotNormalHistogram(residuals(M2Tuk))
plotNormalHistogram(residuals(M2Log))

plot_model(M2, type="pred", terms=c("sample"), show.data = TRUE, )
plot_model(M2Var, type="pred", terms=c("sample"), show.data = TRUE, )
plot_model(M2Tuk, type="pred", terms=c("sample"), show.data = TRUE, )
plot_model(M2TukVar, type="pred", terms=c("sample"), show.data = TRUE, )
plot_model(M2Log, type="pred", terms=c("sample"), show.data = TRUE, )
plot_model(M2LogVar, type="pred", terms=c("sample"), show.data = TRUE, )

nds.data$sample<-as.factor(nds.data$sample)
DunnettTest(x=nds.data$tuktransform, g=nds.data$sample)
DunnettTest(x=nds.data$logtransform, g=nds.data$sample)
DunnettTest(x=nds.data$chla_corrected, g=nds.data$sample)

ght
emmeans(M2, specs=pairwise~sample)
emmeans(M2Var, specs=pairwise~sample)
emmeans(M2Tuk, specs=pairwise~sample)
emmeans(M2TukVar, specs=pairwise~sample)
emmeans(M2Log, specs=pairwise~sample)
emmeans(M2LogVar, specs=pairwise~sample)



#############################
#TO DO: Insert Facundo's script to statistically remove outliers and compare the results -> done. AIC determines tukeys ladder of transformation is the best model in 5/6 depths
#but log is better for answering my question (differences)

#TO DO: Boxplots after data removal
#TO DO: Standardize units for comparison of my results to great lakes, castle lit, (Ugcm-2 is great -> keep it as is and sunmmarize comparisons
# as some papers use log)

## N O T E S ##
#calculated mean and SD in excel per sample -> removed only one value for control from 5m (3SDs higher).  All other chl-a data fits between 2 SD
#colors in r
nds.data = read.csv("NDS_Castle_Sep07.26.21_NR_7.5m.csv", header = T)
adjnds.data <- cbind(nds.data[,2-3])


library(ggplot2)
library(dplyr)

                         
stax <- adjnds.data %>%
group_by(sample) %>%
summarise( 
  n=n(),
  mean=mean(chla_corrected),
  sd=sd(chla_corrected),
  sd=(sd*1),
  sd2=(sd*2),
  sdmax=(sd2+mean),
  smin=(mean-sd2)
) %>%
mutate( se=sd/sqrt(n))  %>%
mutate( ic=se * qt((1-0.05)/2 + .5, n-1))
stax


                         #reorder
                         
stax$sample <- factor(stax$sample,levels = c("Control", "P", "Nitrate", "Nitrate_P", "NH4", "NH4_P"))
nds.data_bp <- nds.data
nds.data_bp$sample <- factor(nds.data_bp$sample,levels = c("Control", "P", "Nitrate", "Nitrate_P", "NH4", "NH4_P"))

nds.data$logtransform
                         # Standard deviation
ggplot(stax) +
geom_bar( aes(x=sample, y=mean), stat="identity", color = "black", fill=c("gray", "springgreen3", "springgreen4", "darkslategray2", "darkslategray4", "slategray4"), alpha=1) +
geom_errorbar( aes(x=sample, ymax=mean+sd, ymin=mean-sd), width=0.15, colour="black", alpha=1, size=.5) +
labs(title="Castle 5 Meter NDS Sep7-26", 
x = "sample", y = expression('chlorphyll-a'~(ugcm)^-2)) +
ylim(0, .85) + #yaxis limits
coord_fixed(10) +   #constrain image
#scale_fill_manual(values='black') + 
theme(panel.background = element_rect(fill = "white",  #adapted from theme_bw()
colour = NA), panel.border = element_rect(fill = NA, 
colour = "grey20"),  plot.title = element_text(hjust = 0.5),
panel.grid.minor = element_line(size = rel(0.5)), 
strip.background = element_rect(fill = "grey85", 
colour = "grey20"), legend.key = element_rect(fill = "white", 
colour = NA), complete = TRUE)

                      
                         # Standard Error
ggplot(stax) +
geom_bar( aes(x=sample, y=mean), stat="identity", fill="forestgreen", alpha=0.5) +
geom_errorbar( aes(x=sample, ymin=mean-se, ymax=mean+se), width=0.4, colour="orange", alpha=0.9, size=1.5) +
ggtitle("using standard error")
                         
                         # Confidence Interval
ggplot(stax) +
geom_bar( aes(x=sample, y=mean), stat="identity", color = "black", fill=c("gray", "springgreen3", "springgreen4", "darkslategray2", "darkslategray4", "slategray4"), alpha=1) +
geom_errorbar( aes(x=sample, ymin=mean-ic, ymax=mean+ic), width=0.4, colour="black", alpha=0.9, size=.5) +
ggtitle("using confidence interval")
                        
                         
                        

M2Log
ggplot(nds.data_bp, aes(x=sample, y=chla_corrected)) + 
geom_boxplot(color = "black", fill=c("gray", "slategray4", "darkslategray2", "darkslategray4","springgreen3", "springgreen4"), alpha=1)
                          
ggplot(nds.data_bp, aes(x=sample, y=logtransform)) + 
geom_boxplot(color = "black", fill=c("gray", "slategray4", "darkslategray2", "darkslategray4","springgreen3", "springgreen4"), alpha=1)

plot_model(M2Log, 
           type="pred", 
           terms=c("sample"), 
           show.data = TRUE, )





