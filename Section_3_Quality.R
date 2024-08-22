# Package loading and Working Directory -----------------------------------
library(lme4) # building mixed effects models
library(lmerTest) # statistical tests for MEMs 
setwd("C:/Users/de656/OneDrive - University of York/Research/VerticalFarm_observationalStudy/")

##########################
#  The aim of this script is:
##########################

####Load data
Microgreen_Biomass_Data <- read.table("Farm_Biomass_Data.txt", 
                                      header = T, sep = "")


# Formatting spectral peak data ------------------------------------------------

Spectral_Data <- read.table("Perfect_Sunrise_MRes_Aggregated_Spectral_Data", 
                            header = T, sep = "")

Microgreen_Biomass_Data$Unique_Identifier <- paste(Microgreen_Biomass_Data$X,
                                                   Microgreen_Biomass_Data$Y)

Spectral_Data$Unique_Identifier <- paste(Spectral_Data$X, Spectral_Data$Y)


Reordered_X <- match(unique(Spectral_Data$Unique_Identifier), 
                     Microgreen_Biomass_Data$Unique_Identifier)

Reordered_Biomass_Data <- Microgreen_Biomass_Data[Reordered_X,] # To get x and y to match

# Add peak red and blue to the biomass data

Reordered_Biomass_Data$Red_Light_Intensity_umols <- Spectral_Data[which(Spectral_Data[,5] == "Red"),4]
Reordered_Biomass_Data$Blue_Light_Intensity_umols <- Spectral_Data[which(Spectral_Data[,5] == "Blue"),4]
Reordered_Biomass_Data$Far_Red_Light_Intensity_umols <- Spectral_Data[which(Spectral_Data[,5] == "Far Red"),4]

write.csv(Reordered_Biomass_Data, file='RBFR_aggregatedData.txt')


Reduced_All_Model <- lmer(Biomass_g ~ 1 + Microgreen+ (1 | Bed : Microgreen) + (1 | Position_in_Bed : Microgreen), data=Reordered_Biomass_Data)


summary(Reduced_All_Model)

ranova(Reduced_All_Model)

plot(Reduced_All_Model)
plot(predict(Reduced_All_Model), Reordered_Biomass_Data$Biomass_g)

PAR_All_Model <- lmer(Biomass_g ~ 1 + Light_Intensity_umols * Microgreen + (1 | Bed : Microgreen) + (1 | Position_in_Bed : Microgreen), data=Reordered_Biomass_Data)

summary(PAR_All_Model)

ranova(PAR_All_Model)

anova(Reduced_All_Model, PAR_All_Model)

temp=anova(Reduced_All_Model, PAR_All_Model)

Red_All_Model <- lmer(Biomass_g ~ 1 + Red_Light_Intensity_umols * Microgreen + (1 | Bed : Microgreen) + (1 | Position_in_Bed : Microgreen), data=Reordered_Biomass_Data)

summary(Red_All_Model)

ranova(PAR_All_Model)

anova(Reduced_All_Model, PAR_All_Model)

FR_All_Model <- lmer(Biomass_g ~ 1 + Far_Red_Light_Intensity_umols * Microgreen + (1 | Bed : Microgreen) + (1 | Position_in_Bed : Microgreen), data=Reordered_Biomass_Data)

summary(FR_All_Model)

ranova(PAR_All_Model)

Blue_All_Model <- lmer(Biomass_g ~ 1 + Blue_Light_Intensity_umols * Microgreen + (1 | Bed : Microgreen) + (1 | Position_in_Bed : Microgreen), data=Reordered_Biomass_Data)


Quality_All_Model <- lmer(Biomass_g ~ 1 + Blue_Light_Intensity_umols * Microgreen + Red_Light_Intensity_umols * Microgreen + Far_Red_Light_Intensity_umols * Microgreen +
                         (1 | Bed : Microgreen) + (1 | Position_in_Bed : Microgreen), data=Reordered_Biomass_Data)



summary(Blue_All_Model)

anova(Reduced_All_Model, Quality_All_Model)

anova(PAR_All_Model, Quality_All_Model)

anova(PAR_All_Model, Reduced_All_Model)

plot(predict(Quality_All_Model), Reordered_Biomass_Data$Biomass_g)

plot(predict(Quality_All_Model)[which(Reordered_Biomass_Data$Microgreen=='Sunflower')], Reordered_Biomass_Data$Biomass_g[which(Reordered_Biomass_Data$Microgreen=='Sunflower')])
plot(predict(Quality_All_Model)[which(Reordered_Biomass_Data$Microgreen=='Kale(CN)')], Reordered_Biomass_Data$Biomass_g[which(Reordered_Biomass_Data$Microgreen=='Kale(CN)')])
plot(predict(Quality_All_Model)[which(Reordered_Biomass_Data$Microgreen=='Kale(Tozer)')], Reordered_Biomass_Data$Biomass_g[which(Reordered_Biomass_Data$Microgreen=='Kale(Tozer)')])
plot(predict(Quality_All_Model)[which(Reordered_Biomass_Data$Microgreen=='Radish')], Reordered_Biomass_Data$Biomass_g[which(Reordered_Biomass_Data$Microgreen=='Radish')])



plot(predict(Reduced_All_Model)[which(Reordered_Biomass_Data$Microgreen=='Sunflower')], Reordered_Biomass_Data$Biomass_g[which(Reordered_Biomass_Data$Microgreen=='Sunflower')])
plot(predict(Reduced_All_Model)[which(Reordered_Biomass_Data$Microgreen=='Kale(CN)')], Reordered_Biomass_Data$Biomass_g[which(Reordered_Biomass_Data$Microgreen=='Kale(CN)')])
plot(predict(Reduced_All_Model)[which(Reordered_Biomass_Data$Microgreen=='Kale(Tozer)')], Reordered_Biomass_Data$Biomass_g[which(Reordered_Biomass_Data$Microgreen=='Kale(Tozer)')])
plot(predict(Reduced_All_Model)[which(Reordered_Biomass_Data$Microgreen=='Radish')], Reordered_Biomass_Data$Biomass_g[which(Reordered_Biomass_Data$Microgreen=='Radish')])



temp=anova(Reduced_All_Model, PAR_All_Model)


Reordered_Biomass_Data$R_B=log10(Reordered_Biomass_Data$Red_Light_Intensity_umols/Reordered_Biomass_Data$Blue_Light_Intensity_umols)
Reordered_Biomass_Data$R_FR=log10(Reordered_Biomass_Data$Red_Light_Intensity_umols/Reordered_Biomass_Data$Far_Red_Light_Intensity_umols)
Reordered_Biomass_Data$FR_B=log10(Reordered_Biomass_Data$Far_Red_Light_Intensity_umols/Reordered_Biomass_Data$Blue_Light_Intensity_umols)


Ratio_All_Model <- lmer(Biomass_g ~ 1 + R_B * Microgreen + R_FR * Microgreen + FR_B * Microgreen +
                            (1 | Bed : Microgreen) + (1 | Position_in_Bed : Microgreen), data=Reordered_Biomass_Data)

anova(Reduced_All_Model, Ratio_All_Model)

################### TO DELETE
# Unique Identifiers match
# 
# Spectral_model <- lmer(Biomass_g ~ Red_Light_Intensity_umols + 
#                          Blue_Light_Intensity_umols + 
#                          Far_Red_Light_Intensity_umols + (1 | Microgreen) + 
#                     (1 | Bed) + (1 | Position_in_Bed) + 
#                       (1 | Sowing_Density) , data = Reordered_Biomass_Data)
# 
# plot(Reordered_Biomass_Data$Biomass_g, predict(Spectral_model))
# abline(a = 0, b = 1)
# 
# summary(Spectral_model)
# 
# rtemp=anova(Spectral_model) # Drop sowing density and bed position
# 
# Spectral_model_No_sow_or_bed <- lmer(Biomass_g ~ Red_Light_Intensity_umols + 
#                          Blue_Light_Intensity_umols + 
#                          Far_Red_Light_Intensity_umols + (1 | Microgreen) + 
#                          (1 | Bed) + (1 | Position_in_Bed), 
#                        data = Reordered_Biomass_Data)
# 
# summary(Spectral_model_No_sow_or_bed) # Remove blue and far-red data
# rtemp=anova(Spectral_model_No_sow_or_bed)

# Reduced_Spectral_model <- lmer(Biomass_g ~ 1 + (1 | Microgreen) + 
#                                  (1 | Bed) + (1 | Position_in_Bed), 
#                                data = Reordered_Biomass_Data)

#summary(Reduced_Spectral_model)

#temp=anova(Reduced_Spectral_model ,Spectral_model_No_sow_or_bed) # Removing spectral data significantly changes the model

# Spectral_model_No_sow_bed_blue_or_far_red <- lmer(Biomass_g ~ 
#                                                     Red_Light_Intensity_umols + 
#                                                     (1 | Microgreen) + (1 | Bed) 
#                                                   + (1 | Position_in_Bed), 
#                                      data = Reordered_Biomass_Data) 
# 
# plot(Reordered_Biomass_Data$Biomass_g, predict(Spectral_model_No_sow_bed_blue_or_far_red))
# abline(a = 0, b = 1)
# 
# temp=anova(Reduced_Spectral_model, Spectral_model_No_sow_bed_blue_or_far_red)
# 
# # Blue or far_red are important
# 
# Spectral_model_No_sow_bed_or_blue <- lmer(Biomass_g ~ Red_Light_Intensity_umols + 
#                                        Far_Red_Light_Intensity_umols + (1 | Microgreen) + 
#                                        (1 | Bed) + (1 | Position_in_Bed), 
#                                      data = Reordered_Biomass_Data)
# 
# temp=anova(Reduced_Spectral_model, Spectral_model_No_sow_bed_or_blue)
# 
# Spectral_model_No_sow_bed_or_far_red <- lmer(Biomass_g ~ Red_Light_Intensity_umols + 
#                                                Blue_Light_Intensity_umols + (1 | Microgreen) + 
#                                             (1 | Bed) + (1 | Position_in_Bed), 
#                                           data = Reordered_Biomass_Data)
# 
# temp=anova(Reduced_Spectral_model, Spectral_model_No_sow_bed_or_far_red)
# 
# 
# 
# 
# 
# 
# Spectral_model_No_sow_bed_or_far_red_red <- lmer(Biomass_g ~ 1 + 
#                                                Blue_Light_Intensity_umols + (1 | Microgreen) + 
#                                                (1 | Bed) + (1 | Position_in_Bed), 
#                                              data = Reordered_Biomass_Data)
# 
# temp=anova(Reduced_Spectral_model ,Spectral_model_No_sow_bed_or_far_red_red)
# 
# Spectral_model_No_sow_bed_or_blue_red <- lmer(Biomass_g ~ 1 + 
#                                                    Far_Red_Light_Intensity_umols + (1 | Microgreen) + 
#                                                    (1 | Bed) + (1 | Position_in_Bed), 
#                                                  data = Reordered_Biomass_Data)
# 
# 
# temp=anova(Reduced_Spectral_model ,Spectral_model_No_sow_bed_or_blue_red)
# 

# MEMs on a species basis -------------------------------------------------

# Subset spectral data in to species 

Sunflower_Data <- Reordered_Biomass_Data[which(
  Reordered_Biomass_Data[,] == "Sunflower"),]

Radish_Data <- Reordered_Biomass_Data[which(
  Reordered_Biomass_Data[,] == "Radish"),]

Kale_CN_Data <- Reordered_Biomass_Data[which(
  Reordered_Biomass_Data[,] == "Kale(CN)"),]

Kale_Tozer_Data <- Reordered_Biomass_Data[which(
  Reordered_Biomass_Data[,] == "Kale(Tozer)"),]

# Start with sunflower MEMs

# Sunflower_model <- lmer(Biomass_g ~ Red_Light_Intensity_umols + 
#                          Blue_Light_Intensity_umols + 
#                          Far_Red_Light_Intensity_umols + (1 | Bed) + 
#                         (1 | Position_in_Bed), 
#                        data = Sunflower_Data)
# 
# summary(Sunflower_model)
# 
# rtemp=anova(Sunflower_model)

pvalues=matrix(1, nrow=3, ncol=4)
rownames(pvalues)=c('C-PAR', 'C-Qua', 'PAR-Qua') #c('PAR', 'R', 'B', 'FR')
colnames(pvalues)=c('Sunflower', 'Radish', 'Kale_CN', 'Kale_Tozer')

Reduced_Sunflower_Model <- lmer(Biomass_g ~ 1 + (1 | Bed) + 
                                  (1 | Position_in_Bed), 
                                data = Sunflower_Data)

#temp=anova(Reduced_Sunflower_Model, Sunflower_model)

Sunflower_model_PAR <- lmer(Biomass_g ~ Light_Intensity_umols + (1 | Bed) + 
                              (1 | Position_in_Bed), 
                            data = Sunflower_Data)

temp=anova(Reduced_Sunflower_Model, Sunflower_model_PAR)
pvalues['C-PAR', 'Sunflower']=temp$`Pr(>Chisq)`[2]

Sunflower_model_QUALITY <- lmer(Biomass_g ~ Red_Light_Intensity_umols + Blue_Light_Intensity_umols + Far_Red_Light_Intensity_umols + (1 | Bed) + 
                                  (1 | Position_in_Bed), 
                                data = Sunflower_Data)

temp=anova(Reduced_Sunflower_Model, Sunflower_model_QUALITY)
pvalues['C-Qua', 'Sunflower']=temp$`Pr(>Chisq)`[2]

temp=anova(Sunflower_model_PAR, Sunflower_model_QUALITY)
pvalues['PAR-Qua', 'Sunflower']=temp$`Pr(>Chisq)`[2]



Reduced_Radish_Model <- lmer(Biomass_g ~ 1 + (1 | Bed) + 
                               (1 | Position_in_Bed), 
                             data = Radish_Data)


Radish_model_PAR <- lmer(Biomass_g ~ Light_Intensity_umols + (1 | Bed) + 
                           (1 | Position_in_Bed), 
                         data = Radish_Data)

temp=anova(Reduced_Radish_Model, Radish_model_PAR)
pvalues['C-PAR', 'Radish']=temp$`Pr(>Chisq)`[2]

Radish_model_QUALITY <- lmer(Biomass_g ~ Red_Light_Intensity_umols + Blue_Light_Intensity_umols + Far_Red_Light_Intensity_umols + (1 | Bed) + 
                               (1 | Position_in_Bed), 
                             data = Radish_Data)

temp=anova(Reduced_Radish_Model, Radish_model_QUALITY)
pvalues['C-Qua', 'Radish']=temp$`Pr(>Chisq)`[2]

temp=anova(Radish_model_PAR, Radish_model_QUALITY)
pvalues['PAR-Qua', 'Radish']=temp$`Pr(>Chisq)`[2]


Reduced_Kale_CN_Model <- lmer(Biomass_g ~ 1 + (1 | Bed) + 
                                (1 | Position_in_Bed), 
                              data = Kale_CN_Data)


Kale_CN_model_PAR <- lmer(Biomass_g ~ Light_Intensity_umols + (1 | Bed) + 
                            (1 | Position_in_Bed), 
                          data = Kale_CN_Data)

temp=anova(Reduced_Kale_CN_Model, Kale_CN_model_PAR)
pvalues['C-PAR', 'Kale_CN']=temp$`Pr(>Chisq)`[2]

Kale_CN_model_QUALITY <- lmer(Biomass_g ~ Red_Light_Intensity_umols + Blue_Light_Intensity_umols + Far_Red_Light_Intensity_umols + (1 | Bed) + 
                                (1 | Position_in_Bed), 
                              data = Kale_CN_Data)

temp=anova(Reduced_Kale_CN_Model, Kale_CN_model_QUALITY)
pvalues['C-Qua', 'Kale_CN']=temp$`Pr(>Chisq)`[2]

temp=anova(Kale_CN_model_PAR, Kale_CN_model_QUALITY)
pvalues['PAR-Qua', 'Kale_CN']=temp$`Pr(>Chisq)`[2]


Reduced_Kale_Tozer_Model <- lmer(Biomass_g ~ 1 + (1 | Bed) + 
                                   (1 | Position_in_Bed), 
                                 data = Kale_Tozer_Data)


Kale_Tozer_model_PAR <- lmer(Biomass_g ~ Light_Intensity_umols + (1 | Bed) + 
                               (1 | Position_in_Bed), 
                             data = Kale_Tozer_Data)

temp=anova(Reduced_Kale_Tozer_Model, Kale_Tozer_model_PAR)
pvalues['C-PAR', 'Kale_Tozer']=temp$`Pr(>Chisq)`[2]

Kale_Tozer_model_QUALITY <- lmer(Biomass_g ~ Red_Light_Intensity_umols + Blue_Light_Intensity_umols + Far_Red_Light_Intensity_umols + (1 | Bed) + 
                                   (1 | Position_in_Bed), 
                                 data = Kale_Tozer_Data)

temp=anova(Reduced_Kale_Tozer_Model, Kale_Tozer_model_QUALITY)
pvalues['C-Qua', 'Kale_Tozer']=temp$`Pr(>Chisq)`[2]

#temp=anova(Kale_Tozer_model_PAR, Kale_Tozer_model_QUALITY)
#pvalues['PAR-Qua', 'Kale_Tozer']=temp$`Pr(>Chisq)`[2]

###adjust each row
apply(pvalues, 1, function(i){
  p.adjust(i, method='hochberg')
})

# 

###################################OLD


Sunflower_model_Blue <- lmer(Biomass_g ~ Blue_Light_Intensity_umols + (1 | Bed) + 
                               (1 | Position_in_Bed), 
                             data = Sunflower_Data)

temp=anova(Reduced_Sunflower_Model, Sunflower_model_Blue)
pvalues['B', 'Sunflower']=temp$`Pr(>Chisq)`[2]

Sunflower_model_Far_red <- lmer(Biomass_g ~ Far_Red_Light_Intensity_umols + (1 | Bed) + 
                                  (1 | Position_in_Bed), 
                                data = Sunflower_Data)

temp=anova(Reduced_Sunflower_Model, Sunflower_model_Far_red)
pvalues['FR', 'Sunflower']=temp$`Pr(>Chisq)`[2]
# Radish

# Radish_model <- lmer(Biomass_g ~ Red_Light_Intensity_umols + 
#                           Blue_Light_Intensity_umols + 
#                           Far_Red_Light_Intensity_umols + (1 | Bed) + 
#                           (1 | Position_in_Bed), 
#                         data = Radish_Data)
# 
# summary(Radish_model)
# 
# rtemp=anova(Radish_model) # Remove non significant fixed effects

# Radish_model_No_bed_or_bed_position <- lmer(Biomass_g ~ Red_Light_Intensity_umols + 
#                        Blue_Light_Intensity_umols + 
#                        Far_Red_Light_Intensity_umols + (1 | Position_in_Bed), 
#                      data = Radish_Data)

Reduced_Radish_Model <- lmer(Biomass_g ~ 1 + (1 | Bed) + (1 | Position_in_Bed), 
                             data = Radish_Data)

#temp=anova(Reduced_Radish_Model, Radish_model_No_bed_or_bed_position)

Radish_model_PAR <- lmer(Biomass_g ~ Light_Intensity_umols + (1 | Bed) + (1 | Position_in_Bed), 
                         data = Radish_Data)

temp=anova(Reduced_Radish_Model, Radish_model_PAR) # PAR is not important (but borderline)
pvalues['PAR', 'Radish']=temp$`Pr(>Chisq)`[2]

Radish_model_RED <- lmer(Biomass_g ~ Red_Light_Intensity_umols + (1 | Bed) + 
                           (1 | Position_in_Bed), 
                         data = Radish_Data)

temp=anova(Reduced_Radish_Model, Radish_model_RED) #red is important
pvalues['R', 'Radish']=temp$`Pr(>Chisq)`[2]

Radish_model_Blue <- lmer(Biomass_g ~ Blue_Light_Intensity_umols + (1 | Bed) + 
                            (1 | Position_in_Bed), 
                          data = Radish_Data)


temp=anova(Reduced_Radish_Model, Radish_model_Blue)
pvalues['B', 'Radish']=temp$`Pr(>Chisq)`[2]

Radish_model_Far_red <- lmer(Biomass_g ~ Far_Red_Light_Intensity_umols + (1 | Bed) + 
                               (1 | Position_in_Bed), 
                             data = Radish_Data)

temp=anova(Reduced_Radish_Model, Radish_model_Far_red)
pvalues['FR', 'Radish']=temp$`Pr(>Chisq)`[2]

# Tozer Kale

# Kale_Tozer_model <- lmer(Biomass_g ~ Red_Light_Intensity_umols + 
#                            Blue_Light_Intensity_umols + 
#                            Far_Red_Light_Intensity_umols + (1 | Bed) + 
#                            (1 | Position_in_Bed), 
#                          data = Kale_Tozer_Data)
# 
# summary(Kale_Tozer_model)
# 
# rtemp=anova(Kale_Tozer_model) # All  fixed effectss are significant 

Reduced_Kale_Tozer_Model <- lmer(Biomass_g ~ 1 + (1 | Bed) + 
                                   (1 | Position_in_Bed), 
                                 data = Kale_Tozer_Data)

#temp=anova(Reduced_Kale_Tozer_Model, Kale_Tozer_model)

Kale_Tozer_model_PAR <- lmer(Biomass_g ~ Light_Intensity_umols + (1 | Bed) + (1 | Position_in_Bed), 
                             data = Kale_Tozer_Data)

temp=anova(Reduced_Kale_Tozer_Model, Kale_Tozer_model_PAR) # PAR is important
pvalues['PAR', 'KaleT']=temp$`Pr(>Chisq)`[2]

Kale_Tozer_model_RED <- lmer(Biomass_g ~ Red_Light_Intensity_umols + (1 | Bed) + 
                               (1 | Position_in_Bed), 
                             data = Kale_Tozer_Data)

temp=anova(Reduced_Kale_Tozer_Model, Kale_Tozer_model_RED) #red is important
pvalues['R', 'KaleT']=temp$`Pr(>Chisq)`[2]

Kale_Tozer_model_Blue <- lmer(Biomass_g ~ Blue_Light_Intensity_umols + (1 | Bed) + 
                                (1 | Position_in_Bed), 
                              data = Kale_Tozer_Data)

temp=anova(Reduced_Kale_Tozer_Model, Kale_Tozer_model_Blue) 
pvalues['B', 'KaleT']=temp$`Pr(>Chisq)`[2]

Kale_Tozer_model_Far_red <- lmer(Biomass_g ~ Far_Red_Light_Intensity_umols + (1 | Bed) + 
                                   (1 | Position_in_Bed), 
                                 data = Kale_Tozer_Data)

temp=anova(Reduced_Kale_Tozer_Model, Kale_Tozer_model_Far_red) #far red is important
pvalues['FR', 'KaleT']=temp$`Pr(>Chisq)`[2]
# Kale CN

# Kale_CN_model <- lmer(Biomass_g ~ Red_Light_Intensity_umols + 
#                         Blue_Light_Intensity_umols + 
#                         Far_Red_Light_Intensity_umols + (1 | Bed) + 
#                         (1 | Position_in_Bed), 
#                       data = Kale_CN_Data)
# 
# summary(Kale_CN_model)
# 
# rtemp=anova(Kale_CN_model) # All fixed effects are not significant

# Kale_CN_model <- lmer(Biomass_g ~ Red_Light_Intensity_umols + 
#                                                Blue_Light_Intensity_umols + 
#                                                Far_Red_Light_Intensity_umols + 
#                                                (1 | Position_in_Bed) + (1 | Bed) + 
#                                                (1 | Bed_Position), 
#                                              data = Kale_CN_Data)

Reduced_Kale_CN_Model <- lmer(Biomass_g ~ 1 + (1 | Bed) + (1 | Position_in_Bed), 
                              data = Kale_CN_Data)

#temp=anova(Reduced_Kale_CN_Model, Kale_CN_model)

Kale_CN_model_PAR <- lmer(Biomass_g ~ Light_Intensity_umols + 
                            (1 | Position_in_Bed) + (1 | Bed), 
                          data = Kale_CN_Data)

temp=anova(Reduced_Kale_CN_Model, Kale_CN_model_PAR) # PAR is not important
pvalues['PAR', 'KaleC']=temp$`Pr(>Chisq)`[2]

Kale_CN_model_RED <- lmer(Biomass_g ~ Red_Light_Intensity_umols + (1 | Bed) + 
                            (1 | Position_in_Bed), 
                          data = Kale_CN_Data)

temp=anova(Reduced_Kale_CN_Model, Kale_CN_model_RED) # Red is important, but near singularity with bed and position in bed
pvalues['R', 'KaleC']=temp$`Pr(>Chisq)`[2]


Kale_CN_model_Blue <- lmer(Biomass_g ~ Blue_Light_Intensity_umols + (1 | Bed) + 
                             (1 | Position_in_Bed), 
                           data = Kale_CN_Data)

temp=anova(Reduced_Kale_CN_Model, Kale_CN_model_Blue) 
pvalues['B', 'KaleC']=temp$`Pr(>Chisq)`[2]
