# Package loading and Working Directory -----------------------------------
library(lme4) # building mixed effects models
library(lmerTest) # statistical tests for MEMs 
setwd("C:/Users/de656/OneDrive - University of York/Research/VerticalFarm_observationalStudy/")

##########################
#  In this script we evaluate whether PAR is significantly associated with biomass for any of our species, with position as a random effect.
##########################

####Load data
Microgreen_Biomass_Data <- read.table("Farm_Biomass_Data.txt", 
                                      header = T, sep = "")


Reduced_All_Model <- lmer(Biomass_g ~ 1 + Microgreen+ (1 | Bed : Microgreen) + (1 | Position_in_Bed : Microgreen), data=Microgreen_Biomass_Data)


summary(Reduced_All_Model)

ranova(Reduced_All_Model)

plot(Reduced_All_Model)
plot(predict(Reduced_All_Model), Microgreen_Biomass_Data$Biomass_g)

PAR_All_Model <- lmer(Biomass_g ~ 1 + Light_Intensity_umols * Microgreen + (1 | Bed : Microgreen) + (1 | Position_in_Bed : Microgreen), data=Microgreen_Biomass_Data)

summary(PAR_All_Model)

ranova(PAR_All_Model)

anova(Reduced_All_Model, PAR_All_Model)

temp=anova(Reduced_All_Model, PAR_All_Model)











##################################################################################

####################Possibly remove
# MEMs on a species basis -------------------------------------------------

# Subset spectral data in to species 

Sunflower_Data <- Microgreen_Biomass_Data[which(
  Microgreen_Biomass_Data[,] == "Sunflower"),]

Radish_Data <- Microgreen_Biomass_Data[which(
  Microgreen_Biomass_Data[,] == "Radish"),]

Kale_CN_Data <- Microgreen_Biomass_Data[which(
  Microgreen_Biomass_Data[,] == "Kale(CN)"),]

Kale_Tozer_Data <- Microgreen_Biomass_Data[which(
  Microgreen_Biomass_Data[,] == "Kale(Tozer)"),]


############### Compare PAR with baseline model


pvals=c()
bics=matrix(NA, nrow=4, ncol=2)
rownames(bics)=c('Sunflower', 'Radish', 'Kale_CN', 'Kale_Tozer')

# Start with sunflower MEMs

Reduced_Sunflower_Model <- lmer(Biomass_g ~ 1 + (1 | Bed) + 
                                  (1 | Position_in_Bed), 
                                data = Sunflower_Data)

Sunflower_model_PAR <- lmer(Biomass_g ~ Light_Intensity_umols + (1 | Bed) + 
                              (1 | Position_in_Bed), 
                            data = Sunflower_Data)

temp=anova(Reduced_Sunflower_Model, Sunflower_model_PAR)
pvals['Sunflower']=temp$`Pr(>Chisq)`[2]
bics['Sunflower',]=temp$BIC

Reduced_Radish_Model <- lmer(Biomass_g ~ 1 + (1 | Bed) + 
                               (1 | Position_in_Bed), 
                             data = Radish_Data)

Radish_model_PAR <- lmer(Biomass_g ~ Light_Intensity_umols + (1 | Bed) + 
                           (1 | Position_in_Bed), 
                         data = Radish_Data)

temp=anova(Reduced_Radish_Model, Radish_model_PAR)
pvals['Radish']=temp$`Pr(>Chisq)`[2]
bics['Radish',]=temp$BIC

Reduced_Kale_CN_Model <- lmer(Biomass_g ~ 1 + (1 | Bed) + 
                                (1 | Position_in_Bed), 
                              data = Kale_CN_Data)

Kale_CN_model_PAR <- lmer(Biomass_g ~ Light_Intensity_umols + (1 | Bed) + 
                            (1 | Position_in_Bed), 
                          data = Kale_CN_Data)

temp=anova(Reduced_Kale_CN_Model, Kale_CN_model_PAR)
pvals['Kale_CN']=temp$`Pr(>Chisq)`[2]
bics['Kale_CN',]=temp$BIC

Reduced_Kale_Tozer_Model <- lmer(Biomass_g ~ 1 + (1 | Bed) + 
                                   (1 | Position_in_Bed), 
                                 data = Kale_Tozer_Data)

Kale_Tozer_model_PAR <- lmer(Biomass_g ~ Light_Intensity_umols + (1 | Bed) + 
                               (1 | Position_in_Bed), 
                             data = Kale_Tozer_Data)

temp=anova(Reduced_Kale_Tozer_Model, Kale_Tozer_model_PAR)
pvals['Kale_Tozer']=temp$`Pr(>Chisq)`[2]
bics['Kale_Tozer',]=temp$BIC

pvals

p.adjust(pvals)
bics