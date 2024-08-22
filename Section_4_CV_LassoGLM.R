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

write.table(Reordered_Biomass_Data, file='Table_S2.txt')

# Lasso analysis using glmnet ------------------------------------------------


# Format the data

# Select only relevant columns

Biomass_Data <- Reordered_Biomass_Data[,-c(3,4,5,8, 10)]

# Remaining Factors need to be converted into dummies 

Biomass_Data_With_Dummies <- makeX(Biomass_Data)

# Response variable
y <- Biomass_Data$Biomass_g

# Matrix of predictors
x <- data.matrix(Biomass_Data_With_Dummies[,-5])  ###Note: need to add all the interaction terms!
#lightRatios
x=cbind(x, 'R_to_FR'=log10(x[,'Red_Light_Intensity_umols']/x[,'Far_Red_Light_Intensity_umols']))
x=cbind(x, 'R_to_B'=log10(x[,'Red_Light_Intensity_umols']/x[,'Blue_Light_Intensity_umols']))
x=cbind(x, 'FR_to_B'=log10(x[,'Far_Red_Light_Intensity_umols']/x[,'Blue_Light_Intensity_umols']))

id_micro=1:4
id_bed=5:8
id_intensity=9
id_pos=9:11+1
id_col=12:14+1
id_ratios=15:17+1
interactionColumns<- function(id1, id2, data){
  temp=lapply(id1, function(i){
    sapply(id2, function(j){
      data[,i]*data[,j]
    })
  })
  namesTemp=lapply(id1, function(i){
    sapply(id2, function(j){
      paste(colnames(data)[i], colnames(data)[j], sep='-')
    })
  })
  temp=do.call(cbind, temp)
  colnames(temp)=unlist(namesTemp)
  temp
}


#microgreen/position
interactionTerms=interactionColumns(id_micro, id_pos, x)
#microgreen/bed
interactionTerms=cbind(interactionTerms, interactionColumns(id_micro, id_bed, x))
#microgreen/colour
interactionTerms=cbind(interactionTerms, interactionColumns(id_micro, id_col, x))
#lightRatiosMicrogreen
interactionTerms=cbind(interactionTerms, interactionColumns(id_micro, id_ratios, x))
#light intensity-microgreen
interactionTerms=cbind(interactionTerms, interactionColumns(id_micro, id_intensity, x))

x_with_inter=cbind(x, interactionTerms)

CV_model <- cv.glmnet(x_with_inter, y, nfold=length(x[,1]), grouped=FALSE) #, alpha = 1)
lambdaB <- CV_model$lambda.min
png('sampleLambdaPlot.png')
plot(CV_model) 
dev.off()

Improved_model <- glmnet(x_with_inter,y, lambda = lambdaB) #alpha=1
round(coef(Improved_model) * 100)

coef(Improved_model)

# Lasso - plotting Coefs----------------------------------------------------------

improved_coefs <- coef(Improved_model)
improved_coefs <- data.frame(improved_coefs[-1,])
improved_coefs <- cbind(rownames(improved_coefs), improved_coefs)
colnames(improved_coefs)=c('coef_names', 'coefs')

#plot categories of coefficients

#positional (for supplement)
ids=improved_coefs[which(improved_coefs[,2]!=0 & grepl('Bed',improved_coefs[,1])),1]
ggplot(improved_coefs[ids,], aes(x = coef_names, y = coefs, fill = coef_names)) +
  geom_col() +
  theme_classic() +
  labs(x = "Coefficients", y = "Value", title = "Positional")+
  theme(axis.text.x = element_text(angle = 90)) +
  guides(fill = "none") +
  geom_hline( yintercept = 0 )


#light
ids=improved_coefs[which(improved_coefs[,2]!=0 & grepl('Light',improved_coefs[,1])),1]
ggplot(improved_coefs[ids,], aes(x = coef_names, y = coefs, fill = coef_names)) +
  geom_col() +
  theme_classic() +
  labs(x = "Coefficients", y = "Value", title = "Light")+
  theme(axis.text.x = element_text(angle = 90)) +
  guides(fill = "none") +
  geom_hline( yintercept = 0 )


#light ratios
ids=improved_coefs[which(improved_coefs[,2]!=0 & grepl('_to_',improved_coefs[,1])),1]
ggplot(improved_coefs[ids,], aes(x = coef_names, y = coefs, fill = coef_names)) +
  geom_col() +
  theme_classic() +
  labs(x = "Coefficients", y = "Value", title = "Ratios")+
  theme(axis.text.x = element_text(angle = 90)) +
  guides(fill = "none") +
  geom_hline( yintercept = 0 )











######################################################
####Now, we need to do leave one out cross validation:


lambdas=sapply(c(1:length(x[,1])), function(i){
  # make training set:
  train_x=x_with_inter[-i, ]
  train_y=y[-i]
  
  #find optimal lambda:
  CV_model <- cv.glmnet(train_x, train_y, nfold=length(x[,1])-1, grouped=FALSE) #, alpha = 1)
  lambdaB <- CV_model$lambda.min
})
png('histogramOfLambdas.png')
hist(lambdas)
dev.off()
models=lapply(c(1:length(x[,1])), function(i){
  # make training set:
  train_x=x_with_inter[-i, ]
  train_y=y[-i]
  
  #find optimal lambda:
  glmnet(train_x, train_y, lambda=lambdas[i]) #, alpha = 1)
})

coefs_combined=sapply(c(1:length(x[,1])), function(i){
  
  coef(models[[i]])[,1]
  
})

write.table(coefs_combined, file='allBeds_coefficients.txt', sep='\t')

#How consistent are the coefficients?2
#how many times does each coef appear
freq_coef=apply(coefs_combined, 1, function(i){length(which(i!=0))})
plot(sort(freq_coef))
abline(v=44)

pdf('Figure_S_variableSelection.pdf', width=3, height=4.5)
plot(sort(freq_coef)/256*100, xlab='variables (sorted)', ylab='% of models with selected variable')
abline(v=44)
abline(h=80)
dev.off()

png('Figure_S_variableSelection.png')
plot(sort(freq_coef)/256*100, xlab='variables (sorted)', ylab='% of models with selected variable')
abline(v=44)
abline(h=80)
dev.off()

used_at_least_once=names(which(freq_coef>=sort(freq_coef)[44]))
#Turn this into tidyverse-y version 
temp=lapply(used_at_least_once, function(i){
  data.frame('coef_names'=rep(i, length(coefs_combined[1,])), 'coefs'=coefs_combined[i,])
})
tidy_coefs <- do.call(rbind, temp)

#make nicer labels
words=c('MicrogreenKale',
        'Tozer',
        'CN',
        'MicrogreenRadish',
        'MicrogreenSunflower',
        ',Intensity_umols',
        'Position_in_Bed',
        'Red_Light_Intensity_umols',
        'Far_Red_Light_Intensity_umols',
        'Blue_Light_Intensity_umols',
        '_to_',
        '_')

replaceWith=c('Kale',
                'T',
                'C',
                'Radish',
                'Sunflower',
                '',
                '',
                'Red',
                'Far red',
                'Blue',
                ':',
                ',')

oldNames=tidy_coefs[,1]
for(i in c(1:length(replaceWith))){
  tidy_coefs[,1]=sub(words[i], replaceWith[i], tidy_coefs[,1])
}

tidy_coefs_all=tidy_coefs
#vioplot for ratios

cols=c('darkcyan','darksalmon', 'mediumorchid', 'hotpink2', 'grey')
names(cols)=unique(Reordered_Biomass_Data[,1], 'other')
names(cols)[1]='Kale\\(C\\)'
names(cols)[4]='Kale\\(T\\)'

tidy_coefs$col=rep('grey', length(tidy_coefs[,1]))
for(i in names(cols)){
  tidy_coefs$col[grep(i, tidy_coefs[,1])]=cols[i]
}

cols_gg=cols
names(cols_gg)=cols

ids=grep('R:FR', tidy_coefs[,1])
a=ggplot(tidy_coefs[ids,], aes(x = coef_names, y = coefs, colour=col)) +
  geom_boxplot() +
  scale_fill_manual(values=cols_gg)+
  scale_colour_manual(values=cols_gg)+
  theme_classic() +
  labs(x = "coefficients", y = "value (g/)", title = "R:FR")+
  theme(axis.text.x = element_text(angle = 90)) +
  guides(fill = "none") +
  geom_hline( yintercept = 0 )

ggsave('Figure_obs_2A.pdf', a, width=5, height=4)

#histograms of red, PAR, blue
ids=which(tidy_coefs[,1]=="Light,Intensity_umols")
par_hist=ggplot(tidy_coefs[ids,], aes(x = coefs, fill='darkmagenta')) +
  geom_histogram() +
  scale_fill_manual(values='darkmagenta', 'darkmagenta')+
  theme_classic() +
  guides(fill = "none") + xlim(NA, 0)+
  labs(x = "coefficient value (g/umols/area)", y = "count", title = "PAR reading")

ggsave('Figure_obs_2Bi.pdf', par_hist, width=5, height=2)

 
ids=which(tidy_coefs[,1]=="Blue")
blue_hist=ggplot(tidy_coefs[ids,], aes(x = coefs, fill='cornflowerblue')) +
  geom_histogram() +
  scale_fill_manual(values='cornflowerblue', 'cornflowerblue')+
  theme_classic() + xlim(0, NA)+
  guides(fill = "none") +
  labs(x = "coefficient value (g/(uW/cm2/nm))", y = "count", title = "Blue")

ggsave('Figure_obs_2Bii.pdf', blue_hist, width=5, height=2)


ids=grep('Red', tidy_coefs[,1])
red_hist=ggplot(tidy_coefs[ids,], aes(x = coefs, fill='red3')) +
  geom_histogram() +
  scale_fill_manual(values='red3', 'red3')+
  theme_classic() + xlim(NA, 0) +
  guides(fill = "none") +
  labs(x = "coefficient value (g/(uW/cm2/nm))", y = "count", title = "Red")

ggsave('Figure_obs_2Biii.pdf', red_hist, width=5, height=2)


#other ratios:
#c("Sunflower-R:B", 'Sunflower-FR:B", "Kale(C)-FR:B", "Kale(T)-Light,Intensity_umols")

#vioplot for lights
ids=grep('Light',oldNames)
#ids=ids[which(tidy_coefs[ids,1]!='Kale(C)-Red')]
ids=ids[which(tidy_coefs[ids,1]!='Red')]
ggplot(tidy_coefs[ids,], aes(x = coef_names, y = coefs, fill = coef_names)) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Coefficients", y = "Value (g/umol)", title = "Light Quality")+
  theme(axis.text.x = element_text(angle = 90)) +
  guides(fill = "none") +
  geom_hline( yintercept = 0 )

#hist(tidy_coefs[which(tidy_coefs[,1]=='Red'),2])
#hist(tidy_coefs[which(tidy_coefs[,1]=='Kale(C)-Red'),2])


predictions_combined=sapply(c(1:length(x[,1])), function(i){
  
   testSet=data.frame(t(x_with_inter[i,]))
  
   predict.glmnet(object = models[[i]], newx = testSet, s = models[[i]]$lambda, type = "response")
  
  
})

plot(y, predictions_combined, xlab='actual biomass (g)', ylab='predicted biomass (g)')
abline(c(0,1))



plot(log10(y), log10(predictions_combined))




#colour-code by species:
cols=c('darkcyan','darksalmon', 'mediumorchid', 'hotpink2')
names(cols)=unique(Reordered_Biomass_Data[,1])

pdf('Figure_obs_2C.pdf', width=5, height=5)
plot(y, predictions_combined, xlab='actual biomass (g)', pch=19, col=cols[Reordered_Biomass_Data[,1]],
     ylab='predicted biomass (g)', xlim=c(min(y), max(y)), ylim=c(min(y), max(y)))
abline(c(0,1))
legend(0, 38, names(cols), col=cols, bty='n', pch=19)
dev.off()

cor(y, predictions_combined)
cor.test(y, predictions_combined)$p.value

cors_by_species=sapply(unique(Reordered_Biomass_Data[,1]), function(i){
  ids=which(Reordered_Biomass_Data[,1]==i)
  plot(y[ids], predictions_combined[ids], xlab='log10(actual biomass (g))', pch=19, col=cols[Reordered_Biomass_Data[ids,1]],
       ylab='log10(predicted biomass (g))')
  abline(c(0,1))
  c(cor(y[ids], predictions_combined[ids]),
    cor.test(y[ids], predictions_combined[ids])$p.value,
    var(y[ids]))
  
})

plot(y, predictions_combined, xlab='log10(actual biomass (g))', pch=19, col=cols[Reordered_Biomass_Data[,1]],
     ylab='log10(predicted biomass (g))', xlim=c(min(y), max(y)), ylim=c(min(y), max(y)))
abline(c(0,1))
cor(y, predictions_combined)
cor.test(y, predictions_combined)$p.value
















# We need to check if red light still comes up if we look at the top two beds and 
#bottom two beds separately

####Now, we need to do leave one out cross validation:
y <- Biomass_Data$Biomass_g
topBeds=which((x[,'BedA']+x[,'BedB'])>0)


lambdas=sapply(1:length(topBeds), function(i){
  # make training set:
  train_x=x_with_inter[topBeds,][-i,]
  train_y=y[topBeds][-i]
  
  #find optimal lambda:
  CV_model <- cv.glmnet(train_x, train_y, nfold=length(x[,1])-1, grouped=FALSE) #, alpha = 1)
  lambdaB <- CV_model$lambda.min
})
hist(lambdas)
models=lapply(c(1:length(topBeds)), function(i){
  # make training set:
  train_x=x_with_inter[topBeds,][-i, ]
  train_y=y[topBeds][-i]
  
  #find optimal lambda:
  glmnet(train_x, train_y, lambda=lambdas[i]) #, alpha = 1)
})

coefs_combined=sapply(c(1:length(topBeds)), function(i){
  
  coef(models[[i]])[,1]
  
})


#How consistent are the coefficients?2
#how many times does each coef appear
freq_coef=apply(coefs_combined, 1, function(i){length(which(i!=0))})
plot(sort(freq_coef))
abline(v=57)

used_at_least_once=names(which(freq_coef>=sort(freq_coef)[57]))
#Turn this into tidyverse-y version 
temp=lapply(used_at_least_once, function(i){
  data.frame('coef_names'=rep(i, length(coefs_combined[1,])), 'coefs'=coefs_combined[i,])
})
tidy_coefs <- do.call(rbind, temp)


###try on bottom bed
tidy_coefs_top=tidy_coefs


####Now, we need to do leave one out cross validation:
y <- Biomass_Data$Biomass_g
topBeds=which((x[,'BedA']+x[,'BedB'])==0)


lambdas=sapply(1:length(topBeds), function(i){
  # make training set:
  train_x=x_with_inter[topBeds,][-i,]
  train_y=y[topBeds][-i]
  
  #find optimal lambda:
  CV_model <- cv.glmnet(train_x, train_y, nfold=length(x[,1])-1, grouped=FALSE) #, alpha = 1)
  lambdaB <- CV_model$lambda.min
})
hist(lambdas)
models=lapply(c(1:length(topBeds)), function(i){
  # make training set:
  train_x=x_with_inter[topBeds,][-i, ]
  train_y=y[topBeds][-i]
  
  #find optimal lambda:
  glmnet(train_x, train_y, lambda=lambdas[i]) #, alpha = 1)
})

coefs_combined=sapply(c(1:length(topBeds)), function(i){
  
  coef(models[[i]])[,1]
  
})


#How consistent are the coefficients?2
#how many times does each coef appear
freq_coef=apply(coefs_combined, 1, function(i){length(which(i!=0))})
plot(sort(freq_coef))
abline(v=52)

used_at_least_once=names(which(freq_coef>=sort(freq_coef)[57]))
#Turn this into tidyverse-y version 
temp=lapply(used_at_least_once, function(i){
  data.frame('coef_names'=rep(i, length(coefs_combined[1,])), 'coefs'=coefs_combined[i,])
})
tidy_coefs <- do.call(rbind, temp)


tidy_coefs_bottom=tidy_coefs


#compare  top bottom all

for(i in c(1:length(replaceWith))){
  tidy_coefs_top[,1]=sub(words[i], replaceWith[i], tidy_coefs_top[,1])
  tidy_coefs_bottom[,1]=sub(words[i], replaceWith[i], tidy_coefs_bottom[,1])
}

unique(tidy_coefs_all[,1])[which(unique(tidy_coefs_all[,1]) %in% unique(tidy_coefs_top[,1]))]
unique(tidy_coefs_all[,1])[which(unique(tidy_coefs_all[,1]) %in% unique(tidy_coefs_bottom[,1]))]

#are the values consistent too?
allThree=unique(tidy_coefs_all[,1])[which(unique(tidy_coefs_all[,1]) %in% unique(tidy_coefs_bottom[,1]) &
                                          unique(tidy_coefs_all[,1]) %in% unique(tidy_coefs_top[,1]))]

var="Light,Intensity_umols"
tidy_comparison_by_space=data.frame('space'=c(rep('all',256), rep('top', 128), rep('bottom', 128)),
           'coefs'=c(tidy_coefs_all[which(tidy_coefs_all[,1]==var),2],
                     tidy_coefs_top[which(tidy_coefs_top[,1]==var),2],
                     tidy_coefs_bottom[which(tidy_coefs_bottom[,1]==var),2]))

ggplot(tidy_comparison_by_space, aes(x=space, y=coefs))+geom_boxplot()

#"Kale(T)-R:FR"
var="Kale(T)-R:FR"
tidy_comparison_by_space=data.frame('space'=c(rep('all',256), rep('top', 128), rep('bottom', 128)),
                                    'coefs'=c(tidy_coefs_all[which(tidy_coefs_all[,1]==var),2],
                                              tidy_coefs_top[which(tidy_coefs_top[,1]==var),2],
                                              tidy_coefs_bottom[which(tidy_coefs_bottom[,1]==var),2]))

ggplot(tidy_comparison_by_space, aes(x=space, y=coefs))+geom_boxplot()



#"Radish-R:FR"
var="Radish-R:FR"
tidy_comparison_by_space=data.frame('col'=rep('darksalmon', 256*2), 'beds'=c(rep('all',256), rep('top', 128), rep('bottom', 128)),
                                    'coefficients'=c(tidy_coefs_all[which(tidy_coefs_all[,1]==var),2],
                                              tidy_coefs_top[which(tidy_coefs_top[,1]==var),2],
                                              tidy_coefs_bottom[which(tidy_coefs_bottom[,1]==var),2]))

temp1=ggplot(tidy_comparison_by_space, aes(x=beds, y=coefficients, colour=col))+
geom_boxplot()+
scale_fill_manual('white', 'white') +
  scale_colour_manual(values='darksalmon','darksalmon')+
  theme_classic() +
  ylim(c(-12, 12)) +
  labs(x = "beds", y = "coefficient value (g)", title = "Radish R:FR")+
  guides(fill = "none") +
  geom_hline( yintercept = 0 )

var="Kale(T)-R:FR"
tidy_comparison_by_space=data.frame('col'=rep('hotpink2', 256*2), 'beds'=c(rep('all',256), rep('top', 128), rep('bottom', 128)),
                                    'coefficients'=c(tidy_coefs_all[which(tidy_coefs_all[,1]==var),2],
                                                     tidy_coefs_top[which(tidy_coefs_top[,1]==var),2],
                                                     tidy_coefs_bottom[which(tidy_coefs_bottom[,1]==var),2]))

temp2=ggplot(tidy_comparison_by_space, aes(x=beds, y=coefficients, colour=col))+
  geom_boxplot()+
  scale_fill_manual('white', 'white') +
  scale_colour_manual(values='hotpink2','hotpink2')+
  theme_classic() +
  ylim(c(-40, 40)) +
  labs(x = "beds", y = "coefficient value (g)", title = "Kale 2 R:FR")+
  guides(fill = "none") +
  geom_hline( yintercept = 0 )

ggsave('Figure_obs_2Di.pdf', temp1, width=3, height=4)
ggsave('Figure_obs_2Dii.pdf', temp2, width=3, height=4)












#############################################
# Repeat this analysis on log10(biomass) to get rid of heteroscedasticity

####Now, we need to do leave one out cross validation:
y_old=y
y=log10(y)

lambdas=sapply(c(1:length(x[,1])), function(i){
  # make training set:
  train_x=x_with_inter[-i, ]
  train_y=y[-i]
  
  #find optimal lambda:
  CV_model <- cv.glmnet(train_x, train_y, nfold=length(x[,1])-1, grouped=FALSE) #, alpha = 1)
  lambdaB <- CV_model$lambda.min
})
hist(lambdas)
models=lapply(c(1:length(x[,1])), function(i){
  # make training set:
  train_x=x_with_inter[-i, ]
  train_y=y[-i]
  
  #find optimal lambda:
  glmnet(train_x, train_y, lambda=lambdas[i]) #, alpha = 1)
})

coefs_combined=sapply(c(1:length(x[,1])), function(i){
  
  coef(models[[i]])[,1]
  
})

#How consistent are the coefficients?2
#how many times does each coef appear
freq_coef=apply(coefs_combined, 1, function(i){length(which(i!=0))})
plot(sort(freq_coef))
abline(v=52)

used_at_least_once=names(which(freq_coef>=sort(freq_coef)[52]))
#Turn this into tidyverse-y version 
temp=lapply(used_at_least_once, function(i){
  data.frame('coef_names'=rep(i, length(coefs_combined[1,])), 'coefs'=coefs_combined[i,])
})
tidy_coefs <- do.call(rbind, temp)

#make nicer labels
words=c('MicrogreenKale',
        'Tozer',
        'CN',
        'MicrogreenRadish',
        'MicrogreenSunflower',
        ',Intensity_umols',
        'Position_in_Bed',
        'Red_Light_Intensity_umols',
        'Far_Red_Light_Intensity_umols',
        'Blue_Light_Intensity_umols',
        '_to_',
        '_')

replaceWith=c('Kale',
              'T',
              'C',
              'Radish',
              'Sunflower',
              '',
              '',
              'Red',
              'Far red',
              'Blue',
              ':',
              ',')

oldNames=tidy_coefs[,1]
for(i in c(1:length(replaceWith))){
  tidy_coefs[,1]=sub(words[i], replaceWith[i], tidy_coefs[,1])
}

#vioplot for ratios
ids=grep('_to_', oldNames)
ggplot(tidy_coefs[ids,], aes(x = coef_names, y = coefs, fill = coef_names)) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Coefficients", y = "Value (log(g))", title = "Ratios")+
  theme(axis.text.x = element_text(angle = 90)) +
  guides(fill = "none") +
  geom_hline( yintercept = 0 )


#vioplot for lights
ids=grep('Light',oldNames)
#ids=ids[which(tidy_coefs[ids,1]!='Kale(C)-Red')]
#ids=ids[which(tidy_coefs[ids,1]!='Red')]
ggplot(tidy_coefs[ids,], aes(x = coef_names, y = coefs, fill = coef_names)) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Coefficients", y = "Value (mg/umol)", title = "Light Quality")+
  theme(axis.text.x = element_text(angle = 90)) +
  guides(fill = "none") +
  geom_hline( yintercept = 0 )

#hist(tidy_coefs[which(tidy_coefs[,1]=='Red'),2])
#hist(tidy_coefs[which(tidy_coefs[,1]=='Kale(C)-Red'),2])


predictions_combined=sapply(c(1:length(x[,1])), function(i){
  
  testSet=data.frame(t(x_with_inter[i,]))
  
  predict.glmnet(object = models[[i]], newx = testSet, s = models[[i]]$lambda, type = "response")
  
  
})

#colour-code by species:
cols=c('darkcyan','darksalmon', 'mediumorchid', 'hotpink2')
names(cols)=unique(Reordered_Biomass_Data[,1])

plot(y, predictions_combined, xlab='log10(actual biomass (g))', pch=19, col=cols[Reordered_Biomass_Data[,1]],
     ylab='log10(predicted biomass (g))', xlim=c(min(y), max(y)), ylim=c(min(y), max(y)))
abline(c(0,1))

cor(y, predictions_combined)
cor.test(y, predictions_combined)$p.value

cors_by_species=sapply(unique(Reordered_Biomass_Data[,1]), function(i){
  ids=which(Reordered_Biomass_Data[,1]==i)
  plot(y[ids], predictions_combined[ids], xlab='log10(actual biomass (g))', pch=19, col=cols[Reordered_Biomass_Data[ids,1]],
       ylab='log10(predicted biomass (g))')
  abline(c(0,1))
  c(cor(y[ids], predictions_combined[ids]),
  cor.test(y[ids], predictions_combined[ids])$p.value,
  var(y[ids]))
  
})

plot(y, predictions_combined, xlab='log10(actual biomass (g))', pch=19, col=cols[Reordered_Biomass_Data[,1]],
     ylab='log10(predicted biomass (g))', xlim=c(min(y), max(y)), ylim=c(min(y), max(y)))
abline(c(0,1))
cor(y, predictions_combined)
cor.test(y, predictions_combined)$p.value

install.packages('lmtest')
library(lmtest)

predictions_combined[which(predictions_combined<0)]=0
a=lm(predictions_combined~y)
bptest(a)
shapiro.test(predictions_combined-y)
hist(predictions_combined-y)










































#separate out for each species
ids=which(x[,"MicrogreenKale(CN)"]==1)
plot(predictions_combined[ids], y[ids])
cor.test(predictions_combined[ids], y[ids])
shapiro.test(predictions_combined[ids]-y[ids])

ids=which(x[,"MicrogreenKale(Tozer)"]==1)
plot(predictions_combined[ids], y[ids])
cor.test(predictions_combined[ids], y[ids])
shapiro.test(predictions_combined[ids]-y[ids])


ids=which(x[,"MicrogreenRadish"]==1)
plot(predictions_combined[ids], y[ids])
cor.test(predictions_combined[ids], y[ids])
shapiro.test(predictions_combined[ids]-y[ids])

ids=which(x[,"MicrogreenSunflower"]==1)
plot(predictions_combined[ids], y[ids])
cor.test(predictions_combined[ids], y[ids])
shapiro.test(predictions_combined[ids]-y[ids])

plot(predictions_combined, y)

qqplot(predictions_combined, y)
hist(predictions_combined-y)
shapiro.test(predictions_combined-y)





Kale_T_Models <- lapply(Kale_Tozer_Data$Position, function(a){
model_data_indecies <- which(Kale_Tozer_Data$Position != a)
model_data <- Kale_Tozer_Data[model_data_indecies,]
y <- model_data$Biomass_g
Kale_Tozer_Data_With_Dummies <- makeX(model_data)
x <- data.matrix(Kale_Tozer_Data_With_Dummies[,5:15])
CV_model_Kale_Tozer <- cv.glmnet(x,y, alpha = 1)
lambdaT <- CV_model_Kale_Tozer$lambda.min
Improved_Kale_Tozer_Model <- glmnet(x,y,alpha = 1, lambda = lambdaT)
})

names(Kale_T_Models) <- Kale_Tozer_Data$Position

T_Kale_Coefs <- sapply(Kale_T_Models, function(a){
  as.matrix(coef(a))
})

T_Kale_Coefs

rownames(T_Kale_Coefs) <- c("intercept", "BedA","BedB","BedC","BedD","Light_Intensity_umols",
                            "Position_in_BedCentre","Position_in_BedCorner","Position_in_BedEdge",
                            "Peak_Red","Peak_Blue","Peak_Far_Red")


Tozer_Kale_Predictions <- sapply(Kale_Tozer_Data$Position, function(a){
  this_model <- Kale_T_Models[[a]]
  predictions <- predict.glmnet(object = this_model, newx = Kale_Tozer_Matrix, s = this_model$lambda, type = "response")
  names(predictions) <- Kale_Tozer_Data$Position
  predictions[a]
})

as.matrix(Tozer_Kale_Predictions)



# Lasso - Prediction plots ------------------------------------------------

#Get predictions


Kale_CN_Data <- Biomass_Data[which(Biomass_Data == "Kale(CN)"),]
Kale_CN_Matrix <- as.matrix(Kale_CN_Data_With_Dummies[,5:15])
Predicted_Kale_CN_Values <- c(predict.glmnet(object = Improved_Kale_CN_Model, newx = Kale_CN_Matrix,  s = lambdaKCN, type = "response"))

Kale_CN_Predicted_Data <- data.frame(Kale_CN_Data$Biomass_g , Predicted_Kale_CN_Values)
colnames(Kale_CN_Predicted_Data) <- c("Actual", "Predicted")

CN_plot <- ggplot(Kale_CN_Predicted_Data,aes(x = Predicted, y = Actual)) +
  geom_point() +
  theme_classic() +
  geom_abline(colour = "#7F00FF", size = 2) +
  labs(title = "Kale (CN)") +
  theme(text = element_text(size = 32),axis.line = element_line(colour = 'black', size = 2))

# Tozer

Kale_Tozer_Data <- Biomass_Data[which(Biomass_Data == "Kale(Tozer)"),]
Kale_Tozer_Matrix <- as.matrix(Kale_Tozer_Data_With_Dummies[,5:15])
Predicted_Kale_Tozer_Values <- c(predict.glmnet(object = Improved_Kale_Tozer_Model, newx = Kale_Tozer_Matrix,  s = lambdaT, type = "response"))

Kale_Tozer_Predicted_Data <- data.frame(Kale_Tozer_Data$Biomass_g , Predicted_Kale_Tozer_Values)
colnames(Kale_Tozer_Predicted_Data) <- c("Actual", "Predicted")

T_Plot <- ggplot(Kale_Tozer_Predicted_Data,aes(x = Predicted, y = Actual)) +
  geom_point() +
  theme_classic() +
  geom_abline(colour = "#7F00FF", size = 2) +
  labs(title = "Kale (Tozer)") +
  theme(text = element_text(size = 32),axis.line = element_line(colour = 'black', size = 2))

T_Plot

#Sunflower
Sunflower_Data <- Biomass_Data[which(Biomass_Data == "Sunflower"),]
Sunflower_Matrix <- as.matrix(Sunflower_Data_With_Dummies[,5:15])
Predicted_Sunflower_Values <- c(predict.glmnet(object = Improved_Sunflower_Model, newx = Sunflower_Matrix,  s = lambdaS, type = "response"))

Sunflower_Predicted_Data <- data.frame(Sunflower_Data$Biomass_g , Predicted_Sunflower_Values)
colnames(Sunflower_Predicted_Data) <- c("Actual", "Predicted")

S_Plot <- ggplot(Sunflower_Predicted_Data,aes(x = Predicted, y = Actual)) +
  geom_point() +
  theme_classic() +
  geom_abline(colour = "#7F00FF", size = 2) +
  labs(title = "Sunflower") +
  theme(text = element_text(size = 32),axis.line = element_line(colour = 'black', size = 2))

#Radish

Radish_Data <- Biomass_Data[which(Biomass_Data == "Radish"),]
Radish_Matrix <- as.matrix(Radish_Data_With_Dummies[,5:15])
Predicted_Radish_Values <- c(predict.glmnet(object = Improved_Radish_Model, newx = Radish_Matrix,  s = lambdaR, type = "response"))

Radish_Predicted_Data <- data.frame(Radish_Data$Biomass_g , Predicted_Radish_Values)
colnames(Radish_Predicted_Data) <- c("Actual", "Predicted")

R_plot <- ggplot(Radish_Predicted_Data,aes(x = Predicted, y = Actual)) +
  geom_point() +
  theme_classic() +
  geom_abline(colour = "#7F00FF", size = 2) +
  labs(title = "Radish") +
  theme(text = element_text(size = 32),axis.line = element_line(colour = 'black', size = 2))

# Bring it all together 

Radish_Predicted_Data$Microgreen <- rep("Radish" , nrow(Radish_Predicted_Data))
Sunflower_Predicted_Data$Microgreen <- rep("Sunflower" , nrow(Sunflower_Predicted_Data))
Kale_Tozer_Predicted_Data$Microgreen <- rep("Kale (Tozer)" , nrow(Kale_Tozer_Predicted_Data))
Kale_CN_Predicted_Data$Microgreen <- rep("Kale (CN)" , nrow(Kale_CN_Predicted_Data))

Predicted_data <- rbind(Radish_Predicted_Data, Sunflower_Predicted_Data, Kale_Tozer_Predicted_Data, Kale_CN_Predicted_Data)

ggplot(Predicted_data,aes(x = Predicted, y = Actual)) +
  geom_point() +
  theme_classic() +
  geom_abline(colour = "#7F00FF", size = 2) +
  facet_wrap(~Microgreen, scales = "free") +
  theme(text = element_text(size = 32),axis.line = element_line(colour = 'black', size = 2)) +
  labs(x = "Predicted biomass (g)", y = "Actual biomass (g)")

multiplot(S_Plot,R_plot,CN_plot,T_Plot, ncol = 2)

# Removing overall intensity ----------------------------------------------

Biomass_Data_No_overall_light_intensity <- Biomass_Data[,c(1,2,3,4,5,6,8,9,10,11,12)]

########################################################Look at this.
# Convert from watts to umols ---------------------------------------------

as_quantum_mol(s.e.irrad = (Biomass_Data$Peak_Red/ 100), w.length =  ) 

# Leave one out cross validation ------------------------------------------

Sunflower_Data <- Biomass_Data[which(Biomass_Data[,1] == "Sunflower"),-c(8,1)]
Radish_Data <- Biomass_Data[which(Biomass_Data[,1] == "Radish"),-c(8,1)]
Kale_Tozer_Data <- Biomass_Data[which(Biomass_Data[,1] == "Kale(Tozer)"),-c(8,1)]
Kale_CN_Data <- Biomass_Data[which(Biomass_Data[,1] == "Kale(CN)"),-c(8,1)]

Kale_Tozer_Data$Position <- paste(Kale_Tozer_Data$X, Kale_Tozer_Data$Y, Kale_Tozer_Data$Bed)

Kale_T_Models <- lapply(Kale_Tozer_Data$Position, function(a){
  model_data_indecies <- which(Kale_Tozer_Data$Position != a)
  model_data <- Kale_Tozer_Data[model_data_indecies,]
  y <- model_data$Biomass_g
  Kale_Tozer_Data_With_Dummies <- makeX(model_data)
  x <- data.matrix(Kale_Tozer_Data_With_Dummies[,5:15])
  CV_model_Kale_Tozer <- cv.glmnet(x,y, alpha = 1)
  lambdaT <- CV_model_Kale_Tozer$lambda.min
  Improved_Kale_Tozer_Model <- glmnet(x,y,alpha = 1, lambda = lambdaT)
})

names(Kale_T_Models) <- Kale_Tozer_Data$Position

T_Kale_Coefs <- sapply(Kale_T_Models, function(a){
  as.matrix(coef(a))
})

T_Kale_Coefs

rownames(T_Kale_Coefs) <- c("intercept", "BedA","BedB","BedC","BedD","Light_Intensity_umols",
                            "Position_in_BedCentre","Position_in_BedCorner","Position_in_BedEdge",
                            "Peak_Red","Peak_Blue","Peak_Far_Red")


Tozer_Kale_Predictions <- sapply(Kale_Tozer_Data$Position, function(a){
  this_model <- Kale_T_Models[[a]]
  predictions <- predict.glmnet(object = this_model, newx = Kale_Tozer_Matrix, s = this_model$lambda, type = "response")
  names(predictions) <- Kale_Tozer_Data$Position
  predictions[a]
})

as.matrix(Tozer_Kale_Predictions)

# Normalising light intensities -------------------------------------------
#NOt necessary, based on how method works.

Biomass_Data_No_overall_light_intensity$Peak_Red <- (Biomass_Data$Peak_Red - min(Biomass_Data$Peak_Red)) / max(Biomass_Data$Peak_Red - min(Biomass_Data$Peak_Red))

Biomass_Data_No_overall_light_intensity$Peak_Blue <-  (Biomass_Data$Peak_Blue - min(Biomass_Data$Peak_Blue)) / max(Biomass_Data$Peak_Blue - min(Biomass_Data$Peak_Blue))

Biomass_Data_No_overall_light_intensity$Peak_Far_Red <- (Biomass_Data$Peak_Far_Red - min(Biomass_Data$Peak_Far_Red)) / max(Biomass_Data$Peak_Far_Red - min(Biomass_Data$Peak_Far_Red))

Sunflower_Data <- Biomass_Data_No_overall_light_intensity[which(Biomass_Data_No_overall_light_intensity[,1] == "Sunflower"),-c(7,1)]
Radish_Data <- Biomass_Data_No_overall_light_intensity[which(Biomass_Data_No_overall_light_intensity[,1] == "Radish"),-c(7,1)]
Kale_Tozer_Data <- Biomass_Data_No_overall_light_intensity[which(Biomass_Data_No_overall_light_intensity[,1] == "Kale(Tozer)"),-c(7,1)]
Kale_CN_Data <- Biomass_Data_No_overall_light_intensity[which(Biomass_Data_No_overall_light_intensity[,1] == "Kale(CN)"),-c(7,1)]

# Start with sunflower

y <- Sunflower_Data$Biomass_g
Sunflower_Data_With_Dummies <- makeX(Sunflower_Data)
x <- data.matrix(Sunflower_Data_With_Dummies[,5:14])

CV_model_Sunflower <- cv.glmnet(x,y, alpha = 1)
lambdaS <- CV_model_Sunflower$lambda.min
plot(CV_model_Sunflower)

Improved_Sunflower_Model <- glmnet(x,y,alpha = 1, lambda = lambdaS)

Sunflower_coefs <- coef(Improved_Sunflower_Model)

rm(Sunflower_Data, y, x, 
   CV_model_Sunflower, CV_model)

# Kale CN

y <- Kale_CN_Data$Biomass_g
Kale_CN_Data_With_Dummies <- makeX(Kale_CN_Data)
x <- data.matrix(Kale_CN_Data_With_Dummies[,5:14])

CV_model_Kale_CN <- cv.glmnet(x,y, alpha = 1)
lambdaKCN <- CV_model_Kale_CN$lambda.min
plot(CV_model_Kale_CN)

Improved_Kale_CN_Model <- glmnet(x,y,alpha = 1, lambda = lambdaKCN)

Kale_CN_coefs <- coef(Improved_Kale_CN_Model)

rm( Kale_CN_Data, y, x, 
    CV_model_Kale_CN, CV_model_Kale_CN)

# Kale(Tozer)

y <- Kale_Tozer_Data$Biomass_g
Kale_Tozer_Data_With_Dummies <- makeX(Kale_Tozer_Data)
x <- data.matrix(Kale_Tozer_Data_With_Dummies[,5:14])

CV_model_Kale_Tozer <- cv.glmnet(x,y, alpha = 1)
lambdaT <- CV_model_Kale_Tozer$lambda.min
plot(CV_model_Kale_Tozer)

Improved_Kale_Tozer_Model <- glmnet(x,y,alpha = 1, lambda = lambdaT)

Kale_Tozer_coefs <-  coef(Improved_Kale_Tozer_Model)

rm(Kale_Tozer_Data, y, x, 
   CV_model_Kale_Tozer, CV_model_Kale_Tozer)

# Finally Radish

y <- Radish_Data$Biomass_g
Radish_Data_With_Dummies <- makeX(Radish_Data)
x <- data.matrix(Radish_Data_With_Dummies[,5:14])

CV_model_Radish <- cv.glmnet(x,y, alpha = 1)
lambdaR <- CV_model_Radish$lambda.min
plot(CV_model_Radish)

Improved_Radish_Model <- glmnet(x,y,alpha = 1, lambda = lambdaR)

Radish_Coefs <- matrix(coef(Improved_Radish_Model))
coef(Improved_Radish_Model)

rm(Radish_Data, y, x, 
   CV_model_Radish, CV_model_Radish)

coef_names <- c("BedA","BedB","BedC","BedD",
                "Position_in_BedCentre","Position_in_BedCorner","Position_in_BedEdge",
                "Peak_Red","Peak_Blue","Peak_Far_Red")

Radish_Coefs <- matrix(coef(Improved_Radish_Model))
Radish_Coefs <- Radish_Coefs[-1,]
print(Radish_Coefs)
Radish_coefs <- data.frame(coef_names,Radish_Coefs)
Radish_Plot <- ggplot(Radish_coefs, aes(x = coef_names, y = Radish_Coefs, fill = coef_names)) +
  geom_col() +
  theme_classic() +
  labs(x = "Coefficients", y = "Value", title = "Radish")+
  theme(axis.text.x = element_text(angle = 90)) +
  guides(fill = "none") +
  geom_hline( yintercept = 0 )

Sunflower_Coefs <- matrix(coef(Improved_Sunflower_Model))
Sunflower_Coefs <- Sunflower_Coefs[-1,]
Sunflower_coefs <- data.frame(coef_names,Sunflower_Coefs)
Sunflower_Plot <- ggplot(Sunflower_coefs, aes(x = coef_names, y = Sunflower_Coefs, fill = coef_names)) +
  geom_col() +
  theme_classic() +
  labs(x = "Coefficients", y = "Value", title = "Sunflower")+
  theme(axis.text.x = element_text(angle = 90)) +
  guides(fill = "none") +
  geom_hline( yintercept = 0 )


Kale_Tozer_Coefs <- matrix(coef(Improved_Kale_Tozer_Model))
Kale_Tozer_Coefs <- Kale_Tozer_Coefs[-1,]
Kale_Tozer_coefs <- data.frame(coef_names,Kale_Tozer_Coefs)
Kale_Tozer_Plot <- ggplot(Kale_Tozer_coefs, aes(x = coef_names, y = Kale_Tozer_Coefs, fill = coef_names)) +
  geom_col() +
  theme_classic() +
  labs(x = "Coefficients", y = "Value", title = "Kale (Tozer)") +
  theme(axis.text.x = element_text(angle = 90)) +
  guides(fill = "none") +
  geom_hline( yintercept = 0 )


Kale_CN_Coefs <- matrix(coef(Improved_Kale_CN_Model))
Kale_CN_Coefs <- Kale_CN_Coefs[-1,]
Kale_CN_coefs <- data.frame(coef_names,Kale_CN_Coefs)
Kale_CN_Plot <- ggplot(Kale_CN_coefs, aes(x = coef_names, y = Kale_CN_Coefs, fill = coef_names)) +
  geom_col() +
  theme_classic() +
  labs(x = "Coefficients", y = "Value", title = "Kale (CN)", fill = "Coefficients") +
  theme(axis.text.x = element_text(angle = 90)) +
  guides(fill = "none") +
  geom_hline( yintercept = 0 )


multiplot(Sunflower_Plot, Radish_Plot, Kale_CN_Plot, Kale_Tozer_Plot, cols = 2)

colnames(Kale_CN_coefs) <- c("Coefficient", "Value")
Kale_CN_coefs$Microgreen <- rep("Kale (CN)", nrow(Kale_CN_coefs))
colnames(Kale_Tozer_coefs) <- c("Coefficient", "Value")
Kale_Tozer_coefs$Microgreen <- rep("Kale (Tozer)", nrow(Kale_Tozer_coefs))
colnames(Sunflower_coefs) <- c("Coefficient", "Value")
Sunflower_coefs$Microgreen <- rep("Sunflower", nrow(Sunflower_coefs))
colnames(Radish_coefs) <- c("Coefficient", "Value")
Radish_coefs$Microgreen <- rep("Radish", nrow(Radish_coefs))

Coefficients <- rbind(Kale_CN_coefs, Kale_Tozer_coefs, Radish_coefs, Sunflower_coefs)

colours <- c("#5D3FD3","#BF40BF","pink","purple","#DA70D6","#7F00FF","#CF9FFF", "#E0B0FF","#CBC3E3","#7C2278","#AA98A9")

Coef_Plot <- ggplot(Coefficients, aes(x = Coefficient, y = Value, fill = Coefficient)) +
  geom_col() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline( yintercept = 0, size = 1.5 ) +
  facet_wrap(~Microgreen) +
  scale_fill_manual(values = colours, labels = c("Bed A", "Bed B", "Bed C", "Bed D", "Blue light intensity", "Far red light intensity", "Red light intensity","Centre of bed", "Corner of bed", "Edge of bed")) +
  scale_x_discrete(labels = c("Bed A", "Bed B", "Bed C", "Bed D", "Blue light intensity", "Far red light intensity", "Red light intensity","Centre of bed", "Corner of bed", "Edge of bed")) + 
  theme(axis.text.x = element_text(vjust = 0.5, size = 20),text = element_text(size = 32),axis.line = element_line(colour = 'black', size = 2))

Coef_Plot

# CN

Kale_CN_Data <- Biomass_Data_No_overall_light_intensity[which(Biomass_Data_No_overall_light_intensity == "Kale(CN)"),]
Kale_Matrix <- as.matrix(Kale_CN_Data_With_Dummies[,5:14])
Predicted_Kale_CN_Values <- c(predict.glmnet(object = Improved_Kale_CN_Model, newx = Kale_Matrix,  s = lambdaKCN, type = "response"))

Kale_CN_Predicted_Data <- data.frame(Kale_CN_Data$Biomass_g , Predicted_Kale_CN_Values)
colnames(Kale_CN_Predicted_Data) <- c("Actual", "Predicted")

CN_plot <- ggplot(Kale_CN_Predicted_Data,aes(x = Predicted, y = Actual)) +
  geom_point() +
  theme_classic() +
  geom_abline(colour = "#7F00FF", size = 2) +
  labs(title = "Kale (CN)") +
  theme(text = element_text(size = 32),axis.line = element_line(colour = 'black', size = 2))

# Tozer

Kale_Tozer_Data <- Biomass_Data_No_overall_light_intensity[which(Biomass_Data_No_overall_light_intensity == "Kale(Tozer)"),]
Kale_Matrix <- as.matrix(Kale_Tozer_Data_With_Dummies[,5:14])
Predicted_Kale_Tozer_Values <- c(predict.glmnet(object = Improved_Kale_Tozer_Model, newx = Kale_Matrix,  s = lambdaT, type = "response"))

Kale_Tozer_Predicted_Data <- data.frame(Kale_Tozer_Data$Biomass_g , Predicted_Kale_Tozer_Values)
colnames(Kale_Tozer_Predicted_Data) <- c("Actual", "Predicted")

T_Plot <- ggplot(Kale_Tozer_Predicted_Data,aes(x = Predicted, y = Actual)) +
  geom_point() +
  theme_classic() +
  geom_abline(colour = "#7F00FF", size = 2) +
  labs(title = "Kale (Tozer)") +
  theme(text = element_text(size = 32),axis.line = element_line(colour = 'black', size = 2))

T_Plot

#Sunflower
Sunflower_Data <- Biomass_Data_No_overall_light_intensity[which(Biomass_Data_No_overall_light_intensity == "Sunflower"),]
Sunflower_Matrix <- as.matrix(Sunflower_Data_With_Dummies[,5:14])
Predicted_Sunflower_Values <- c(predict.glmnet(object = Improved_Sunflower_Model, newx = Sunflower_Matrix,  s = lambdaS, type = "response"))

Sunflower_Predicted_Data <- data.frame(Sunflower_Data$Biomass_g , Predicted_Sunflower_Values)
colnames(Sunflower_Predicted_Data) <- c("Actual", "Predicted")

S_Plot <- ggplot(Sunflower_Predicted_Data,aes(x = Predicted, y = Actual)) +
  geom_point() +
  theme_classic() +
  geom_abline(colour = "#7F00FF", size = 2) +
  labs(title = "Sunflower") +
  theme(text = element_text(size = 32),axis.line = element_line(colour = 'black', size = 2))

#Radish

Radish_Data <- Biomass_Data_No_overall_light_intensity[which(Biomass_Data_No_overall_light_intensity == "Radish"),]
Radish_Matrix <- as.matrix(Radish_Data_With_Dummies[,5:14])
Predicted_Radish_Values <- c(predict.glmnet(object = Improved_Radish_Model, newx = Radish_Matrix,  s = lambdaR, type = "response"))

Radish_Predicted_Data <- data.frame(Radish_Data$Biomass_g , Predicted_Radish_Values)
colnames(Radish_Predicted_Data) <- c("Actual", "Predicted")

R_plot <- ggplot(Radish_Predicted_Data,aes(x = Predicted, y = Actual)) +
  geom_point() +
  theme_classic() +
  geom_abline(colour = "#7F00FF", size = 2) +
  labs(title = "Radish") +
  theme(text = element_text(size = 32),axis.line = element_line(colour = 'black', size = 2))

# Bring it all together 

Radish_Predicted_Data$Microgreen <- rep("Radish" , nrow(Radish_Predicted_Data))
Sunflower_Predicted_Data$Microgreen <- rep("Sunflower" , nrow(Sunflower_Predicted_Data))
Kale_Tozer_Predicted_Data$Microgreen <- rep("Kale (Tozer)" , nrow(Kale_Tozer_Predicted_Data))
Kale_CN_Predicted_Data$Microgreen <- rep("Kale (CN)" , nrow(Kale_CN_Predicted_Data))

Predicted_data <- rbind(Radish_Predicted_Data, Sunflower_Predicted_Data, Kale_Tozer_Predicted_Data, Kale_CN_Predicted_Data)

ggplot(Predicted_data,aes(x = Predicted, y = Actual)) +
  geom_point() +
  theme_classic() +
  geom_abline(colour = "#7F00FF", size = 2) +
  facet_wrap(~Microgreen, scales = "free") +
  theme(text = element_text(size = 32),axis.line = element_line(colour = 'black', size = 2)) +
  labs(x = "Predicted biomass (g)", y = "Actual biomass (g)")

multiplot(S_Plot,R_plot,CN_plot,T_Plot, ncol = 2)

# Determining R squared values

rm(CN_plot,Coef_Plot,Coefficients)


# First for sunflower

sst <- sum((Sunflower_Data$Biomass_g - mean(Sunflower_Data$Biomass_g))^2)
sse <- sum((Predicted_Sunflower_Values - Sunflower_Data$Biomass_g)^2)

rsq_Sunflower <- 1 - sse/sst
rsq_Sunflower

# Radish Data

sst <- sum((Radish_Data$Biomass_g - mean(Radish_Data$Biomass_g))^2)
sse <- sum((Predicted_Radish_Values - Radish_Data$Biomass_g)^2)

rsq_Radish <- 1 - sse/sst
rsq_Radish

# Kale CN

sst <- sum((Kale_CN_Data$Biomass_g - mean(Kale_CN_Data$Biomass_g))^2)
sse <- sum((Predicted_Kale_CN_Values - Kale_CN_Data$Biomass_g)^2)

rsq_Kale_CN <- 1 - sse/sst
rsq_Kale_CN

# Kale Tozer 

sst <- sum((Kale_Tozer_Data$Biomass_g - mean(Kale_Tozer_Data$Biomass_g))^2)
sse <- sum((Predicted_Kale_Tozer_Values - Kale_Tozer_Data$Biomass_g)^2)

rsq_Kale_Tozer <- 1 - sse/sst
rsq_Kale_Tozer


