# STAT 309 Final Project
# Yildiz Derinkok
# 11 May 2023


# Read in Data
karolinska <- read.table("karolinska.txt", header = TRUE)
head(karolinska)
dim(karolinska)
dim(karolinska[karolinska$HighVolDiagHosp == 1, ]) #number of units in active


# Clean Variables (date)
typeof(karolinska$NumericDateOfDiagnosis)
karolinska$DateOfDiagnosis <- as.Date(karolinska$DateOfDiagnosis, format = '%Y-%m-%d') #make date into a date type variable
karolinska$NumericDateOfDiagnosis <- as.numeric(karolinska$DateOfDiagnosis, format = '%Y-%m-%d') #create new variable with date as a numeric variable
#install.packages("lubridate")
library("lubridate")


#Create new variables:
karolinska$MonthOfDiagnosis <- month(karolinska$DateOfDiagnosis) #create new variable for month

karolinska$GoodWeather <- 0 #variable for if weather is good (higher temp, more sun) based on month variable (may thru sept)
karolinska$GoodWeather[karolinska$MonthOfDiagnosis >= 5 & karolinska$MonthOfDiagnosis <= 9] <- 1

karolinska$OldRuralMale <- 0 #variable for old rural male
karolinska$OldRuralMale[karolinska$AgeAtDiagnosis >= 65 & karolinska$FromRuralArea == 1 & karolinska$Male == 1] <- 1



# 1. Considering high- v. low-volume diagnosing hospital as the patient-selected treatment, create subclasses or matches to establish covariate balance.
#do observational study design on diagnosing hospital
# z: diagnosing hospital
# w: treating hospital
# design study on z


#Removing the outcome variable, the w, and character date variable:
obs <- subset(karolinska, select = -c(YearsSurvivingAfterDiagnosis, HighVolTreatHosp, DateOfDiagnosis)) 
head(obs)

#Covariate balance with no subclasses:
round(apply(obs[obs$HighVolDiagHosp == 1, ], 2, "mean"), 2) #means for active treatment
round(apply(X = obs[obs$HighVolDiagHosp == 0, ], MARGIN = 2, FUN = "mean"), 2) #means for control


#Trying different subclassification designs:

###### Propensity scores using regression on all covariates ######
logregmodel1 <- glm(obs$HighVolDiagHosp ~ ., data = obs, family = binomial)
#summary(logregmodel1)
l1.pscores <- logregmodel1$fitted.values #vector of estimated propensity scores
summary(l1.pscores)
sort(l1.pscores)
boxplot(l1.pscores~obs$HighVolDiagHosp, main = "Propensity Scores by Treatment Group")


#Discarding units:
summary(l1.pscores[obs$HighVolDiagHosp == 1])
obs[l1.pscores < 0.17, ] #all of these people are old rural men who were diagnosed in good weather. none of them have equivalents in active, so i am discarding
obs$l1.pscores <- l1.pscores
obs <- obs[obs$l1.pscores >= 0.17, ]
dim(obs)


# Create subclasses by propensity score for l1
obs$l1.pclass <- 1 # >= .17 $ < .27
obs$l1.pclass[obs$l1.pscores > 0.27 & obs$l1.pscores <= 0.50] <- 2
obs$l1.pclass[obs$l1.pscores > 0.50 & obs$l1.pscores <= 0.61] <- 3
obs$l1.pclass[obs$l1.pscores > 0.61] <- 4

table(obs$l1.pclass[obs$HighVolDiagHosp ==0]) #num units in each subclass for control
table(obs$l1.pclass[obs$HighVolDiagHosp ==1]) #num units in each subclass for active


#For each subclass created in 5, calculate the mean within each treatment group for each covariate:
l1.s1.a <- round(apply(obs[obs$HighVolDiagHosp == 1 & obs$l1.pclass == 1, ], 2, "mean"), 2)
l1.s1.c <- round(apply(obs[obs$HighVolDiagHosp == 0 & obs$l1.pclass == 1, ], 2, "mean"), 2)
l1.s2.a <- round(apply(obs[obs$HighVolDiagHosp == 1 & obs$l1.pclass == 2, ], 2, "mean"), 2)
l1.s2.c <- round(apply(obs[obs$HighVolDiagHosp == 0 & obs$l1.pclass == 2, ], 2, "mean"), 2)
l1.s3.a <- round(apply(obs[obs$HighVolDiagHosp == 1 & obs$l1.pclass == 3, ], 2, "mean"), 2)
l1.s3.c <- round(apply(obs[obs$HighVolDiagHosp == 0 & obs$l1.pclass == 3, ], 2, "mean"), 2)
l1.s4.a <- round(apply(obs[obs$HighVolDiagHosp == 1 & obs$l1.pclass == 4, ], 2, "mean"), 2)
l1.s4.c <- round(apply(obs[obs$HighVolDiagHosp == 0 & obs$l1.pclass == 4, ], 2, "mean"), 2)

l1.weights <- table(obs$l1.pclass)/dim(obs)[1] #vector of weights by num units in each subclass

#Get weighted means for each covariate:
l1.wtd.ave.a <- rep(NA, 7) #num of covariates
l1.wtd.ave.c <- rep(NA, 7)

for(ii in 2:8){
	l1.wtd.ave.a[ii-1] <- l1.s1.a[ii]*l1.weights[1] + l1.s2.a[ii]*l1.weights[2] + 
	l1.s3.a[ii]*l1.weights[3] + l1.s4.a[ii]*l1.weights[4] #multiply the means of each 	covariate by the proportion of people in active (weights)
	
	l1.wtd.ave.c[ii-1] <- l1.s1.c[ii]*l1.weights[1] + l1.s2.c[ii]*l1.weights[2] + 		l1.s3.c[ii]*l1.weights[3] + l1.s4.c[ii]*l1.weights[4]
}

#The weighted covariate means (use to check balance):
round(l1.wtd.ave.a, 2) 
round(l1.wtd.ave.c, 2)




###### Propensity scores using regression on some covariates ######
obs <- subset(karolinska, select = -c(YearsSurvivingAfterDiagnosis, HighVolTreatHosp, DateOfDiagnosis)) #remove old subclassification and add back removed data
head(obs)
dim(obs)

logregmodel2 <- glm(obs$HighVolDiagHosp ~ AgeAtDiagnosis + FromRuralArea + Male + NumericDateOfDiagnosis, data = obs, family = binomial)
l2.pscores <- logregmodel2$fitted.values #vector of estimated propensity scores
summary(l2.pscores)
sort(l2.pscores)
boxplot(l2.pscores~obs$HighVolDiagHosp, main = "Propensity Scores by Treatment Group")


#Discarding units:
summary(l2.pscores[obs$HighVolDiagHosp == 1])
obs[l2.pscores < 0.20, ] #all of these people are old rural men. none of them have equivalents in active, so i am discarding
obs$l2.pscores <- l2.pscores
obs <- obs[obs$l2.pscores >= 0.20, ]
dim(obs)


# Create subclasses by propensity score for l2
obs$l2.pclass <- 1
obs$l2.pclass[obs$l2.pscores > 0.30 & obs$l2.pscores <= 0.50] <- 2
obs$l2.pclass[obs$l2.pscores > 0.50 & obs$l2.pscores <= 0.66] <- 3
obs$l2.pclass[obs$l2.pscores > 0.66] <- 4

table(obs$l2.pclass[obs$HighVolDiagHosp ==0]) #num units in each subclass for control
table(obs$l2.pclass[obs$HighVolDiagHosp ==1]) #num units in each subclass for active


#For each subclass created in 5, calculate the mean within each treatment group for each covariate:
l2.s1.a <- round(apply(obs[obs$HighVolDiagHosp == 1 & obs$l2.pclass == 1, ], 2, "mean"), 2)
l2.s1.c <- round(apply(obs[obs$HighVolDiagHosp == 0 & obs$l2.pclass == 1, ], 2, "mean"), 2)
l2.s2.a <- round(apply(obs[obs$HighVolDiagHosp == 1 & obs$l2.pclass == 2, ], 2, "mean"), 2)
l2.s2.c <- round(apply(obs[obs$HighVolDiagHosp == 0 & obs$l2.pclass == 2, ], 2, "mean"), 2)
l2.s3.a <- round(apply(obs[obs$HighVolDiagHosp == 1 & obs$l2.pclass == 3, ], 2, "mean"), 2)
l2.s3.c <- round(apply(obs[obs$HighVolDiagHosp == 0 & obs$l2.pclass == 3, ], 2, "mean"), 2)
l2.s4.a <- round(apply(obs[obs$HighVolDiagHosp == 1 & obs$l2.pclass == 4, ], 2, "mean"), 2)
l2.s4.c <- round(apply(obs[obs$HighVolDiagHosp == 0 & obs$l2.pclass == 4, ], 2, "mean"), 2)

l2.weights <- table(obs$l2.pclass)/dim(obs)[1] #vector of weights by num units in each subclass

#Get weighted means for each covariate:
l2.wtd.ave.a <- rep(NA, 7) #num of covariates
l2.wtd.ave.c <- rep(NA, 7)

for(ii in 2:8){
	l2.wtd.ave.a[ii-1] <- l2.s1.a[ii]*l2.weights[1] + l2.s2.a[ii]*l2.weights[2] + 
	l2.s3.a[ii]*l2.weights[3] + l2.s4.a[ii]*l2.weights[4] #multiply the means of each 	covariate by the proportion of people in active (weights)
	
	l2.wtd.ave.c[ii-1] <- l2.s1.c[ii]*l2.weights[1] + l2.s2.c[ii]*l2.weights[2] + 		l2.s3.c[ii]*l2.weights[3] + l2.s4.c[ii]*l2.weights[4]
}

#The weighted covariate means (use to check balance):
round(l2.wtd.ave.a, 2) 
round(l2.wtd.ave.c, 2)




###### Propensity scores using tree on all covariates ######
obs <- subset(karolinska, select = -c(YearsSurvivingAfterDiagnosis, HighVolTreatHosp, DateOfDiagnosis)) #remove old subclassification and add back removed datahead(obs)

library(party)
t1.output <- ctree(obs$HighVolDiagHosp ~., data = obs)
plot(t1.output)
t1.pscores <- predict(output) #vector of estimated propensity scores

# Create subclasses by propensity score for t1
obs$t1.pclass <- 1
obs$t1.pclass[t1.pscores > 0.50] <- 2

table(obs$t1.pclass[obs$HighVolDiagHosp ==0]) #num units in each subclass for control
table(obs$t1.pclass[obs$HighVolDiagHosp ==1]) #num units in each subclass for active


#For each subclass created in 5, calculate the mean within each treatment group for each covariate:
t1.s1.a <- round(apply(obs[obs$HighVolDiagHosp == 1 & obs$t1.pclass == 1, ], 2, "mean"), 2)
t1.s1.c <- round(apply(obs[obs$HighVolDiagHosp == 0 & obs$t1.pclass == 1, ], 2, "mean"), 2)
t1.s2.a <- round(apply(obs[obs$HighVolDiagHosp == 1 & obs$t1.pclass == 2, ], 2, "mean"), 2)
t1.s2.c <- round(apply(obs[obs$HighVolDiagHosp == 0 & obs$t1.pclass == 2, ], 2, "mean"), 2)

t1.weights <- table(obs$t1.pclass)/dim(obs)[1] #vector of weights by num units in each subclass

#Get weighted means for each covariate:
t1.wtd.ave.a <- rep(NA, 7) #num of covariates
t1.wtd.ave.c <- rep(NA, 7)

for(ii in 2:8){
	t1.wtd.ave.a[ii-1] <- t1.s1.a[ii]*t1.weights[1] + t1.s2.a[ii]*t1.weights[2] 		#multiply the means of each covariate by the proportion of people in active 		(weights)
	
	t1.wtd.ave.c[ii-1] <- t1.s1.c[ii]*t1.weights[1] + t1.s2.c[ii]*t1.weights[2]
}

#The weighted covariate means (use to check balance):
round(t1.wtd.ave.a, 2) 
round(t1.wtd.ave.c, 2)



#Chose subclass design L1:
head(karolinska)
karolinska$pscores <- l1.pscores
karolinska <- karolinska[karolinska$pscores >= 0.17, ]
dim(karolinska)

karolinska$pclass <- 1
karolinska$pclass[karolinska$pscores > 0.27 & karolinska$pscores <= 0.50] <- 2
karolinska$pclass[karolinska$pscores > 0.50 & karolinska$pscores <= 0.61] <- 3
karolinska$pclass[karolinska$pscores > 0.61] <- 4

table(karolinska$pclass[karolinska$HighVolDiagHosp ==0]) #num units in each subclass
table(karolinska$pclass[karolinska$HighVolDiagHosp ==1])

head(karolinska)

#2. Using diagnosing hospital as an instrument, estimate the average causal effect of high- v. low-volume treating hospital on survival and calculate an interval for this estimate.
#itt = neyman estimate for diagnosing hospital


#Clean outcome variable YearsSurvivingAfterDiagnosis by creating 3 binary variables
typeof(karolinska$YearsSurvivingAfterDiagnosis)
karolinska$Survive1 <- 1
karolinska$Survive1[karolinska$YearsSurvivingAfterDiagnosis == "1"] <- 0 #survived longer than 1 year or no
karolinska$Survive2_4 <- 0
karolinska$Survive2_4[karolinska$YearsSurvivingAfterDiagnosis == "2-4"] <- 1 #survived 2-4 years
karolinska$Survive5 <- 0
karolinska$Survive5[karolinska$YearsSurvivingAfterDiagnosis == "5+"] <- 1 #survived more than 5 years or no

weights <- table(karolinska$pclass)/dim(karolinska)[1]


#Find itt and pc for each subclass to calculate ittc for each subclass. Then calculate overall ittc by weighing by subclass:
itt.s1 <- mean(karolinska$Survive1[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 1]) - mean(karolinska$Survive1[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 1])
itt.s2 <- mean(karolinska$Survive1[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 2]) - mean(karolinska$Survive1[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 2])
itt.s3 <- mean(karolinska$Survive1[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 3]) - mean(karolinska$Survive1[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 3])
itt.s4 <- mean(karolinska$Survive1[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 4]) - mean(karolinska$Survive1[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 4])


#Get proportion of people who are always-takers. So (assuming no defiers) proportion of people who were diagnosed in low and treated in high (divide by people in subclass diagnosed at low) (out of the people diagnosed in low, the ones who got treated in high)
pa.s1 <- dim(karolinska[karolinska$HighVolDiagHosp == 0 & karolinska$HighVolTreatHosp == 1 & karolinska$pclass == 1, ])[1]/dim(karolinska[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 1, ])[1]
pa.s2 <- dim(karolinska[karolinska$HighVolDiagHosp == 0 & karolinska$HighVolTreatHosp == 1 & karolinska$pclass == 2, ])[1]/dim(karolinska[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 2, ])[1]
pa.s3 <- dim(karolinska[karolinska$HighVolDiagHosp == 0 & karolinska$HighVolTreatHosp == 1 & karolinska$pclass == 3, ])[1]/dim(karolinska[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 3, ])[1]
pa.s4 <- dim(karolinska[karolinska$HighVolDiagHosp == 0 & karolinska$HighVolTreatHosp == 1 & karolinska$pclass == 4, ])[1]/dim(karolinska[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 4, ])[1]


#Proportion of people who are compliers or always-takers. So proportion of people who were diagnosed and treated in high
pca.s1 <- dim(karolinska[karolinska$HighVolDiagHosp == 1 & karolinska$HighVolTreatHosp == 1 & karolinska$pclass == 1, ])[1]/dim(karolinska[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 1, ])[1]
pca.s2 <- dim(karolinska[karolinska$HighVolDiagHosp == 1 & karolinska$HighVolTreatHosp == 1 & karolinska$pclass == 2, ])[1]/dim(karolinska[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 2, ])[1]
pca.s3 <- dim(karolinska[karolinska$HighVolDiagHosp == 1 & karolinska$HighVolTreatHosp == 1 & karolinska$pclass == 3, ])[1]/dim(karolinska[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 3, ])[1]
pca.s4 <- dim(karolinska[karolinska$HighVolDiagHosp == 1 & karolinska$HighVolTreatHosp == 1 & karolinska$pclass == 4, ])[1]/dim(karolinska[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 4, ])[1]


#Get proportion who are compliers (diagnosed and treated in high-volume) by subtracting proportion of always-takers from proportion of comppliers or always-takers:
pc.s1 <- pca.s1 - pa.s1
pc.s2 <- pca.s2 - pa.s2
pc.s3 <- pca.s3 - pa.s3
pc.s4 <- pca.s4 - pa.s4

ittc.s1 <- itt.s1/pc.s1
ittc.s2 <- itt.s2/pc.s2
ittc.s3 <- itt.s3/pc.s3
ittc.s4 <- itt.s4/pc.s4

ittc.s.vector <- c(ittc.s1, ittc.s2, ittc.s3, ittc.s4)
ittc.total <- sum(ittc.s.vector*weights) #Causal effect of being treated at a high vs low volume treating hospital


#Calculate variances for ittc in each subclass
# z: diagnosing hospital
# w: treating hospital
# y: outcome
# var(ittc) = (1/pc)^2*var(itt) + (itt/pc)^2*var(pc) - 2(1/pc)(itt/pc)cov(itt, pc) 
# var(itt) = (Syt^2/Nt)(Syc^2/Nc)
# var(pc) = (Swt^2/Nt)(Swc^2/Nc)
#Syt^2 = sample variance of y among z = 1
#Swt^2 = sample variance of w among z = 1

itt.var.s1 <- (var(karolinska$Survive1[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 1])/dim(karolinska[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 1, ])[1] + var(karolinska$Survive1[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 1])/dim(karolinska[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 1, ])[1])
pc.var.s1 <- (var(karolinska$HighVolTreatHosp[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 1])/dim(karolinska[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 1, ])[1] + var(karolinska$HighVolTreatHosp[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 1])/dim(karolinska[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 1, ])[1])
cov.s1 <- (1/dim(karolinska[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 1, ])[1])*cov(karolinska$Survive1[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 1], karolinska$HighVolTreatHosp[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 1]) + (1/dim(karolinska[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 1, ])[1])*cov(karolinska$Survive1[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 1], karolinska$HighVolTreatHosp[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 1])

ittc.var.s1 <- (1/pc.s1)^2 * itt.var.s1 + (itt.s1/pc.s1)^2 * pc.var.s1 - 2*(1/pc.s1)*(itt.s1/pc.s1)*cov.s1 #not surprised if its big. small sample size in subclasses. also used method of moments


itt.var.s2 <- (var(karolinska$Survive1[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 2])/dim(karolinska[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 2, ])[1] + var(karolinska$Survive1[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 2])/dim(karolinska[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 2, ])[1])
pc.var.s2 <- (var(karolinska$HighVolTreatHosp[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 2])/dim(karolinska[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 2, ])[1] + var(karolinska$HighVolTreatHosp[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 2])/dim(karolinska[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 2, ])[1])
cov.s2 <- (1/dim(karolinska[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 2, ])[1])*cov(karolinska$Survive1[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 2], karolinska$HighVolTreatHosp[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 2]) + (1/dim(karolinska[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 2, ])[1])*cov(karolinska$Survive1[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 2], karolinska$HighVolTreatHosp[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 2])

ittc.var.s2 <- (1/pc.s2)^2 * itt.var.s2 + (itt.s2/pc.s2)^2 * pc.var.s2 - 2*(1/pc.s2)*(itt.s2/pc.s2)*cov.s2 


itt.var.s3 <- (var(karolinska$Survive1[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 3])/dim(karolinska[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 3, ])[1] + var(karolinska$Survive1[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 3])/dim(karolinska[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 3, ])[1])
pc.var.s3 <- (var(karolinska$HighVolTreatHosp[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 3])/dim(karolinska[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 3, ])[1] + var(karolinska$HighVolTreatHosp[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 3])/dim(karolinska[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 3, ])[1])
cov.s3 <- (1/dim(karolinska[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 3, ])[1])*cov(karolinska$Survive1[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 3], karolinska$HighVolTreatHosp[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 3]) + (1/dim(karolinska[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 3, ])[1])*cov(karolinska$Survive1[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 3], karolinska$HighVolTreatHosp[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 3])

ittc.var.s3 <- (1/pc.s3)^2 * itt.var.s3 + (itt.s3/pc.s3)^2 * pc.var.s3 - 2*(1/pc.s3)*(itt.s3/pc.s3)*cov.s3


itt.var.s4 <- (var(karolinska$Survive1[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 4])/dim(karolinska[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 4, ])[1] + var(karolinska$Survive1[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 4])/dim(karolinska[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 4, ])[1])
pc.var.s4 <- (var(karolinska$HighVolTreatHosp[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 4])/dim(karolinska[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 4, ])[1] + var(karolinska$HighVolTreatHosp[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 4])/dim(karolinska[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 4, ])[1])
cov.s4 <- (1/dim(karolinska[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 4, ])[1])*cov(karolinska$Survive1[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 4], karolinska$HighVolTreatHosp[karolinska$HighVolDiagHosp == 1 & karolinska$pclass == 4]) + (1/dim(karolinska[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 4, ])[1])*cov(karolinska$Survive1[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 4], karolinska$HighVolTreatHosp[karolinska$HighVolDiagHosp == 0 & karolinska$pclass == 4])

ittc.var.s4 <- (1/pc.s4)^2 * itt.var.s4 + (itt.s4/pc.s4)^2 * pc.var.s4 - 2*(1/pc.s4)*(itt.s4/pc.s4)*cov.s4


# Get total variance
var.s.vector <- c(ittc.var.s1, ittc.var.s2, ittc.var.s3, ittc.var.s4)
ittc.var.total <- sum(var.s.vector*weights^2)

# 95% CI:
ittc.total - 2*sqrt(ittc.var.total) #lower bound
ittc.total + 2*sqrt(ittc.var.total) #upper bound
#ittc has to be btw -1 and 1. so if interval bigger than that not helpful

