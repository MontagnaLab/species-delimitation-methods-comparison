##################################################################################
##### Factors affecting the efficiency of molecular species delimitation in  #####
##### a species-rich insect family ###############################################

## Giulia Magoga, Diego Fontaneto & Matteo Montagna

##################################################################################

## R code to generate the analyses, using R 3.6.3
## Last update: Milano, 4 January 2021

##################################################################################

# Loading R package --------------------------------------------------------------

library('lme4')
library('car')
library('MuMIn')

# Working directory --------------------------------------------------------------

# check the folder

getwd()

# change it if needed
setwd("") #<- Change me

# Loading the database -----------------------------------------------------------

giulia <- read.csv("table2_plus_revision.csv", header=T, sep="\t")
dim(giulia)
names(giulia)
summary(giulia)

# Analyses -----------------------------------------------------------------------

# 1. overall assessment of predictors on efficiency, without interactions --------

model1 <- lme4::glmer(cbind(match, failure) ~ as.factor(Method) +
							 scale(Diff_id_2e3) + 
							 scale(Diff_id_3) + 
							 scale(Dist_median_medianL) + 
							 as.factor(Liv_tass) + 
							 scale(N_hap_meanL) + 
							 (1|Dataset), 
							 data=giulia, 
							 family=binomial)

# check model fit

plot(model1)

# print results in a table

GLMER1 <- capture.output(car::Anova(model1))
write("Results of the GLMER:", "results 1 global GLMER.csv")
write(GLMER1, "results 1 global GLMER.csv", append=T)

# obtain relative importance value and update the table
options(na.action="na.fail")
summary(MuMIn::model.avg(get.models(dredge(model1), cumsum(weight)<=.999)))
sw(MuMIn::model.avg(get.models(dredge(model1), cumsum(weight)<=.999)))
GLMER2 <- capture.output(sw(MuMIn::model.avg(get.models(dredge(model1),
				cumsum(weight)<=.999))))
options(na.action="na.omit")
write(paste("","###","Results of the Multimodel Averaging:", sep="\n"),  
		"results 1 global GLMER.csv", 
		append=T)
write(GLMER2, "results 1 global GLMER.csv", 
		append=T)


# 2. overall assessment of predictors on efficiency, interacting with methods ----

model3 <- lme4::glmer(cbind(match, failure) ~ as.factor(Method) + 
							scale(Diff_id_2e3) +  
							scale(Diff_id_3) + 
							scale(Dist_median_medianL) + 
							as.factor(Liv_tass) + 
							scale(N_hap_meanL) + 
							scale(Diff_id_2e3):as.factor(Method) + 
							scale(Diff_id_3):as.factor(Method) + 
							scale(Dist_median_medianL):as.factor(Method) + 
							as.factor(Liv_tass):as.factor(Method) + 
							scale(N_hap_meanL):as.factor(Method) + 
							(1|Dataset), 
							data=giulia, 
							family=binomial)

# check model fit

plot(model3)

# print results in a table

GLMER3 <- capture.output(car::Anova(model3))
write(paste("","###","Results of the GLMER with interactions:", sep="\n"), 
		"results 1 global GLMER.csv", append=T)
write(GLMER3, "results 1 global GLMER.csv", append=T)
options(na.action="na.fail")
GLMER4 <- capture.output(sw(MuMIn::model.avg(get.models(dredge(model1),
				cumsum(weight)<=.999))))
options(na.action="na.omit")
write(paste("","###","Results of Multimodel Averaging with interactions:",
			sep="\n"),  
		"results 1 global GLMER.csv", 
		append=T)
write(GLMER4, "results 1 global GLMER.csv", 
		append=T)

# interactions are significant: keep them in the final model


# 3. separate analyses for each delimitation method ------------------------------

# overall analyses with results saved in tables

write("Results of the GLM by method:", "results glm by method.csv")
for(i in levels(giulia$Method)) {
    model.0 <- glm(cbind(match, failure) ~ scale(Diff_id_3) + 
    					scale(Diff_id_2e3) + 
    					scale(Dist_median_medianL) + 
    					as.factor(Liv_tass) + 
    					scale(N_hap_meanL), 
    					data=giulia[giulia$Method==as.character(i),], 
    					family=binomial)
#    pdf(paste("2.0 Model plot",as.character(i),".pdf"))
#    par(mfrow=c(2,2))
#	    plot(model.0, main=as.character(i))
#    dev.off()
    output <- capture.output(car::Anova(model.0))
    write(paste("","###",as.character(i), sep="\n"), 
    	"results glm by method.csv", append=T)
    write(output, "results glm by method.csv", append=T)
}

# relative importance values for each model

ss=split(giulia, giulia$Method) 

abgd=ss[[1]]
ASAP=ss[[2]]
gmyc=ss[[3]]
loc_min=ss[[4]]
mptp=ss[[5]]
thr=ss[[6]]

model_a <- glm(cbind(match, failure) ~ scale(Diff_id_3) + 
					scale(Diff_id_2e3) + 
					scale(Dist_median_medianL) + 
					as.factor(Liv_tass) + 
					scale(N_hap_meanL), 
					data=abgd, 
					family=binomial)
options(na.action="na.fail")
sw(MuMIn::model.avg(get.models(dredge(model_a), cumsum(weight)<=.999)))
options(na.action="na.omit")

model_as <- glm(cbind(match, failure) ~ scale(Diff_id_3) + 
					scale(Diff_id_2e3) + 
					scale(Dist_median_medianL) + 
					as.factor(Liv_tass) + 
					scale(N_hap_meanL), 
					data=ASAP, 
					family=binomial)
options(na.action="na.fail")
sw(MuMIn::model.avg(get.models(dredge(model_as), cumsum(weight)<=.999)))
options(na.action="na.omit")

model_g <- glm(cbind(match, failure) ~ scale(Diff_id_3) + 
					scale(Diff_id_2e3) + 
					scale(Dist_median_medianL) + 
					as.factor(Liv_tass) + 
					scale(N_hap_meanL), 
					data=gmyc, 
					family=binomial)
options(na.action="na.fail")
sw(MuMIn::model.avg(get.models(dredge(model_g), cumsum(weight)<=.999)))
options(na.action="na.omit")

model_l <- glm(cbind(match, failure) ~ scale(Diff_id_3) + 
					scale(Diff_id_2e3) + 
					scale(Dist_median_medianL) + 
					as.factor(Liv_tass) + 
					scale(N_hap_meanL), 
					data=loc_min, 
					family=binomial)
options(na.action="na.fail")
sw(MuMIn::model.avg(get.models(dredge(model_l), cumsum(weight)<=.999)))
options(na.action="na.omit")

model_m <- glm(cbind(match, failure) ~ scale(Diff_id_3) + 
					scale(Diff_id_2e3) + 
					scale(Dist_median_medianL) + 
					as.factor(Liv_tass) + 
					scale(N_hap_meanL), 
					data=mptp, 
					family=binomial)
options(na.action="na.fail")
sw(MuMIn::model.avg(get.models(dredge(model_m), cumsum(weight)<=.999)))
options(na.action="na.omit")

model_t <- glm(cbind(match, failure) ~ scale(Diff_id_3) + 
					scale(Diff_id_2e3) + 
					scale(Dist_median_medianL) + 
					as.factor(Liv_tass) + 
					scale(N_hap_meanL), 
					data=thr, 
					family=binomial)
options(na.action="na.fail")
sw(MuMIn::model.avg(get.models(dredge(model_t), cumsum(weight)<=.999)))
options(na.action="na.omit")
