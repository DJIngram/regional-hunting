############## Script to structure the bayesian mixed-effects models on wild meat hunting

# clear environment
rm(list=ls())

# load packages
library(raster)
library(rgdal)
library(ggplot2)
library(sf)
library(car)
library(gdalUtils)
library(maptools)
library(reshape2)
library(stringr)
library(ggpubr)
library(brms)
library(GGally)
library(DHARMa)
library(DHARMa.helpers)
library(emmeans)

# Load in mean hunter site-level data for analyses, which includes site level
# information and spatially extracted environmental variables used in the analyses
studyTotYes<-read.csv("input-local-file-path.csv")


################# Variables include:
# referenceID = a unique identifier for the reference of the source data
# siteID = a unique identifier for the site data was collected in
# studyID = a unique identifier for the study
# country = the country name
# access10kaeaCrop = site-level estimate of travel time to nearest town
# forestCorrect20km = site-level estimate of mean forest cover within 20km of the site
# midYear = the year the data was collected in (mid-year if it spans multiple)
# pop20kCorrect = site-level estimate of mean human population density within 20km of the site
# hunterMonths = number of hunters * number of months monitored
# HDI = site-level estimate of the human development index
# distPA = site-level estimate of the distance from the site to the nearest protected area
# propPrimates = the proportion of animals caught that were primates
# propGun = the proportion of animals killed by gun
# propSold = the proportion of animals that were sold
# propVillage = the proportion of hunters that were village hunters rather than hunter-gatherers
# forInt_20km = site-level estimate of the forest condition index
# bmHunterDay = mean hunter offtake per day in kg at the site level
# monitorDays = number of days hunters were monitored
# numSurveyed = number of hunters surveyed
# numRodents = number of rodents caught
# numCetart = number of cetartiodactyls caught
# adjustedSpeciesRichness = adjusted species richness of animals caught
# duikRatio = the ratio of larger to smaller-bodied duikers caught


# checking data structure and preparing columns as factors
str(studyTotYes)
studyTotYes$referenceID <- as.factor(studyTotYes$referenceID)
studyTotYes$siteID <- as.factor(studyTotYes$siteID)
studyTotYes$studyID <- as.factor(studyTotYes$studyID)
studyTotYes$country <- as.factor(studyTotYes$country)
studyTotYes$country <- as.factor(studyTotYes$country)


# There is only 1 zero value in the proportion gun variable
# so we cannot use a zero-inflated model, therefore add 0.001 to the single zero value and run with Beta
# as described in the manuscript
studyTotYes <- studyTotYes %>% 
  mutate(propGun = ifelse(propGun == 0, 0.001, propGun))


# Add Scaled and centred variables for use in analysis
# scale and centre with this function
stdSC_pred <- function(predictor){
  (predictor - mean(predictor,na.rm = T)) / 
    sd(predictor, na.rm = T)}

studyTotYes$access10kaeaCropSC <- stdSC_pred(studyTotYes$access10kaeaCrop)
studyTotYes$forestCorrect20kmSC <- stdSC_pred(studyTotYes$forestCorrect20km)
studyTotYes$midYearSC <- stdSC_pred(studyTotYes$midYear)
studyTotYes$pop20kCorrectSC <- stdSC_pred(studyTotYes$pop20kCorrect)
studyTotYes$hunterMonthsSC <- stdSC_pred(studyTotYes$hunterMonths)
studyTotYes$HDISC <- stdSC_pred(studyTotYes$HDI)
studyTotYes$distPASC <- stdSC_pred(studyTotYes$distPA)
studyTotYes$propPrimatesSC <- stdSC_pred(studyTotYes$propPrimates)
studyTotYes$propGunSC <- stdSC_pred(studyTotYes$propGun)
studyTotYes$propSoldSC <- stdSC_pred(studyTotYes$propSold)
studyTotYes$propVillageSC <- stdSC_pred(studyTotYes$propVillage)
studyTotYes$forestCorrectProp20km <- studyTotYes$forestCorrect20km/100
studyTotYes$forestCorrectProp20kmSC <- stdSC_pred(studyTotYes$forestCorrectProp20km)
studyTotYes$forInt20SC <- stdSC_pred(studyTotYes$forInt_20km)
studyTotYes$bmHunterDaySC <- stdSC_pred(studyTotYes$bmHunterDay)
studyTotYes$monitorDaysSC <- stdSC_pred(studyTotYes$monitorDays)
studyTotYes$numSurveyedSC <- stdSC_pred(studyTotYes$numSurveyed)



# calculate ratio variable rodent:cetart
studyTotYes$ruRatio <- studyTotYes$numRodents / studyTotYes$numCetart
str(studyTotYes)

# check data
studyTotYes %>% glimpse

# Add 0.001 to the few values of ruRatio that are zero, to run the Gamma model
studyTotYes <- studyTotYes %>% 
  mutate(ruRatio = ifelse(ruRatio == 0, 0.001, ruRatio))


# one value for numRodents = 0, so urRatio becomes Inf. so model fails.
# Add 1 to that study only as described in the paper
studyTotYes <- studyTotYes %>% 
  mutate(numRodents = ifelse(numRodents == 0, 1, numRodents))

# calculate ungulates:rodents ratio
studyTotYes$urRatio <- studyTotYes$numCetart / studyTotYes$numRodents

hist(studyTotYes$urRatio) # check data


### checking correlations of variables used in models

ggpairs(studyTotYes, columns = c("propGunSC","distPASC","pop20kCorrectSC", "propVillageSC",
                                "HDISC", "access10kaeaCropSC", "midYearSC", "forInt20SC", "propSoldSC"))

ggpairs(studyTotYes, columns = c("distPA", "propGun"))
ggpairs(studyTotYes, columns = c("HDISC", "hunterMonthsSC"))
ggpairs(studyTotYes, columns = c("propGun", "access10kaeaCropSC"))
ggpairs(studyTotYes, columns = c("forInt20SC", "access10kaeaCropSC"))
ggpairs(studyTotYes, columns = c("midYear", "monitorDaysSC"))
ggpairs(studyTotYes, columns = c("midYear", "numSurveyedSC"))
ggpairs(studyTotYes, columns = c("propGun", "midYear", "HDI"))



# Example of running a simple linear model to check variance inflation factors of possible predictors
lmSEM <- lm(bmHunterDay ~ monitorDaysSC + numSurveyedSC + propGunSC + propVillageSC
            + distPASC + forInt20SC + HDISC + propSoldSC + pop20kCorrectSC
            + access10kaeaCropSC, data = studyTotYes)
vif(lmSEM) # check results



###### Prepare data for modelling

# Priors to include in the model for each response variable

prior3 = c(prior(normal(0, 5), class = "b", resp = "bmHunterDay") +
             prior(normal(0, 5), class = "b", resp = "propGun") +
             prior(normal(0, 5), class = "b", resp = "propSold") +
             prior(normal(0, 5), class = "b", resp = "adjustedSpeciesRichness") +
             prior(normal(0, 5), class = "b", resp = "numPrim") +
             prior(normal(0, 5), class = "b", resp = "urRatio") +
             prior(normal(0, 5), class = "b", resp = "duikRatio"))


# Running models

bioRPf3 <- bf(bmHunterDay ~ numSurveyedSC + monitorDaysSC + propGunSC + distPASC + pop20kCorrectSC
              + HDISC + propVillageSC + forInt20SC + propSoldSC + access10kaeaCropSC
              + (1 |ID| referenceID), decomp="QR") + Gamma(link="log")

gunf3 <- bf(propGun ~ numSurveyedSC + monitorDaysSC + midYearSC*propVillageSC + pop20kCorrectSC
            + access10kaeaCropSC + propSoldSC
            + (1 |ID| referenceID), decomp="QR") + Beta()

usef3 <- bf(propSold ~ numSurveyedSC + monitorDaysSC + midYearSC*propVillageSC + pop20kCorrectSC
            + HDISC + access10kaeaCropSC
            + (1 |ID| referenceID), decomp="QR") + Beta()

specRichf3 <- bf(adjustedSpeciesRichness ~ numSurveyedSC + monitorDaysSC + propGunSC + pop20kCorrectSC
                 + distPASC + forInt20SC + access10kaeaCropSC
                 + (1 |ID| referenceID), decomp="QR") + negbinomial()

urR3 <- bf(urRatio ~ numSurveyedSC + monitorDaysSC + propGunSC + pop20kCorrectSC
           + distPASC + forInt20SC + access10kaeaCropSC
           + (1 |ID| referenceID), decomp="QR") + Gamma(link="log")

primt3 <- bf(numPrim | trials(totHarvest) ~ numSurveyedSC + monitorDaysSC + propGunSC + pop20kCorrectSC
             + distPASC + forInt20SC + access10kaeaCropSC
             + (1 |ID| referenceID), decomp="QR") + beta_binomial(link = "logit", link_phi = "log")

duikRatf3 <- bf(duikRatio ~ numSurveyedSC + monitorDaysSC + propGunSC + pop20kCorrectSC
                + distPASC + forInt20SC + access10kaeaCropSC
                + (1 |ID| referenceID), decomp="QR") + Gamma(link="log")

final_R3 <- brm(bioRPf3 +
                  specRichf3 +
                  urR3 +
                  primt3 +
                  duikRatf3 +
                  gunf3 +
                  usef3 +
                  set_rescor(FALSE), 
                data=studyTotYes, prior = prior3,
                chains = 4, iter = 6000, warmup = 1500,
                cores = 2, seed = 12, init = 0,
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.95, max_treedepth = 20),
                file = "local-file-path")


summary(final_R3) # see results with 95% uncertainty interval
summary(final_R3, prob = 0.89) # see results with 89% uncertainty interval


### checking model assumptions
# residuals
library(bayesplot)
color_scheme_set("purple")

pp_check(final_R3, cores = 1, ndraws = 200, resp = "adjustedSpeciesRichness")
pp_check(final_R3, cores = 1, ndraws = 200, resp = "propGun")
pp_check(final_R3, cores = 1, ndraws = 200, resp = "propSold")
pp_check(final_R3, cores = 1, ndraws = 200, resp = "bmHunterDay")
pp_check(final_R3, cores = 1, ndraws = 200, resp = "urRatio")
pp_check(final_R3, cores = 1, ndraws = 200, resp = "numPrim")
pp_check(final_R3, cores = 1, ndraws = 200, resp = "duikRatio")
bayes_R2(final_R3)

## plot means
pp_check(final_R3, type = "stat", stat = "mean", resp = "adjustedSpeciesRichness")
pp_check(final_R3, type = "stat", stat = "mean", resp = "propGun")
pp_check(final_R3, type = "stat", stat = "mean", resp = "propSold")
pp_check(final_R3, type = "stat", stat = "mean", resp = "bmHunterDay")
pp_check(final_R3, type = "stat", stat = "mean", resp = "urRatio")
pp_check(final_R3, type = "stat", stat = "mean", resp = "numPrim")
pp_check(final_R3, type = "stat", stat = "mean", resp = "duikRatio")


# Assessing convergence
plot(final_R3)


# Check plot
color_scheme_set("darkgray")
mcmc_scatter(
  as.matrix(final_R3),
  pars = c("b_bmHunterDay_HDISC", "b_bmHunterDay_Intercept"), 
  np = nuts_params(final_R3), 
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
)



#plot posterior draws of predictive errors 
plot(predictive_error(final_R3, resp = "bmHunterDay"))


# more predictive error checks
pp_check(final_R3, resp = "bmHunterDay", type='error_scatter_avg')
pp_check(final_R3, resp = "bmHunterDay", type='intervals')


### Extract the probability of direction (Maximum Probability of Effect - MPE)
library(bayestestR)
p_direction(final_R3, method = "kernel")



# check residuals and test overdispersion using DHARMa package

ASRsim <- dh_check_brms(final_R3, resp = "adjustedSpeciesRichness")
DHARMa::testDispersion(ASRsim)

GUNsim <- dh_check_brms(final_R3, resp = "propGun")
DHARMa::testDispersion(GUNsim)

SOLDsim <- dh_check_brms(final_R3, resp = "propSold")
DHARMa::testDispersion(SOLDsim)

HUNTERsim <- dh_check_brms(final_R3, resp = "bmHunterDay")
DHARMa::testDispersion(HUNTERsim)

URRATIOsim <- dh_check_brms(final_R3, resp = "urRatio")
DHARMa::testDispersion(URRATIOsim)

PRIMsim <- dh_check_brms(final_R3, resp = "numPrim")
DHARMa::testDispersion(PRIMsim)

DUIKsim <- dh_check_brms(final_R3, resp = "duikRatio")
DHARMa::testDispersion(DUIKsim)




########### spatial autocorrelation test with DHARMa

# first need to add site coordinates back to dataframe
stackEXTRACT<-read.csv("local-file-path.csv")

studyTotYes$lat <- stackEXTRACT$latitude[match(studyTotYes$siteID, stackEXTRACT$siteID)]
studyTotYes$long <- stackEXTRACT$longitude[match(studyTotYes$siteID, stackEXTRACT$siteID)]

# need to aggregate residuals per location as some studies were conducted at the same sites

studyTotYes$coords <- paste(studyTotYes$long,", ",studyTotYes$lat)
coords <- c(unique(studyTotYes$coords))
x_unique <- c(str_extract(coords, "^.+(?=,)"))
y_unique <- c(str_extract(coords, "(?<=, ).+$"))

# spatial autocorrelation for bmHunterDay model
HUNTERsimCALC<-recalculateResiduals(HUNTERsim,group = studyTotYes$coords)
testSpatialAutocorrelation(HUNTERsimCALC, x = x_unique, y = y_unique) # not spatially autocorrelated

# spatial autocorrelation for duikRatio model
DUIKsimCALC<-recalculateResiduals(DUIKsim,group = studyTotYes$coords)
testSpatialAutocorrelation(DUIKsimCALC, x = x_unique, y = y_unique) # not spatially autocorrelated

# spatial autocorrelation for propGun model
GUNsimCALC<-recalculateResiduals(GUNsim,group = studyTotYes$coords)
testSpatialAutocorrelation(GUNsimCALC, x = x_unique, y = y_unique) # not spatially autocorrelated

# spatial autocorrelation for propSold model
SOLDsimCALC<-recalculateResiduals(SOLDsim,group = studyTotYes$coords)
testSpatialAutocorrelation(SOLDsimCALC, x = x_unique, y = y_unique) # not spatially autocorrelated

# spatial autocorrelation for adjustedSpeciesRichness model
ASRsimCALC<-recalculateResiduals(ASRsim,group = studyTotYes$coords)
testSpatialAutocorrelation(ASRsimCALC, x = x_unique, y = y_unique) # not spatially autocorrelated

# spatial autocorrelation for urRatio model
URRATIOsimCALC<-recalculateResiduals(URRATIOsim,group = studyTotYes$coords)
testSpatialAutocorrelation(URRATIOsimCALC, x = x_unique, y = y_unique) # not spatially autocorrelated

# spatial autocorrelation for numPrim model
PRIMsimCALC<-recalculateResiduals(PRIMsim,group = studyTotYes$coords)
testSpatialAutocorrelation(PRIMsimCALC, x = x_unique, y = y_unique) # not spatially autocorrelated





