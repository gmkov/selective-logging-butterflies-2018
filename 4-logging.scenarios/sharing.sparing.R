### modeling land-sharing vs land-sparing logging concessions ###

### workflow ###

# 1. Prepare data, remove primary sites, remove sp found in <3sites
# 2. Create individual species abundance-logging intenstity curves & obtain predictions at different intensity levels
# 3. Make hypothetical concessions with predicted species numbers at each intensity level
 

rm(list=ls())
dev.off()
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

###################################
######### 1. PREPARE DATA #########
###################################

master <- read.csv("butterfly.capture.data.csv")
site_data <- read.csv("site_info.csv") #contains site logging intensity, capture rates, sp richness
sp_pres <- read.csv("site_sp_ab.csv") # contains individual species abundances per site

# make presence sp df
# remove riodinidae
sp_pres$presence <- ifelse(sp_pres$ab>0, 1, 0)

#disable scinetific notation
options(scipen=999)

#add intensity
sp_pres$yield <- site_data$int_vol[match(site_sp_pres$site, site_data$site)]

#add trap hours
head(site_data)
sp_pres$trap_hours <- site_data$effort_hours[match(site_sp_pres$site, site_data$site)]
colnames(p_pres)[3] <- "value"
colnames(sp_pres)[2] <- "species"

# to count logged sites with logging intensity 0 as logged
sp_pres$forest_type <- ifelse(grepl("lo", sp_pres$site), "logged", "primary")

## calculate no of singleton sp
sp_by_forest<- summarise(group_by(sp_pres, forest_type, species ),
          ab=sum(value))
sp_by_forest_wide <- spread(sp_by_forest, forest_type, ab)
sp_by_forest_wide$logged_unique <- ifelse(sp_by_forest_wide$logged-sp_by_forest_wide$primary==sp_by_forest_wide$logged, "yes", "no")
sp_by_forest_wide$primary_unique <- ifelse(sp_by_forest_wide$primary-sp_by_forest_wide$logged==sp_by_forest_wide$primary, "yes", "no")
sp_by_forest_wide$unique <- ifelse(sp_by_forest_wide$primary_unique=="no"&sp_by_forest_wide$logged_unique=="no", "no", "yes")

logged_unique <- subset(sp_by_forest_wide, sp_by_forest_wide$logged_unique=="yes")
primary_unique <- subset(sp_by_forest_wide, sp_by_forest_wide$primary_unique=="yes")

singletons<- summarise(group_by(sp_pres,species ),
                         ab=sum(value))
no.singletons <- (subset(singletons,ab==1))
logged.singletons <- subset(no.singletons, no.singletons$species %in% logged_unique$species)
primary.singletons <- subset(no.singletons, no.singletons$species %in% primary_unique$species)

write.csv(logged_unique, "/Users/gabrielamontejokovacevich/Dropbox (Cambridge University)/MRES/PROJECT/PAPER1/REANALYSIS/1.METHOD5/data/logged_unique.csv", row.names = FALSE)
write.csv(primary_unique, "/Users/gabrielamontejokovacevich/Dropbox (Cambridge University)/MRES/PROJECT/PAPER1/REANALYSIS/1.METHOD5/data/primary_unique.csv", row.names = FALSE)


###################################################
###### 2. SPECIES CURVES and PREDICTIONS ##########
###################################################

#### effort adjusted abundance for ten days collection ####

library(pscl)
library(boot)
library(AER)
library(MASS)
library(mgcv)

# remove primary sites for species abundance-yield curves
sp_pres.logged <- subset(sp_pres, sp_pres$forest_type=="logged")
sp_pres.logged$presence <- ifelse(sp_pres.logged$value>0, 1, 0)

# adjust for effort 
sp_pres.logged$value.adj <- as.numeric((sp_pres.logged$value/sp_pres.logged$trap_hours)*1300)

# remove species were less than three occurences (less than three sites)
n.occurrences <- tapply(sp_pres.logged$presence, INDEX= sp_pres.logged$species, FUN = sum)
no.n.deficient <- n.occurrences[n.occurrences >=3 ] 
uniquespecies <- names(no.n.deficient) 
spids <- names(no.n.deficient) 

#creates list of species names (sp ids)
ll<-length(spids)
spids

#collection zone 
dispersion_tests <- list()
inflation_tests <- list()
occurrence_summaries <- list()
GLM_summaries <- list() #creates a list to store your model summaries
GLM.occ_summaries <- list()
POLYModel_summaries <- list() #creates a list to store your model summaries
simple_summaries <- list()

slopes <- rep(NA, length(uniquespecies))
sp <- c(uniquespecies)
results <- data.frame(slopes, sp)
options(scipen=999)
preds_list <- list()
predictions <- list()

AIC <- rep(NA, length(uniquespecies))
GVC <- rep(NA, length(uniquespecies))
RSQ <- rep(NA, length(uniquespecies))
sp <- c(uniquespecies)
table_GAM <- data.frame(AIC, GCV, RSQ, sp)

dev.off()

for (i in 1:length(uniquespecies)) #
{ 
  #Preamble
  cat(paste('species=',uniquespecies[i]), '\n')  #just prints the iteration number, for reference only.
  #subset the data by species
  speciesname_data_simp <- sp_pres.logged[sp_pres.logged$species==uniquespecies[i], ] #subsets data by name of species
  
  ## GAM
  GAM <- gam(value.adj~s(yield, k=3), data=speciesname_data_simp, family='poisson', select=TRUE)
  simple_summaries[[i]]  <- GAM$coefficients #stores model coefficients in a list
  names(simple_summaries)[[i]] <- paste(uniquespecies[i])
  
  #predict at diff yields, obtain CI for plotting and store
  pred_simu <- predict.gam(GAM, data.frame(yield=0:75), type="response", se.fit = T)
  preds_df <- data.frame(value=pred_simu, yield=0:75)
  predictions[[i]] <- uniquespecies
  predictions[[i]] <- pred_simu
  preds_df$up95 <- preds_df$value.fit+(1.96*preds_df$value.se.fit)
  preds_df$low95 <- preds_df$value.fit-(1.96*preds_df$value.se.fit)
  
  # plot and save curves 
  png(filename = paste0("/plots/time.adjusted/", uniquespecies[i], ".png"))
  #remember ggplots need to be names and printed to be saved
  a <- ggplot(speciesname_data_simp, aes(x=yield, y=value.adj)) +
    geom_site() +
    geom_line(data=preds_df,aes(x=yield, y=value.fit))+
    geom_ribbon(data=preds_df, aes(x=yield, y=value.fit, ymax=up95, ymin=low95), alpha=0.5)+
    theme_classic()+
    ylab("abundance (time adjusted)")+
    xlab("logging intensity")
  #ylim(0,50)
  print(a)
  dev.off()
  
  #species abundance predictions for simulations (just preds)
  pred_values <- predict.gam(GAM, data.frame(yield=c(10, 15, 20, 25, 30, 35, 40)), type="response")
  preds_list[[i]] <- pred_values
}


options(scipen=999)
#create df with all predicitons from all sp
predictions_all <- do.call("rbind", preds_list)
#make into dataframe
predictions_all <- data.frame(predictions_all)
#make column with sp names
predictions_all$species <- uniquespecies

#change order columns
predictions_all <- predictions_all[,c(8,1,2,3,4,5,6,7)]

#change column names (i10= intensity 10m3 per site)
colnames(predictions_all) <- c("species", "i10", "i15", "i20", "i25", "i30", "i35", "i40" )
head(predictions_all)

# make prediction of primary sites per species and add it to predictions all
sp_pres.primary <- subset(sp_pres, sp_pres$forest_type=="primary")

# use only those with more than 3 occurences (89 nymph sp)
sp_pres.primary <- subset(sp_pres.primary, sp_pres.primary$species %in% uniquespecies)

# adjust for effort.
sp_pres.primary$value.adj <- as.numeric((sp_pres.primary$value/sp_pres.primary$trap_hours)*1300)

# take mean occurence as primary forest prediction
primary.site.preds <- summarise(group_by(sp_pres.primary, species),
                                mean.ab.adj=mean(value.adj))

predictions_all$i0 <- primary.site.preds$mean.ab.adj
predictions_all <- predictions_all[,c(1,9,2,3,4,5,6,7,8)]

#save the predictions (predicted abundance of species at each logging intensity)
write.csv(predictions_all, "/predictions.timeadj.all.csv", row.names = FALSE)



##################################################
########### 3. Hypothetical concessions ##########
##################################################

# create hypothetical concessions following Fig. 2.
# concession numbering goes from most sharing to most sparing

conc_1 <- data.frame(species=predictions_all$species, 
                     abundance= predictions_all$i20*4,
                     conc="conc_1")

conc_2 <- data.frame(species=predictions_all$species, 
                     abundance= predictions_all$i30 + predictions_all$i20 +predictions_all$i20 +predictions_all$i10,
                     conc="conc_3")

conc_3 <- data.frame(species=predictions_all$species, 
                     abundance= predictions_all$i30 + predictions_all$i25 +predictions_all$i15 +predictions_all$i10,
                     conc="conc_5")

conc_4 <- data.frame(species=predictions_all$species, 
                     abundance= predictions_all$i35 + predictions_all$i30 +predictions_all$i15 +predictions_all$i0,
                     conc="conc_8")

conc_5 <- data.frame(species=predictions_all$species, 
                      abundance= predictions_all$i40 + predictions_all$i40 +predictions_all$i0 +predictions_all$i0,
                      conc="conc_10")

### combine concessions
all_conc <- rbind(primary,
                  conc_1,
                  conc_2,
                  conc_3,
                  conc_4,
                  conc_5)

write.csv(all_conc, "/Users/gabrielamontejokovacevich/Dropbox (Cambridge University)/MRES/PROJECT/PAPER1/REANALYSIS/1.METHOD5/data/all_conc.time.adj.csv", row.names = FALSE)


# hypothetical concessions ready to explore



