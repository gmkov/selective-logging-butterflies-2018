#### GMK 2018 #####
#### impacts of selective logging on butterflies, Biological Conservation ####
#### on the effects of logging intensity on species richness, abundance, and composition ###
### example code, not tested. any queries to mgm49@cam.ac.uk ###
### FIG. 3 ####

rm(list=ls())
dev.off()
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(AER)
library(MASS)
library(vegan)
library(mgcv)
library(cowplot)
library(gridExtra)
library(gridBase)
library(tibble)

### data load ####
master <- read.csv("butterfly.capture.data.csv")
site_rich_master <- read.csv("site_info.csv") #contains site logging intensity, capture rates, sp richness
matrix <- read.csv("matrix_spab.csv", row.names = 1) #matrix of species abundance per site and species

# load colourblind palette
cbbPalette <- c("#000000", "#009E73", "#e79f00", "#9ad0f3", "#0072B2", "#D55E00", 
                "#CC79A7", "#F0E442")


#############################################################################
##### 1. FIGURE 3 A1 boxplot ab (capture_site in 48h) primary vs logged ######

# add column with forest type
site_rich_master$type <- ifelse(grepl("lo", site_rich_master$site), "logged", "primary")

## ab in 48h 
site_rich_master$cap_per_hour <- site_rich_master$total_cap_site/site_rich_master$effort_hours
site_rich_master$cap_48h <- site_rich_master$cap_per_hour*48

## boxplot of butterflies captured in a 48h period (primary vs logged)

a1_48h <- ggplot(site_rich_master, aes(x=type, y=cap_48h))+
  geom_boxplot(aes(fill=type),
               outlier.size = 1)+
  scale_fill_manual(values=c("#CC79A7", "#009E73"))+
  geom_text(size=5,color="black",aes(x=0.65,y=18, label="A1"))+
  theme_classic()+
  theme(
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.text = element_text(size=14), 
    axis.title = element_text(size = 14))+
  theme(legend.position="none")+
  labs(y="Total abundance per site in 48h")+
  labs(x="")+
  scale_x_discrete(labels=c("Logged \n sites", "Primary \n sites"))

a1_48h

####################################################################
###### 2. FIGURE 3, A2 gam just logged sites abundance (48h) vs intensity ########

# just logged
logg_ab_int <- subset(site_rich_master, site_rich_master$type=="logged")

# analyses
pred_simu <- list()

GAM <- gam(cap_48h~s(int_vol, k=2), data=logg_ab_int, family='poisson', select=TRUE)
summary(GAM)
gam.check(GAM)

#predict at diff ints, obtain CI for plotting and store
pred_simu <- predict.gam(GAM, data.frame(int_vol=0:70), type="response", se.fit = T)
preds_df <- data.frame(value=pred_simu, int_vol=0:70)

preds_df$up95 <- preds_df$value.fit+(1.96*preds_df$value.se.fit)
preds_df$low95 <- preds_df$value.fit-(1.96*preds_df$value.se.fit)


a2_48h <- ggplot(logg_ab_int, aes(x=int_vol, y=cap_48h)) +
  geom_site() +
  geom_jitter(width = .1, aes(-2), data = site_rich_master[site_rich_master$type=="primary",], 
              fill="#009E73", colour="black", shape=24, size=2.2, alpha=.8)+
  geom_line(data=preds_df,aes(x=int_vol, y=value.fit))+
  geom_text(size=5,color="black",aes(x=1,y=18, label="A2"))+
  geom_ribbon(data=preds_df, aes(x=int_vol, y=value.fit, ymax=up95, ymin=low95), alpha=0.5)+
  labs(y="Total abundance per site in 48h")+
  xlab(expression(paste("Logging intensity (", m^3/ha^-1, ")", sep="")))+
  theme_classic()+
  theme(
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.text = element_text(size=14), 
    axis.title = element_text(size = 14))
#ylim(0,18)

a2_48h


#########################################################################
###### 3. FIGURE 3, B1 boxplot sp_rich (48h) site prim vs logged #######

site_rich_master$sp_rich48h <- (site_rich_master$sp_rich/site_rich_master$effort_hours)*48

## plot

b1_48h <- ggplot(site_rich_master, aes(x=type, y=sp_rich48h))+
  geom_boxplot(aes(fill=type),
               outlier.size = 1)+
  scale_fill_manual(values=c("#CC79A7", "#009E73"))+
  geom_text(size=5,color="black",aes(x=0.65,y=4, label="B1")) +
  theme_classic() +
  theme(
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.text = element_text(size=14), 
    axis.title = element_text(size = 14))+
  theme(legend.position="none")+
  labs(y="Species richness per site in 48h") +
  labs(x="")+
  scale_x_discrete(labels=c("Logged \n sites", "Primary \n sites"))
b1_48h


####################################################################
###### 4. FIGURE 3, B2 gam just logged sprich in 48h vs int ########

#take logged sites only
logg_sprich_int <- subset(site_rich_master, site_rich_master$type=="logged")

# preambles
pred_simu <- list()

GAM <- gam(sp_rich48h~s(int_vol, k=2), data=logg_sprich_int, family='poisson', select=TRUE)
summary(GAM)
gam.check(GAM)

#predict at diff ints, obtain CI for plotting and store
pred_simu <- predict.gam(GAM, data.frame(int_vol=0:70), type="response", se.fit = T)
preds_df <- data.frame(value=pred_simu, int_vol=0:70)

preds_df$up95 <- preds_df$value.fit+(1.96*preds_df$value.se.fit)
preds_df$low95 <- preds_df$value.fit-(1.96*preds_df$value.se.fit)

# plot

b2_48h <- ggplot(logg_sprich_int, aes(x=int_vol, y=sp_rich48h)) +
  geom_site() +
  geom_jitter(width = .1, aes(-2), data = site_rich_master[site_rich_master$type=="primary",], 
             fill="#009E73", colour="black", shape=24, size=2.2, alpha=.8)+
  geom_line(data=preds_df,aes(x=int_vol, y=value.fit))+
  geom_text(size=5,color="black",aes(x=1,y=4, label="B2"))+
  geom_ribbon(data=preds_df, aes(x=int_vol, y=value.fit, ymax=up95, ymin=low95), alpha=0.5)+
  theme_classic()+
  theme(
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.text = element_text(size=14), 
    axis.title = element_text(size = 14))+
  labs(y="Species richness per site in 48h") +
  xlab(expression(paste("Logging intensity (", m^3/site, ")", sep="")))+
  ylim(0,4)


b2_48h


###################################
###### 5. FIGURE 3, C- NMDS #######

# make matrix
summary <- summarise(group_by(master, site, genus_sp),
                                ab=n())

matrix <- spread(summary, 
                            genus_sp,
                            ab,
                            fill = "0")

matrix <- as.data.frame(matrix)
row.names(matrix) <- matrix[,1]
matrix <- matrix[,-1]
matrix <- data.matrix(matrix)

# order site_rich_master columns
target <- rownames(matrix)
site_order <- site_rich_master[match(target, site_rich_master$site),]

##### NMDS anlyses 
nmds <- metaMDS(matrix, k=2, trymax = 100)
nmds$stress
nmds

## NMDS plot first look
plot(nmds)

# plot 
data.scores <- as.data.frame(scores(nmds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$grp <- factor(site_order$type)  #  add the grp variable created earlier
data.scores

c <- ggplot() + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_site(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=grp,colour=grp),size=5) + # add the site markers
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site),size=3,vjust=0) +  # add the site labels
  scale_colour_manual(values=c("logged" = "#CC79A7","primary"="#009E73"),name="Site type") +
  scale_shape_manual(values=c(16,17), name="Site type") +
  geom_text(size=5,color="black",aes(x=-1,y=1, label="C", fontface=2))+
  #ylim(-1,1)+
  #xlim(-1.5,1) +
  theme_classic()+
  theme(
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.text = element_text(size=14), 
    axis.title = element_text(size = 14),
  #theme(legend.position=c(0.15,0.2))
  legend.position = "none")

c


#################################################################
###### 6. FIGURE 2D- NMDS axis 1 vs logging intensity GAM #######

###### just logged
sc <- scores(nmds, choices=c(1,2))

# grab nmds values from logged sites from nmds in C) 
sc.logged <- sc[c(1:16,37:60),]
site_order.logged <- subset(site_order, site_order$type=="logged")
plotting.logged  <- data_frame(int=site_order.logged$int_vol, NMDS1=sc.logged[,1])

sc.primary <- sc[17:36,]
site_order.primary <- subset(site_order, site_order$type=="primary")
plotting.primary  <- data_frame(int=site_order.primary$int_vol, NMDS1=sc.primary[,1])

preds_list <- list()
predictions <- list()
pred_simu <- list()
preds_df <-list()

GAM.logged <- gam(NMDS1~s(int, k=2), data=plotting.logged, select=TRUE)
summary(GAM.logged)

#predict at diff ints, obtain CI for plotting and store
pred_simu <- predict.gam(GAM.logged, data.frame(int=0:70), type="response", se.fit = T)
preds_df <- data.frame(value=pred_simu, int=0:70)
preds_df$up95 <- preds_df$value.fit+(1.96*preds_df$value.se.fit)
preds_df$low95 <- preds_df$value.fit-(1.96*preds_df$value.se.fit)

#plot
d_logged  <- ggplot(plotting.logged, aes(x=int, y=NMDS1)) +
  geom_site() +
  geom_jitter(width = .1, aes(-2), data = plotting.primary, 
              fill="#009E73", colour="black", shape=24, size=2.2, alpha=.8)+
  geom_line(data=preds_df,aes(x=int, y=value.fit))+
  geom_text(size=5,color="black",aes(x=0.5,y=1.1, label="D"))+
  geom_ribbon(data=preds_df, aes(x=int, y=value.fit, ymax=up95, ymin=low95), alpha=0.5)+
  theme_classic()+
  theme(
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.text = element_text(size=14), 
    axis.title = element_text(size = 14))+
  ylab("NMDS1")+
  xlab(expression(paste("Logging intensity (", m^3 , "/site)", sep="")))

d_logged

############## COMBINE PLOTS #############

setwd("/figs/")
png("fig3_comm.png", bg="transparent", width = 480, height = 280, units = "cm",
    res = 300)
plot_grid(a1_48h, a2_48h, c, b1_48h, b2_48h, d_logged, rel_widths=c(0.6,1.2,1.2,0.6,1.2,1.2), ncol = 3, nrow = 2)
dev.off()
