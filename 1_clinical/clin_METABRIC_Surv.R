# Script: Survival analyses in the METABRIC cohort based on OS, RFI
# Author: Lennart Hohmann
# Date: 07.11.2024
#-------------------
# empty environment
rm(list=ls())
# set working directory to the project directory
setwd("~/PhD_Workspace/Project_Basal/")
# cohort
cohort <- "METABRIC"
#-------------------
# packages
source("./scripts/src/general_functions.R")
source("./scripts/1_clinical/src/clin_functions.R")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, 
               ggplot2,
               ggfortify,
               survival,
               survminer,
               grid,
               gridExtra,
               ggplotify,
               remotes,
               stringi)
#-------------------
# set/create output directories
# for plots
output.path <- "./output/1_clinical/"
dir.create(output.path)
# for data
data.path <- "./data/METABRIC/1_clinical/processed/"
dir.create(data.path)
#-------------------
# input paths
infile.1 <- "./data/METABRIC/0_GroupSamples/ERpHER2n_sampleIDs.RData"
infile.2 <- "./data/Parameters/color_palette.RData"
infile.3 <- "./data/METABRIC/1_clinical/raw/Merged_annotations.RData"
# output paths
plot.file <- paste0(output.path,cohort,"_SA.pdf")
txt.file <- paste0(output.path,cohort,"_SA.txt")
#outfile.1 <- paste0(data.path,"SA.RData")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

# load sampleIDs
sampleIDs <- loadRData(infile.1)[c("ERpHER2n_Basal", "ERpHER2n_LumA", "ERpHER2n_LumB")]
# load palette
color.palette <- loadRData(infile.2)[c("LumA","LumB","Basal")]
names(color.palette) <- paste0("PAM50=", names(color.palette))

# Load data
svdata <- loadRData(infile.3)
# Sample selection
svdata <- svdata[svdata$METABRIC_ID %in% unname(unlist(sampleIDs)), ]
# Handle NA values in Chemotherapy and Endocrine
svdata$Chemotherapy[is.na(svdata$Chemotherapy)] <- 0
svdata$Endocrine[is.na(svdata$Endocrine)] <- 0
# Create Treatment variable
svdata$Treatment <- ifelse(svdata$Chemotherapy == 1 & svdata$Endocrine == 1, "CE",
                           ifelse(svdata$Chemotherapy == 0 & svdata$Endocrine == 1, "E", NA))
# Create LN variable based on lymph_nodes_positive
svdata$LN <- ifelse(svdata$lymph_nodes_positive > 0, "N+", "N0")
# Convert specified columns to numeric using lapply
numeric_columns <- c("OS", "OSbin", "RFI", "RFIbin", "Age", "TumSize")
svdata[numeric_columns] <- lapply(svdata[numeric_columns], as.numeric)
# Set as a factor with specific levels
svdata$Grade <- factor(svdata$Grade, levels = c("1", "2", "3"))
svdata$PAM50 <- factor(svdata$PAM50, levels = c("LumA", "LumB", "Basal"))
svdata$LN <- factor(svdata$LN, levels = c("N0", "N+"))
# Select relevant columns
svdata <- svdata[, c("METABRIC_ID", "Age", "Grade", "TumSize", "LN", 
                     "PAM50", "OS", "OSbin", "RFI", "RFIbin", "Treatment")]

#######################################################################
# Part 2: OS
#######################################################################

OM <- "OS"
OMbin <- "OSbin"

######################################
# Investigate the EC treatment group #
######################################

# data, surv object, and fit
EC.dat <- svdata %>% 
  filter(Treatment == "CE")
EC.surv <- Surv(EC.dat[[OM]], EC.dat[[OMbin]])
EC.fit <- survminer::surv_fit(EC.surv~PAM50, data=EC.dat, 
                              conf.type="log-log") 
# label output
plot.title <- paste0("ERpHER2n; ",cohort,"; CT+ET; ",OM)
txt.out <- append(txt.out,
                  c(plot.title, "\n", 
                    paste0("Analyses with clinical endpoint: ",OM, " in treatment group: CT+ET"),"\n###########################################\n"))
# subtype numbers
txt.out <- append(txt.out,
                  c(capture.output(table(EC.dat$PAM50)),"\n###########################################\n"))
# add the median OM of censored patients
median.EC <- median(EC.dat[which(EC.dat[[OMbin]]==0),][[OM]])
txt.out <- append(txt.out,c(paste0("Median ",OM, " for CT+ET censored patients = ",median.EC),"\n###########################################\n"))

##########################

# uv cox
main.pam50 <- coxph(EC.surv~PAM50, data=EC.dat)
res <- summary(main.pam50)
plot <- ggforest(main.pam50,data=EC.dat,
                 main=plot.title) + theme_bw()

plot.list <- append(plot.list,list(plot))
txt.out <- append(txt.out,c(capture.output(res),"\n###########################################\n"))

##########################

# KM plot
plot <- ggsurvplot(EC.fit, 
                   data = EC.dat,
                   pval = TRUE, conf.int = FALSE,
                   xlab = paste0(OM," (years)"), 
                   break.x.by = 1,
                   break.y.by = 0.1,
                   ylab = paste0(OM," probability"),
                   ylim = c(0,1),
                   title = plot.title,
                   legend = c(0.9,0.9),
                   legend.title = "Subtypes",
                   palette = color.palette)[["plot"]]

plot.list <- append(plot.list,list(plot))

##########################

# mv cox
main.all <- coxph(EC.surv~PAM50+Age+LN+TumSize+Grade, 
                  data=EC.dat) 
res <- summary(main.all)
plot <- ggforest(main.all,
                 data=EC.dat,
                 main=plot.title) + theme_bw()

plot.list <- append(plot.list,list(plot))
txt.out <- append(txt.out,c(capture.output(res),"\n###########################################\n"))

########################################
# Investigate the Endo treatment group
########################################

# data, surv object, and fit
E.dat <- svdata %>% 
  filter(Treatment == "E")
E.surv <- Surv(E.dat[[OM]], E.dat[[OMbin]])
E.fit <- survminer::surv_fit(E.surv~PAM50, data=E.dat, 
                             conf.type="log-log") 
# label output
plot.title <- paste0("ERpHER2n; ",cohort,"; ET; ",OM)
txt.out <- append(txt.out,
                  c(plot.title, "\n", 
                    paste0("Analyses with clinical endpoint: ",OM, " in treatment group: ET"),"\n###########################################\n"))
# subtype numbers
txt.out <- append(txt.out,
                  c(capture.output(table(E.dat$PAM50)),"\n###########################################\n"))
# add the median OM of censored patients
median.E <- median(E.dat[which(E.dat[[OMbin]]==0),][[OM]])
txt.out <- append(txt.out,c(paste0("Median ",OM, " for ET censored patients = ",median.E),"\n###########################################\n"))

##########################

# uv cox
main.pam50 <- coxph(E.surv~PAM50, data=E.dat)
res <- summary(main.pam50)
plot <- ggforest(main.pam50,data=E.dat,
                 main=plot.title) + theme_bw()

plot.list <- append(plot.list,list(plot))
txt.out <- append(txt.out,c(capture.output(res),"\n###########################################\n"))

##########################

# KM plot
plot <- ggsurvplot(E.fit, 
                   data = E.dat,
                   pval = TRUE, conf.int = FALSE,
                   xlab = paste0(OM," (years)"), 
                   break.x.by = 1,
                   break.y.by = 0.1,
                   ylab = paste0(OM," probability"),
                   ylim = c(0,1),
                   title = plot.title,
                   legend = c(0.9,0.96),
                   legend.title = "Subtypes",
                   palette = color.palette)[["plot"]]

plot.list <- append(plot.list,list(plot))

##########################

# mv cox
main.all <- coxph(E.surv~PAM50+Age+LN+TumSize+Grade, 
                  data=E.dat) 
res <- summary(main.all)
plot <- ggforest(main.all,
                 data=E.dat,
                 main=plot.title) + theme_bw()

plot.list <- append(plot.list,list(plot))
txt.out <- append(txt.out,c(capture.output(res),"\n###########################################\n"))


#######################################################################
# Part 3: RFI
#######################################################################

OM <- "RFI"
OMbin <- "RFIbin"

######################################
# Investigate the EC treatment group #
######################################

# data, surv object, and fit
EC.dat <- svdata %>% 
  filter(Treatment == "CE")
EC.surv <- Surv(EC.dat[[OM]], EC.dat[[OMbin]])
EC.fit <- survminer::surv_fit(EC.surv~PAM50, data=EC.dat, 
                              conf.type="log-log") 
# label output
plot.title <- paste0("ERpHER2n; ",cohort,"; CT+ET; ",OM)
txt.out <- append(txt.out,
                  c(plot.title, "\n", 
                    paste0("Analyses with clinical endpoint: ",OM, " in treatment group: CT+ET"),"\n###########################################\n"))
# subtype numbers
txt.out <- append(txt.out,
                  c(capture.output(table(EC.dat$PAM50)),"\n###########################################\n"))
# add the median OM of censored patients
median.EC <- median(EC.dat[which(EC.dat[[OMbin]]==0),][[OM]])
txt.out <- append(txt.out,c(paste0("Median ",OM, " for CT+ET censored patients = ",median.EC),"\n###########################################\n"))

##########################

# uv cox
main.pam50 <- coxph(EC.surv~PAM50, data=EC.dat)
res <- summary(main.pam50)
plot <- ggforest(main.pam50,data=EC.dat,
                 main=plot.title) + theme_bw()

plot.list <- append(plot.list,list(plot))
txt.out <- append(txt.out,c(capture.output(res),"\n###########################################\n"))

##########################

# KM plot
plot <- ggsurvplot(EC.fit, 
                   data = EC.dat,
                   pval = TRUE, conf.int = FALSE,
                   xlab = paste0(OM," (years)"), 
                   break.x.by = 1,
                   break.y.by = 0.1,
                   ylab = paste0(OM," probability"),
                   ylim = c(0,1),
                   title = plot.title,
                   legend = c(0.9,0.3),
                   legend.title = "Subtypes",
                   palette = color.palette)[["plot"]]

plot.list <- append(plot.list,list(plot))

##########################

# mv cox
#exp(confint(main.all)) # the problem is that grade 1 has no event
# exclude and redefine surv object
EC.dat <- EC.dat[which(EC.dat$Grade != 1),]
EC.dat$Grade <- relevel(EC.dat$Grade, ref = "2")
EC.dat$Grade <- droplevels(EC.dat$Grade)
EC.surv <- Surv(EC.dat[[OM]], EC.dat[[OMbin]])
#
main.all <- coxph(EC.surv~PAM50+Age+LN+TumSize+Grade,
                  data=EC.dat) 
res <- summary(main.all)
plot <- ggforest(main.all,
                 data=EC.dat,
                 main=plot.title)

plot.list <- append(plot.list,list(plot))
txt.out <- append(txt.out,c(capture.output(res),"\n###########################################\n"))

########################################
# Investigate the Endo treatment group
########################################

# data, surv object, and fit
E.dat <- svdata %>% 
  filter(Treatment == "E")
E.surv <- Surv(E.dat[[OM]], E.dat[[OMbin]])
E.fit <- survminer::surv_fit(E.surv~PAM50, data=E.dat, 
                             conf.type="log-log") 
# label output
plot.title <- paste0("ERpHER2n; ",cohort,"; ET; ",OM)
txt.out <- append(txt.out,
                  c(plot.title, "\n", 
                    paste0("Analyses with clinical endpoint: ",OM, " in treatment group: ET"),"\n###########################################\n"))
# subtype numbers
txt.out <- append(txt.out,
                  c(capture.output(table(E.dat$PAM50)),"\n###########################################\n"))
# add the median OM of censored patients
median.E <- median(E.dat[which(E.dat[[OMbin]]==0),][[OM]])
txt.out <- append(txt.out,c(paste0("Median ",OM, " for ET censored patients = ",median.E),"\n###########################################\n"))

##########################

# uv cox
main.pam50 <- coxph(E.surv~PAM50, data=E.dat)
res <- summary(main.pam50)
plot <- ggforest(main.pam50,data=E.dat,
                 main=plot.title) + theme_bw()

plot.list <- append(plot.list,list(plot))
txt.out <- append(txt.out,c(capture.output(res),"\n###########################################\n"))

##########################

# KM plot
plot <- ggsurvplot(E.fit, 
                   data = E.dat,
                   pval = TRUE, conf.int = FALSE,
                   xlab = paste0(OM," (years)"), 
                   break.x.by = 1,
                   break.y.by = 0.1,
                   ylab = paste0(OM," probability"),
                   ylim = c(0,1),
                   title = plot.title,
                   legend = c(0.9,0.3),
                   legend.title = "Subtypes",
                   palette = color.palette)[["plot"]]

plot.list <- append(plot.list,list(plot))

##########################

# mv cox
main.all <- coxph(E.surv~PAM50+Age+LN+TumSize+Grade, 
                  data=E.dat) 
res <- summary(main.all)
plot <- ggforest(main.all,
                 data=E.dat,
                 main=plot.title) + theme_bw()

plot.list <- append(plot.list,list(plot))
txt.out <- append(txt.out,c(capture.output(res),"\n###########################################\n"))

#######################################################################
#######################################################################

# save plots
pdf(file = plot.file, onefile = TRUE)#, width = 8.3/8, height = 11.7/8) 
for (i in 1:length(plot.list)) {
  print(plot.list[[i]])
}
dev.off()

# save text
writeLines(txt.out, txt.file)