## Riva anf Fahrig 2022
## Protecting many small patches will maximize biodiversity conservation for most taxa: the SS > SL principle

# prepare data
library(data.table)
library(tidyr)
library(dplyr)
library(purrr)
library(rlist)
library(splitstackshape)

# diversity estimation
library(iNEXT)
library(SpadeR)
library(vegan)
library(paleotree)

# plotting
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(ggeffects)
library(effects)
#font
library(extrafont)
loadfonts(device = "win")

# modeling
library(glmmTMB)
library(AICcmodavg)
library(MuMIn)


# remove scientific notation
options(scipen=999)
set.seed(654)

## open FragSAD dataset available at https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.2861
## open FragSAD dataset
data = fread("Data\\fragSAD_predicts_ewers.csv", header = TRUE) # version from Chase et al. 2020
metadata = fread("Data\\fragSAD_metadata.csv", header = TRUE)

metadata$dataset_id <- as.factor(metadata$dataset_id)
data$dataset_label <- as.factor(data$dataset_label)

# list of paper to be retained based on assessment
papers_reassessment = fread("Data\\papers_reassessment.csv", header = TRUE)
papers_reassessment_filter <- papers_reassessment[,c(2,3)]

## data cleaning 
## remove 'continuous' patches for which we don't have real size
data <- subset(data, data$frag_size_char != 'continuous') 
## remove studies based on a one-by-one assessment of fit for analysis (see methods)
data <- merge(data, papers_reassessment_filter, by = 'dataset_label')
data <- subset(data, data$included == 1) 

## remove patches with patch size assigned arbitrarily or with patch size larger than 50% of the total HA; see "paper_reassessment"
## unless otherwise stated, removed only arbitrary sized large patches
## datasets where inspected by FR
data <- data[!(data$dataset_label == "Almeida-Gomes_2014" & data$frag_size_num > 7000),] 
data <- data[!(data$dataset_label == "Baz_1995" & data$frag_size_num > 1000),] # largest patch bigger than 50%
data <- data[!(data$dataset_label == "Bell_2006_a" & data$frag_size_num > 1000),] 
data <- data[!(data$dataset_label == "Bell_2006_b" & data$frag_size_num > 1000),] 
data <- data[!(data$dataset_label == "Benedick_2006" & data$frag_size_num > 5000),] # two large patches removes in addition to arbitrary control forests
data <- data[!(data$dataset_label == "Berg_1997" & data$frag_size_num > 999),] 
data <- data[!(data$dataset_label == "Bernard_2007" & data$frag_size_num > 350),] # one large patch removed in addition to arbitrary 750 ha controls
data <- data[!(data$dataset_label == "Bolger_1997" & data$frag_size_num > 999),] 
data <- data[!(data$dataset_label == "Bragagnolo_2007" & data$frag_size_num > 999),] 
data <- data[!(data$dataset_label == "Cabrera-Guzman_2012" & data$frag_size_num > 500),] 
data <- data[!(data$dataset_label == "Davies_2003" & data$frag_size_num > 100),] 
data <- data[!(data$dataset_label == "Gavish_2012_a" & data$frag_size_num > 100),] # remove large patch that is larger than 50% HA
data <- data[!(data$dataset_label == "Gavish_2012_b" & data$frag_size_num > 100),] # remove large patch that is larger than 50% HA
data <- data[!(data$dataset_label == "Gavish_2012_c" & data$frag_size_num > 100),] # remove large patch that is larger than 50% HA
data <- data[!(data$dataset_label == "Giladi_2011_a" & data$frag_size_num > 100),] # remove large patch that is larger than 50% HA
data <- data[!(data$dataset_label == "Giladi_2011_c" & data$frag_size_num > 100),] # remove large patch that is larger than 50% HA
data <- data[!(data$dataset_label == "Giladi_2011_b" & data$frag_size_num > 35),] # remove three large patches that are too larger to compare HA 
data <- data[!(data$dataset_label == "Jung_2014" & data$frag_size_num > 999),] 
data <- data[!(data$dataset_label == "Lima_2015" & data$frag_size_num > 800),] # remove two large patches that are too larger to compare HA
data <- data[!(data$dataset_label == "Manu_2007" & data$frag_size_num > 500),] # remove large patch that is larger than 50% HA
data <- data[!(data$dataset_label == "Montgomery_2014" & data$frag_size_num > 3000),] # remove large patch that is larger than 50% HA
data <- data[!(data$dataset_label == "Nogueira_2016" & data$frag_size_num > 200),] 
data <- data[!(data$dataset_label == "Nyeko_2009" & data$frag_size_num > 200),] 
data <- data[!(data$dataset_label == "Owen_2008" & data$frag_size_num > 200),] 
data <- data[!(data$dataset_label == "Raheem_2009" & data$frag_size_num > 15000),] # remove large patch that is larger than 50% HA
data <- data[!(data$dataset_label == "Silva_2009" & data$frag_size_num > 999),] # remove large patch that is larger than 50% HA
data <- data[!(data$dataset_label == "Slade_2013" & data$frag_size_num > 100),] # remove large patch that is larger than 50% HA
data <- data[!(data$dataset_label == "Struebig_2008" & data$frag_size_num > 15000),] 
data <- data[!(data$dataset_label == "Vulinec_2008" & data$frag_size_num > 70),] # remove two large patches that are too larger to compare HA
data <- data[!(data$dataset_label == "Vasconcelos_2006" & data$frag_size_num > 999),] 
data <- data[!(data$dataset_label == "Brosi_2008" & data$frag_size_num > 200),] # remove large patch that is larger than 50% HA
data <- data[!(data$dataset_label == "Lambert_2003" & data$frag_size_num > 300),] # remove large patch that is larger than 50% HA
data <- data[!(data$dataset_label == "Savilaakso_2009" & data$frag_size_num > 250),] 
data <- data[!(data$dataset_label == "Wang_2012" & data$frag_size_num > 999),] # remove large patch that is larger than 50% HA
data <- data[!(data$dataset_label == "Kapoor_2008" & data$frag_size_num > 1999),] 


## remove site filtering column
data <- data[,-c(10)]

## label as factor
data$dataset_label <- as.factor(data$dataset_label)

## create a unique ID for study name and its patches
data$study_and_patch <- as.factor(paste(data$dataset_label, data$frag_id))

## in "patches", each row represent a patch in the original data
patches <- data[, c(1,2,10,5)]
patches<- unique(patches)

## in "assemblages" and "individuals", every row represents one taxon observed in a patch
assemblages <-  data[,c(1,2,10,8,9)] 
individuals <- assemblages[,-c(4)]

### create lists where each study is one element of a list and replace studies where 
### abundance data are reported as non-integers with that value multiplied by 100
list_individuals <- list()
for (i in 1 : 76){list_individuals[[i]] <- subset(individuals,dataset_label == levels(individuals$dataset_label)[i] )}
# Kappes_2009, Lambert_2003 and Dickman_1987 are expressed as density per area or trap; multiply by 100 to obtain integer values
for (i in 1 : 76){if (min(list_individuals[[i]]$abundance) < 1) {list_individuals[[i]]$abundance <- (list_individuals[[i]]$abundance * 100)}}  
individuals <- do.call(rbind.data.frame, list_individuals)

### number of individuals sampled in each study, in each patch, and total patch area sampled in each study
study_individuals <- aggregate(abundance ~ dataset_label, data=individuals, FUN=sum)
patch_individuals <- aggregate(abundance ~ study_and_patch, data=individuals, FUN=sum)
study_total_area <- aggregate(frag_size_num ~ dataset_label, data=patches, FUN=sum)

### add information to the dataset borrowing from original data
patch_individuals$dataset_label <- individuals$dataset_label[match(patch_individuals$study_and_patch,individuals$study_and_patch)]
patch_individuals$frag_id <- individuals$frag_id[match(patch_individuals$study_and_patch,individuals$study_and_patch)]
patch_individuals$total_abundance <- study_individuals$abundance[match(patch_individuals$dataset_label,study_individuals$dataset_label)]
patch_individuals$patch_area <- patches$frag_size_num[match(patch_individuals$study_and_patch,patches$study_and_patch)]
patch_individuals$total_patch_area <- study_total_area$frag_size_num[match(patch_individuals$dataset_label,study_total_area$dataset_label)]

#################################################################
##### CALCULATE SAMPLING EFFORT PER PATCH
#################################################################

# create a new column to split into list based on sample design
# "sample_id" was "plot_id" in FragSAD, defined as "plot_id: Plot-identifier of a sampling plot or a transect. If only one plot per site, 1 is given." in supplementary information
# note that in Owen_2008 the sample_id does not start from 1 in every patch
data$unique_sampling_event <- as.factor(paste(data$study_and_patch, data$sample_id))

# retain the effort (column 6), study_and_patch (column 10), and also the sampling event (column 11, created one line of code above) 
data_sampling_effort <- data[, c(6,10,11)]

# unique, because every sampling (column 11) event has the same effort (column 6)
data_sampling_effort <- unique(data_sampling_effort) 

# sum the total effort put in every patch
data_sampling_effort <- aggregate(sample_eff ~ study_and_patch, data = data_sampling_effort, FUN=sum) # by aggregating, we have a sum of the sampling events

# add a dataset_label from original data
data_sampling_effort$dataset_label <- data$dataset_label[match(data_sampling_effort$study_and_patch, data$study_and_patch)]

# transform sampling effort in values between 0 and 1, where 1 is the maximum sampling effort observed in a dataset
data_sampling_effort <- split(data_sampling_effort, data_sampling_effort$dataset_label)
for (i in 1: length(data_sampling_effort)) {
  data_sampling_effort[[i]]$sample_eff <- data_sampling_effort[[i]]$sample_eff/max(data_sampling_effort[[i]]$sample_eff)
}

data_sampling_effort <- do.call(rbind.data.frame, data_sampling_effort)

# match the relative sampling effort with information on every patch
patch_individuals$sampling_effort <- data_sampling_effort$sample_eff[match(patch_individuals$study_and_patch, data_sampling_effort$study_and_patch)]
# the number of individuals expected if all patches were sampled with the same effort
patch_individuals$abundance_per_effort <- patch_individuals$abundance / patch_individuals$sampling_effort


# calculate normalized # of individuals adjusted (i.e., at equal sampling effort) per patch 
patch_individuals <- split(patch_individuals, patch_individuals$dataset_label)
for (i in 1:length(patch_individuals)) {
  patch_individuals[[i]]$normalized_individuals_per_sampl_unit <- patch_individuals[[i]]$abundance_per_effort / max(patch_individuals[[i]]$abundance_per_effort)
}

patch_individuals <- do.call(rbind.data.frame, patch_individuals)

# normalized_individuals_per_sampl_unit2 is a transformation to avoind values equal to 1 or zero; necessary to fit beta regression
# transformation follows Smithson, M., & Verkuilen, J. (2006). A better lemon squeezer? Maximum-likelihood regression with beta-distributed dependent variables. Psychological Methods, 11(1), 54–71.
patch_individuals$normalized_individuals_per_sampl_unit2 <- (patch_individuals$normalized_individuals_per_sampl_unit * (nrow(patch_individuals)-1) + 0.5) / nrow(patch_individuals)
patch_individuals$log_patch_size <- log(patch_individuals$patch_area)

#  model structure follows Chase et al. 2020 in Nature 
model <- glmmTMB(normalized_individuals_per_sampl_unit2 ~ log_patch_size  +  
                   (log_patch_size|dataset_label), 
                 data = patch_individuals, family = "gaussian") ## family=beta_family(link="logit")) ## 
# gaussian distribution does a better job than beta distribution 
#summary(model)
#AIC(model)

patch_individuals$predicted_normalized_individuals <- predict(model, type = "response")


#################################################################
##### CALCULATE INDIVIDUALS PER RANDOMIZATIONS
#################################################################

# calculate the proportion of individuals per patch
patch_individuals$abundance_per_ha <- patch_individuals$abundance / patch_individuals$patch_area
patch_individuals$relative_patch_size <- patch_individuals$patch_area / patch_individuals$total_patch_area 
# calculate the minimum to define the # of individuals simulated 
min_abundance_per_ha_study <- aggregate(abundance_per_ha ~ dataset_label, data=patch_individuals, FUN=min)
colnames(min_abundance_per_ha_study)[2] <- "min_abundance_per_ha_study"

# merge with patch-level information
patch_individuals <- merge(patch_individuals, min_abundance_per_ha_study, by = 'dataset_label')
patch_individuals$simulation_constant_abundance <- patch_individuals$min_abundance_per_ha_study * patch_individuals$patch_area
#patch_individuals$simulations_adjusted_abundance <- patch_individuals$simulation_constant_abundance * patch_individuals$predicted_normalized_individuals

# increase the number of individuals to be resampled in each dataset by standardizing the predictions within each dataset
patch_individuals <- split(patch_individuals, patch_individuals$dataset_label)
for (i in 1:length(patch_individuals)) {
  patch_individuals[[i]]$predicted_normalized_individuals_adjusted <- patch_individuals[[i]]$predicted_normalized_individuals / max(patch_individuals[[i]]$predicted_normalized_individuals)
}
patch_individuals <- do.call(rbind.data.frame, patch_individuals)

patch_individuals$simulations_adjusted_abundance <- patch_individuals$simulation_constant_abundance * patch_individuals$predicted_normalized_individuals_adjusted



#################################################################
##### SIMULATE 100 TIMES THE NUMBER OF INDIVIDUALS TO BE RESAMPLED IN EACH PATCH
#################################################################

# create lists where each study is one element of a list
list_patches <- list()
for (i in 1 : 76){list_patches[[i]] <- subset(patches,dataset_label == levels(data$dataset_label)[i] )}

list_assemblages <- list()
for (i in 1 : 76){list_assemblages[[i]] <- subset(assemblages,dataset_label == levels(data$dataset_label)[i] )}

# if the smallest value of a list is < 1, multiply all the abundances of the species recorded as relative occurrences for 100 
# so that rows can be created proportionally to their incidence
for (i in 1 : 76){if (min(list_assemblages[[i]]$abundance) < 1) {list_assemblages[[i]]$abundance <- (list_assemblages[[i]]$abundance * 100)}}  

# extend so that each row has nrow = nindividuals
for (i in 1 : 76) {list_assemblages[[i]] <- splitstackshape::expandRows(list_assemblages[[i]], "abundance") }
list_assemblages2 <- rbindlist(list_assemblages)

# list of all species seen in all patches in the original datasets
list_assemblages_patches <- list()
for (i in 1 : nlevels(list_assemblages2$study_and_patch)) {list_assemblages_patches[[i]] <-  subset(list_assemblages2,study_and_patch == levels(list_assemblages2$study_and_patch)[i] )}


### create 100 number of individuals sampled per patch
list_sim_individuals <- list()

for (i in 1 : nrow(patch_individuals)){
  list_sim_individuals[[i]] <- replicate(100, (floor(patch_individuals$simulations_adjusted_abundance[[i]]) + rbinom(1, 1,(patch_individuals$simulations_adjusted_abundance[[i]] - floor(patch_individuals$simulations_adjusted_abundance[[i]])))))
}
## these are the numbers of individuals that will be selected in each patch; 100 randomizations



################################################################
####### CREATE RANDOMIZED TABLES
################################################################

sim_ind <- do.call(rbind.data.frame, list_sim_individuals)
hist(log10(rowSums(sim_ind)))

### SAMPLE ROWS ON A PER-PATCH BASIS TO HAVE # OF INDIVIDUALS PROPORTIONAL TO PATCH SIZE
get_rows <- function(df, rows) df[rows, , drop = FALSE]


# followed https://jennybc.github.io/purrr-tutorial/ls12_different-sized-samples.html
list_sampled_data <- list()
split_data <- split(list_assemblages2[,1:4], list_assemblages2$study_and_patch)
group_sizes <- vapply(split_data, nrow, integer(1))

for (i in 1:100){
  n <-sim_ind[,i]
  sampled_obs <- mapply(sample, group_sizes, n)
  sampled_data <- mapply(get_rows, split_data, sampled_obs, SIMPLIFY = FALSE)
  sampled_data2 <- do.call(rbind, sampled_data) 
  
  list_sampled_data[[i]] <- sampled_data2
}


# ## if needed, break in lists
for (i in 1:100){
  list_sampled_data[[i]] <- split(list_sampled_data[[i]], list_sampled_data[[i]][,3]) #
}

# convert from a series of rows each representing a sampled individuals, into a vectors of abundances per patch
# when a patch has zero simulated individuals, this still create a column that I will remove later to avoid losing the patch from the dataset
for (i in 1:length(list_sampled_data)){
  for (j in 1:length(list_sampled_data[[i]])){
    if(nrow(list_sampled_data[[i]][j][[1]]) == 0) {
      mydf <- (data.frame(Var1 = "no_species_present", Freq = 0))
      list_sampled_data[[i]][j][[1]]<- setNames(data.frame(t(mydf[,-1])), mydf[,1])
    } else {
      list_sampled_data[[i]][j][[1]] <- setNames(data.frame(t(as.data.frame(table(list_sampled_data[[i]][j][[1]]$species))[,-1])), as.data.frame(table(list_sampled_data[[i]][j][[1]]$species))[,1])
    }
  }
}


# split the 100 lists into 76 datasets
study_index <- patches$dataset_label
for (i in 1: length(list_sampled_data)){
  list_sampled_data[[i]] <- split(list_sampled_data[[i]], study_index)
}

# function to convert list of vectors of species abundances into table
to_table <- function(list_to_table) {
  table_prova <- rbindlist(lapply(list_to_table, function(x) as.data.frame.list(x)), fill=TRUE)
  table_prova[is.na(table_prova)] <- 0
  table_prova$study_and_patch <- names(list_to_table)
  list_to_table <- table_prova
}

# make each list (100 simulations per 76 studies) into a table of species per sites
for (i in 1:length(list_sampled_data)){
  for (j in 1: length(list_sampled_data[[i]])){
    list_sampled_data[[i]][[j]] <- to_table(list_sampled_data[[i]][[j]]) #rbindlist(lapply(list_sampled_data[[i]][[j]], function(x) as.data.frame.list(x)), fill=TRUE)
  }
}

# remove empty column that was necessary to retain patches with zero simulated individuals
for (i in 1:length(list_sampled_data)){
  for (j in 1: length(list_sampled_data[[i]])){
    
    if("no_species_present" %in% colnames(list_sampled_data[[i]][j][[1]])){
      list_sampled_data[[i]][j][[1]] <- select(list_sampled_data[[i]][j][[1]], -no_species_present)
    } 
    
  }
}


# add patch size to each list
list_patches_unlisted <- do.call(rbind.data.frame, list_patches)
list_patches_unlisted <- list_patches_unlisted[,3:4]

for (i in 1:length(list_sampled_data)){
  for (j in 1: length(list_sampled_data[[i]])){
    
    list_sampled_data[[i]][j][[1]] <- merge(list_sampled_data[[i]][j][[1]], list_patches_unlisted, by = "study_and_patch")
    
  }
}


################################################################
####### SLOSS COMPARISONS
################################################################

#### sloss function adapted from package Lexiguel https://github.com/kamapu/Lexiguel
sloss <- function(table, env=data.frame(), area) {
  if(!is.matrix(table)) table <- as.matrix(table)
  area <- substitute(area)
  area <- eval(area, env, parent.frame())
  SLOSS <- list(SL=list(), LS=list())
  # First the calculation from small to large
  SLOSS$SL$area <- c(0, cumsum(area[order(area)]))
  Flor <- apply(table[order(area),], 2, cumsum)
  Flor[Flor > 0] <- 1
  SLOSS$SL$species <- c(0, apply(Flor, 1, sum))
  # Now the calculation from large to small
  SLOSS$LS$area <- c(0, cumsum(area[order(area, decreasing=TRUE)]))
  Flor <- apply(table[order(area, decreasing=TRUE),], 2, cumsum)
  Flor[Flor > 0] <- 1
  SLOSS$LS$species <- c(0, apply(Flor, 1, sum))
  # Calculation of SLOSS index
  SLOSS$Index <- with(SLOSS$SL, curve_area(area,
                                           species))/with(SLOSS$LS, curve_area(area, species))
  # Final object
  class(SLOSS) <- c("SLOSS","list")
  return(SLOSS)
}

curve_area <- function(x, y, bottom=0) {
  D1 <- c(diff(x))
  D2 <- c(diff(y))
  Area <- sum(D1*((y - bottom)[-length(y)]) + D1*D2/2, na.rm=TRUE)
  return(Area)
}


# edited original function in Lexiguel to include points
plot.SLOSS <- function(x, y=NULL, sl.lty=2, sl.lwd=1, sl.col="black", ls.lty=1,
                       ls.lwd=1, ls.col="black", show.index=TRUE, digits.index=2, cex.index=1,
                       pos.index=c(0.05,0.95), show.legend=FALSE, pos.legend="bottomright",
                       bty.legend="o", main="SLOSS curves",...) {
  with(x$SL, plot(area, species, type="l", lty=sl.lty, lwd=sl.lwd, col=sl.col,
                  main=main, ...))
  with(x$LS, lines(area, species, lty=ls.lty, lwd=ls.lwd, col=ls.col))
  
  with(x$SL, points(area, species, pch = 18, cex = 1))
  with(x$LS, points(area, species, pch = 18, cex = 1))
  
  if(show.legend) {
    legend(pos.legend, lty=c(sl.lty,ls.lty), lwd=c(sl.lwd,ls.lwd),
           legend=c("small to large","large to small"), bty=bty.legend)
  }
  if(show.index) {
    with(x$SL, text(max(area)*pos.index[1], max(species)*pos.index[2],
                    labels=paste("SLOSS-index =",
                                 round(x$Index, digits.index)),
                    cex=cex.index, pos=4))
  }
}




# 
INPUT <- list_sampled_data[[10]][74][[1]]

SL_OR_SS <- function(INPUT) {
  
  table_sloss <- as.matrix(INPUT)
  area_sloss <-   table_sloss[,ncol(table_sloss)]
  
  table_sloss <- table_sloss[, - c(1, ncol(table_sloss))]
  sloss_comparison <- sloss(table_sloss, area = area_sloss)
  
  # sl <- approx(sloss_comparison$SL$area, 
  #              sloss_comparison$SL$species, 
  #              xout = seq(from = as.numeric(min(area_sloss)), to = (max(sloss_comparison$SL$area) - as.numeric(max(area_sloss))), length.out = 100), 
  #              method = "linear")
  # 
  # ls <- approx(sloss_comparison$LS$area, 
  #              sloss_comparison$LS$species, 
  #              xout = seq(from = as.numeric(min(area_sloss)), to = (max(sloss_comparison$SL$area) - as.numeric(max(area_sloss))), length.out = 100), 
  #              method = "linear")
  
  sl <- approx(sloss_comparison$SL$area,
               sloss_comparison$SL$species,
               xout = seq(from = as.numeric(max(area_sloss)), to = (max(sloss_comparison$SL$area) - as.numeric(max(area_sloss))), length.out = 100),
               method = "linear")
  
  ls <- approx(sloss_comparison$LS$area,
               sloss_comparison$LS$species,
               xout = seq(from = as.numeric(max(area_sloss)), to = (max(sloss_comparison$SL$area) - as.numeric(max(area_sloss))), length.out = 100),
               method = "linear")
  
  plot.SLOSS(sloss_comparison, show.index = F)
  points(ls, col = "blue")
  points(sl, col = "red")
  
  comparison_sloss <- as.data.frame(t(rbind(ls$y, sl$y)))
  
  comparison_sloss$outcome <- comparison_sloss$V1 - comparison_sloss$V2 # V1 is LS
  
  outcome <- if (all(comparison_sloss$outcome > 0)){
    print("SL > SS")
  } else {
    if (all(comparison_sloss$outcome < 0)){print("SS > SL")
    } else {
      print("SS = SL")
    }
  }
  
  #return(outcome)
}

# test function
SL_OR_SS(list_sampled_data[[88]][7][[1]]) # first number simulations from 1 to 100, second number represent dataset from 1 to 76

## create data for analysis
ELE_data <- list_sampled_data

for (i in 1:length(ELE_data)){
  for (j in 1: length(ELE_data[[i]])){
    
    ELE_data[[i]][j][[1]] <- SL_OR_SS(ELE_data[[i]][j][[1]])
    
  }
}

study_names <- names(ELE_data[[1]])


for (i in 1:length(ELE_data)){
  ELE_data[[i]] <- do.call(rbind.data.frame, ELE_data[[i]])
}

for (i in 1:length(ELE_data)){
  colnames(ELE_data[[i]]) <- "comparison"
}

ele_data_analysis <- do.call(rbind.data.frame, ELE_data)
ele_data_analysis$dataset_id <- rep(study_names, 100)

## correct Edwards 2010, which was incorrectly classified, from plants to birds
metadata[35,14] <- metadata[9,14]

ele_data_analysis <- merge(ele_data_analysis, metadata, by = "dataset_id")

# check proportion of SLOSS comparisons
n_dataset <- table(ele_data_analysis$taxa)/100
n_comparisons <- table(ele_data_analysis$taxa, ele_data_analysis$comparison)/100
# SL > SS, SS = SL, SS > SL
n_comparisons[1:5]/n_dataset
n_comparisons[6:10]/n_dataset
n_comparisons[11:15]/n_dataset

# plot raw outcome of SLOSS comparison
ggplot_data <- as.data.frame(table(ele_data_analysis$comparison, ele_data_analysis$dataset_id))

# create a dataframe with 3 columns, one per SLOSS outcome
ggplot_data2 <- reshape(ggplot_data, idvar = "Var2", timevar = "Var1", direction = "wide")
# order the dataframe

ggplot_data2 <- ggplot_data2[order(ggplot_data2[,4], -ggplot_data2[,2] ),]
names(ggplot_data2)[1] <- "dataset_id"

# set the factors in the same order as the ordered dataframe
ggplot_data2$dataset_id <- factor(ggplot_data2$dataset_id,ggplot_data2$dataset_id)

ggplot_data <- melt(ggplot_data2, id.vars="dataset_id")
ggplot_data$variable <- gsub("Freq.", "", ggplot_data$variable)
ggplot_data <- subset(ggplot_data, ggplot_data$variable != "order")

ggplot_data_taxa <- metadata[,c(1, 14)]
ggplot_data <- merge(ggplot_data, ggplot_data_taxa)

ggplot_data_label <- unique(ggplot_data[, c(1,4)])


## FIGURE 2
p1 <- ggplot(ggplot_data, aes(fill=variable, y=value, x= dataset_id )) + 
  geom_bar(position="fill", stat="identity")+ 
  scale_fill_manual(values = c("blue", "grey", "red"))+
  theme_few()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(text=element_text(size=18,  family="serif"))+
  theme(legend.title=element_blank())+
  ylab("Proportion of SLOSS comparisons in a dataset")+ theme(legend.position="top")#+ 
  #ggtitle("Full dataset")+ theme(plot.title = element_text(hjust = 0.5))
  #facet_wrap(~ taxa)#+
  #theme(axis.text.x = element_text(angle = 90,hjust =0 ))

p2 <- ggplot(ggplot_data[ggplot_data$taxa == "amphibians & reptiles",], aes(fill=variable, y=value, x= dataset_id )) + 
  geom_bar(position="fill", stat="identity")+ 
  scale_fill_manual(values = c("blue", "grey", "red"))+
  theme_few()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(text=element_text(size=14,  family="serif"))+
  theme(legend.title=element_blank())+
  ylab("Proportion of SLOSS comparisons in a dataset")+ 
  ggtitle("amphibians & reptiles") + theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="none", axis.title.y = element_blank())

p3 <- ggplot(ggplot_data[ggplot_data$taxa == "mammals",], aes(fill=variable, y=value, x= dataset_id )) + 
  geom_bar(position="fill", stat="identity")+ 
  scale_fill_manual(values = c("blue", "grey", "red"))+
  theme_few()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(text=element_text(size=14,  family="serif"))+
  theme(legend.title=element_blank())+
  ylab("Proportion of SLOSS comparisons in a dataset")+ 
  ggtitle("mammals") + theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="none", axis.title.y = element_blank())
#f

p5 <- ggplot(ggplot_data[ggplot_data$taxa == "invertebrates",], aes(fill=variable, y=value, x= dataset_id )) + 
  geom_bar(position="fill", stat="identity")+ 
  scale_fill_manual(values = c("blue", "grey", "red"))+
  theme_few()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(text=element_text(size=14,  family="serif"))+
  theme(legend.title=element_blank())+
  ylab("Proportion of SLOSS comparisons in a dataset")+ 
  ggtitle("invertebrates") + theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="none", axis.title.y = element_blank())
#f

p4 <- ggplot(ggplot_data[ggplot_data$taxa == "birds",], aes(fill=variable, y=value, x= dataset_id )) + 
  geom_bar(position="fill", stat="identity")+ 
  scale_fill_manual(values = c("blue", "grey", "red"))+
  theme_few()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(text=element_text(size=14,  family="serif"))+
  theme(legend.title=element_blank())+
  ylab("Proportion of SLOSS comparisons in a dataset")+ 
  ggtitle("birds") + theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="none", axis.title.y = element_blank())
#f

p6 <- ggplot(ggplot_data[ggplot_data$taxa == "plants",], aes(fill=variable, y=value, x= dataset_id )) + 
  geom_bar(position="fill", stat="identity")+ 
  scale_fill_manual(values = c("blue", "grey", "red"))+
  theme_few()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(text=element_text(size=14,  family="serif"))+
  theme(legend.title=element_blank())+
  ylab("Proportion of SLOSS comparisons in a dataset")+ 
  ggtitle("plants") + theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="none", axis.title.y = element_blank())


# 
ggarrange(p1, 
          ggarrange (p2, p3, p4, p5, p6, nrow = 1),
          ncol = 1,
          heights = c(0.75, 0.25))

ggsave("fig2.jpg", path = "Figures", width = 3600, height = 2200, units = "px",  device='jpeg', dpi=300)



#### MODELING ANALYSIS
ele_data_analysis$SS <- as.factor(ele_data_analysis$comparison)
#levels(ele_data_analysis$SS)
levels(ele_data_analysis$SS)[match("SL > SS",levels(ele_data_analysis$SS))] <- "0"
levels(ele_data_analysis$SS)[match("SS = SL",levels(ele_data_analysis$SS))] <- "0"
levels(ele_data_analysis$SS)[match("SS > SL",levels(ele_data_analysis$SS))] <- "1"

ele_data_analysis$SL <- as.factor(ele_data_analysis$comparison)
levels(ele_data_analysis$SL)[match("SL > SS",levels(ele_data_analysis$SL))] <- "1"
levels(ele_data_analysis$SL)[match("SS = SL",levels(ele_data_analysis$SL))] <- "0"
levels(ele_data_analysis$SL)[match("SS > SL",levels(ele_data_analysis$SL))] <- "0"

ele_data_analysis$SSSL <- as.factor(ele_data_analysis$comparison)
levels(ele_data_analysis$SSSL)[match("SL > SS",levels(ele_data_analysis$SSSL))] <- "0"
levels(ele_data_analysis$SSSL)[match("SS = SL",levels(ele_data_analysis$SSSL))] <- "1"
levels(ele_data_analysis$SSSL)[match("SS > SL",levels(ele_data_analysis$SSSL))] <- "0"

# set up factors
ele_data_analysis$dataset_id <- as.factor(ele_data_analysis$dataset_id)
ele_data_analysis$taxa <- as.factor(ele_data_analysis$taxa)


# calculate patch size evenness for every study
evenness <- list()
for (i in 1: length(list_patches)){
  evenness[[i]] <- (vegan::diversity(list_patches[[i]]$frag_size_num))/log(length(list_patches[[i]]$frag_size_num))
}

evenness <- do.call(rbind.data.frame, evenness)
evenness$dataset_id <- study_names
colnames(evenness) <- c("evenness", "dataset_id")

ele_data_analysis <- merge(ele_data_analysis, evenness, by = "dataset_id")


# calculate species richness for every study
richness <- list()
for (i in 1: length(list_patches)){
  richness[[i]] <- length(table(list_assemblages[[i]]$species))
}

richness <- do.call(rbind.data.frame, richness)
richness$dataset_id <- study_names
colnames(richness) <- c("richness", "dataset_id")
ele_data_analysis <- merge(ele_data_analysis, richness, by = "dataset_id")


#
data_m <- ele_data_analysis
data_m$comparison <- as.factor(data_m$comparison)

library(brms)
# models
fit <- brm(comparison ~ taxa + evenness + (1|dataset_id), data = data_m, 
           family = categorical(link = "logit"), #multinomial(), 
           inits = 0,
           seed   = 123,
           chains = 4, cores = 10)

# save.image("SLOSS3.RData")
# # Load workspace back to RStudio
#load("SLOSS3.RData")


conditions <- make_conditions(fit, vars = c("taxa", "evenness"))
conditions_herp <- conditions[1:3,]

ce <- conditional_effects(
  fit, 
 conditions = conditions_herp,
  prob = 0.89,
  categorical = TRUE,
  spaghetti = FALSE)#,
  #spaghetti = TRUE,
  #,
  #effects = "taxa:evenness")

plot(ce, plot = FALSE)[[1]] + 
  scale_color_manual(values=c("blue", "grey", "red"))+
  theme_few()+
  theme(text=element_text(size=18,  family="serif"))+
  theme(axis.title.x=element_blank(),
        strip.text.x = element_blank())+
  theme(legend.title=element_blank())+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  #scale_x_discrete(guide = guide_axis(n.dodge=2))

ggsave("fig35.jpg", path = "C:\\Users\\feder\\OneDrive\\Desktop\\Riva Fahrig SLOSS #3\\figures", 
       width = 4100, height = 1000, units = "px",  device='jpeg', dpi=300)


## figure 4
# mean patch size
patches <- data[, c(1,2,10,5)]
patches<- unique(patches)
names(patches)[1] <- "dataset_id"
patches <- merge(patches, metadata[,c(1,14)], by = "dataset_id")
patches$log10area <- log10(patches$frag_size_num)

patches$taxa <- factor(patches$taxa, levels = c("amphibians & reptiles", "mammals", "birds", "invertebrates", "plants"))

ggplot(patches, aes(x=log10area, fill=taxa)) +
  geom_histogram( color="olivedrab",fill = "chartreuse4", position = 'identity', col=I("olivedrab"))+
  facet_wrap(~taxa, ncol = 5)+
  theme_few() +
  scale_y_continuous(expand = c(0,0),limits = c(0,60))+
  geom_rect(aes(xmin=2, xmax=Inf, ymin=0,ymax=Inf), alpha=0.002, fill="gold")+
  geom_vline(xintercept = 2, linetype="longdash", color = "gold", size=0.5)+
    geom_vline(xintercept = 3, linetype="longdash", color = "gold", size=0.5)+
  geom_vline(xintercept = 3.47, linetype="longdash", color = "gold", size=0.5)+
    geom_vline(xintercept = 4, linetype="longdash", color = "gold", size=0.5)+
  geom_vline(xintercept = 4.69, linetype="longdash", color = "gold", size=0.5)+
  scale_fill_manual(values = c("black", "black", "black", "black", "black"))+
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2, 3, 4), limits =c(-3, 5),
                     labels = c("0.01","0.1","1","10","100","1.000","10.000"),
                     guide = guide_axis(n.dodge=2))+
  xlab("Patch area (ha)") + ylab("Number of patches")+
  theme(legend.position="none")+
  theme(text=element_text(size=18,  family="serif"))

ggsave("fig4.jpg", path = "C:\\Users\\feder\\OneDrive\\Desktop\\Riva Fahrig SLOSS #3\\figures", 
       width = 4500, height = 1500, units = "px",  device='jpeg', dpi=300)


# how many patches for each taxa?
patches[patches$taxa == "amphibians & reptiles",]
patches[patches$taxa == "mammals",]
patches[patches$taxa == "birds",]
patches[patches$taxa == "invertebrates",]
patches[patches$taxa == "plants",]


# how many patches are larger than the thresholds in fig 4?
sum(patches$frag_size_num > 1000)/nrow(patches) #27 patches
sum(patches$frag_size_num > 10000)/nrow(patches) # 2 patches
sum(patches$frag_size_num > 50000)/nrow(patches) # 0 patches

# examples of evenness for figure 3 (see patches of habitat in the top row)
set_patches <- c(1,1,2,10,16,30)
vegan::diversity(set_patches)/log(length(set_patches))

set_patches <- c(1,2,2,5,10)
vegan::diversity(set_patches)/log(length(set_patches))

set_patches <- c(2,5,7)
(vegan::diversity(set_patches))/log(length(set_patches))
