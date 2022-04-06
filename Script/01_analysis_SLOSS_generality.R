## Riva anf Fahrig 2022
## Generalities in SLOSS

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

# modeling
library(glmmTMB)
library(DHARMa)

#font
library(extrafont)
loadfonts(device = "win")


# remove scientific notation
options(scipen=999)
set.seed(654)
#memory.limit(size=30000) # increase default memory use to 32 gb = 32000 (1 gb = 1000 mb)

## open FragSAD dataset available at https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.2861
## .csv data and metadata provided for convenience
## open FragSAD dataset
data = fread("C:\\Users\\feder\\OneDrive\\Desktop\\Riva Fahrig SLOSS #2\\data\\fragSAD_predicts_ewers.csv", header = TRUE)
metadata = fread("C:\\Users\\feder\\OneDrive\\Desktop\\Riva Fahrig SLOSS #2\\data\\new_meta_2_merge.csv", header = TRUE)
metadata_original = fread("C:\\Users\\feder\\OneDrive\\Desktop\\Riva Fahrig SLOSS #2\\data\\FragSAD_metadata_original.csv", header = TRUE)

metadata$dataset_id <- as.factor(metadata$dataset_id)
data$dataset_label <- as.factor(data$dataset_label)

papers_reassessment = fread("C:\\Users\\feder\\OneDrive\\Desktop\\Riva Fahrig SLOSS #2\\data\\papers_reassessment.csv", header = TRUE)
papers_reassessment_filter <- papers_reassessment[,c(2,3)]
#data$frag_id <- as.factor(data$frag_id)

## data cleaning 
## remove 'continuous' patches for which we don't have real size
data <- subset(data, data$frag_size_char != 'continuous') 
## remove studies that I assessed on a one-by-one basis
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
# note that in Owen_2008 the sample_id do not start from 1 in every patch
data$unique_sampling_event <- as.factor(paste(data$study_and_patch, data$sample_id))

# retain the effort (column 6), study_and_patch (column 10), and also the sampling event (column 11, created one line of code above) 
data_sampling_effort <- data[, c(6,10,11)]

# unique, because every sampling (column 11) event has the same effort (column 6)
data_sampling_effort <- unique(data_sampling_effort) # 1475 sampling events for 1222 patches, some have been sampled multiple times

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

# c.lfs + (c.lfs | dataset_label)  we followed Chase et al. 2020 model structure
model <- glmmTMB(normalized_individuals_per_sampl_unit2 ~ log_patch_size  +  
                   (log_patch_size|dataset_label), 
                 data = patch_individuals, family = "gaussian")## family=beta_family(link="logit")) ## 
# gaussian distributio ndoes a better job than beta distribution 
summary(model)
AIC(model)


patch_individuals$predicted_normalized_individuals <- predict(model, type = "response")

# strong effect of dataset random effect
plot(patch_individuals$predicted_normalized_individuals ~ patch_individuals$dataset_label)

# some datasets have a positive area effect, some have a negative patch area effect
plot(predict(model, type = "response"), patch_individuals$log_patch_size)

# relationship between model predictions and observed number of individuals per sampling unit (normalized to a value between 0 and 1)
plot(predict(model, type = "response"), patch_individuals$normalized_individuals_per_sampl_unit2)
abline(0, 1, col = "red")
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
####### create randomized tables
################################################################

sim_ind <- do.call(rbind.data.frame, list_sim_individuals)
hist(log10(rowSums(sim_ind)))

## TEST IF SIMULATIONS ARE CORRECT: check if simulations are bigger than the number of individuals observed in a patch
# list_test_sim <- rep( list(list()), 1186 )
# for (i in 1:100){
#   for (k in 1: 1186){
#     list_test_sim[[k]][i] <- isTRUE(sim_ind[k,i] > patch_individuals[k, 3])
#  }
# }
# 
# test_sim <- do.call(rbind.data.frame, list_test_sim)
# any(test_sim==TRUE)
## none of the simulated # of individuals in each patch is larger than the original number of individuals observed


### SAMPLE ROWS ON A PER-PATCH BASIS TO HAVE # OF INDIVIDUALS PROPORTIONAL TO PATCH SIZE
get_rows <- function(df, rows) df[rows, , drop = FALSE]

# split_data <- split(list_assemblages2, list_assemblages2$study_and_patch)
# group_sizes <- vapply(split_data, nrow, integer(1))
# n <-sim_ind[,3]
# sampled_obs <- mapply(sample, group_sizes, n)
# sampled_data <- mapply(get_rows, split_data, sampled_obs, SIMPLIFY = FALSE)
# sampled_data2 <- do.call(rbind, sampled_data)


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

# TEST: did we remove the "no_species_present" column added to retain empty patches? Note that empty patches must be retained for the simulations 
# list_sampled_data[[1]][74][[1]] # does not have the "no_species_present" column anymore

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


# edited from original function in Lexiguel to include points
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

SL_OR_SS(list_sampled_data[[88]][7][[1]]) # first number simulations from 1 to 100, second number dataset from 1 to 75

## ELE study


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


ele_data_analysis <- merge(ele_data_analysis, metadata, by = "dataset_id")

# check proportion of SLOSS comparisons
n_dataset <- table(ele_data_analysis$taxa)/100
n_comparisons <- table(ele_data_analysis$taxa, ele_data_analysis$comparison)/100
# SL > SS, SS + SL, SS > SL
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

## plot
ggplot(ggplot_data, aes(fill=variable, y=value, x= dataset_id )) + 
  geom_bar(position="fill", stat="identity")+ 
  scale_fill_manual(values = c("blue", "grey", "red"))+
  theme_few()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(text=element_text(size=16,  family="serif"))+
  theme(legend.title=element_blank())+
  ylab("Proportion of SLOSS comparisons in a dataset")
  

# ggplot_data$taxa2 <- ggplot_data$taxa
# ggplot_data$taxa2[ggplot_data$taxa2 == "invertebrates"] <- 1
# ggplot_data$taxa2[ggplot_data$taxa2 == "mammals"] <- 2
# ggplot_data$taxa2[ggplot_data$taxa2 == "plants"] <- 3
# ggplot_data$taxa2[ggplot_data$taxa2 == "amphibians & reptiles"] <- 4
# ggplot_data$taxa2[ggplot_data$taxa2 == "birds"] <- 5
# ggplot_data$taxa2 <- as.numeric(ggplot_data$taxa2)


#### analysis
ele_data_analysis$SS <- as.factor(ele_data_analysis$comparison)
#levels(ele_data_analysis$SS)
levels(ele_data_analysis$SS)[match("SL > SS",levels(ele_data_analysis$SS))] <- "0"
levels(ele_data_analysis$SS)[match("SS = SL",levels(ele_data_analysis$SS))] <- "0"
levels(ele_data_analysis$SS)[match("SS > SL",levels(ele_data_analysis$SS))] <- "1"

ele_data_analysis$SL <- as.factor(ele_data_analysis$comparison)
levels(ele_data_analysis$SL)[match("SL > SS",levels(ele_data_analysis$SL))] <- "1"
levels(ele_data_analysis$SL)[match("SS = SL",levels(ele_data_analysis$SL))] <- "0"
levels(ele_data_analysis$SL)[match("SS > SL",levels(ele_data_analysis$SL))] <- "0"

# set up factors
ele_data_analysis$dataset_id <- as.factor(ele_data_analysis$dataset_id)
ele_data_analysis$sphere.fragment <- as.factor(ele_data_analysis$sphere.fragment)
ele_data_analysis$sphere.matrix <- as.factor(ele_data_analysis$sphere.matrix)
ele_data_analysis$biome <- as.factor(ele_data_analysis$biome)
ele_data_analysis$taxa <- as.factor(ele_data_analysis$taxa)
ele_data_analysis$time.since.fragmentation <- as.factor(ele_data_analysis$time.since.fragmentation)
ele_data_analysis$Matrix.category <- as.factor(ele_data_analysis$Matrix.category)


# calculate patch size evenness for every study
evenness <- list()
for (i in 1: length(list_patches)){
  evenness[[i]] <- (vegan::diversity(list_patches[[i]]$frag_size_num))/log(length(list_patches[[i]]$frag_size_num))
}

evenness <- do.call(rbind.data.frame, evenness)
evenness$dataset_id <- study_names
colnames(evenness) <- c("evenness", "dataset_id")

ele_data_analysis <- merge(ele_data_analysis, evenness, by = "dataset_id")

# mean patch size
mean_patch_size <- list()
for (i in 1: length(list_patches)){
  mean_patch_size[[i]] <- mean(list_patches[[i]]$frag_size_num)
}
mean_patch_size <- do.call(rbind.data.frame, mean_patch_size)
mean_patch_size$dataset_id <- study_names
colnames(mean_patch_size) <- c("mean_patch_size", "dataset_id")


# calculate species richness for every study
richness <- list()
for (i in 1: length(list_patches)){
  richness[[i]] <- length(table(list_assemblages[[i]]$species))
}

richness <- do.call(rbind.data.frame, richness)
richness$dataset_id <- study_names
colnames(richness) <- c("richness", "dataset_id")

ele_data_analysis <- merge(ele_data_analysis, richness, by = "dataset_id")

# add ectothermic vs endothermic
ele_data_analysis$temperature <- ele_data_analysis$taxa

levels(ele_data_analysis$temperature)[match("amphibians & reptiles",levels(ele_data_analysis$temperature))] <- "ectotherm"
levels(ele_data_analysis$temperature)[match("birds",levels(ele_data_analysis$temperature))] <- "endotherm"
levels(ele_data_analysis$temperature)[match("invertebrates",levels(ele_data_analysis$temperature))] <- "ectotherm"
levels(ele_data_analysis$temperature)[match("mammals",levels(ele_data_analysis$temperature))] <- "endotherm"
levels(ele_data_analysis$temperature)[match("plants",levels(ele_data_analysis$temperature))] <- "ectotherm"


# add verts vs inv
ele_data_analysis$vert <- ele_data_analysis$taxa

levels(ele_data_analysis$vert)[match("amphibians & reptiles",levels(ele_data_analysis$vert))] <- "vert"
levels(ele_data_analysis$vert)[match("birds",levels(ele_data_analysis$vert))] <- "vert"
levels(ele_data_analysis$vert)[match("invertebrates",levels(ele_data_analysis$vert))] <- "invert"
levels(ele_data_analysis$vert)[match("mammals",levels(ele_data_analysis$vert))] <- "vert"
levels(ele_data_analysis$vert)[match("plants",levels(ele_data_analysis$vert))] <- "plant"

# add whether a taxon flies or not
ele_data_analysis$group <- c("bees", "bees", "lizards","orthoptera", "lepidoptera",#
                             "frogs","lizards","lepidoptera", "birds", "bats", #
                             "mammals", "spiders", "bees", "bees", "amphibians",#
                             "plants", "plants", "birds", "ants", "termites",#
                             "mammals", "birds", "birds", "mammals", "beetles",#
                             "beetles","mammals", "spiders", "spiders","spiders",#
                             "mammals", "plants", "plants", "plants", "birds",#
                             "bats", "spiders", "bees", "beetles", "birds",#
                             "spiders", "snails", "beetles", "spiders", "mammals",#
                             "amphibians","amphibians","reptiles","birds","birds",#
                             "birds", "bats", "beetles", "bees", "bees", "spiders",#added Nemesio 2010; 
                             "orthoptera","beetles", "beetles","frogs", "snails",
                             "plants","bees", "plants","plants","bees",
                             "lepidoptera","mammals","lepidoptera","bats","birds",
                             "lepidoptera", "ants", "beetles", "reptiles", "lepidoptera")


ele_data_analysis$fly <- c("fly", "fly", "no_fly","fly", "fly",
                           "no_fly","no_fly","fly", "fly", "fly", 
                           "no_fly", "no_fly", "fly", "fly", "no_fly",
                           "plants", "plants", "fly", "fly", "fly",
                           "no_fly", "fly", "fly", "no_fly", "fly", #filgueiras, dung beetles fly
                           "no_fly","no_fly", "no_fly", "no_fly","no_fly", # Fujita, ground beetles, no_fly
                           "no_fly", "plants", "plants", "plants", "fly",
                           "fly", "no_fly", "fly", "no_fly", "fly",#jung, ground beetles
                           "no_fly", "no_fly", "no_fly", "no_fly", "no_fly", #knapp, ground beetles
                           "no_fly","no_fly","no_fly","fly","fly",
                           "fly", "fly", "fly", "fly", "fly", "no_fly", #added Nemesio 2010; Montgomery "fly" because we don't know the types of beetles
                           "fly","fly", "fly","no_fly", "no_fly", #Nyeko 2009 and Owen 2008 we don't know whether beetles were flying or not.
                           "plants","fly", "plants","plants","fly",
                           "fly","no_fly","fly","fly","fly",
                           "fly", "fly", "fly", "no_fly", "fly")

ele_data_analysis$animal <- c("animal", "animal", "animal","animal", "animal",
                              "animal","animal","animal", "animal", "animal", 
                              "animal", "animal", "animal", "animal", "animal",
                              "plants", "plants", "animal", "animal", "animal",
                              "animal", "animal", "animal", "animal", "animal",
                              "animal","animal", "animal", "animal","animal",
                              "animal", "plants", "plants", "plants", "animal",
                              "animal", "animal", "animal", "animal", "animal",
                              "animal", "animal", "animal", "animal", "animal",
                              "animal","animal","animal","animal","animal",
                              "animal", "animal", "animal", "animal", "animal", "animal", #added Nemesio 2010
                              "animal","animal", "animal","animal", "animal",
                              "plants","animal", "plants","plants","animal",
                              "animal","animal","animal","animal","animal",
                              "animal", "animal", "animal", "animal", "animal")
# Edwards 2010 is a study on, not on plants
ele_data_analysis$logrich <- log2(ele_data_analysis$richness)

#
data_m <- ele_data_analysis[, c(18, 19, 1, 7, 10, 15, 16, 17, 20:27)]


# models
library(glmmTMB)
library(effects)
library(AICcmodavg)


# models for SS > SL
# univariate models
model1 <- glmmTMB(SS ~ 1 + (1|dataset_id), data = data_m, family = "binomial")
model2 <- glmmTMB(SS ~ taxa + (1|dataset_id), data = data_m, family = "binomial")
model3 <- glmmTMB(SS ~ logrich + (1|dataset_id), data = data_m, family = "binomial")
model4 <- glmmTMB(SS ~ temperature + (1|dataset_id), data = data_m, family = "binomial")
model5 <- glmmTMB(SS ~ fly + (1|dataset_id), data = data_m, family = "binomial")
model6 <- glmmTMB(SS ~ evenness + (1|dataset_id), data = data_m, family = "binomial")

# bivariate, additive
model7 <- glmmTMB(SS ~ taxa + logrich + (1|dataset_id), data = data_m, family = "binomial")
model8 <- glmmTMB(SS ~ taxa + temperature + (1|dataset_id), data = data_m, family = "binomial")
model9 <- glmmTMB(SS ~ taxa + fly + (1|dataset_id), data = data_m, family = "binomial")
model10 <- glmmTMB(SS ~ taxa + evenness + (1|dataset_id), data = data_m, family = "binomial")
model11 <- glmmTMB(SS ~ logrich + temperature + (1|dataset_id), data = data_m, family = "binomial")
model12 <- glmmTMB(SS ~ logrich + fly + (1|dataset_id), data = data_m, family = "binomial")
model13 <- glmmTMB(SS ~ logrich + evenness + (1|dataset_id), data = data_m, family = "binomial")
model14 <- glmmTMB(SS ~ temperature + fly + (1|dataset_id), data = data_m, family = "binomial")
model15 <- glmmTMB(SS ~ temperature + evenness + (1|dataset_id), data = data_m, family = "binomial")
model16 <- glmmTMB(SS ~ fly + evenness + (1|dataset_id), data = data_m, family = "binomial")

#bivariate, interactive
model17 <- glmmTMB(SS ~ taxa * logrich + (1|dataset_id), data = data_m, family = "binomial")
model18 <- glmmTMB(SS ~ taxa * temperature + (1|dataset_id), data = data_m, family = "binomial")
# model 18 does not converge; hack by hand-coding the categories
model19 <- glmmTMB(SS ~ taxa * fly + (1|dataset_id), data = data_m, family = "binomial")
model20 <- glmmTMB(SS ~ taxa * evenness + (1|dataset_id), data = data_m, family = "binomial")
model21 <- glmmTMB(SS ~ logrich * temperature + (1|dataset_id), data = data_m, family = "binomial")
model22 <- glmmTMB(SS ~ logrich * fly + (1|dataset_id), data = data_m, family = "binomial")
model23 <- glmmTMB(SS ~ logrich * evenness + (1|dataset_id), data = data_m, family = "binomial")
model24 <- glmmTMB(SS ~ temperature * fly + (1|dataset_id), data = data_m, family = "binomial")
model25 <- glmmTMB(SS ~ temperature * evenness + (1|dataset_id), data = data_m, family = "binomial")
model26 <- glmmTMB(SS ~ fly * evenness + (1|dataset_id), data = data_m, family = "binomial")

model_id <- seq(1,26,1)

model_df <- c(attributes(logLik(model1))$df, attributes(logLik(model2))$df, attributes(logLik(model3))$df, 
              attributes(logLik(model4))$df, attributes(logLik(model5))$df, attributes(logLik(model6))$df,
              attributes(logLik(model7))$df, attributes(logLik(model8))$df, attributes(logLik(model9))$df, 
              attributes(logLik(model10))$df, attributes(logLik(model11))$df, attributes(logLik(model12))$df,
              attributes(logLik(model13))$df, attributes(logLik(model14))$df, attributes(logLik(model15))$df, 
              attributes(logLik(model16))$df, attributes(logLik(model17))$df, attributes(logLik(model18))$df,
              attributes(logLik(model19))$df, attributes(logLik(model20))$df, attributes(logLik(model21))$df, 
              attributes(logLik(model22))$df, attributes(logLik(model23))$df, attributes(logLik(model24))$df,
              attributes(logLik(model25))$df, attributes(logLik(model26))$df)

model_aic <- c(AIC(model1), AIC(model2), AIC(model3), AIC(model4),AIC(model5),
               AIC(model6), AIC(model7), AIC(model8), AIC(model9),AIC(model10),
               AIC(model11), AIC(model12), AIC(model13), AIC(model14),AIC(model15),
               AIC(model16), AIC(model17), AIC(model18), AIC(model19),AIC(model20),
               AIC(model21), AIC(model22), AIC(model23), AIC(model24),AIC(model25),AIC(model26))


model_ranking <- data.frame(model_id, model_df, model_aic)
model_ranking$model_aic[18] <- 4822.34 #
model_ranking$delta_aic <- - (min(model_ranking$model_aic) - model_ranking$model_aic)
model_ranking$wi_step <- exp( -0.5 * model_ranking$delta_aic)
model_ranking$wi <- model_ranking$wi_step / sum(model_ranking$wi_step)
model_ranking <- model_ranking[, -5]

#plot(allEffects(model), type = "response")
model10 <- glmmTMB(SS ~ 0 + taxa + evenness + (1|dataset_id), data = data_m, family = "binomial")
summary(model10)

library(pscl)
pR2(model10)


plot_model <- ggpredict(model10, 
                        c("evenness [all]", "taxa" ),
                        type = "re") 


# ggplot(plot_model, aes(x = x, y = predicted, fill = group, col = group))+
#   #geom_ribbon(aes(ymin = conf.low, ymax = conf.high, alpha = 0.99))+
#   geom_line(size = 3)+
#   scale_color_manual(values = c("blue", "brown", "gold", "purple", "darkgreen"))+
#   #scale_fill_manual(values = c("blue", "brown", "gold", "purple", "darkgreen"))+
#   theme_few()+
#   labs(x = "Patch size evenness", y = "Probability of SS > SL")
#   

levels(plot_model$group)
plot_model$group <- factor(plot_model$group, levels = c("amphibians & reptiles", "mammals", "invertebrates", "birds", "plants"))

plot(plot_model, #colors = "social", 
     ci = TRUE, 
     ci.style = "ribbon", 
     facets = TRUE,
     add.data = FALSE,
     collapse.group = FALSE,
     line.size = 1,
     colors = c("blue", "slateblue4", "darkmagenta", "violetred", "red")) + 
  labs(
    x = "Patch size evenness", 
    y = "Probability of SS > SL", 
    title = ""  ) + 
  labs(colour = "Taxa") +
  facet_wrap(~group, ncol = 5)+
  theme_few() +
  theme(text=element_text(size=16,  family="serif"))




## patches
patches <- data[, c(1,2,10,5)]
patches<- unique(patches)
names(patches)[1] <- "dataset_id"
patches <- merge(patches, metadata[,c(1,14)], by = "dataset_id")
patches$log10area <- log10(patches$frag_size_num)

patches$taxa <- factor(patches$taxa, levels = c("amphibians & reptiles", "mammals", "invertebrates", "birds", "plants"))



evenness_mean <- merge(evenness, mean_patch_size)
evenness_mean$log10_mean <- log10(evenness_mean$mean_patch_size)
evenness_mean <- merge(evenness_mean, metadata[,c(1,14)])

ggplot(evenness_mean, aes(x = log10_mean, y = evenness, col = taxa))+
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_few() +
  scale_x_continuous(breaks = c( -1, 0, 1, 2, 3, 4), limits =c(-1, 4),
                     labels = c("0.1","1","10","100","1.000","10.000"),
                     guide = guide_axis(n.dodge=2))+
  xlab("Log10(Patch area)") + ylab("Patch size evenness")+
  theme(text=element_text(size=20,  family="serif"))

summary(glmmTMB(evenness_mean$evenness ~ evenness_mean$log10_mean + (1|evenness_mean$taxa)))
summary(glmmTMB(evenness_mean$evenness ~ evenness_mean$log10_mean + evenness_mean$taxa)) # AIC random effect model much better
summary(glmmTMB(evenness_mean$evenness ~ evenness_mean$log10_mean * evenness_mean$taxa)) # AIC random effect model much better

ggplot(patches, aes(x=log10area, fill=taxa)) +
  geom_histogram( color="#e9ecef", position = 'identity')+
  facet_wrap(~taxa, ncol = 5)+
  theme_few() +
  scale_y_continuous(expand = c(0,0),limits = c(0,60))+
  geom_vline(xintercept = 3, linetype="longdash", color = "grey", size=1)+
  geom_vline(xintercept = 4, linetype="longdash", color = "darkgrey", size=1)+
  geom_vline(xintercept = 4.69, linetype="longdash", color = "black", size=1)+
  scale_fill_manual(values = c("blue", "slateblue4", "darkmagenta", "violetred", "red"))+
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2, 3, 4), limits =c(-3, 5),
                     labels = c("0.01","0.1","1","10","100","1.000","10.000"),
                     guide = guide_axis(n.dodge=2))+
  xlab("Log10(Patch area)") + ylab("Number of patches")+
  theme(legend.position="none")+
  theme(text=element_text(size=20,  family="serif"))

# library(ggpubr)
# ggarrange(p1, p2)

sum(patches$frag_size_num > 1000)/nrow(patches) #27 patches
sum(patches$frag_size_num > 10000)/nrow(patches) # 2 patches
sum(patches$frag_size_num > 50000)/nrow(patches) # 0 patches


prova1 <- c(10,10,10, 100)
(vegan::diversity(list_patches[[i]]$frag_size_num))/log(length(list_patches[[i]]$frag_size_num))

prova1 <- c(50,50,100)
prova1 <- c(25,25,25,25, 100)
prova1 <- c(10,10,10,10,10,10,10,10,10,10,
            100) # 0.76 evenness

prova1 <- c(rep(3,33),100)
prova1 <- c(rep(10,10),100)
prova1 <- c(rep(5,20),100)
prova1 <- c(rep(20,5),100)
prova1 <- c(rep(50,2),100)


(vegan::diversity(prova1))/log(length(prova1)) # https://cran.r-project.org/web/packages/vegan/vignettes/diversity-vegan.pdf







# ## for SL
# # models for SS > SL
# # univariate models
# model1 <- glmmTMB(SL ~ 1 + (1|dataset_id), data = data_m, family = "binomial")
# model2 <- glmmTMB(SL ~ taxa + (1|dataset_id), data = data_m, family = "binomial")
# model3 <- glmmTMB(SL ~ logrich + (1|dataset_id), data = data_m, family = "binomial")
# model4 <- glmmTMB(SL ~ temperature + (1|dataset_id), data = data_m, family = "binomial")
# model5 <- glmmTMB(SL ~ fly + (1|dataset_id), data = data_m, family = "binomial")
# model6 <- glmmTMB(SL ~ evenness + (1|dataset_id), data = data_m, family = "binomial")
# 
# # bivariate, additive
# model7 <- glmmTMB(SL ~ taxa + logrich + (1|dataset_id), data = data_m, family = "binomial")
# model8 <- glmmTMB(SL ~ taxa + temperature + (1|dataset_id), data = data_m, family = "binomial")
# model9 <- glmmTMB(SL ~ taxa + fly + (1|dataset_id), data = data_m, family = "binomial")
# model10 <- glmmTMB(SL ~ taxa + evenness + (1|dataset_id), data = data_m, family = "binomial")
# model11 <- glmmTMB(SL ~ logrich + temperature + (1|dataset_id), data = data_m, family = "binomial")
# model12 <- glmmTMB(SL ~ logrich + fly + (1|dataset_id), data = data_m, family = "binomial")
# model13 <- glmmTMB(SL ~ logrich + evenness + (1|dataset_id), data = data_m, family = "binomial")
# model14 <- glmmTMB(SL ~ temperature + fly + (1|dataset_id), data = data_m, family = "binomial")
# model15 <- glmmTMB(SL ~ temperature + evenness + (1|dataset_id), data = data_m, family = "binomial")
# model16 <- glmmTMB(SL ~ fly + evenness + (1|dataset_id), data = data_m, family = "binomial")
# 
# #bivariate, interactive
# model17 <- glmmTMB(SL ~ taxa * logrich + (1|dataset_id), data = data_m, family = "binomial")
# model18 <- glmmTMB(SL ~ taxa * temperature + (1|dataset_id), data = data_m, family = "binomial")
# model19 <- glmmTMB(SL ~ taxa * fly + (1|dataset_id), data = data_m, family = "binomial")
# model20 <- glmmTMB(SL ~ taxa * evenness + (1|dataset_id), data = data_m, family = "binomial")
# model21 <- glmmTMB(SL ~ logrich * temperature + (1|dataset_id), data = data_m, family = "binomial")
# model22 <- glmmTMB(SL ~ logrich * fly + (1|dataset_id), data = data_m, family = "binomial")
# model23 <- glmmTMB(SL ~ logrich * evenness + (1|dataset_id), data = data_m, family = "binomial")
# model24 <- glmmTMB(SL ~ temperature * fly + (1|dataset_id), data = data_m, family = "binomial")
# model25 <- glmmTMB(SL ~ temperature * evenness + (1|dataset_id), data = data_m, family = "binomial")
# model26 <- glmmTMB(SL ~ fly * evenness + (1|dataset_id), data = data_m, family = "binomial")
# 
# model_id <- seq(1,26,1)
# 
# model_df <- c(attributes(logLik(model1))$df, attributes(logLik(model2))$df, attributes(logLik(model3))$df,
#               attributes(logLik(model4))$df, attributes(logLik(model5))$df, attributes(logLik(model6))$df,
#               attributes(logLik(model7))$df, attributes(logLik(model8))$df, attributes(logLik(model9))$df,
#               attributes(logLik(model10))$df, attributes(logLik(model11))$df, attributes(logLik(model12))$df,
#               attributes(logLik(model13))$df, attributes(logLik(model14))$df, attributes(logLik(model15))$df,
#               attributes(logLik(model16))$df, attributes(logLik(model17))$df, attributes(logLik(model18))$df,
#               attributes(logLik(model19))$df, attributes(logLik(model20))$df, attributes(logLik(model21))$df,
#               attributes(logLik(model22))$df, attributes(logLik(model23))$df, attributes(logLik(model24))$df,
#               attributes(logLik(model25))$df, attributes(logLik(model26))$df)
# 
# model_aic <- c(AIC(model1), AIC(model2), AIC(model3), AIC(model4),AIC(model5),
#                AIC(model6), AIC(model7), AIC(model8), AIC(model9),AIC(model10),
#                AIC(model11), AIC(model12), AIC(model13), AIC(model14),AIC(model15),
#                AIC(model16), AIC(model17), AIC(model18), AIC(model19),AIC(model20),
#                AIC(model21), AIC(model22), AIC(model23), AIC(model24),AIC(model25),AIC(model26))
# 
# 
# model_ranking2 <- data.frame(model_id, model_df, model_aic)
# model_ranking2$model_aic[18] <- 3021.336 #
# model_ranking2$delta_aic <- - (min(model_ranking2$model_aic) - model_ranking2$model_aic)
# model_ranking2$wi_step <- exp( -0.5 * model_ranking2$delta_aic)
# model_ranking2$wi <- model_ranking2$wi_step / sum(model_ranking2$wi_step)
# model_ranking2 <- model_ranking2[, -5]
# 
# #plot(allEffects(model), type = "response")
# plot_model <- ggpredict(model20,
#                         c("evenness [all]", "taxa" ),
#                         type = "re")
# 
# plot(plot_model) +
#   labs(
#     x = "Patch size evenness",
#     y = "Probability of SS > SL",
#     title = ""  ) +
#   labs(colour = "Log10(species richness)")
