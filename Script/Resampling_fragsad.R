
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
data = fread("C:\\Users\\feder\\OneDrive\\Desktop\\GITHUB\\SLOSS_3\\Data\\fragSAD_predicts_ewers.csv", header = TRUE) # version from Chase et al. 2020
metadata = fread("C:\\Users\\feder\\OneDrive\\Desktop\\GITHUB\\SLOSS_3\\Data\\fragSAD_metadata.csv", header = TRUE)

metadata$dataset_id <- as.factor(metadata$dataset_id)
data$dataset_label <- as.factor(data$dataset_label)

# list of paper to be retained based on assessment
papers_reassessment = fread("C:\\Users\\feder\\OneDrive\\Desktop\\GITHUB\\SLOSS_3\\Data\\papers_reassessment.csv", header = TRUE)
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

###
###
###

## example of resampled metacommunity
list_sampled_data[[1]][["Aguiar_2012"]]
list_sampled_data[[2]][["Aguiar_2012"]]
list_sampled_data[[3]][["Aguiar_2012"]]


