library("dplyr")
library("tidyr")
library("reshape2")
library("stringr")
library("ggplot2")

hplc_msms_analysis <- function(stdfile,
                               mgfile,
                               pythonfile,
                               plot,
                               output){
  # Load the concentration values file, the mg_sample file, & the metabolite data from Python. Change 
  std_concentrations <- read.csv(file = stdfile)
  mg_sample <- read.csv(file = mgfile)
  data <- read.csv(file = pythonfile)
  
  # Reorganize python output
  data[1] = NULL # Delete first column from dataframe
  
  # Invert data to make each sample name the 1st column and each corresponding area value the 2nd column
  data <- gather(data, "sample", "area")
  
  # Delete X in front, .CSV, and r1 or r2 in sample names
  data$sample = gsub("^X|.CSV|final_coding_project.rawdata.rawdata.|_r1|_r2", "", data$sample)
  
  # Average every area value with identical sample names
  data_avg <- aggregate(area~sample,data,mean)
  data_avg
  
  # Extract the metabolite label into a new column
  data_split <- str_split_fixed(data_avg$sample, "_Pos_", 2)
  colnames(data_split) <- c("sample1", "metabolite")
  data_split2 <- as.data.frame(data_split)
  # When split occurs, the result is a matrix, but needs to be a dataframe
  #class(data_split)
  
  # Delete the "ID" column from both dataframes
  data_avg <- tibble::rowid_to_column(data_avg, "ID")
  data_split2 <- tibble::rowid_to_column(data_split2, "ID")
  
  # Merge the metabolite column into the average data frame
  data_by_metab <- merge(data_avg, data_split2)
  data_by_metab <- subset(data_by_metab, select = -c(ID,sample))
  
  # Subset the water values, save as a separate dataframe
  water_metab <- subset(data_by_metab, grepl("Water", sample1))
  water_metab <- subset(water_metab, select = -c(sample1))
  names(water_metab)[names(water_metab) == 'area'] <- 'water_area'
  
  # Delete Water rows, subtract water area from samples area, make new column with normalized area
  data_norm <- merge(data_by_metab,water_metab,sort=FALSE)
  data_norm <- data_norm[data_norm$sample1 != "Water", ]
  data_norm$area_norm <- (data_norm$area - data_norm$water_area)
  data_norm <- subset(data_norm, select = -c(area,water_area))
  
  # Add column with concentrations for each standard and metabilite
  data_conc <- full_join(data_norm, std_concentrations, by = c("sample1" = "sample1", "metabolite" = "metabolite"))
  
  # Extract only Std information
  Std <- subset(data_conc, grepl("PosStd", sample1))
  #Std_cre <- subset(Std, grepl("creatine", metabolite))
  
  # Calculate linear regression for each set of standards, pull out the slope values for each metabolite and save as dataframe
  slopes <- c()
  for(i in unique(Std$metabolite)) {
    temp <- subset(Std, Std$metabolite == i)
    temp_lm = lm(area_norm ~ concentration, data = temp)
    result <- temp_lm$coefficients[2]
    names(result) <- i
    slopes <- c(slopes, result)
  }
  
  # Convert vector to list, then to data frame
  slopes <- as.list(slopes)
  slopes <- as.data.frame(slopes)
  slopes <- t(slopes)
  slopes <- as.data.frame(slopes)
  colnames(slopes)[colnames(slopes) == 'V1'] <- 'slope'
  slopes <- tibble::rownames_to_column(slopes, "metabolite")
  
  # Add slope values to original data frame, remove standard rows
  data_slopes <- merge(data_conc, slopes)
  data_slopes <- data_slopes[!grepl("PosStd", data_slopes$sample1),]
  
  # Normalize area values to slope value from standard curve linear regression
  data_slopes$pg_uL_extract <- (data_slopes$area_norm / data_slopes$slope)
  
  # Normalize pg/uL extract values to 25uL extract, and from picograms to nanograms
  data_slopes$total_ng_extract <- ((data_slopes$pg_uL_extract*25)/1000)
  
  # Extract NMS values
  NMS <- subset(data_slopes, grepl("nms", metabolite))
  NMS <- select(NMS, sample1, total_ng_extract)
  NMS$ng_added <- 2.5
  NMS$extraction_efficiency <- (NMS$total_ng_extract / NMS$ng_added)
  NMS <- select(NMS, sample1, extraction_efficiency)
  #colnames(NMS)[colnames(NMS) == 'total_ng_extract'] <- 'extraction_efficiency'
  
  data_NMS <- merge(data_slopes, NMS)
  
  # Normalize total ng extracted to extraction efficiency
  data_NMS$total_ng_sample <- (data_NMS$total_ng_extract / data_NMS$extraction_efficiency)
  
  # Calculate final value of ng/mg of each metabolite in each sample; convert negative values to '0'
  data_final <- merge(data_NMS, mg_sample)
  data_final$ng_mg_sample <- (data_final$total_ng_sample / data_final$mg_sample)
  data_final <- select(data_final, sample1, metabolite, ng_mg_sample)
  data_final[data_final < 0] <- 0
  
  # Plot data
  bargraph <- ggplot(data_final,aes(x=metabolite,y=ng_mg_sample,fill=factor(sample1)))+
    geom_bar(stat="identity",position="dodge")+
    xlab("Metabolites")+ylab("ng/mg in Sample")
  pdf(plot)
  print(bargraph)
  dev.off()
  
  # Save data_final as .csv file. Change name of output file for your specifi data!
  write.csv(data_final, file = output)

}