# Set working directory to where *20bq.20mq*.txt 
# from CoverageOfKeyGenes output and coverage csv file 
# are on computer
# *20bq.20mq*.txt has the per base coverage for regions of interest

setwd("~/OneDrive-UPMC/Germline-Validation/ACMG/WGS-underperforming-regions/data")

# First create the normalization table
# Import mean coverage for each sample measured
# by Picard CollectWgsMetrics 
# (taken from GoogleDrive Sentieon QC tracking spreadsheet)

all_sample_cov<-read.csv("sample_cov.csv", header=FALSE)
colnames(all_sample_cov)<-c("Sample", "Mean_Coverage")
head(all_sample_cov)

# Calculate normailzation factor based on min coverage cutoff (20)
all_sample_cov$norm_cov<-20/all_sample_cov$Mean_Coverage

# For each samaple 20bq.20mq*.txt file, multiple the actual depth
# by the normalization factor calculated in the step above

# Get list of sample files
myFiles<-list.files(".", "*.txt")

# Loop through files, normalize depth, and output new normalized depth file
for (x in 1:length(myFiles)){
  
  #Get Sample Name
  sample_name<-tools::file_path_sans_ext(myFiles[[x]])
  sample_name<-sub("20bq.20mq.txt_", "", sample_name)
  
  #Import coverage data from CoverageOfKeyGenes Output
  normalized_cov<-read.table(myFiles[x], header=TRUE, stringsAsFactors = FALSE)
  
  #Find normalization factor by sample name
  norm_factor<-all_sample_cov[all_sample_cov==sample_name,]$norm_cov
  
  #Normalize depth of coverage at each base
  normalized_cov$Total_Depth<-round(norm_factor*normalized_cov$Total_Depth)
  normalized_cov$Average_Depth_sample<-round(norm_factor*normalized_cov$Average_Depth_sample)
  normalized_cov[,4]<-round(norm_factor*normalized_cov[,4])
  
  #Write new normalized file
  write.table(normalized_cov, paste(sample_name, "_normalized.txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE)
}

##########################################################################
# Reprocess using partial Coverage of Key Genes Script before next step 
# using the new 20bq.20mq.*normalized.txt file -
# run coverage_or_key_genes_normalize_cov.sh
##########################################################################


# Import the reprocessed normalized data to determine overall 
# Underperforming Regions for WGS test

# Set working directory to location where the new *normalized_coverage_metrics.txt
# files are from CoverageOfKeyGenes

setwd("~/OneDrive-UPMC/Germline-Validation/ACMG/WGS-underperforming-regions/data/normalized/coverage")

# Get list of files in directory
myFiles<-list.files(".", "coverage_metrics.txt")

# Set counter to 0
x=0

# Loop through summary coverage_metrics.txt files and create
# One dataframe that includes the fractionMissing column 
# for each sample

for (x in 1:length(myFiles)){
  
  #Get Sample Name
  sample_name<-tools::file_path_sans_ext(myFiles[[x]])
  sample_name<-sub(".10cov.coverage_metrics", "", sample_name)
  
  
  #Import normalized coverage data from CoverageOfKeyGenes Output
  dat<-read.table(myFiles[x], header=TRUE)
  
  # initiate the dataframe with first sample
  if (x==1){ 
    all_samples_fracMissing<-sample_bed[,c(1,3)]
    colnames(all_samples_fracMissing)<-c("Exon", sample_name)}
  
  # Bind additional sample columns to dataframe
  else{
    all_samples_fracMissing<-cbind(all_samples_fracMissing, dat[,3])
    colnames(all_samples_fracMissing)[x+1]<-sample_name
  }
  
}

# Transpose dataframe so columns are positions and rows are samples

dat_transpose<-as.data.frame(t(all_samples_fracMissing[,-1]))
colnames(dat_transpose)<-all_samples_fracMissing$Exon
summary(dat_transpose)$mean

# Can't find a good way to store the summary stats data
# as a dataframe so just export and open in excel
test<-as.data.frame(summary(test_transpose))

# Save summary stats for fracMisssing per region
write.csv(test, "all_samples_fracMissing_summary_data.csv")

# Save the merged dataframe containing fraction mission for each position in every samples
write.csv(all_samples_fracMissing, "all_samples_fracMissing_matrix.csv", quote=FALSE, row.names = FALSE)

