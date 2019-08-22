# underperforming-regions_analysis

# Determining Underperforming Regions

To do the analysis
1. Run CoverageOfKeyGenes (noplot) - using the sample bam files and regions of interest bed file on DNAnexus
2. The output will contain *<bq>bq.<mq>mq.<sample_name>.txt 

   *This file contains the read depth at each position in the regions of interest

   *We will use this to create a normalized read depth file

3. Using the file above and a coverage csv file run the first part of 'determining-underperforming-regions-used-for-analysis.R' 

   *The coverage csv file is a file with 2 columns. 

     *Column 1 = Sample Name

     *Column 2 = Mean Coverage calculated using Picard CollectWgsMetrics (and found in the GoogleDrive Sentieon QC Tracking Spreadsheet)

   * The first part of the Rscript will output a new normalized *bq.mq_normalized.txt file 

4. Use the *normalized.txt file to run covverage_of_key_genes_normalized_cov.sh

   * This script will run the rest of the steps from the CoverageOfKeyGenes applet using the normalized read depth files produced in step 3 and output a summary *_coverage_metrics.txt file containing the fracMissing for each region of interest
   * fracMissing is the percent of bases that were not covered for a particular region of interest
   * This column can be used to determine the underperforming regions

5. Run the rest of the 'determining-underperforming-regions-used-for-analyis.R' script sing the *_coverage_metrics.txt file

  *Output will include 

    1. A table containing the fracMissing for each region in each sample being looked at
    2. A table containing the summary statistics for each region of interest (Mean, Min, Max, etc.)

6. I'm still trying to automate the R script to just output the underperforming regions - but for sake of time I just used excel after this point to filter the summary statistics table by fracMissing Mean >= [Cutoff (0.05)]
