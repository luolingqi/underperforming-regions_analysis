#!/bin/bash
# process the normalized files where read depth at each position was mutliplied by 20/mean_sample_cov

gene_list_path="/Users/samanthamentch/OneDrive-UPMC/Germline-Validation/ACMG/WGS-underperforming-regions/ACMG59_preferred_transcripts.txt.format.mod.RefSeqCurated.exon.sum.no-alt.bed"

files=$(ls *.txt)
for file in $files; do
	bam_prefix="${file%.*}"


		################################################
		# 2. Turn GATK output text file into a bed file
		################################################

		awk -v OFS="\t" -v c=10 'NR>1{if($2 >= c){split($1, a, ":"); print a[1],a[2]-1,a[2]}}' $file   >  coverage/"$bam_prefix".bed


		################################################
		# 3. Intersect output bed with orig gene list
		################################################
		words=$(wc -l < coverage/"$bam_prefix".bed)

		if [ $words -eq 0 ];then 
			awk -v OFS="\t" '{arr[$4]+=0} END {for (i in arr) {print i, arr[i]}}' $gene_list_path > coverage/"$bam_prefix".genes.sum.bed  
			 # echo "None of the bases in any of your genes pass your thresholds. Make sure your gene coordinates are correct or try more lenient thresholds."
			  #  exit 1
			   
			  else 
			  	bedtools intersect  -a $gene_list_path -b coverage/"$bam_prefix".bed -wao > coverage/"$bam_prefix".genes.bed


		################################################
		# 4. Sum By Gene
		################################################

		awk -v OFS="\t" '{arr[$4]+=$8} END {for (i in arr) {print i, arr[i]}}' coverage/"$bam_prefix".genes.bed > coverage/"$bam_prefix".genes.sum.bed
		fi

		sort coverage/"$bam_prefix".genes.sum.bed >  coverage/"$bam_prefix".genes.sum.sorted.bed
		awk -v OFS="\t" '{arr[$4]+=$3-$2} END {for (i in arr) {print i,arr[i]}}' $gene_list_path  > coverage/totByGene.bed
		sort coverage/totByGene.bed > coverage/sorted.totByGene.bed
				  
		join -1 1 -2 1 -t $'\t' coverage/"$bam_prefix".genes.sum.sorted.bed coverage/sorted.totByGene.bed > coverage/"$bam_prefix".summary.txt
				    
				      
		awk -v OFS="\t" 'BEGIN {print "Gene\tnumMissing\tfracMissing"} {tot_fail = $3-$2; frac_fail = tot_fail / $3; print $1,tot_fail,frac_fail }' coverage/"$bam_prefix".summary.txt > coverage/"$bam_prefix".10cov.coverage_metrics.txt
done

