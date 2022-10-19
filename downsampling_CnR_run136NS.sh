
#!bin/bash
#activate BEDOPS environment

source activate BEDOPS_env

# shuffle with BEDOPS

		## calculate the value to insert in SampleSIze for downsampling
			# loop over all WT samples
			# split name according to _ separator and keep second field saving in x variable
			# split x variable according to - separator and keep just the fields from 4 on containing the ab name and seller
			# find all tIMEC-ANP-H1 samples and count lines keeping just the nr
			# find all tIMEC-ANP samples and count lines keeping just the nr
			# find all tIMEC samples and count lines keeping just the nr
			
			# subtract wt and c6
			# calculate which is the sample with lower reads number between WT and C6 (lowerSample) and define the nr of reads to use for downsampling the larger file (sampleSize)
			# define which of the two file is the one to be downsampled and save it in bedelemt variable
for i in PAIRED*tIMEC-ANP32e-H1*merged_sorted.bed; do
	x=$(echo $i| cut -f2 -d_);
	ab=$(echo $x | cut -f5,6 -d-);
	tIAH1count=$(wc -l $(find PAIRED*tIMEC-ANP32e-H1-$ab*) | cut -f1 -d ' ');
	tIAcount=$(wc -l $(find PAIRED*tIMEC-ANP32e-$ab*) | cut -f1 -d ' ');
	tIcount=$(wc -l $(find PAIRED*tIMEC-$ab*) | cut -f1 -d ' ');
	# create array with counts for each cell line
	count_array=($tIAH1count, $tIAcount, $tIcount)
	bed_array=("tIMEC-ANP32e-H1-" "tIMEC-ANP32e-" "tIMEC-")
	# define minimum as first values in the array
	min=${count_array[0]}
	# loop to find smaller values
	for e in ${!count_array[@]}; do 
		echo " Processing: ${bed_array[$e]}: ${count_array[$e]}";
		if (("${count_array[$e]}" <= $min));
		then min=${count_array[$e]};
		lowerSample=$(echo ${bed_array[$e]} $ab ${count_array[$e]});
		index_lowerSample=$e;
		#echo "lower sample is: $lowerSample -- Index lower sample is: $index_lowerSample";
		else min=$min;
		#echo "lower sample is: $lowerSample -- Index lower sample is: $index_lowerSample";
		fi;
		done
	
	# Remove from the initial sample list the one corresponding to the smaller library based on the index
	unset bed_array[$index_lowerSample]

	echo "Results after looping are the following:\n LOWER SAMPLE: $lowerSample \n BED ELEMENT TO DOWNSAMPLE: ${bed_array[@]}"
	
	
	# Define files to be downsampled:
	bedelements=()
	for sample in ${bed_array[@]}; do
	bedelements+=($(find PAIRED*$sample$ab*));
	echo "${bedelements[@]}"; done
	
	# extract the size of the smaller library
	sampleSize=$(echo $lowerSample | cut -f3 -d' '); 
	echo -e $lowerSample "read nr to use:" $sampleSize "\nbedelement is:" $bedelement;
	
	# create one subfolder for every antibody
	mkdir $ab;
	
	# Calculate how reads are distributed along each chromosome and save the results in count.txt file
	# Loops along bedelements to be downsampled
	for bed in ${bedelements[@]}; do
		for chr in `bedextract --list-chr $bed`; do
		bedextract $chr $bed > $ab/$bed.$chr.bed;
		count=`wc -l $ab/$bed.$chr.bed | cut -d' ' -f1`;
		echo -e "$chr\t$count"; done > $ab/$bed.counts.txt;
	
		echo "proportions calculated for:" $ab;
	
		# Calculate the proportions of reads to keep on each chromosome depending on the final number of reads to keep (sampleSize) and save it in the proportion.txt file
		sum=`cut -f2 $ab/$bed.counts.txt | perl -nle '$sum += $_ } END { print $sum' -`;
		awk -v sum=$sum -v sampleSize=$sampleSize '{print $0"\t"($2/sum)"\t"int(sampleSize*$2/sum)}' $ab/$bed.counts.txt > $ab/$bed.proportions.txt;
	
		# do a random sampling from each chromosome based on the calculated proportions and save it into a new bed file
		awk -v samplename=$ab/$bed '{
		perChrName=$1;
		perChrN=$4;
		cmd="shuf "samplename"."perChrName".bed | head -n "perChrN;
			system(cmd);
		}' $ab/$bed.proportions.txt | tee $ab/$bed.downsampled.bed;
	
		sort-bed $ab/$bed.downsampled.bed > $ab/$bed.downsampled.sorted.bed;
	
		echo "downsampling done for:" $ab; done
	
	
