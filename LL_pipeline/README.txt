Hello!
To set this pipeline up:
0) Please back up your data. Its my first time working with the docker environment.
1) Download our pipeline and cd in to LL_pipeline directory.
	git clone https://github.com/mortunco/project_dockerv2.git and cd LL_pipeline
2) Build your image with. 
	$docker build -t ll_pipeline . 
3) Create following directoryies under /project/Batch1 with:
	mkdir -p project/Batch1/input/run1
	mkdir -p project/Batch1/final 
4) Move your input.vcf.gz in to /project/Batch1/input/run1/
5) Run the container. 
	$docker run --user=$UID -it  -v `pwd`\/project/:/project ll_pipeline


Notes:
“`pwd`\/project/” attaches the project its child directories to the docker container as a volume.
—user=$UID was solved my read/write permission issues.

Method Description:
###Analysis Part 1###
With the given VCF file and BED file input. Pipeline first filters the dbsnp, non-PASS mutations and homozygous events. This filtrated version is saved as annotated_full.myspecial.vcf. Later on, the intersecting mutation with respect to bed is extracted and saved to final.myspecial.vcf. To test the significance of this event, we generate random bed files with same number of peaks and peak widths across genome and test the mutation numbers.
-The result of the analysis is kept under snv_number.txt with the false discovery rate.
-Mutation numbers of the randomly generated beds are kept under snv_randomized_mutation_counts.txt
-snv_variant_classification_files* fare for ICGC data. Please do not mind.

### Analysis Part 2 ###
In calculation of mutation ratios of each TF (in ~/project/bedfiles/). Mutation numbers are found for each TF. After that they are ranked. 

###Analysis Part 3###
In this part, distribution of the mutations across TF binding regions are calculated. To view the broad picture, each peak of TF is extended by 5000bp from both sides. Then mutations are mapped on these regions. Finally, all of the peaks are comprised and general distribution of mutation on these peaks are visualised. 

