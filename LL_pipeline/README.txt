Hello!
To set this pipeline up:

Please back up your data. Its my first time working with singularity environement. 

## setup ##
In your home directory or desired directory. Download SIF file. Then create a project directory.

```$mkdir project && cd project && git clone -b singularity_adoptation https://github.com/mortunco/project_dockerv2.git```

Create neccasary file structure on Batch1. This is very important because otherwie the software will not work. Go into ```project_dockerv2/LL_pipeline/project/Batch1/```. Generate ```input/run1``` and ```final``` directories under Batch1. Finall move your compressed vcf file in to input/run1 directory.

```$mkdir -p input/run1 && mkdir final```

## Running the pipeline ##
Go back to the location where ll_pipeline_environment.sif is located and initiate shell promt with ll_pipeline_environement.

```$singularity shell ll_pipeline_environment.sif```

While in the environment. Go to the Batch1 directory. firestarter.sh bash script will initiate the series of scripts. However, because of the singularity adoptation I couldnt figure out the pathing very well. So I need your help with inputting the correct path. As stated in the firestarter.sh, absolute path to /some/path/to/project/directory/project_dockerv2/LL_pipeline/ is required. Please prove this and initiate the script.

```$sh firestarter.sh path``` 

Notes:
(Current build of the pipeline is grch38 also vcf format is set to strelka. If these are not your parameters, please contact to me.)
(I totally agree that firestarter.sh is not the best solution for the case but it is originally built for docker where it directly locates image at ~. To be able to run singularity shell without sudo, i came up with this work around. So please ask if you need assitance on this issue.)


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


