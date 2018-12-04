echo "Please run the script with the absolute path to project_dockerv2/LL_pipeline/project directory"
echo "For example, if this is located at your home directory than it should be"
echo "/home/favorite_user_name/some/directories/project_dockerv2/LL_pipeline/"
echo "Script will add this prefix to the relative paths"
echo ">To input run it like sh firestarter.sh /home/favorite_user_name/some/directories/project_dockerv2/LL_pipeline/ \n\n\n\n\n"
mypath=$1

echo '!!! --- Welcome to our initial pipeline --- !!!'
echo '!! First stage is starting !!'
### Phase 1: Mutation Number Calculation & Randomsation ###
python -u $mypath/project/codes/sfsdf_short.py -p "$mypath/project/Batch1" -l "$mypath/project/bedfiles/Clinical_AR_raw.bed" -it 5 -vf strelka -ref $mypath
echo '!! First stage is done !!'

echo '!! Second stage is starting. TF mutation density analysis initiated !!'
### Phase x: Mutation numbers of each TF ###
mkdir mutation_numbers_per_TF
for i in $mypath/project/bedfiles/*.bed;
do
	mybedfilename=`basename $i`
	for b in ./final/run1/annotated_full*.vcf;
	do
		myvcffilename=`basename $b`
		echo "Intersecting $mybedfilename and $myvcffilename"
		bedtools intersect -a $i -b $b -c > ./mutation_numbers_per_TF/$mybedfilename\_$myvcffilename\.txt 
	done
done

Rscript $mypath/project/codes/codetobuildthistable_peakOverlapWithAR.R
Rscript $mypath/project/codes/COLLABORATION_Mutation_TF_dotplot.R
echo '!! TF mutation density analysis finished !!\n\n\n'

echo '!! Long distance mutation density analysis initiated !! '
mkdir  -p long_range_density_analysis/extended_beds
mkdir -p long_range_density_analysis/mutation_files_in_txt
for i in $mypath/project/bedfiles/*.bed;
do 
	mybedfilename=`basename $i`
	echo "Extending $mybedfilename 5000bp from both sides"
	bedtools slop -i $i -b 5000 -g $mypath/project/reference_files/hg37.chrom.sizes.txt > long_range_density_analysis/extended_beds/$mybedfilename\_5kbothsides.bed
	
	for b in ./final/run1/annotated_full*.vcf;
        do
        myvcffilename=`basename $b`
        echo "Intersecting $mybedfilename with $myvcffilename"
	bedtools intersect -a ./long_range_density_analysis/extended_beds/$mybedfilename\_5kbothsides.bed -b $b >long_range_density_analysis/mutation_files_in_txt/$mybedfilename\_$myvcffilename\.txt
	done
	
	Rscript $mypath/project/codes/COLLABORATION_mutation_location_density.R long_range_density_analysis/mutation_files_in_txt/$mybedfilename\_$myvcffilename\.txt \
	long_range_density_analysis/extended_beds/$mybedfilename\_5kbothsides.bed $mybedfilename\_$myvcffilename\_analysis long_range_density_analysis/
done

echo '!! All three analyses finished. Bye bye! :D !!'
