

bcftools query -l ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz > all.poplabels



bcftools query -l chr3_460_540_phased.vcf > chr3_460_540_phased.ids
for POP in YRI LWK GWD MSL ESN; do
    cut -f 1 -d ' ' ${POP}_inds.txt > ${POP}_IDs.txt
    grep -v -f ${POP}_IDs.txt chr3_460_540_phased_inds.txt > all_except_${POP}.txt
    grep -v -f all_except_${POP}.txt all_inds.txt > ${POP}.poplabels
done

cat MSL_IDs.txt YRI_IDs.txt GWD_IDs.txt ESN_IDs.txt LWK_IDs.txt > AFR_IDs.txt

bcftools view -S AFR_IDs.txt chr3_460_540_phased.vcf > chr3_460_540_phased_AFR.vcf

for POP in YRI LWK GWD MSL ESN; do
    cut -f 1 -d ' ' ${POP}_inds.txt > ${POP}_IDs.txt
    grep -v -f ${POP}_IDs.txt AFR_IDs.txt > all_except_${POP}.txt
    grep -v -f all_except_${POP}.txt all_inds.txt > ${POP}.poplabels
done

grep -f AFR_IDs.txt all_inds.txt > AFR.poplabels



bcftools view -S AFR_IDs.txt --regions 3 -O z -o 1000g_chr3_AFR.vcf.gz ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
mv 1000g_chr3_AFR.vcf.gz 1000g_chr3_AFR.vcf
gzip 1000g_chr3_AFR.vcf


for POP in YRI LWK GWD MSL ESN; do
    bcftools view -S ${POP}_IDs.txt > 1000g_chr3_${POP}.vcf
done

for POP in YRI LWK GWD MSL ESN; do
    bcftools view -S ${POP}_IDs.txt --regions 3:46000000-54000000 ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz > 1000g_chr3_46_54_${POP}.vcf
done
bcftools view -S $AFR_IDs.txt --regions 3:46000000-54000000 ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz > 1000g_chr3_46_54_AFR.vcf


ln -s ~/populationgenomics/project_data/chr3region/1000g_chr3_AFR.vcf

ln -s ~/populationgenomics/project_data/chr3region/1000g_chr3_46_54_AFR.vcf
ln -s ~/populationgenomics/project_data/chr3region/1000g_chr3_46_54_YRI.vcf
ln -s ~/populationgenomics/project_data/chr3region/1000g_chr3_46_54_LWK.vcf
ln -s ~/populationgenomics/project_data/chr3region/1000g_chr3_46_54_GWD.vcf
ln -s ~/populationgenomics/project_data/chr3region/1000g_chr3_46_54_MSL.vcf
ln -s ~/populationgenomics/project_data/chr3region/1000g_chr3_46_54_ESN.vcf


ln -s ~/populationgenomics/project_data/chr3region/AFR.poplabels
ln -s ~/populationgenomics/project_data/chr3region/YRI.poplabels
ln -s ~/populationgenomics/project_data/chr3region/LWK.poplabels
ln -s ~/populationgenomics/project_data/chr3region/GWD.poplabels
ln -s ~/populationgenomics/project_data/chr3region/MSL.poplabels
ln -s ~/populationgenomics/project_data/chr3region/ESN.poplabels
ln -s ~/populationgenomics/project_data/chr3region/all_except_YRI.txt
ln -s ~/populationgenomics/project_data/chr3region/all_except_LWK.txt
ln -s ~/populationgenomics/project_data/chr3region/all_except_GWD.txt
ln -s ~/populationgenomics/project_data/chr3region/all_except_MSL.txt
ln -s ~/populationgenomics/project_data/chr3region/all_except_ESN.txt
ln -s ~/populationgenomics/project_data/chr3region/human_ancestor_3.fa
ln -s ~/populationgenomics/project_data/chr3region/20140520.chr3.strict_mask.fasta.gz
ln -s ~/populationgenomics/project_data/chr3region/genetic_map_chr3_combined_b37.txt


~/populationgenomics/software/relate/bin/RelateFileFormats --mode ConvertFromVcf --haps 1000g_chr3_46_54_AFR.haps --sample 1000g_chr3_46_54_AFR.sample -i 1000g_chr3_46_54_AFR --poplabels AFR.poplabels

~/populationgenomics/software/relate/scripts/PrepareInputFiles/PrepareInputFiles.sh --haps 1000g_chr3_46_54_AFR.haps --sample 1000g_chr3_46_54_AFR.sample --ancestor human_ancestor_3.fa --mask 20140520.chr3.strict_mask.fasta.gz -o 1000g_chr3_46_54_AFR_ALL
~/populationgenomics/software/relate/bin/Relate --mode All -m 1.25e-8 -N 20000 --sample 1000g_chr3_46_54_AFR_ALL.sample.gz --haps 1000g_chr3_46_54_AFR_ALL.haps.gz --map genetic_map_chr3_combined_b37.txt --annot 1000g_chr3_46_54_AFR_ALL.annot --dist 1000g_chr3_46_54_AFR_ALL.dist.gz --memory 20 -o 1000g_chr3_46_54_AFR_ALL
~/populationgenomics/software/relate/scripts/EstimatePopulationSize/EstimatePopulationSize.sh -m 1.25e-8 -N 20000 -i 1000g_chr3_46_54_AFR_ALL --poplabels LWK.poplabels -o 1000g_chr3_46_54_AFR_ALL_popsize --threshold 0 --num_iter 5 --years_per_gen 29 --threads 14 --threshhold 0
~/populationgenomics/software/relate/scripts/DetectSelection/DetectSelection.sh -i 1000g_chr3_46_54_AFR_ALL -m 1.25e-8 --poplabels LWK.poplabels -o 1000g_chr3_46_54_AFR_ALL_selection

# for POP in YRI LWK GWD MSL ESN; do
for POP in ESN; do
    ~/populationgenomics/software/relate/scripts/PrepareInputFiles/PrepareInputFiles.sh --haps 1000g_chr3_46_54_AFR.haps --sample 1000g_chr3_46_54_AFR.sample --ancestor human_ancestor_3.fa --mask 20140520.chr3.strict_mask.fasta.gz --remove_ids all_except_${POP}.txt -o 1000g_chr3_46_54_AFR_${POP} 
    ~/populationgenomics/software/relate/bin/Relate --mode All -m 1.25e-8 -N 20000 --sample 1000g_chr3_46_54_AFR_${POP}.sample.gz --haps 1000g_chr3_46_54_AFR_${POP}.haps.gz --map genetic_map_chr3_combined_b37.txt --annot 1000g_chr3_46_54_AFR_${POP}.annot --dist 1000g_chr3_46_54_AFR_${POP}.dist.gz --memory 20 -o 1000g_chr3_46_54_AFR_${POP} 
    ~/populationgenomics/software/relate/scripts/EstimatePopulationSize/EstimatePopulationSize.sh -m 1.25e-8 -N 20000 -i 1000g_chr3_46_54_AFR_${POP} --poplabels ${POP}.poplabels -o 1000g_chr3_46_54_AFR_${POP}_popsize --threshold 0 --num_iter 5 --years_per_gen 29 --threads 14 --threshhold 0 
    /populationgenomics/software/relate/scripts/DetectSelection/DetectSelection.sh -i 1000g_chr3_46_54_AFR_${POP} -m 1.25e-8 --poplabels ${POP}.poplabels -o 1000g_chr3_46_54_AFR_${POP}_selection
done



# 
~/populationgenomics/software/relate/bin/RelateFileFormats --mode ConvertFromVcf --haps 1000g_chr3_AFR.haps --sample 1000g_chr3_AFR.sample -i 1000g_chr3_AFR --poplabels AFR.poplabels

~/populationgenomics/software/relate/scripts/PrepareInputFiles/PrepareInputFiles.sh --haps 1000g_chr3_AFR.haps --sample 1000g_chr3_AFR.sample --ancestor human_ancestor_3.fa --mask 20140520.chr3.strict_mask.fasta.gz -o 1000g_chr3_AFR_ALL
~/populationgenomics/software/relate/bin/Relate --mode All -m 1.25e-8 -N 20000 --sample 1000g_chr3_AFR_ALL.sample.gz --haps 1000g_chr3_AFR_ALL.haps.gz --map genetic_map_chr3_combined_b37.txt --annot 1000g_chr3_AFR_ALL.annot --dist 1000g_chr3_AFR_ALL.dist.gz --memory 20 -o 1000g_chr3_AFR_ALL
~/populationgenomics/software/relate/scripts/EstimatePopulationSize/EstimatePopulationSize.sh -m 1.25e-8 -N 20000 -i 1000g_chr3_AFR_ALL --poplabels LWK.poplabels -o 1000g_chr3_AFR_ALL_popsize --threshold 0 --num_iter 5 --years_per_gen 29 --threads 14 --threshhold 0
~/populationgenomics/software/relate/scripts/DetectSelection/DetectSelection.sh -i 1000g_chr3_AFR_ALL -m 1.25e-8 --poplabels LWK.poplabels -o 1000g_chr3_AFR_ALL_selection

for POP in YRI LWK GWD MSL ESN; do
    ~/populationgenomics/software/relate/scripts/PrepareInputFiles/PrepareInputFiles.sh --haps 1000g_chr3_AFR.haps --sample 1000g_chr3_AFR.sample --ancestor human_ancestor_3.fa --mask 20140520.chr3.strict_mask.fasta.gz --remove_ids all_except_${POP}.txt -o 1000g_chr3_AFR_${POP}
    ~/populationgenomics/software/relate/bin/Relate --mode All -m 1.25e-8 -N 20000 --sample 1000g_chr3_AFR_${POP}.sample.gz --haps 1000g_chr3_AFR_${POP}.haps.gz --map genetic_map_chr3_combined_b37.txt --annot 1000g_chr3_AFR_${POP}.annot --dist 1000g_chr3_AFR_${POP}.dist.gz --memory 20 -o 1000g_chr3_AFR_${POP}
    ~/populationgenomics/software/relate/scripts/EstimatePopulationSize/EstimatePopulationSize.sh -m 1.25e-8 -N 20000 -i 1000g_chr3_AFR_${POP} --poplabels ${POP}.poplabels -o 1000g_chr3_AFR_${POP}_popsize --threshold 0 --num_iter 5 --years_per_gen 29 --threads 14 --threshhold 0
    ~/populationgenomics/software/relate/scripts/DetectSelection/DetectSelection.sh -i 1000g_chr3_AFR_${POP} -m 1.25e-8 --poplabels ${POP}.poplabels -o 1000g_chr3_AFR_${POP}_selection
done



# ln -s ~/populationgenomics/project_data/chr3region/chr3_460_540_phased_AFR.vcf

# ~/populationgenomics/software/relate/bin/RelateFileFormats --mode ConvertFromVcf --haps chr3_460_540_phased_AFR.haps --sample chr3_460_540_phased_AFR.sample -i chr3_460_540_phased_AFR --poplabels AFR.poplabels
# ~/populationgenomics/software/relate/scripts/PrepareInputFiles/PrepareInputFiles.sh --haps chr3_460_540_phased_AFR.haps --sample chr3_460_540_phased_AFR.sample --ancestor human_ancestor_3.fa --mask 20140520.chr3.strict_mask.fasta.gz --remove_ids all_except_LWK.txt -o chr3_460_540_phased_AFR_LWK
# ~/populationgenomics/software/relate/bin/Relate --mode All -m 1.25e-8 -N 20000 --sample chr3_460_540_phased_AFR_LWK.sample.gz --haps chr3_460_540_phased_AFR_LWK.haps.gz --map genetic_map_chr3_combined_b37.txt --annot chr3_460_540_phased_AFR_LWK.annot --dist chr3_460_540_phased_AFR_LWK.dist.gz --memory 20 -o chr3_460_540_phased_AFR_LWK
# ~/populationgenomics/software/relate/scripts/EstimatePopulationSize/EstimatePopulationSize.sh -m 1.25e-8 -N 20000 -i chr3_460_540_phased_AFR_LWK --poplabels LWK.poplabels -o chr3_460_540_phased_AFR_LWK_popsize --threshold 0 --num_iter 5 --years_per_gen 29 --threads 14 --threshhold 0
# ~/populationgenomics/software/relate/scripts/DetectSelection/DetectSelection.sh -i chr3_460_540_phased_AFR_LWK -m 1.25e-8 --poplabels LWK.poplabels -o chr3_460_540_phased_AFR_LWK_selection


