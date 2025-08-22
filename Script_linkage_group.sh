VCF=/Users/marina/Library/CloudStorage/OneDrive-Menntaský/QTL_mapping/Salpinus_align_LinkageMaps/populations.snps.vcf
PEDIGREE=/Users/marina/Library/CloudStorage/OneDrive-Menntaský/QTL_mapping/Salpinus_align_LinkageMaps/pedigree_file.txt
mkdir LM-Run-$(date '+%y%m%d@%H%M')_SeparateModules
cd LM-Run-250514\@1548_SeparateModules

# Appeler le genotype des grands parents

java -cp /Users/marina/Desktop/lepmap3/bin ParentCall2 data=$PEDIGREE vcfFile=$VCF outputParentPosterior=1 halfSibs=1 > p.call

# Filtering2 (Seulement la ligne utilisee pour les analyses)

java -cp /Users/marina/Desktop/lepmap3/bin Filtering2 data=/Users/marina/Library/CloudStorage/OneDrive-Menntaský/QTL_mapping/Salpinus_align_LinkageMaps/H19_101/LM-Run-250514@1548_SeparateModules/p.call  removeNonInformative=1 dataTolerance=0.0001 missingLimit=0.3 MAFLimit=0.05 > filtering

awk 'BEGIN { line=0; print "marker_id contig pos contig_pos"} !/#java|CHR/ { line++; print line, $1, $2, $1"_"$2}' filtering > markerlist.txt
NUMBR=$(expr $(wc -l < markerlist.txt) - 1)
printf "\n%s\n" "After filtering, there are $NUMBR markers"
java -cp /Users/marina/Desktop/lepmap3/bin SeparateChromosomes2 data=filtering numThreads=4 lodLimit=10 sizeLimit=7 distortionLod=1 > map.txt

awk '{print $1}' map.txt | sort -n | uniq -c

java -cp /Users/marina/Desktop/lepmap3/bin JoinSingles2All data=filtering lodLimit=5 sizeLimit=7 numThreads=4 iterate=1 map=map.txt > map_2.txt

awk '(NR>=7)' filtering| cut -f 1,2 > snps.txt

for chr in {1..38}; do java -cp /Users/marina/Desktop/lepmap3/bin OrderMarkers2 data=filtering map=map_2.txt numThreads=4 chromosome=${chr} outputPhasedData=1 useKosambi=1 usePhysical=1 recombination1=0.039 recombination2=0.096 calculateIntervals=interval_1.txt > order_phased_${chr}.txt; done

cut -f1,2 p.call  > cut_pcall.txt # remplacer p.call par son chemin entier s'il n'est pas dans le meme dossier

awk '(NR>=7)' filtering | cut -f 1,2 > snps.txt

awk -vfullData=1 -f /Users/marina/Library/CloudStorage/OneDrive-Menntaský/QTL_mapping/Salpinus_align_LinkageMaps/H19_101/LM-Run-250514@1548_SeparateModules/scripts/map2genotypes.awk snps.txt > genotypes.txt

for chr in {1..38}; do; awk -vfullData=1 -f /Users/marina/Library/CloudStorage/OneDrive-Menntaský/QTL_mapping/Salpinus_align_LinkageMaps/H19_101/LM-Run-250514@1548_SeparateModules/scripts/map2genotypes.awk order_phased_${chr}.txt > genotypes_${chr}.txt; done

