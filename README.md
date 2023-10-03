# Project Abstract
Overfishing remains a threat to coral reef fishes worldwide, with large carnivores often disproportionately vulnerable. Marine protected areas (MPAs) can restore fish populations and biodiversity, but their impact has been understudied in mesophotic coral ecosystems (MCEs – extensions of shallow reefs from 30 to 150 m) particularly in countries of the Coral Triangle – the global center of marine biodiversity. We analyzed videos from baited remote underwater video systems deployed in 2016 to investigate the assemblage structure of large carnivorous fishes at shallow (4-12 m) and mesophotic (45-96 m) depths in two of the largest and most isolated MPAs in the Philippines: Tubbataha Reefs Natural Park (TRNP), a fully no-take MPA enacted in 1988; and Cagayancillo, an archipelagic municipality surrounded by an extensive but not fully no-take MPA declared in 2016. We focused on groupers (Serranidae), snappers (Lutjanidae), emperors (Lethrinidae), jacks (Carangidae); and the endangered Cheilinus undulatus (Labridae). Surprisingly, mean fish abundance and species richness were not greater in TRNP than Cagayancillo regardless of depth despite long-term protection in the former. Limited impacts of fishing in Cagayancillo may largely explain this result. Differentiation of fish assemblages was evident between TRNP and Cagayancillo but more obvious between depths at each location, probably due more to habitat than MPA effects. In Cagayancillo, overall carnivorous reef fish, grouper, and jack mean abundance were 2, 2, and 10 times higher, respectively, at mesophotic depths, suggesting that MCEs can serve as deep refugia from fishing. The findings emphasize that MCEs are distinct from shallow reefs, serve as important habitat for species susceptible to overfishing, and thus, must be explicitly included in the design of MPAs. This study also highlights the value of maintaining strict protection of MPAs like TRNP for the Coral Triangle and an opportunity to safeguard intact fish assemblages in Cagayancillo by expanding its no-take zones.
Keywords: BRUVs, fish assemblages, deep reefs, mesophotic reefs, predatory fish


---
# Description of Data Files and Subdirectories

## **`sandbox`** 

subdirectory contains R scripts from tutorials that were used to construct the scripts for the manuscript. It also contains first drafts of R scripts. 

## **`Misc. Figures`** 

subdirectory contains drafts of supplementary figures and first drafts of figures that were not included in the main manuscript. 

Every `.png` file was used as a figure for the final version of the manuscript. 

Data files used in R scripts to create figures in the manuscript are **meso_euphotic_carniv_fish_videobaitstations_all.rds**, **WorkingData_CLEANED_TUB,CAG.xls**x, and **PHIRES_MetaData.xlsx**. We used the file **model_fitting_functions.R **as a function path. 

## **`data_wrangling_vis_salvador_ordination_rarefaction.R`** 

Input Files 
	* data: `./WorkingData_CLEANED_TUB,CAG.xlsx`
	* `./PHIRES_MetaData.xlsx`

R script used to create the nMDS ordination to look at the differences in fish assemblage between shallow and mesophotic reefs and between Tubbataha Reefs National Park (TRNP) and Cagayancillo. Similarities of fish assemblage structure across the two study locations and depths was quantified using the Bray-Curtis Dissimilarity Index. From that matrix and the vegan::metaMDS command, we generated a non-metric multidimensional scaling (nMDS) to visualize the assemblage structure. We did have to remove one of the sampling sites from Cagayancillo (CAG_017 listed in the data files) as an outlier. We tested for the differences in fish assemblages between the study locations and depths using PERMANOVA. In this R script, we also generated species rarefaction curve based on the abundance of large carnivorous fishes. We used non-parametric Chao1 estimator for overall species richness and the estaccumR command within the vegan package to generate the rarefaction curves. Individual rarefaction curves were created for each treatment combination: shallow reefs at TRNP, shallow reefs at Cagayancillo, MCEs at TRNP, and MCEs at Cagayancillo. Individual rarefaction curves were also generated for each taxonomic grouping. Since Cheilinus undulatus is only one species, there was not an individual rarefaction curve generated for it. 

Within both files hypothesistesting_MaxN.R and hypothesistesting_speciesrichness.R, the mixed command in the afex R package was used to test the effects of the study location and depth category on species richness and abundance (Singmann et al., 2017). We used the following statistical formula: y ~ depth_category * study_location + (1|study_location: bait_type). Bait type was used as random blocking factor nested within the MPA because no bait types were shared between the two study locations. After each statistical distribution was chosen to satisfy the assumptions of the model, we calculated the estimated marginal means and confidence intervals with the emmeans R package (Lenth & Love, 2017). The command emmeans::contrast was used to test for differences between treatment combinations, and the false discovery rate was controlled at .05. The multcomp:cld command was used to label and separate significantly different treatment combinations from the estimated marginal means. 

## **`hypothesistesting_MaxN.R`** 

is the R script used to make the barplot of the estimated marginal Means of MaxN between TRNP and Cagayancillo and the faceted version of that barplot among our six taxonomic groupings: Cheilinus undulatus, Serranidae, Lutjanidae, Lethrinidae, Carangidae, and Galeomorphii. Since MaxN is count data, we fit the overall MaxN data and the MaxN data between groupings to the Poisson distribution. 

## **`hypothesistesting_speciesrichness.R`** 

is the R script used to make the barplot of the estimated marginal means of species richness between TRNP and Cagayancillo and the facted version of that barplot among five taxonomic groupings: Serranidae, Lutjanidae, Lethrinidae, Carangidae, and Galeomorphii. We used non-parametric Chao1 estimator for overall species richness. However, when dividing the species richness data per taxonomic grouping, there did not seem to be enough data to use the Chao1 estimator. Therefore, instead of using the Chao1 estimator, we used the observed number of species. We fit the overall Chao1 estimate of species richness to the Gamma distribution and the number of species observations for each taxonomic grouping to the Poisson distribution.


---

# Work Completed

* [Chris's Log of Changes Made on 9/14/2022](log_ceb_2022-09-14.md)
* [Chris's Log of Changes Made on 7/26/2022](log_ceb_2022-07-26.md)
