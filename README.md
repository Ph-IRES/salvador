# Project Abstract
In the Philippines, overfishing has caused large declines in abundance of many fish species, with large carnivorous fish often being the most vulnerable. In response, marine protected areas (MPAs) have been implemented throughout the country, but the impact of long-term MPA protection has been understudied in mesophotic coral ecosystems (MCEs). MCEs may be critical in providing refuge from overfishing and storms, but this hypothesis has rarely been tested. We analyzed baited remote underwater videos (n = 65) from 2016 to investigate fish assemblage structure between mesophotic and shallow coral reefs at two different study locations: 1) Tubbataha Reefs Natural Park (TRNP), which has been a strictly protected, fully no-take MPA since 1988 and 2) Cagayancillo, an island group which was recently declared as an MPA in 2016. We identified and recorded the abundance of large carnivorous fish species from the families Serranidae, Lutjanidae, Lethrinidae, and Carangidae; the critically endangered species Cheilinus undulatus; and sharks from the superorder Galeomorphii. Species richness was higher in MCEs than shallow reefs at both locations, and the relative abundance was highest in Cagayancillo MCEs. Our results indicate that MCEs can serve as a depth refuge because Cagayancillo exhibited equal or greater differences in richness, abundance, and composition between MCEs and shallow reefs compared to TRNP. At TRNP, there was little difference in abundance between depths, suggesting that abundance may have recovered after twenty-eight years of protection. Our study provides a glimpse into how long-term MPA protection can positively impact large carnivorous fish assemblages and preserve biodiversity.
Keywords: BRUVs, fish assemblages, deep reefs, mesophotic reefs, predatory fish


---

##Description of Data Files and Subdirectories ##

sandbox subdirectory contains R scripts from tutorials that were used to construct the scripts for the manuscript. It also contains first drafts of R scripts. 

Misc. Figures subdirectory contains drafts of supplementary figures and first drafts of figures that were not included in the main manuscript. 

Every .png file was used as a figure for the final version of the manuscript. 

Files used in R scripts to create figures in the manuscript are meso_euphotic_carniv_fish_videobaitstations_all.rds, WorkingData_CLEANED_TUB,CAG.xlsx, and PHIRES_MetaData.xlsx. We used the file model_fitting_functions.R as a function path. 

data_wrangling_vis_salvador_ordination_rarefaction.R is the R script used to create the nMDS ordination to look at the differences in fish assemblage between shallow and mesophotic reefs and between Tubbataha Reefs National Park (TRNP) and Cagayancillo. Similarities of fish assemblage structure across the two study locations and depths was quantified using the Bray-Curtis Dissimilarity Index. From that matrix and the vegan::metaMDS command, we generated a non-metric multidimensional scaling (nMDS) to visualize the assemblage structure. We did have to remove one of the sampling sites from Cagayancillo (CAG_017 listed in the data files) as an outlier. We tested for the differences in fish assemblages between the study locations and depths using PERMANOVA. In this R script, we also generated species rarefaction curve based on the abundance of large carnivorous fishes. We used non-parametric Chao1 estimator for overall species richness and the estaccumR command within the vegan package to generate the rarefaction curves. Individual rarefaction curves were created for each treatment combination: shallow reefs at TRNP, shallow reefs at Cagayancillo, MCEs at TRNP, and MCEs at Cagayancillo. Individual rarefaction curves were also generated for each taxonomic grouping. Since Cheilinus undulatus is only one species, there was not an individual rarefaction curve generated for it. 

Within both files hypothesistesting_MaxN.R and hypothesistesting_speciesrichness.R, the mixed command in the afex R package was used to test the effects of the study location and depth category on species richness and abundance (Singmann et al., 2017). We used the following statistical formula: y ~ depth_category * study_location + (1|study_location: bait_type). Bait type was used as random blocking factor nested within the MPA because no bait types were shared between the two study locations. After each statistical distribution was chosen to satisfy the assumptions of the model, we calculated the estimated marginal means and confidence intervals with the emmeans R package (Lenth & Love, 2017). The command emmeans::contrast was used to test for differences between treatment combinations, and the false discovery rate was controlled at .05. The multcomp:cld command was used to label and separate significantly different treatment combinations from the estimated marginal means. 

hypothesistesting_MaxN.R is the R script used to make the barplot of the estimated marginal Means of MaxN between TRNP and Cagayancillo and the faceted version of that barplot among our six taxonomic groupings: Cheilinus undulatus, Serranidae, Lutjanidae, Lethrinidae, Carangidae, and Galeomorphii. Since MaxN is count data, we fit the overall MaxN data and the MaxN data between groupings to the Poisson distribution. 

hypothesistesting_speciesrichness.R is the R script used to make the barplot of the estimated marginal means of species richness between TRNP and Cagayancillo and the facted version of that barplot among five taxonomic groupings: Serranidae, Lutjanidae, Lethrinidae, Carangidae, and Galeomorphii. We used non-parametric Chao1 estimator for overall species richness. However, when dividing the species richness data per taxonomic grouping, there did not seem to be enough data to use the Chao1 estimator. Therefore, instead of using the Chao1 estimator, we used the observed number of species. We fit the overall Chao1 estimate of species richness to the Gamma distribution and the number of species observations for each taxonomic grouping to the Poisson distribution.


---

# Work Completed

* [Chris's Log of Changes Made on 7/26/2022](log_ceb_2022-07-26.md)
