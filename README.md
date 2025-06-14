# General Description
This repository contains our BRUV data, the associated metadata, and all of the R scripts used to generate the figures in our published manuscript at the peer-reviewed journal: Aquatic Conservation: Marine and Freshwater Ecosystems. Here is the citation for our manuscript: 

Salvador, M. , Utzurrum, J., Murray, R., Delijero, K., Conales, S., Bird, C., Gauthier,,D., Abesamis, R (2024). Intact shallow and mesophotic assemblages of large carnivorous reef fishes underscore the importance of large and remote protected areas in the Coral Triangle. *Aquatic Conservation: Marine and Freshwater Ecosystems*, *34*(2), e4108. https://doi.org/10.1002/aqc.4108

A short presentation of the study's main results, which was used to present at the annual American Geophysical Union Conference, can be downloaded from the file AGUpresentation_Salvador.mp4. 

# Project Abstract
Overfishing remains a threat to coral reef fishes worldwide, with large carnivores often disproportionately vulnerable. Marine protected areas (MPAs) can restore fish populations and biodiversity, but their impact has been understudied in mesophotic coral ecosystems (MCEs – extensions of shallow reefs from 30 to 150 m) particularly in countries of the Coral Triangle – the global center of marine biodiversity. We analyzed videos from baited remote underwater video systems deployed in 2016 to investigate the assemblage structure of large carnivorous fishes at shallow (4-12 m) and mesophotic (45-96 m) depths in two of the largest and most isolated MPAs in the Philippines: Tubbataha Reefs Natural Park (TRNP), a fully no-take MPA enacted in 1988; and Cagayancillo, an archipelagic municipality surrounded by an extensive but not fully no-take MPA declared in 2016. We focused on groupers (Serranidae), snappers (Lutjanidae), emperors (Lethrinidae), jacks (Carangidae); and the endangered *Cheilinus undulatus* (Labridae). Surprisingly, mean fish abundance and species richness were not greater in TRNP than Cagayancillo regardless of depth despite long-term protection in the former. Limited impacts of fishing in Cagayancillo may largely explain this result. Differentiation of fish assemblages was evident between TRNP and Cagayancillo but more obvious between depths at each location, probably due more to habitat than MPA effects. In Cagayancillo, overall carnivorous reef fish, grouper, and jack mean abundance were 2, 2, and 10 times higher, respectively, at mesophotic depths, suggesting that MCEs can serve as deep refugia from fishing. The findings emphasize that MCEs are distinct from shallow reefs, serve as important habitat for species susceptible to overfishing, and thus, must be explicitly included in the design of MPAs. This study also highlights the value of maintaining strict protection of MPAs like TRNP for the Coral Triangle and an opportunity to safeguard intact fish assemblages in Cagayancillo by expanding its no-take zones.
Keywords: BRUVs, fish assemblages, deep reefs, mesophotic reefs, predatory fish


---
# Description of Data Files and Subdirectories

Overall, data files used in R scripts to create figures in the manuscript are **WorkingData_CLEANED_TUB,CAG.xlsx** and **PHIRES_MetaData.xlsx**. We also used the file **model_fitting_functions.R** as a function path and **meso_euphotic_carniv_fish_videobaitstations_all.rds** to clean the original data by removing rows with missing data or multiple bands. 

Within both files hypothesistesting_MaxN.R and hypothesistesting_speciesrichness.R, the mixed command in the afex R package was used to test the effects of the study location and depth category on species richness and abundance. We used the following statistical formula: y ~ depth_category * study_location + (1|study_location: bait_type). Bait type was used as random blocking factor nested within the MPA because no bait types were shared between the two study locations. After each statistical distribution was chosen to satisfy the assumptions of the model, we calculated the estimated marginal means and confidence intervals with the emmeans R package. The command emmeans::contrast was used to test for differences between treatment combinations, and the false discovery rate was controlled at .05. The multcomp:cld command was used to label and separate significantly different treatment combinations from the estimated marginal means. 

## **`hypothesistesting_MaxN.R`** 

Inputs:

R Files:
   * `./meso_euphotic_carniv_fish_videobaitstations_all.rds`
   * `./model_fitting_functions.R`

Data Files:
   * `./WorkingData_CLEANED_TUB,CAG.xlsx`
   * `./PHIRES_MetaData.xlsx`

R script used to make the barplot of the overall estimated marginal Means of MaxN between TRNP and Cagayancillo and the faceted version of that barplot among our five taxonomic groupings: Serranidae, Lutjanidae, Lethrinidae, Carangidae, and *Cheilinus undulatus*. Since MaxN is count data, we fit the overall MaxN data and the MaxN data between groupings to the Poisson distribution. 

## **`hypothesistesting_speciesrichness.R`** 

R script used to make the barplot of the overall estimated marginal means of species richness between TRNP and Cagayancillo and the faceted version of that barplot among four taxonomic groupings: Serranidae, Lutjanidae, Lethrinidae, and Carangidae. We used non-parametric Chao1 estimator for overall species richness. However, when dividing the species richness data per taxonomic grouping, there did not seem to be enough data to use the Chao1 estimator. Therefore, instead of using the Chao1 estimator, we used the observed number of species. We fit the overall Chao1 estimate of species richness to the Gamma distribution and the number of species observations for each taxonomic grouping to the Poisson distribution.

## **`data_wrangling_vis_salvador_permanova.R`** 
Input Files:
- `./WorkingData_CLEANED_TUB,CAG.xlsx`
- `./PHIRES_MetaData.xlsx`

R script used to apply the two-way PERMANOVAs to test differences in fish assemblages with the vegan::adonis2 package. The first tested the effects of depth category and study locations. The second tested the effects of depth category and bait type because bait type was not consistent across study locations and depths. Pairwise post-hoc tests were conducted using pairwiseAdonis::pairwise.adonis while accounting for a false discovery rate. 

## **`data_wrangling_vis_salvador_ordination.R`** 

Input Files:
- `./WorkingData_CLEANED_TUB,CAG.xlsx`
- `./PHIRES_MetaData.xlsx`

R script used to create the nMDS and dbRDA ordinations to look at the differences in fish assemblage between shallow and mesophotic reefs and between Tubbataha Reefs National Park (TRNP) and Cagayancillo. Similarities of fish assemblage structure across the two study locations and depths was quantified using the Bray-Curtis Dissimilarity Index. From that matrix and the vegan::metaMDS command, we generated a non-metric multidimensional scaling (nMDS) to visualize the assemblage structure. We did have to remove one of the sampling sites from Cagayancillo (CAG_017 listed in the data files) as an outlier. 

Distance-based ReDundancy Analysis (dbRDA) was also performed with vegan::dbrda to visualize the same data as used to generate the nMDS plot, constraining variance to the location and depth factors, as well as their interaction. Confidence ellipses (95%) were plotted to highlight treatment groupings.The same outlier from Cagayncillo (CAG_017) was also removed from the dbRDA plots to better discern patterns of similarity
amongst the remaining stations.

This R script was also used to conduct the SIMPER and indicspecies analyses. SIMPER was used to identify species that were influential in distinguishing fish assemblages (i.e., accounted for 75% of the total
dissimilarity) between shallow and mesophotic depths within TRNP and Cagayancillo, whereas indicspecies was used to identify indicator species associated with each of the station groups defined by the four treatment groupings (shallow-TRNP, mesophotic-TRNP,shallow-Cagayancillo and mesophotic-Cagayancillo) or combinations of these groups. In indicspecies, species were identified based on an indicator value index (IndVal), which is the
product of two quantities: A (a measure of specificity)—the mean abundance of a species in a group, or combinations of groups, compared to all other groups; B (a measure of fidelity)—the relative frequency of occurrence of a species in the stations within a group or combinations of groups.

## **`data_wrangling_vis_salvador_rarefaction.R`**

Input Files:
- `./WorkingData_CLEANED_TUB,CAG.xlsx`
- `./PHIRES_MetaData.xlsx`
  
R script used to generate individual, sample size-based rarefaction and extrapolation curves from the iNEXT package. Separate rarefaction and extrapolation curves for overall species richness were created using pooled data from all stations and data from stations within each treatment combination (shallow-TRNP, mesophotic-TRNP, shallow-Cagayancillo and mesophotic=Cagayancillo). Rarefaction and extrapolation curves were also generated for the taxonomic groups Serranidae, Lutjanidae,Lethrinidae and Carangidae using pooled data from all stations. All pooled abundance data from the stations was transformed into incidence data. Incidence data was then input into iNEXT to derive estimates of species richness (q = 0). Extrapolation models were applied to examine species richness beyond the number of BRUVs that were deployed and to identify the asymptote across treatment combinations and taxonomic groups.

Since Cheilinus undulatus is only one species, there was not an individual rarefaction curve generated for it. 

## **`supplementary_information`** 
subdirectory that contains supplementary figures and tables that were published as Supporting Information. 

## **`final_ figures`** 
subdirectory contains the drafts of figures generated from the R scripts that were included in the main manuscript. Also contains another subdirectory `.manuscript_figures` that has the final versions of the figures with the added fish silhouettes and manually updated labels and legends used in the main manuscript. 

## **`test_outputs`** 
contains outputs from indicspecies, unique taxa, and SIMPER analyses as well as post-hoc tests on bait type and depth category. 

## **`sandbox`** 

subdirectory that contains R scripts from tutorials that were used to construct the scripts for the manuscript. It also contains first drafts of R scripts. 

## **`other_R_files`** 
subdirectory that contains old drafts of R scripts that were used but eventually scrapped. Not necessary to generate final figures or do any model testing. 

## **`misc_figures`** 

subdirectory that contains drafts of figures that were NOT included in the main manuscript. 

## **`log_changes`** 
subdirectory that contains coauthor's Chris Bird's descriptions of logged changes in the Github Repository:
* [Chris's Log of Changes Made on 9/14/2022](log_ceb_2022-09-14.md)
* [Chris's Log of Changes Made on 7/26/2022](log_ceb_2022-07-26.md)


