# human_variation
This repo contains links to processed data and code for reproducing analyses in the [manuscript on bioRxiv](https://www.biorxiv.org/content/10.1101/2022.10.07.511366v1.full): "Inter-individual genomic and transcriptomic variation in human cortical cell types" (Johansen, Somasundaram, et al).  The repo is under construction, but will continue to be updated with scripts used to generate the figures in this manuscript in the near future.

Processed data will be hosted on [CELLxGENE](https://cellxgene.cziscience.com/collections/35928d1c-36fc-4f93-9a8d-0b921ab41745)

## Key links

* Allen Brain Map project page [(LINK)](https://knowledge.brain-map.org/data/Z4KJ7FC4QZWPP8TV4OK/summary) 
* bioRxiv preprint [(LINK)](https://www.biorxiv.org/content/10.1101/2022.10.07.511366v1) 
* Manuscript GitHub repository [(LINK)](https://github.com/AllenInstitute/human_variation) 
* Entry point for raw data (controlled access) [(LINK)](https://nemoarchive.org/resources/accessing-controlled-access-data.php) 
* CELLxGENE tool exploration, visualization, and access to processed RNA-seq data [(LINK)](https://cellxgene.cziscience.com/collections/35928d1c-36fc-4f93-9a8d-0b921ab41745) 

## Code book

Please follow the directions within the code documents below for reproducing analyses from this study.  Prior to running this analysis, you will need to download all the files from `input_files` and `scripts`, as well as the processed data files from [CELLxGENE](https://cellxgene.cziscience.com/collections/35928d1c-36fc-4f93-9a8d-0b921ab41745).

1. **R workbook for abundance and some variable gene analysis** [(LINK TO SCRIPT)](http://htmlpreview.github.io/?https://github.com/AllenInstitute/human_variation/blob/master/scripts/abundance_variableGenes_randomForest.nb.html)  This notebook includes scripts for generating some of the figures and supplemental figures (1E, 2A-B, 2D-F, S2, S3, S4, S5, and S7), with a few points not related to figure panels.  I ran it on RStudio on a Linux box with 1TB GB of memory (but I think 256GB would be sufficient).  The random forest step took ~3 days, but the rest of the code can be run in a few hours.  
2. **R workbook for comparison with SEA-AD data set** [(LINK TO SCRIPT)](http://htmlpreview.github.io/?https://github.com/AllenInstitute/human_variation/blob/master/scripts/SEAAD_comparison.nb.html)  This notebook includes scripts for generating supplemental figures S13 and S14 that that compare data from the adult donors in this study with aged donors with and without Alzheimer’s disease from the Seattle Alzheimer’s Disease Brain Cell Atlas (SEA-AD) at sea-ad.org.  
3. **Additional scripts**  to be added soon

## License

The license for this package is available on Github at: https://github.com/AllenInstitute/human_variation/blob/master/LICENSE

## Level of Support

This code will be updated for completeness after submission, and then only if figures change during review.

## Contribution Agreement

If you contribute code to this repository through pull requests or other mechanisms, you are subject to the Allen Institute Contribution Agreement, which is available in full at: https://github.com/AllenInstitute/human_variation/blob/master/CONTRIBUTION
