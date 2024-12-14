Project: Ecological analysis of EBPR samples; 
Author: Alexander Wang

This is a sub-project from within Mphil Project

# [ Taxaheatmap]
Uses results provided by UQ's ACE facility for 16s rRNA amplicon sequencing, to process the taxonomy to the deepest level, then creates heatmaps to allow for visual comparison between different samples. 

# [ Amplicon sample analysis]
(Based on Patrick Schloss (aka Riffomonas) codeclub tutorials: https://www.youtube.com/c/RiffomonasProject)

Tests significance for individual taxa between two groups using the Wilcoxon test with Benjamini-Hochberg correction. 
Performs Non-metric Dimension scaling (NMDS) ordination with Bray curtis distance to investigate the variation in the data.
Uses avgdist function in Vegan to control for uneven sampling effort between groups using the rarefaction method. 

# [ FCS analysis]
Uses the R- flowCore package to access key parameters and analyse Flow cytometry data 
