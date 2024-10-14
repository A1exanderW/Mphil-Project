# Alexander Wang

# Code was adapted from Riffomonas tutorials : CC-122, CC186

################################################################################
## Testing for significance of class ##

library(tidyverse)
library(vegan)
library(readr)
library(broom)
library(ggtext)

# Clear environment
rm(list=ls())

## Read in taxa data
data<-read_csv("raw data/Albie/level-3.csv")%>%
  pivot_longer(-index,names_to="taxonomy", values_to= "count")

## Use regex to extract class
## gets rid of Size
data_taxa<-data %>%
  mutate(taxa=str_replace_all(string=taxonomy, pattern="\\w_+", replacement=" " ))%>%
  mutate(taxa2=str_replace_all(string=taxa, pattern="_", replacement="Unassigned"))%>%
  separate(taxa2, into=c("kingdom","phylum","class","order","family","genus"),sep=";")

data_taxa<-data_taxa%>%
  mutate(phylum=if_else(phylum==" " , "Unassigned", phylum))%>%
  mutate(class=if_else(class==" " , "Unassigned", class))

## Read in metadata
metadata<-read_csv("raw data/Albie/metadata.csv")%>%
  rename_all(tolower)%>% 
  rename(index=group)

composite<-inner_join(data_taxa,metadata,by="index")%>%
  group_by(index, taxonomy)%>%
  summarize(count=sum(count), .groups="drop")%>%
  group_by(index)%>%
  mutate(rel_abund=count/sum(count))%>%
  ungroup()%>%
  select(-count) %>%
  inner_join(., metadata, by="index")

signif_genera<-composite %>%
  nest(data=-taxonomy)%>%
  mutate(test = map(.x=data, ~wilcox.test(rel_abund~before_status, data=.x) %>% tidy))
  
  ## Account for replicates
  mutate(p.adjust = p.adjust(p.value, method="BH")) %>%
  filter(p.adjust < 0.05) %>%
  select(taxonomy, p.adjust)
  
  
################################################################################
  ## Performing NMDS ##
  
  rownames(data)<- data$Group
  data <- data[,-1] #removes the first column, samples, from the table
  data2<-as.matrix(data)
  
  str(data)
  
  ## so we have can convert to matrix form, where our rows are each sample and the columns are the different OTUS, 
  ##  and data in cells are the counts
  
  dist<- vegdist(data2, method="bray")
  
  nmds <-metaMDS(dist)
  
  ## Convert to tibble or dataframe to plot
  scores(nmds) %>%
    as_tibble(rownames)
  
  ## Perform NMDS
  
  nmds<- metaMDS(data2)
  
  ## Rarefy with avgdist
  set.seed(112024)
  dist <- avgdist(data2, dmethod="bray", iterations=100, sample=1800)
  # where 100 is the default iterations
  
  ggplot(dist)
  
################################################################################
  
view<-signif_genera$data[[1]]

str(view$rel_abund)

wilcox.test(view$rel_abund~view$before_status, data=view)
  

rownames(data)<- data$Group
data <- data[,-1] #removes the first column, samples, from the table
data2<-as.matrix(data)

str(data)

## so we have can convert to matrix form, where our rows are each sample and the columns are the different OTUS, 
##  and data in cells are the counts

dist<- vegdist(data2, method="bray")

nmds <-metaMDS(dist)

## Convert to tibble or dataframe to plot
scores(nmds) %>%
  as_tibble(rownames)

## Perform NMDS
nmds<- metaMDS(data2)

## Rarefy with avgdist
set.seed(112024)
dist <- avgdist(data2, dmethod="bray", iterations=100, sample=1800)
# where 100 is the default iterations

ggplot(dist)
####################################################################################
