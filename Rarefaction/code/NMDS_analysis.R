library(tidyverse)
library(vegan)
library(readr)
library(broom)
library(ggtext)

# Clear environment
rm(list=ls())

## Read in taxa data
data<-read_csv("Ace results/Albie/J5863/P323/qiime_barplot/level-3.csv")%>%
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
metadata<-read_csv("Ace results/Albie/J5863/P323/qiime_barplot/metadata.csv")%>%
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
  
  ## Account for
  mutate(p.adjust = p.adjust(p.value, method="BH")) %>%
  filter(p.adjust < 0.05) %>%
  select(taxonomy, p.adjust)
  
  
############################################################################################
  
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
  
  ###############################################################################################3
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
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

## Extract the taxa from level 3 data ##
data<- read_csv("Ace results/Albie/J5863/P323/qiime_barplot/level-3 _index.csv")
data_matrix <- as_tibble(data)
samples <-data_matrix$index

# convert columns into rows (pivot_wider) and into tibble structure
data_2<- cbind(samples, data_matrix)%>%
  select(-index)

data_initial <- data_2 %>%
  pivot_longer(cols = -samples, names_to = "taxa", values_to = "count")
# data_initial<-as_tibble(pivot_longer(data_2,names_to= "taxa", values_to="count"))

data_temp<-data_initial%>% 
  select(-samples,taxa,count)%>%
  mutate(taxa2=str_replace_all(string=taxa, pattern="\\w_+", replacement=" " ))%>%
  mutate(taxa3=str_replace_all(string=taxa2, pattern="_", replacement="Unassigned"))%>%
  separate(taxa3, into=c("kingdom","phylum","class","order","family","genus"),sep=";")%>%
  select(class,count)


##############################################################  

data_init <- data_2 %>%
  pivot_longer(cols = -samples, names_to = "taxa", values_to = "count")
# data_initial<-as_tibble(pivot_longer(data_2,names_to= "taxa", values_to="count"))

data_temp<-data_init%>% 
  select(samples,taxa,count)%>%
  mutate(taxa2=str_replace_all(string=taxa, pattern="\\w_+", replacement=" " ))%>%
  mutate(taxa3=str_replace_all(string=taxa2, pattern="_", replacement="Unassigned"))%>%
  separate(taxa3, into=c("kingdom","phylum","class","order","family","genus"),sep=";")%>%
  ## Need to find away to change first row to unassigned each column

data_temp<- data_temp%>%
  mutate(class=str_trim(class, side="both"))
  
## Testing for significance of microbial taxa

data_OTU<-read_csv("Ace results/Albie/J5863/P323/qiime_barplot/level-7_OTU.csv")%>%
  select(-Group, starts_with("Otu"))%>%
  as.data.frame()

data_OTU_2 <- data_OTU %>%
  pivot_longer(cols = Group, names_to = "OTU", values_to = "count")%>%
  select(-Group)

data_OTU_mat <- as_tibble(data_OTU_2)%>%
  samples <-data_OTU$samples

# convert columns into rows (pivot_wider) and into tibble structure
data_2<- cbind(samples, data_OTU_mat)%>%
  select(-index)

new_table<-cbind(data_temp,data_OTU_2)

## Testing for signif (cheat method, prepped excel sheet )
data_OTU_cur<-read_csv("Ace results/Albie/J5863/P323/qiime_barplot/level-3_curated.csv")%>%
  as.data.frame(Group)%>%
  mutate(., rel_abund=Count/sum(Count))

sig_class<-data_OTU_cur%>%
  nest(data=-Taxa)%>%
  wilcox.test(., paired=TRUE)


## map lets us iterate our function over each row 
## ~ tells r we are referrign to a function
## rel abund explained by srn
## data=x, whre data is a variable for the wilcox test 

mutate(test=map(.x=data, ~wilcox.test(rel_abund~srn, data=.x) %>% tidy)) %>%
  
# Wrangle data, and join,
# Sort by sample names
# add a relative abundance
# perform wilcoxin test 





