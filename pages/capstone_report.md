---
layout: page
title: Capstone Project Report
permalink: "/capstone_report/"
image: 
---

# Onset of Type II Diabetes
by Marina Witherell and Vini Karumuru

<a id="table-of-contents"></a>
## 0 - Table of Contents
0 [Table of Contents](#table-of-contents)
1. [Introduction](#introduction)
2. [Hypothesis](#hypothesis)
    - 2.1 [Hypothesis 1](#hypothesis-1)
    - 2.2 [Hypothesis 2](#hypothesis-2)
3. [Methods](#methods)
    - 3.0 [Common Methods](#common-methods)
    - 3.1 [Methods for Hypothesis 1](#methods-1)
    - 3.2 [Methods for Hypothesis 2](#methods-2)
    - 3.3 [Methods for Further Exploration of Top Taxa](#methods-3)
4. [Results](#results)
    - 4.1. [Results for Hypothesis 1](#results-1)
    - 4.2. [Results for Hypothesis 2](#results-2)
    - 4.3. [Results for Further Exploration of Top Taxa](#results-3)
5. [Discussion](#discussion)
    - 5.1 [Further Directions](#further-directions)
6. [References](#references)

<a id="introduction"></a>
## 1 - Introduction
Type 2 diabetes (T2D) is a significant health issue in the US. In 2022, the CDC reported 10.7% of adults had type 2 diabetes. There are also an estimated 8.5 million undiagnosed cases, and 79 million cases of pre-diabetes with 50% of these cases converting to a lifetime of diabetes. There has been notice of differences of gut microbiome between diabetic and non-diabetic people, which has been observed to have an effect on glucose levels (Stanford Medicine). This project hopes to better understand the changes of the microbiome during diabetic onset and progression. The microbiomes of the nasal samples and those of the fecal samples and the relationship between the two sample types are of prime interest. \
The dataset used in this project was taken from the NIH’s Integrated Human Genome Project. Due to limited computational resources, a randomly sampled subset of 1000 samples was taken from the original dataset for our analysis. Our dataset contains samples from 77 pre-diabetic patients. A patient is considered to be pre-diabetic when their fasting blood sugar level is higher than normal but is not high enough to be considered diabetic (this would mean a fasting blood sugar level of 100-125 mg/dL) (Maruthur, 2023). Samples were collected from the patients at two body sites: the feces and the nasal cavity. The dataset contained 465 fecal samples and 535 nasal samples. Most patients had samples taken during more than one timepoint, hence why there are many more samples than there are patients. The microbial DNA was extracted from these samples and underwent 16S amplicon sequencing (Stanford Medicine). As our dataset only contained data from pre-diabetic patients, with no set of control, healthy patients, the focus of our analysis was to compare the fecal and nasal microbiomes of these pre-diabetic patients.

***

<a id="hypothesis"></a>
## 2 - Hypothesis
<a id="hypothesis-1"></a>
### 2.1 Hypothesis 1
The fecal microbiome is more diverse than the nasal microbiome in pre-diabetic patients. This is because the fecal microbiome comes in more direct contact with the internal metabolic processes of the body and therefore is more frequently altered.

<a id="hypothesis-2"></a>
### 2.2 Hypothesis 2
The microbial profile of the skin and fecal samples will be significantly different when using the beta diversity measure. If it is found significant, I will focus on the most abundant taxa and identify the group that explains the differences.

***

<a id="methods"></a>
## 3 - Methods
<a id="common-methods"></a>
### 3.0 Common Methods
We were provided with a phyloseq object created by the Snyder Lab under the iPOP project at Stanford Medicine (Zhou et al, 2019). According to their research paper, the phyloseq object was built using the QIIME pipeline, along with custom scripts in Python and R (Caporaso et al, 2010).The phyloseq object contains sample metadata, the number of reads per sample for each OTU (operational taxonomic unit), and taxonomic information for each OTU.
We began our exploration of the dataset by looking at the number of reads per sample. We summed the number of reads per column in the OTU table and performed calculations with conditional statements on these sums. We found that approximately 3% of our samples (31) contained less than 5,000 reads and the remaining samples (969) contained more than 5,000 reads. The largest sample contained 234,415 and the smallest sample contained 2 reads. The average number of reads per sample was 24,934.
The otu_table in the phyloseq object was normalized using the relative abundance method with the transform_sample_counts() function. The dataset was then split into two separate phyloseq objects, one for each body site (feces and nasal cavity). To accomplish this, we used the subset_samples() function, as the body site information was contained in the sample_data table within the phyloseq object. Most of our methods involved us performing computations and manipulations on these two phyloseq objects separately (but simultaneously), then comparing the two body sites.

<a id="methods-1"></a>
### 3.1 Methods for Hypothesis 1 (Vini)
I tested my hypothesis by employing two measures of alpha diversity: species richness and Shannon index. Species richness was calculated per sample (column in the otu_table) by summing the number of OTUs that had an abundance that was non-zero. Then, as a summary statistic, I determined the median species richness.
I then went on to compute the Shannon index of each sample. The Shannon index (H) is an alpha diversity measure that takes into account both species richness and species evenness. A higher Shannon index would indicate a greater amount of diversity. It can be computed using the equation H=-i=1Spilnpi, where pi is the proportion of the sample’s reads that are found to be the ith OTU, ln is the natural logarithm, and s is the total number of OTUs in the community. The phyloseq package already has an inbuilt function that can calculate various measures of diversity, so I did not have to manually code this calculation. Instead, I utilized the estimate_richness() function with the argument measures = “Shannon”. This gave me a dataframe containing a list of Shannon indexes, one corresponding to each sample. I was then able to plot these indexes at the two body sites as side-by-side boxplots with the ggplot package. I conducted a two-sided Wilcoxon rank-sum test between the two body sites to test the null hypothesis that there is no difference between the median Shannon index of the fecal and nasal microbial communities.

<a id="methods-2"></a>
### 3.2 Methods for Hypothesis 2 (Marina)
Using the normalized data, I conducted a permutational multivariate analysis of variation (PERMANOVA) test on the data. This method compares the variation between the two sample sites to the variation within the two sample sites. I calculated this using the Bray-Curtis distance method, which calculates the similarity of each instance in the phyloseq object. The Bray-Curtis dissimilarity is a measure of beta diversity, showing how different the nasal and fecal microbiomes are from each other. I then went on to identify the most abundant taxa in each body site. 
I decided that the best way to do this would be to find the top ten OTUs from the nasal cavity and the top ten from the feces, using the two separate phyloseq objects of data from each body site (in section 3.0). Note that the following steps were done for each of the phyloseq objects separately until stated otherwise. First, I found the row sums of the OTU table and organized the table based on those sums in an increasing to decreasing fashion. Then, I made a function to filter out OTUs that were not present in at least three samples. This was implemented so that an OTU that was highly abundant in only one sample would not make it into the top ten (since it is supposed to represent the most abundant across the body site). From these remaining OTUs, the top ten were subset into a new phyloseq object. I created another filter function that returns the samples that have at least one of the top 10 OTUs present. This was required for the ordination plot. I then updated the phyloseq object to only include those samples, which left me with two filtered phyloseq objects - nasal and feces. I combined these two objects into one and created an ordination plot of the taxa (based on order) and the samples. I also wanted to inspect the top 10% of the OTUs and how that compared to just the top ten from each body site, so I repeated all of the steps above using the top 60 OTUs. This left me with four ordination plots. I also performed another PERMANOVA test on these to find if the top results were significantly different, or if the most abundant from each site were similar.

<a id="methods-3"></a>
### 3.3 Methods for Further Exploration of Top Taxa (Vini)
To further expand our analysis, I decided that it might be interesting to dive deeper into the top 20 OTUs we had identified to be the most abundant in each body site. More specifically, I wanted to see if there were any trends in the abundances of these top taxa among many patients. 
I began by removing any patients that did not have samples taken from both body sites. 65 patients remained after doing so. Then, because the dataset contained multiple samples per patient, I averaged each patient’s reads at each body site. For each of the top 10 fecal OTUs, I conducted a one-sided paired Wilcoxon rank-sum test where the mean abundance of this OTU in each patient’s feces was paired with the mean abundance of this OTU at their nasal cavity. I was testing the null hypothesis that the median abundance of the OTU in the feces is either less than or no different from the median abundance of the OTU in the nasal cavity. I repeated this process for the top 10 nasal OTUs (but flipping the body sites in the null hypothesis, to test whether the abundance of these nasal OTUs was significantly greater in the nasal cavity compared to the feces). For any of these one-sided rank-sum tests that resulted in a non-significant p-value (p > 0.05), I conducted a paired two-sided rank-sum test to see if there was any difference between the body sites. The ggpaired() function within the ggpubr package was utilized to plot paired boxplots for each of the top 20 taxa. For the purpose of better visualization, I log-scaled the mean abundance so that more of the data points would be discernible. I used the plot_list() function within the cowplot package to combine multiple of the individual OTU plots into a single figure so that these plots could be more easily compared to one another.

***

<a id="results"></a>
## 4 - Results
### 4.1 Results for Hypothesis 1
<figure>
    <img src="{{site.baseurl}}/assets/images/fig_4.1.1.png">
    <figcaption>
    <strong>Figure 4.1.1 Table of Median Alpha Diversity Measures by Body Site.</strong> Each column is a body site.The first row contains median species richness (the number of OTUs detected in a sample) and second row contains the median Shannon Index (which takes into account both species richness and species evenness).</figcaption>
</figure>
<p>
<img src="{{site.baseurl}}/assets/images/fig_4.1.2.png">
<em>
<strong>Figure 4.1.2 Side-by-Side Boxplots of Shannon Indexes At Each Body Site.</strong> Each datapoint represents a sample, with the left boxplot displaying samples from the feces and the right boxplot displaying samples from the nasal cavity. A two-sided, unpaired Wilcoxon rank-sum test yielded a p-value many orders of magnitude lower than 0.05.</em>
</p>
The fecal samples had higher alpha diversity values overall. This is supported by the fact that the median species richness in the fecal samples is more than twice the species richness in the nasal samples (Figure 4.1.1), meaning each sample tended to have a greater number of OTUs detected. Additionally, the fecal samples had a set of significantly higher Shannon indexes when compared to the nasal samples (p < 2.2e-16), which is indicative of greater diversity. 