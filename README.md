# Introduction
Clustrator is a clustering tool in Bash. It can groups isolates into clusters based on the genomic distances derived from single nucleotide polymorphisms (SNPs), kmers, genes, or alleles. It's intended for use by public health professionals to conduct epidemiological surveillance and identify outbreaks effectively and researchers to retrospectively identify putative epidemiological clusters. It requires no prior knowledge of Bash scripting.

## Citation

If you use Clustrator, please cite (link to paper in JMM)

Read more about Clustrator and how we used to cluster isolates based on Single Nucleotide Polymorphisms (SNPs): (Link to our publication)

# Description
- Clustrator requires two arguments: name of the input file and the distance threshold (t). 
  - The input is a tab-separated file containing pairwise comparisons with genetic distances, one per line. 
- Clustrator starts by separating the pairwise comparisons based on the distance according to the provided threshold. 
- For each pairwise comparison with a distance below the threshold, the script searches for an existing cluster containing either of the two genomes. 
  - If there is a match, that pairwise comparison is added to the existing cluster. If there is no match, a new cluster is created. Clusters are named by incrementing the cluster number variable. 
- Next, clustrator finds any clusters that contain overlapping genomes and merges those clusters.
- Clustrator also checks to make sure that each cluster contains pairwise comparisons for each genome pair. 
  - Missing pairwise comparisons are found in one of initially created files and added to the appropriate cluster. 
- The final output includes a tsv of the pairwise comparisons and distances for each cluster, a tsv that lists each genome and what cluster it was assigned to, and a statistics file.



# Installation
To install the Clustrator, simply clone this repository to your local machine:

``` git clone (link)
```

# Usage
To use Clustrator, navigate to the directory containing the tool, and run the script with the desired flags:
```
cd /path/to/directory
bash clustrator.sh -s <source> -n <min number of -s> pairwise_distance.tsv <SNP> 
```
The flags are as follows:
-s: source of the isolates. Source can be anything based on your metadata file
-n: denotes minimum number of isolates/samples of a given source 
If there is no need for a source, then the above command can be edited to: 
``` 
bash clustrator.sh pairwise_distance.tsv <SNP>  
```


# Compatibility
Clustrator is compatible with the following operating systems:
  - Ubuntu version 20.04.5 LTS (Focal Fossa)
  - Red Hat Enterprise Linux version 8.4 (Ootpa) 

Please note that while clustrator may work on other Linux distributions, we have only tested it on the above versions of Ubuntu and Red Hat. If you encounter any issues while using the tool on these or other operating systems, please let us know by opening an issue on GitHub.

## Utilities for command line on Mac
[Homebrew](https://brew.sh/) can be used to install the GNU versions of core utilities for Mac.


# Contributing
We welcome contributions! Please create a new branch for your feature or bugfix, then submit a pull request.


# License
This is licensed under MIT License- see the LICENSE.md file (https://github.com/Denes-Lab/GenomeCluster/blob/main/LICENSE.md) for details. 
