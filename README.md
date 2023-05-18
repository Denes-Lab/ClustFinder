# Introduction
ClustFinder is a clustering tool in Bash designed to automate clustering of genomes based on the genetic distances derived from single nucleotide polymorphisms (SNPs), kmers, genes, or alleles. This tool will aid researchers and public health professionals in efficient identification of epidemiological clusters. It requires no prior knowledge of Bash scripting.

## Citation

If you use ClustFinder, please cite (link to paper in JMM)

Read more about ClustFinder and how we used to cluster isolates based on Single Nucleotide Polymorphisms (SNPs): (Link to our publication)

# Description
- The ClustFinder script requires two arguments: the name of the input file and the distance threshold. 
   - The input is a file containing pairwise comparisons with genetic distances, one per line, in TSV (tab-separated values) format. 
- The script starts by separating the pairwise comparisons based on the distance according to the provided threshold. 
   - For each pairwise comparison with a distance below the threshold, the script searches for an existing cluster containing either of the two isolates.      - If there is a match, that pairwise comparison is added to the existing cluster. 
   - If there is no match, a new cluster is created. 
 - Next, the script finds any clusters that contain overlapping isolates and merges those clusters. 
 - Subsequently, the script checks to make sure that each cluster contains pairwise comparisons for each isolate pair. 
 - Missing pairwise comparisons are found in one of initially created files and added to the appropriate cluster. 
 - The final output includes a file for each cluster that contains the pairwise comparisons and distances, a file that lists isolates and the clusters they were assigned to, and a file with cluster statistics.




# Installation
To install the ClustFinder, simply clone this repository to your local machine:

``` git clone (link)
```

# Usage
To use ClustFinder, navigate to the directory containing the tool, and run the script with the desired flags:
```
cd /path/to/directory
bash clustfinder.sh -s <category> -n <min number of -s> pairwise_distance.tsv <SNP> 
```
The flags are as follows:
-s: user defined category of the isolates. It can be anything like- source, location, sampling time etc., based on your requirements
-n: denotes minimum number of genomes of a given category
If there is no need for a category, then the above command can be edited to: 
``` 
bash clustfinder.sh pairwise_distance.tsv <SNP>  
```


# Compatibility
ClustFinder is compatible with the following operating systems:
  - Ubuntu version 20.04.5 LTS (Focal Fossa)
  - Red Hat Enterprise Linux version 8.4 (Ootpa) 

Please note that while ClustFinder may work on other Linux distributions, we have only tested it on the above versions of Ubuntu and Red Hat. If you encounter any issues while using the tool on these or other operating systems, please let us know by opening an issue on GitHub.

## Utilities for command line on Mac
[Homebrew](https://brew.sh/) can be used to install the GNU versions of core utilities for Mac.


# Contributing
We welcome contributions! Please create a new branch for your feature or bugfix, then submit a pull request.


# License
This is licensed under MIT License- see the LICENSE.md file (https://github.com/Denes-Lab/GenomeCluster/blob/main/LICENSE.md) for details. 
