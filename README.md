###### Description
Clustrator is a clustering tool in Bash. It can groups isolates into clusters based on the genomic distances derived from single nucleotide polymorphisms (SNPs), kmers, genes, or alleles. It's intended for use by public health professionals to conduct epidemiological surveillance and identify outbreaks effectively and researchers to retrospectively identify putative epidemiological clusters. It requires no prior knowledge of Bash scripting.

Read more about Clustrator and how we used to cluster isolates based on Single Nucleotide Polymorphisms (SNP) here :
(Link our publication)

###### Installation
To install the Clustrator, simply clone this repository to your local machine:

git clone (link)

##### Utilities for command line in Linux and Mac



###### Usage
To use Clustrator, navigate to the directory containing the tool, and run the script with the desired flags:

cd /path to script
bash clustrator.sh -s <source> -n <min number of -s> pairwise_distance.tsv <SNP>
 
The flags are as follows:
-s: source of the isolates. It can be anything based on the metadata file
-n: denotes minimum number of isolates of a particular source 

  
###### Acknowlegdements
Clustrator would not have been without

###### Contributing
We welcome contributions! Please create a new branch for your feature or bugfix, then submit a pull request.


###### License
This is licensed under MIT License- see the LICENSE.md file (https://github.com/Denes-Lab/GenomeCluster/blob/main/LICENSE.md) for details. 
