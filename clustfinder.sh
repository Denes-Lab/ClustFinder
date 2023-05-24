#!/bin/bash
#Harleen K. Chaggar, Lauren K. Hudson, Ryan Kuster, Margaret E. Staton, Katie N. Garman, John R. Dunn, and Thomas G. Denes

#usage: bash script.sh -s <category> -n <min number of -s> pairwise_distance.tsv <SNP>   #user defined category -s and min number -n provided
#usage: bash script.sh pairwise_distance.tsv <SNP>   #no flags 

# Default values for flags
source_option=""
min_number=""

# Parse optional flags
while getopts "hts:n:" opt; do
  case $opt in
    h)
      # Display help message and exit
      echo "Usage: $0 [-h] [-s <source>] [-n <min_number>]"
      echo "  -h: Show help message"
      echo "  -s: Set source option"
      echo "  -n: Set minimum number"
      exit 0
      ;;
    s)
      # Set the source option
      source_option=$OPTARG
      echo "Set source_option to $source_option"
      ;;
    n)
      # Set the minimum number
      min_number=$OPTARG
      echo "Set min_number to $min_number"
      ;;
    \?)
      # Invalid option
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      # Missing argument for option
      echo "Option -$OPTARG requires an argument" >&2
      exit 1
      ;;
  esac
done

# Shift off the processed options
shift "$((OPTIND-1))"


DATE=$(date +"%F_%H-%M")
LOG=$DATE.log
THRESH=$2
THRESH2=$(( 2*$2 ))

echo "**************************************************" >> $LOG
echo "starting clustering script" >> $LOG
date +"%F %T" >> $LOG
START=$(date +"%s")
echo "--------------------------------------------------" >> $LOG
echo >> $LOG
echo "	script name/version:				        ${0##*/}
	path to script:						${0}
	current directory:					$PWD
	input pairwise distances file:		                $1
	lower distance threshold:	  	                $THRESH
	upper distance threshold:			            $THRESH2" >> $LOG
	
echo >> $LOG
echo >> $LOG
echo >> $LOG
echo >> $LOG
echo >> $LOG
echo "**************************************************" >> $LOG
echo "starting sort using a threshold of " $2 >> $LOG
date +"%F %T">> $LOG
echo "--------------------------------------------------" >> $LOG
echo >> $LOG


#add tests:
#test that input file exists
#test that input is tsv
#test that input has at least 2 lines
#test that input has numbers in the third column
#test that input has blank line at the bottom; if not, add that line
#test that threshold is a number



#reading input pairwise distance file and filters the distances as per specific SNP threshold
while IFS=$'\t' read genome1 genome2 distance; do
    if [[ $distance -le $THRESH ]] ; then
        echo -e "$genome1""\t""$genome2""\t""$distance" >> "$THRESH"dist_less.tsv
    elif [[ $distance -gt $THRESH2 ]] ; then
    	echo -e "$genome1""\t""$genome2""\t""$distance" >> "$THRESH2"dist_more.tsv
    else
        echo -e "$genome1""\t""$genome2""\t""$distance" >> "$THRESH-$THRESH2"dist.tsv
    fi
done < "$1"

echo >> $LOG
echo >> $LOG
echo >> $LOG
echo >> $LOG
echo >> $LOG
echo "**************************************************" >> $LOG
echo "making cluster files from $2"dist_less.tsv >> $LOG
date +"%F %T">> $LOG
echo "--------------------------------------------------" >> $LOG
echo >> $LOG

CLno=0  #counter for the number of clusters
while IFS=$'\t' read genome1 genome2 distance; do
     match="$(grep -l  -e $genome1 -e $genome2 CL*.tsv | head -n 1)"    #if less than, searches for existing cluster containing either genome1 or genome2
	if [[ ! -f "$match" ]] ; then   
     		let CLno=CLno+1   #If no matching file, increment CLno counter
     		CLno_ed=$(printf "%03d" $CLno)
     		echo $genome1 $genome2 "match: n; CL: " $CLno_ed " (new)">> $LOG
        	echo -e "$genome1\t$genome2\t$distance" >> CL"$CLno_ed".tsv
	else
        	echo $genome1 $genome2 "match: y; CL: " $match >> $LOG
        	echo -e "$genome1\t$genome2\t$distance" >> $match
        fi 		   
done < "$THRESH"dist_less.tsv

echo >> $LOG
echo >> $LOG
echo >> $LOG
echo >> $LOG
echo >> $LOG
echo "**************************************************" >> $LOG
echo "started fetching unique genomes from CL.tsv" >> $LOG
date +"%F %T">> $LOG
echo "--------------------------------------------------" >> $LOG
echo >> $LOG

#keeping the only unique genomes from every CL*.tsv file
for f in CL*.tsv; do 
	cut -d $'\t' -f 1 $f >> ids-"$f"
	cut -d $'\t' -f 2 $f >> ids-"$f"
	sort -u ids-"$f" >> ids-unq-"$f"  
done 

#replacing ids-unq- with nothing in ids-unq-clusters-comb_a.tsv; appending contents in ids-unq-clusters-comb_b.tsv; appending contents in ids-unq-clusters-comb.tsv 
for f in ids-unq-CL*.tsv; do 
    sed "s/^/$f\t/" "$f" >> ids-unq-clusters-comb_a.tsv  
done
 
sed "s/ids-unq-//" ids-unq-clusters-comb_a.tsv >> ids-unq-clusters-comb_b.tsv 
sed "s/.tsv//" ids-unq-clusters-comb_b.tsv >> ids-unq-clusters-comb.tsv  
rm ids-unq-clusters-comb_a.tsv
rm ids-unq-clusters-comb_b.tsv  #deleting the intermediate a.tsv and b.tsv files

#sorting, counting and finding duplicates 
sort -k2 ids-unq-clusters-comb.tsv >> ids-unq-clusters-comb-sort.tsv 
uniq -c -f1 ids-unq-clusters-comb-sort.tsv >> ids-unq-clusters-comb-sort-count.tsv

#search clusters-comb-sort-count.tsv for any >1 
#if any, search those ids in ids-unq-clusters-comb-sort.tsv, returning clusters containing duplicates
awk '{if($1>=2) print $2, $3}' ids-unq-clusters-comb-sort-count.tsv >> ids-unq-clusters-comb-sort-dups.tsv #$2 $3 is the column number 

#it stores column2 (IDs) in dups file into an array 'end', then 'next' moves to sort file and prints all duplicate values in a new file all-dups.
awk 'NR==FNR {end[$2];next} ($2 in end)' ids-unq-clusters-comb-sort-dups.tsv ids-unq-clusters-comb-sort.tsv >> ids-unq-clusters_sort-dups-in_columns.tsv

echo >> $LOG
echo >> $LOG
echo >> $LOG
echo >> $LOG
echo >> $LOG
echo "**************************************************" >> $LOG
echo "started merging cluster files based on ids-unq-clusters_in_rows.tsv" >> $LOG
date +"%F %T">> $LOG
echo "--------------------------------------------------" >> $LOG
echo >> $LOG

awk -v OFS="\t" '{a[$2]=a[$2] FS $1} END{for(i in a) print i a[i]}' ids-unq-clusters_sort-dups-in_columns.tsv >> ids-unq-clusters_in_rows.tsv

#making backup of original CL* 
mkdir backup_original-clusters/
cp CL*.tsv ids-unq-clusters_in_rows.tsv backup_original-clusters/

#deleting column1 (ID) in ids-unq-clusters_in_rows.tsv
 awk 'BEGIN {OFS="\t"}; {print $2,$3}' ids-unq-clusters_in_rows.tsv | sort -k1 >> clusters_formatted_in_rows.tsv

cp clusters_formatted_in_rows.tsv clusters_formatted_in_rows-backup.tsv


echo >> $LOG
echo >> $LOG
echo >> $LOG
echo >> $LOG
echo >> $LOG
echo "**************************************************" >> $LOG
echo "started merging cluster* files based on common string in clusters_formatted_in_rows.tsv " >> $LOG
date +"%F %T">> $LOG
echo "--------------------------------------------------" >> $LOG
echo >> $LOG

#merging clusters listed on rows based on a common string/word in clusters_formatted_in_rows.tsv and simultaneously deleting rows in clusters_formatted_in_rows.tsv after merging.
while [ "$(head -n 1 clusters_formatted_in_rows.tsv)" ]; do #expression or variable is not empty (it has a value) 
	c1=$(head -n 1 clusters_formatted_in_rows.tsv | cut -d$'\t' -f1)
	c2=$(head -n 1 clusters_formatted_in_rows.tsv | cut -d$'\t' -f2)
	echo "*****************************"  >> $LOG
	echo $c1 " is being added to " $c1"_to-merge.tsv" >> $LOG
	echo $c1 >> "$c1"_to-merge.tsv
	echo $c2 " is being added to " $c1"_to-merge.tsv" >> $LOG
    echo $c2 >> "$c1"_to-merge.tsv
    echo "removing " $c1 $c2 " from clusters_formatted_in_rows.tsv" >> $LOG
    sed -i.bak '/'$c1''$'\t'$c2'/d' clusters_formatted_in_rows.tsv
    echo "clusters_formatted_in_rows.tsv now contains:" >> $LOG 
	cat clusters_formatted_in_rows.tsv >> $LOG
    
        while read -r lineA; do
		grepno="$(grep $lineA "clusters_formatted_in_rows.tsv" | wc -l)"
		grepmatchall="$(grep $lineA "clusters_formatted_in_rows.tsv")"
		echo "searching for matches for " $lineA " in clusters_formatted_in_rows.tsv" >> $LOG
		grep $lineA clusters_formatted_in_rows.tsv > grepmatchalltmp.txt
		echo "matches for " $lineA " : " $grepno >> $LOG   
		echo "all matches for " $lineA ":" >> $LOG  
		cat grepmatchalltmp.txt  >> $LOG

			if [[ $grepno -ge 1 ]]
				then
					echo $lineA "had 1 or more matches" >> $LOG
					while IFS=$'\t' read c1b c2b; do
						echo "processing line of " $grepmatchall ": "  >> $LOG
						echo $c1b >> $LOG
						echo $c2b >> $LOG
						echo $c1b " is being added to " $c1"_to-merge.tsv" >> $LOG
						echo $c1b >> "$c1"_to-merge.tsv
						echo $c2b " is being added to " $c1"_to-merge.tsv" >> $LOG
						echo $c2b >> "$c1"_to-merge.tsv
						echo "removing " $c1b " " $c2b " from clusters_formatted_in_rows.tsv" >> $LOG
						sed -i.bak '/'$c1b''$'\t'$c2b'/d' clusters_formatted_in_rows.tsv
						echo "clusters_formatted_in_rows.tsv now contains:"  >> $LOG  
						cat clusters_formatted_in_rows.tsv  >> $LOG
					done < grepmatchalltmp.txt
				else	
					echo "no matches for " $lineA    >> $LOG
        	fi
        rm grepmatchalltmp.txt       
        done < "$c1"_to-merge.tsv
	echo "clusters_formatted_in_rows.tsv now contains:"    >> $LOG
	cat clusters_formatted_in_rows.tsv  >> $LOG
done
 
#sorting  CL*_to-merge*.tsv; keeping unique genomes
mkdir backup_pre-merged-clusters/
for f in CL*_to-merge.tsv; do 
	sort -u "$f" >> ids-unq-"$f"
	mv $f backup_pre-merged-clusters/
done

#awk command to remove extra tab character in column1 in ids-unq-CL*_to-merge.tsv
for file in ids-unq-CL*_to-merge.tsv
do                           
   awk -v OFS="\t" '{for(i=1; i<=NF; i++) a[i]=a[i] OFS $i} END{for(i=1; i<=NF; i++) print a[i]}' "$file" | awk '{sub(/^\t/,""); print}' >> clusters_to-merge_rows.tsv
   mv $file backup_pre-merged-clusters/
done

#based on ids-unq-clusters_in_rows.tsv- renaming CL* to add .tsv as the extensions to be able to recognize CL*. 
sed 's/\(CL[[:digit:]]*\)/\1.tsv/g' clusters_to-merge_rows.tsv >> clusters_to-merge_rows_ext.tsv
#here "_added-ext" is adding .tsv extension after the cluster*

echo >> $LOG
echo >> $LOG
echo >> $LOG
echo >> $LOG
echo >> $LOG
echo "**************************************************" >> $LOG
echo "started merging cluster files based on clusters_to-merge_rows_ext.tsv" >> $LOG
date +"%F %T">> $LOG
echo "--------------------------------------------------" >> $LOG
echo >> $LOG

#reading clusters_to-merge_rows_ext.tsv file line by line to merge only those CL*.tsv together that are listed on one row
while IFS=$'\t' read -r line; do
        echo "line: "$line  >> $LOG
        outfile="$(echo $line | sed 's/.tsv[[:space:]]/_/g')"
        files_to_merge="$(echo $line | sed 's/\t/ /g')"
        echo "outfile: "$outfile  >> $LOG
        echo "files_to_merge: "$files_to_merge  >> $LOG
        echo "Merging "$files_to_merge" into "$outfile  >> $LOG
        cat $files_to_merge >> "$outfile"
        for file in $files_to_merge; do
        	echo "moving to backup_pre-merged-clusters/: $file"  
        	mv $file backup_pre-merged-clusters/
        	mv ids-$file backup_pre-merged-clusters/
        	mv ids-unq-$file backup_pre-merged-clusters/
    	done
done < clusters_to-merge_rows_ext.tsv

#moving files (based on clusters_formatted_in_rows-backup.tsv clusters_formatted_in_rows.tsv clusters_to-merge_rows_ext.tsv clusters_to-merge_rows.tsv into backup_pre-merged-clusters/
mv clusters_formatted_in_rows-backup.tsv clusters_formatted_in_rows.tsv clusters_formatted_in_rows.tsv.bak clusters_to-merge_rows_ext.tsv clusters_to-merge_rows.tsv backup_pre-merged-clusters/
mv ids-unq-clusters-comb.tsv ids-unq-clusters-comb-sort.tsv ids-unq-clusters-comb-sort-count.tsv ids-unq-clusters-comb-sort-dups.tsv ids-unq-clusters_in_rows.tsv ids-unq-clusters_sort-dups-in_columns.tsv backup_pre-merged-clusters/

#for merged CL*.tsv files, keeping only unique genome IDs
for f in CL*_CL*.tsv; do 
	cut -d $'\t' -f 1 $f >> ids-"$f"
	cut -d $'\t' -f 2 $f >> ids-"$f"
	sort -u ids-"$f" >> ids-unq-"$f"
done 

# combining all ids-unq-* files into one tsv file 
for f in ids-unq-CL*.tsv; do 
    sed "s/^/$f\t/" "$f" >> ids-unq-clusters-comb_a.tsv
done

sed "s/ids-unq-//" ids-unq-clusters-comb_a.tsv >> ids-unq-clusters-comb_b.tsv
sed "s/.tsv//" ids-unq-clusters-comb_b.tsv >> ids-unq-clusters-comb.tsv
rm ids-unq-clusters-comb_a.tsv
rm ids-unq-clusters-comb_b.tsv

mkdir backup-check-dups-CL
cp CL*.tsv backup-check-dups-CL/
cd backup-check-dups-CL

#making stats file to check for number of observations in these CL*.tsv files
for file in CL*.tsv; do
    count=$(wc -l < "$file")  >> $LOG
    echo "$file: $count observations"  >> stats-before-miss-pairs.tsv  #count is shell variable that stores no. of genomes found in each .tsv file 
	echo "$file: $count observations"  >> $LOG
done
cd ../

echo >> $LOG
echo >> $LOG
echo >> $LOG
echo >> $LOG
echo >> $LOG
echo "**************************************************" >> $LOG
echo "started finding missing pairwise comparisons from files " >> $LOG
date +"%F %T">> $LOG
echo "--------------------------------------------------" >> $LOG
echo >> $LOG

#reading CL* into memory (in form of hash table), stores it in the processed_pairs array, 
for c in CL*.tsv; do
  FILES=("$THRESH""dist_less.tsv" "$THRESH2""dist_more.tsv" "$THRESH""-""$THRESH2""dist.tsv")

  awk -v clust="$c" -v LOG="$LOG" -v f1="${FILES[0]}" -v f2="${FILES[1]}" -v f3="${FILES[2]}" '
  function load_contents(files) {    #load_contents function reads contents of input files and stores in hashtable (contents_hash) for easy and faster lookup
    for (file_idx in files) {
      file = files[file_idx]
      close(file)

      while ((getline line < file) > 0) {
        split(line, f, "\t")
        pair = f[1] "\t" f[2]
        contents_hash[pair] = line
      }
      close(file)
    }
  }

  function find_matching_line(g1, g2, clust) {   #find_matching_line checks if a genome pair is present in contents_hash
    pair = g1 "\t" g2
    pair_r = g2 "\t" g1

    if (!(pair in processed_pairs) && !(pair_r in processed_pairs)) {     #If match, it appends to the cluster file and updates the processed_pairs hash table. 
      if (pair in contents_hash) {
        print contents_hash[pair] >> clust
        processed_pairs[pair] = 1
        missing_pairs_found++
        return 1
      } else if (pair_r in contents_hash) {
        print contents_hash[pair_r] >> clust
        processed_pairs[pair_r] = 1
        missing_pairs_found++        
        return 1
      }
    }

    return 0
  }

  BEGIN {    #initilizes array with filenames, loads contents into contents_hash
    delete genomes
    FILES[1] = f1
    FILES[2] = f2
    FILES[3] = f3

    load_contents(FILES)

    while ((getline genome_id < ("ids-unq-" clust)) > 0) {   #reads unique genome IDs from each cluster into genomes hashtable
      genomes[genome_id] = 1
    }
    missing_pairs_found = 0
  }
  {
    pair = $1 "\t" $2
    pair_r = $2 "\t" $1

    if (!(pair in processed_pairs) && !(pair_r in processed_pairs)) {    #processes each line of the current cluster file and checks if the genome pair is already in the processed_pairs hash table. If not, the genome pair is added to the hash table and the line is written to a temporary file.
      processed_pairs[pair] = 1
      print $0 > clust"_temp"
    }
  }
  END {    #renames the temporary file to the original cluster file name, iterates through all possible genome pairs within the cluster, and calls the find_matching_line function to search for missing pairwise comparisons.
    close(clust"_temp")

    system("mv " clust"_temp " clust)

    for (g1 in genomes) {
      for (g2 in genomes) {
        if (g1 == g2) continue

        find_matching_line(g1, g2, clust)
      }
    }
    print "Cluster:", clust >> LOG
    print "Unique genome pairs processed:", length(processed_pairs) >> LOG
    print "Missing pairwise comparisons found:", missing_pairs_found >> LOG
    print "---------------------------------------" >> LOG
  }' "$c"
done


#counting genomes in CL*.tsv files to check updates in the number of genomes after finding missing pairwise comparisons
for file in CL*.tsv; do
    count=$(wc -l < "$file")  >> $LOG
    echo "$file: $count observations"  >> stats-miss-pairs.tsv  #count is shell variable that stores no. of genomes found in each .tsv file 
	echo "$file: $count observations"  >> $LOG
done

no_clusters=$(ls CL*.tsv | wc -l)
echo "total number of clusters: " $no_clusters  >> $LOG

for ids_unq_c in ids-unq-CL*.tsv; do
    file_name="${ids_unq_c#ids-unq-}"  # removes 'ids-unq-' prefix
    no_ids=$(wc -l < "$ids_unq_c" | awk '{print $1}')  # print only the first field (line count)
    echo -e "$file_name\t$no_ids" >> ids_unq-clusters_stats.tsv
done


#metadata file should be in directory running this whole script
if [ -n "$source_option" ]; then
echo "The source option was set to: $source_option"
echo "started filtering clusters based on source" >> $LOG
		
mkdir filtered_clusters/
cp ids-unq-clusters-comb.tsv filtered_clusters/
mv sources.tsv filtered_clusters/ #can just move the sources.tsv file in the directory where this whole script will run. 
cd filtered_clusters/

genomes_and_sources="sources.tsv"   #script recognizes user defined category file as "source.tsv" but filename can be edited as per needs. This file should be tab separated and contains genome IDs in first column and metadata category in second column.
genomes_in_clusters="ids-unq-clusters-comb.tsv"
output_file="filtered_clusters_with_genomes.tsv"

# Extract genomes from source
grep "$source_option" "$genomes_and_sources" | awk '{print $1}' > filtered_genomes.tsv

# Sort the input files
sort -k 1,1 filtered_genomes.tsv > sorted_filtered_genomes.tsv
sort -k 2,2 "$genomes_in_clusters" > sorted_genomes_in_clusters.tsv

# Combine filtered genomes and cluster information
join -1 1 -2 2 -o 1.1,2.1 sorted_filtered_genomes.tsv sorted_genomes_in_clusters.tsv > filtered_genomes_in_clusters.tsv

# Count genomes per cluster and filter clusters with at least a minimum number of genomes of "provided source"
awk -v min="$min_number" 'BEGIN {OFS="\t"} {count[$2]++; genomes[$2] = genomes[$2] ? genomes[$2] OFS $1 : $1} END {for (cluster in count) if (count[cluster] >= min) print cluster ": " genomes[cluster]}' filtered_genomes_in_clusters.tsv > "$output_file"

# Cleanup intermediate files
rm filtered_genomes.tsv sorted_filtered_genomes.tsv sorted_genomes_in_clusters.tsv filtered_genomes_in_clusters.tsv

else
	echo "No source option was provided; no filtering will be done"
fi

echo >> $LOG
echo >> $LOG
echo >> $LOG
echo >> $LOG
echo >> $LOG
echo "**************************************************" >> $LOG
echo "finished" >> $LOG
date +"%F %T">> $LOG
FINISH=$(date +"%s")
DURATION=$(($FINISH-$START))
echo "Duration: $(($DURATION / 3600 )) hours $((($DURATION % 3600) / 60)) minutes $(($DURATION % 60)) seconds"
echo "Duration: $(($DURATION / 3600 )) hours $((($DURATION % 3600) / 60)) minutes $(($DURATION % 60)) seconds" >> $LOG
echo "--------------------------------------------------" >> $LOG
echo >> $LOG
echo >> $LOG
echo >> $LOG
echo >> $LOG
echo >> $LOG
echo "**************************************************" >> $LOG
