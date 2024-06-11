## Antonio Perez-Castro

# Epigenomics Task 1

The present exercise is built over steps 3.1 and 3.2 of the class tutorial. To run the tutorial 3.1, make sure the docker desktop app is running before running the container in the WSL2.

```bash

git clone https://github.com/bborsari/epigenomics_uvic
cd epigenomics_uvic
sudo docker run -v $PWD:$PWD -w $PWD --rm -it dgarrimar/epigenomics_course

```

## 4. EN‐TEx ATAC‐seq data: downstream analyses

#### 4.1 Move to folder ATAC-seq, and create folders to store bigBed data files and peaks analyses files. Make sure the files are organized in a consistent way as done for ChIP-seq.

```bash

cd ATAC-seq
mkdir data anlyses annotation
mkdir data/bigBed.files data/bigWig.files data/bed.files
mkdir analyses/peaks.analysis

```
#### 4.2 Retrieve from a newly generated metadata file ATAC-seq peaks (bigBed narrow, pseudoreplicated peaks, assembly GRCh38) for stomach and sigmoid_colon for the same donor used in the previous sections. Hint: have a look at what we did here. Make sure your md5sum values coincide with the ones provided by ENCODE.

In encode portal select DNA accessibility and tissues sigmoidal_colon and stomach for sample ENCDO451RUA. Download the text file that contains the URL for the metadata. The metadata will be downloaded using that URL inserted in the code below:

```bash
../bin/download.metadata.sh "https://www.encodeproject.org/metadata/?replicates.library.biosample.donor.uuid=d370683e-81e7-473f-8475-7716d027849b&status=released&status=submitted&status=in+progress&biosample_ontology.term_name=sigmoid+colon&biosample_ontology.term_name=stomach&assay_slims=DNA+accessibility&type=Experiment"
```

The next step consist on parsing the metadata file to search the lines that contain the string "ATAC-seq" with a grep command.
The second grep command looks for lines containing "bigBed_narrowPeak" while the third and fourth look for "pesudoreplicated_peaks" and "GRCh38".
The relevant inforamtion is stored in the file bigBed.peaks.ids.txt. This information is used in a loop to retrive the selected bigBed files from encode and store them in data/bigBed.files

```bash
grep -F ATAC-seq metadata.tsv |\
grep -F "bigBed_narrowPeak" |\
grep -F "pseudoreplicated_peaks" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $11, $23}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u > analyses/bigBed.peaks.ids.txt

cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do
  wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done
```

Similarly, we can retrieve the bigwig files for the same samples

```bash
grep -F ATAC-seq metadata.tsv |\
grep -F "bigWig" |\
grep -F "fold_change_over_control" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $11, $23}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u > analyses/bigWig.FC.ids.txt

cut -f1 analyses/bigWig.FC.ids.txt |\
while read filename; do
  wget -P data/bigWig.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigWig"
done
```

It's a good practice to check the quality of the files just downloaded. We can verify their MD5 hash with the next command.

```bash
file_type="bigBed"; ../bin/selectRows.sh <(cut -f1 analyses/"$file_type".*.ids.txt) metadata.tsv | cut -f1,46 > data/"$file_type".files/md5sum.txt; cat data/"$file_type".files/md5sum.txt | while read filename original_md5sum; do md5sum data/"$file_type".files/"$filename"."$file_type" | awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}'; done > tmp; mv tmp data/"$file_type".files/md5sum.txt; awk '$2 != $3' data/"$file_type".files/md5sum.txt
```
No output was obtained, we can proceed.


#### 4.3 For each tissue, run an intersection analysis using BEDTools: report 1) the number of peaks that intersect promoter regions, 2) the number of peaks that fall outside gene coordinates (whole gene body, not just the promoter regions). Hint: have a look at what we did here and here.

The bedtools intersect requires .bed files. The first next step will transform the bigBed files into .bed files and store them into data/bed.files directory.

```bash
cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do
  bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
done
```
we also need reference .bed files with the non redundant TSS and the protein coding gene body from step 3.2 in the tutorial, we can bring them to the ATAC-seq branch of the directory tree.

```bash
cp ../ChIP-seq/annotation/gencode.v24.protein.coding.non.redundant.TSS.bed annotation
cp ../ChIP-seq/annotation/gencode.v24.protein.coding.gene.body.bed annotation
```
The next step will intercept the start and end regions of TSS with the ATAC-seq peaks. The result of the next code will be a .txt file with the name of all genes that contain an ATAC-seq peak that intercepts a promoter region. After we can count how many we have identified.

```bash
cut -f-2 analyses/bigBed.peaks.ids.txt | \
while read filename tissue; do 
  bedtools intersect -a annotation/gencode.v24.protein.coding.non.redundant.TSS.bed -b data/bed.files/"$filename".bed -u | \
  cut -f7 | \
  sort -u > analyses/peaks.analysis/genes.with.peaks."$tissue".ATAC-seq.bed
done
```

The next commands will count the number of genes that showed an ATAC-seq peak overlapping with the TSS in sigmoid_colon. The result is 14830.

```bash
wc -l analyses/peaks.analysis/genes.with.peaks.sigmoid_colon.ATAC-seq.bed
```

Similarly, we can count how many in stomach. The result is 15029

```bash
wc -l analyses/peaks.analysis/genes.with.peaks.stomach.ATAC-seq.bed
```

In order to count how many peaks fall outside gene bodies we can use the bedtool intersect function with the ATAC-peaks in the term -a and the reference gencode.v24.protein.coding.gene.body.bed in term -b, using option -v. The resulting .bed files are saved in the directory analyses/peaks.analysis to be counted

```bash
cut -f-2 analyses/bigBed.peaks.ids.txt |\
while read filename tissue; do 
  bedtools intersect -a data/bed.files/"$filename".bed -b annotation/gencode.v24.protein.coding.gene.body.bed -v |\
    sort -u > analyses/peaks.analysis/peaks.outside.genes."$tissue".ATAC-seq.bed	
done
```

Now we are ready to count how many peaks per tissue
There are 34537 peaks in stomach

```bash
wc -l analyses/peaks.analysis/peaks.outside.genes.stomach.ATAC-seq.bed
```

There are 37035 peaks in sigmoid_colon

```bash
wc -l analyses/peaks.analysis/peaks.outside.genes.sigmoid_colon.ATAC-seq.txt
```

## 5. Distal regulatory activity.
From the previous section, we can count on the information about ATAC-seq peaks that fall outside any gene. They may contain distal regulatory elements.

#### 5.1. Create a folder regulatory_elements inside epigenomics_uvic. This will be the folder where you store all your subsequent results

```bash
cd ..
mkdir regulatory_elements
cd regulatory_elements
mkdir analyses data
mkdir data/bed.files
mkdir data/bed.files/k27 data/bed.files/k4me

```

#### 5.2.  Distal regulatory regions are usually found to be flanked by both H3K27ac and H3K4me1. From your starting catalogue of open regions in each tissue, select those that overlap peaks of H3K27ac AND H3K4me1 in the corresponding tissue. You will get a list of candidate distal regulatory elements for each tissue. How many are they?

In order to complete this task we need to download the files corresponding to the H3K27ac data and H3K4me1 for the same sample and tissues. The procedure we followed was identical and step 3.1 in the tutorial except for specifying the right histone mark.
Once we obtained the text file, the first line contains the URL to download the metadata file for the selected samples.
We used the next command to create the metadata file for the selected samples.
```bash
../bin/download.metadata.sh "https://www.encodeproject.org/metadata/?replicates.library.biosample.donor.uuid=d370683e-81e7-473f-8475-7716d027849b&status=released&status=submitted&status=in+progress&assay_title=Histone+ChIP-seq&target.label=H3K27ac&target.label=H3K4me1&biosample_ontology.term_name=sigmoid+colon&biosample_ontology.term_name=stomach&type=Experiment"
```
we have to parse the metadata file just created in search for the filenames we need, looking for H3k27ac, bigBed_narrowPeak, pseudoreplicated_peaks and GRCh38 key words. The filenames identified are used by the command wget to retrieve the bigbed files needed.

```bash
grep -F H3K27ac metadata.tsv |\
grep -F "bigBed_narrowPeak" |\
grep -F "pseudoreplicated_peaks" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $11, $23}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u > analyses/bigBed.peaks.ids.k27.txt

cut -f1 analyses/bigBed.peaks.ids.k27.txt |\
while read filename; do
  wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done
```

And similarly, we can retrieve the files for H3K4me1 histone mark.

```bash
grep -F H3K4me1 metadata.tsv |\
grep -F "bigBed_narrowPeak" |\
grep -F "pseudoreplicated_peaks" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $11, $23}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u > analyses/bigBed.peaks.ids.k4me1.txt

cut -f1 analyses/bigBed.peaks.ids.k4me1.txt |\
while read filename; do
  wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done
```

Before running bedtool intercept we have to transform the bigBed files into .bed files. The resulting files were stored in different files for k27 or k4me inside the directory data/bed.files

```bash
cut -f1 analyses/bigBed.peaks.ids.k27.txt |\
while read filename; do
  bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/k27/"$filename".bed
done

cut -f1 analyses/bigBed.peaks.ids.k4me1.txt |\
while read filename; do
  bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/k4me/"$filename".bed
done
```

The next step looks for the ATAC-seq peaks that falls outside gene bodies in stomach that overlap with H3k27ac and with H3k4me1 at the same time. For that, the result of the first bedtool intercept is piped into the second bedtool intercept, and the result is counted or saved as data/bed.files/stomach.peaks.all3.bed
8022 such peaks were found in stomach. Similarly, the commands for sigmoid colon are presented, fiding 14215 of such regions.

```bash
bedtools intersect -a ../ATAC-seq/analyses/peaks.analysis/peaks.outside.genes.stomach.ATAC-seq.bed -b data/bed.files/k27/ENCFF977LBD.bed -u | bedtools intersect -a - -b data/bed.files/k4me/ENCFF844XRN.bed -u | wc -l
bedtools intersect -a ../ATAC-seq/analyses/peaks.analysis/peaks.outside.genes.stomach.ATAC-seq.bed -b data/bed.files/k27/ENCFF977LBD.bed -u | bedtools intersect -a - -b data/bed.files/k4me/ENCFF844XRN.bed -u > data/bed.files/stomach.peaks.all3.bed
bedtools intersect -a ../ATAC-seq/analyses/peaks.analysis/peaks.outside.genes.sigmoid_colon.ATAC-seq.bed -b data/bed.files/k27/ENCFF872UHN.bed -u | bedtools intersect -a - -b data/bed.files/k4me/ENCFF724ZOF.bed -u | wc -l
bedtools intersect -a ../ATAC-seq/analyses/peaks.analysis/peaks.outside.genes.sigmoid_colon.ATAC-seq.bed -b data/bed.files/k27/ENCFF872UHN.bed -u | bedtools intersect -a - -b data/bed.files/k4me/ENCFF724ZOF.bed -u > data/bed.files/sigmoid_colon.peaks.all3.bed
```

#### 5.3  Focus on regulatory elements that are located on chromosome 1 (hint: to parse a file based on the value of a specific column, have a look at what we did here), and generate a file regulatory.elements.starts.tsv that contains the name of the regulatory region (i.e. the name of the original ATAC-seq peak) and the start (5') coordinate of the region.

To solve this task we are filtering the files tissue.peaks.all.bed for lines containing chr1 and we are printing the columns 4 and 2 separated by tab into a .tsv file per tissue.

```bash
grep "chr1" data/bed.files/stomach.peaks.all3.bed | awk '$1 == "chr1" {print $4"\t"$2}' > analysis/regulatory.elements.starts.stomach.tsv
grep "chr1" data/bed.files/sigmoid_colon.peaks.all3.bed | awk '$1 == "chr1" {print $4"\t"$2}' > analysis/regulatory.elements.starts.sigmoid_colon.tsv
```

#### 5.4  Focus on protein-coding genes located on chromosome 1. From the BED file of gene body coordinates that you generated here, prepare a tab-separated file called gene.starts.tsv which will store the name of the gene in the first column, and the start coordinate of the gene on the second column (REMEMBER: for genes located on the minus strand, the start coordinate will be at the 3'). Use the command below as a starting point:

This task is done in a single line of command with grep for chr1 and pipping the result to an awk command written with a conditional that will check the strand of the gene before taking the accurate starting point. The result is stored in data/gene.starts.tsv

```bash
grep "chr1"  ../ChIP-seq/annotation/gencode.v24.protein.coding.gene.body.bed | awk 'BEGIN{FS=OFS="\t"}{if ($6=="+"){start=$2} else {start=$3}; print $4, start}' > data/gene.starts.tsv
```
#### 5.5  Task 5: Download or copy this python script inside the epigenomics_uvic/bin folder. Have a look at the help page of this script to understand how it works:

The link offered downloaded a .py file into the host machine. To transfer it inside the docker container we are using I opened another WSL2 terminal to run the next codes. The first one makes a copy of the .py file into the root, the second uploads it into the bin directory into the container.
The name of the container will be different every time we run the container. 

```bash
cp /mnt/c/Users/anton/Downloads/get.distance.py ~/
docker cp ~/get.distance.py 223255cb0803:/home/antoniojperez/epigenomics_uvic/bin
```

Inside the WSL2 terminal that is running the container we can check the presence of the file get.distance.py into the bin directory, and we can navigate to that directory to be able to modify the file with nano.

```bash
ls ../bin
nano ../bin/get.distance.py
```

The file was modified to the next script:

#!/usr/bin/env python


#************
# LIBRARIES *
#************

import sys
from optparse import OptionParser


#*****************
# OPTION PARSING *
#*****************

parser = OptionParser()
parser.add_option("-i", "--input", dest="input")
parser.add_option("-s", "--start", dest="start")
options, args = parser.parse_args()

open_input = open(options.input)
enhancer_start = int(options.start)


#********
# BEGIN *
#********

x=1000000 # set maximum distance to 1 Mb
selectedGene="" # initialize the gene as empty
selectedGeneStart=0 # initialize the start coordinate of the gene as empty

for line in open_input.readlines(): # for each line in the input file
	gene, y = line.strip().split('\t') # split the line into two columns based on a tab 
	position = int(y) # define a variable called position that correspond to the integer of the start of the gene
	diff = abs(position – enhancer_start) # compute the absolute value of the difference between position and enhancer_start

	# if diff < x:
		x = diff  # this value will now be your current x
		selectedGene = gene # save gene as selectedGene
		selectedGeneStart = position # save position as selectedGeneStart

print "\t".join([selectedGene, str(selectedGeneStart), str(x)])

At this point the exercice offers a test to check the functioning of the .py script. I can't obtain the same result and I can't find where the problem may be. Since the file gene.starts.tsv is generated from provided reference files I'm assuming the reference files may have changed since this tutorial was prepared years ago.

#### 5.6  For each regulatory element contained in the file regulatory.elements.starts.tsv, retrieve the closest gene and the distance to the closest gene using the python script you created above. Use the command below as a starting point:
We use the script get.distance.py recursively, element by element, looking for the regulatory region identified, that means, the region that is an ATAC-seq peak flanked by H3K27ac and H3K4me1, that seats closer to the starting point of each gene in chr1. The results are storaged in .tsv files and contain just two columns, the gene ID and the distance to the closest regulatory element.

```bash
cat analysis/regulatory.elements.starts.stomach.tsv | while read element start; do 
   python ../bin/get.distance.py --input data/gene.starts.tsv --start $start 
done  | awk -F'\t' '{print $1"\t"$3}' > analyses/regulatoryElements.genes.distances.stomach.tsv

cat analysis/regulatory.elements.starts.sigmoid_colon.tsv | while read element start; do 
   python ../bin/get.distance.py --input data/gene.starts.tsv --start $start 
done  | awk -F'\t' '{print $1"\t"$3}' > analyses/regulatoryElements.genes.distances.sigmoid_colon.tsv
```

#### 5.7  Use R to compute the mean and the median of the distances stored in regulatoryElements.genes.distances.tsv.

```bash
cd ../bin
nano
```

we wrote the .R script:

#!/usr/bin/env Rscript

# Load necessary library
library(optparse)

# Define the command line arguments
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Path to the input file", metavar="character")
)

# Parse the command line arguments
parser = OptionParser(option_list=option_list)
arguments = parse_args(parser)

# Read the data
data <- read.table(arguments$options$input, sep="\t", header=FALSE)

# Extract the second column
column2 <- data$V2

# Calculate the mean and median
mean_val <- mean(column2, na.rm = TRUE)
median_val <- median(column2, na.rm = TRUE)

# Print the results
print(paste("Mean: ", mean_val))
print(paste("Median: ", median_val))

The .R script was storaged with the name calculate_stats.R. Now we can apply it to measure mean and median distances between regulatory elements and the start point of the genes.

```bash
cd analyses
Rscript ../../bin/calculate_stats.R --input regulatoryElements.genes.distances.stomach.tsv
Rscript ../../bin/calculate_stats.R --input regulatoryElements.genes.distances.sigmoid_colon.tsv
```
The results for stomach were:
Mean:  27550.5
Median:  10226

The results for the sigmoid colon were:
Mean:  40714.49
Median:  12714

Looks like the median value is similar in both tissues while some long distance regulatory elements are pulling the mean higher in sigmoid colon.
