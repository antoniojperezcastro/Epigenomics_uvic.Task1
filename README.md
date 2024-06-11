# Epigenomics_uvic.Task1
## Antonio Perez-Castro

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
