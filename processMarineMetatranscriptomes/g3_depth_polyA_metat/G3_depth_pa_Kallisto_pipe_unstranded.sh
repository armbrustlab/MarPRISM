## Adapted from Zinka  read_counts_G3PM_pipe.sh
## step 3 only: Kallisto counts - UNSTRANDED! index already done.


# 1. For G3PM assemblies: combine and cluster all contigs in aa space
#    from the bf100 group with linlcust. 99% identity clusters.
# 2. Use these clustered aa reads to select for the same ones in nt space. Get
#    a collection of nt clustered contigs for the whole cruise.
# 3. Map raw pair-end reads for each sample onto the clustered contigs to get
#    raw counts… use Kallisto.

## Build the kallisto index from the clustered assemblies:

# This command will take the input fasta NPac.${PROJECT}.bf100.id99.nt.fasta
# and create the kallisto index NPac.${PROJECT}.bf100.id99.nt.idx

cd /scratch/Sacha/Kallisto/G3depth

#mkdir kallisto
screen -S depth

#set limit for shell script at ~800G 
ulimit -v 830000000
kallisto index -i NPac.G3PA_depth.id99.nt.idx /mnt/nfs/projects/armbrust-metat/gradients3/g3_depth_pa_metat/assemblies/clustered/nucleotideClusterSequences.fasta

############################# 3. Get counts with Kalisto ############################

##### Running kallisto in scratch directory:

#make a sample list - note, I did not change the sample names provided by Elaina
ls /mnt/nfs/projects/armbrust-metat/gradients3/g3_depth_pa_metat/merged/*.fastq.gz | sed 's|/mnt/nfs/projects/armbrust-metat/gradients3/g3_depth_pa_metat/merged/\(.*\).flash.extendedFrags.fastq.gz|\1|' > sample_names.txt

# create the kallisto alignment function:
function run_kallisto {
# create output directory for this sample
mkdir ${SAMPLE}

# declare paths and filenames of short reads
LEFT_READS="/mnt/nfs/projects/armbrust-metat/gradients3/g3_depth_pa_metat/paired_end/${SAMPLE}.R1.fastq.gz"
RIGHT_READS="/mnt/nfs/projects/armbrust-metat/gradients3/g3_depth_pa_metat/paired_end/${SAMPLE}.R2.fastq.gz"

# run kallisto using the index, sample and short reads (writing out a log)
echo "Running kallisto on ${SAMPLE} against ${INDEX}..."
time kallisto quant -i ${INDEX} -o ${SAMPLE} --threads=${N_THREADS} <(zcat ${LEFT_READS}) <(zcat ${RIGHT_READS}) >> G3PA_depth.${SAMPLE}.kallisto.log

# gzip, rename the output to something project and sample specific:
gzip ${SAMPLE}/abundance.tsv
gzip ${SAMPLE}/run_info.json
mv ${SAMPLE}/abundance.tsv.gz ${PROJECT}/${SAMPLE}.unstranded.abundance.tsv.gz
# rename run_info.json to something project and sample specific:
mv ${SAMPLE}/run_info.json.gz ${PROJECT}/${SAMPLE}.unstranded.run_info.json.gz
}

# Declare variables and local paths:
SAMPLE_LIST="sample_names.txt"
PROJECT="/mnt/nfs/projects/armbrust-metat/gradients3/g3_depth_pa_metat/kallisto/unstranded"
INDEX="NPac.G3PA_depth.id99.nt.idx"
N_THREADS=48 

#######################
# Iterate through the samples in the sample list:

for SAMPLE in $(cat ${SAMPLE_LIST}); do
echo $SAMPLE
time run_kallisto
done


