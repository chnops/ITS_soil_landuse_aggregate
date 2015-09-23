## Author: Fan Yang ##
## change pathways before use ##
## need RDPTools (https://github.com/rdpstaff/RDPTools), usearch (http://www.drive5.com/usearch/), python 2.7, biopython module, extended panda-seq (http://rdp.cme.msu.edu/download/RDP_Assembler.tgz) ##
 
NAME="ITS"
LOCATION="/PATH/TO/OUTPUT/DIRECTORY"
INDEX_FILE="/PATH/TO/INDEX/FASTQ"
READ_1="/PATH/TO/READ1/FASTQ"
READ_2="/PATH/TO/READ2/FASTQ"
BINNER="/PATH/TO/bin_reads.py"
UPARSE_PY="/PATH/TO/UPARSE/python_scripts"
REF_ITS="/PATH/TO/fungalits_warcup_trainingdata1/Warcup.fungalITS.fasta"

SEQFILTER="/PATH/TO/RDPTools/SeqFilters.jar"
PANDASEQ="/PATH/TO/RDP_Assembler/pandaseq/pandaseq"
USEARCH="/PATH/TO/usearch70"

cd $LOCATION

echo "assembling pair-ended reads ..."
mkdir 1_"$NAME"_demultiplex 1_"$NAME"_demultiplex/assembled
$PANDASEQ -N -o 10 -e 25 -F -d rbfkms -l 250 -L 280 -f $READ_1 -r $READ_2 1> 1_"$NAME"_demultiplex/assembled/assembled_reads.fastq 2>  1_"$NAME"_demultiplex/assembled/assembled_reads_stats.txt
echo "done."

echo "parsing index file ..."
java -jar $SEQFILTER --seq-file $INDEX_FILE --tag-file $MAP --outdir 1_"$NAME"_demultiplex/parse_index
echo "moving quality trimmed index reads to directory 'trimmed_tags'"
mkdir 1_"$NAME"_demultiplex/parse_index/trimmed_tags
mv 1_"$NAME"_demultiplex/parse_index/result_dir/*/*_trimmed.fasta 1_"$NAME"_demultiplex/parse_index/trimmed_tags/
echo "done."

echo "demultiplexing assembled reads ..."
cd 1_"$NAME"_demultiplex/parse_index/trimmed_tags
python $BINNER $LOCATION/1_"$NAME"_demultiplex/assembled/assembled_reads.fastq
cd $LOCATION
mkdir 1_"$NAME"_demultiplex/demultiplexed_assembled_reads
mv 1_"$NAME"_demultiplex/parse_index/trimmed_tags/*_assem.fastq 1_"$NAME"_demultiplex/demultiplexed_assembled_reads
echo "done."

echo "prune reads to maxee=0.5 ..."
mkdir 2_"$NAME"_uparse 2_"$NAME"_uparse/quality_filtered
cd 1_"$NAME"_demultiplex/demultiplexed_assembled_reads
for i in *.fastq; do $USEARCH -fastq_filter $i -fastq_maxee 0.5 -fastaout $LOCATION/2_"$NAME"_uparse/quality_filtered/"$i"_maxee_0.5.fasta; done
echo "done."

echo "dereplicate reads ..."
mkdir $LOCATION/2_"$NAME"_uparse/derep
cd $LOCATION/2_"$NAME"_uparse/quality_filtered
for i in *_0.5.fasta; do $USEARCH -derep_fulllength $i -output $LOCATION/2_"$NAME"_uparse/derep/"$i"_unique.fasta -sizeout; done
echo "done."

echo "sorting by cluster size and removing singltons ..."
mkdir $LOCATION/2_"$NAME"_uparse/sorted
cd $LOCATION/2_"$NAME"_uparse/derep
for i in *_unique.fasta; do $USEARCH -sortbysize $i -output $LOCATION/2_"$NAME"_uparse/sorted/"$i"_sorted.fa -minsize 2; done
echo "done."

echo "removing chimeras using de novo ..."
mkdir $LOCATION/2_"$NAME"_uparse/chime_denovo
cd $LOCATION/2_"$NAME"_uparse/sorted
for i in *.fa; do $USEARCH -cluster_otus $i -otuid 0.985 -otus $LOCATION/2_"$NAME"_uparse/chime_denovo/"$i"_otu1.fa; done
echo "done."

echo "removing chimeras using references ..."
mkdir $LOCATION/2_"$NAME"_uparse/chime_ref $LOCATION/2_"$NAME"_uparse/chime_ref/stats $LOCATION/2_"$NAME"_uparse/chime_ref/chimeras $LOCATION/2_"$NAME"_uparse/chime_ref/good_otus
cd $LOCATION/2_"$NAME"_uparse/chime_denovo
for i in *_otu1.fa; do $USEARCH -uchime_ref $i -db $REF_ITS -uchimeout $LOCATION/2_"$NAME"_uparse/chime_ref/stats/"$i".uchime -strand plus -selfid -mindiv 1.5 -mindiffs 5 -chimeras $LOCATION/2_"$NAME"_uparse/chime_ref/chimeras/"$i"_chimera.fa -nonchimeras $LOCATION/2_"$NAME"_uparse/chime_ref/good_otus/"$i"_good.fa; done
echo "done."

echo "consolidating all good otus into 1 file ..."
mkdir $LOCATION/2_"$NAME"_uparse/consoilidate_otus
cd $LOCATION/2_"$NAME"_uparse/chime_ref/good_otus
cat *_good.fa > $LOCATION/2_"$NAME"_uparse/consoilidate_otus/cat_otu_good.fa
cd $LOCATION/2_"$NAME"_uparse/consoilidate_otus
$USEARCH -derep_fulllength cat_otu_good.fa -output cat_otu_good_derep.fa -sizeout
$USEARCH -sortbysize cat_otu_good_derep.fa -output cat_otu_good_derep_sorted.fa -minsize 1
$USEARCH -cluster_otus cat_otu_good_derep_sorted.fa -otuid 0.985 -uparse_break -100.0 -otus all_otus.fa

echo "renaming otus that are rid of chimeras to 'OTU_XX' ..."
cd $LOCATION/2_"$NAME"_uparse/consoilidate_otus
python $UPARSE_PY/fasta_number.py all_otus.fa OTU_ > all_otus_renamed.fa
echo "done"

echo "classify each otu ..."
mkdir $LOCATION/final
java -Xmx4g -jar $CLASSIFIER classify -g fungalits_warcup -c 0.5 -f fixrank -o $LOCATION/final/all_otus_renamed_ITS_classified_0.5.txt -h $LOCATION/final/all_otus_renamed_ITS_hier.txt all_otus_renamed.fa

echo "mapping all sequences (including singletons) back to final otu's ..."
mkdir $LOCATION/2_"$NAME"_uparse/mapped $LOCATION/2_"$NAME"_uparse/mapped/uc $LOCATION/2_"$NAME"_uparse/mapped/seqs
cd 2_"$NAME"_uparse/quality_filtered
for i in *.5.fasta; do $USEARCH  -usearch_global $i -db $LOCATION/2_"$NAME"_uparse/consoilidate_otus/all_otus_renamed.fa -strand plus -id 0.985 -uc $LOCATION/2_"$NAME"_uparse/mapped/uc/"$i"_map.uc -matched $LOCATION/2_"$NAME"_uparse/mapped/seqs/"$i"_matched.fa; done
echo "done"

echo "writing out otu table ..."
cd $LOCATION/2_"$NAME"_uparse/mapped/uc
python $UPARSE_PY/uc2otutab.py *.uc > $LOCATION/final/otu_table.txt
