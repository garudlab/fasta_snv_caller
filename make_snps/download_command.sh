species=Neisseria_gonorrhoeae
species_snvs_dir=../${species}_aligned/snps
genome_accession=PRJEB2999 # Grad et al 2014
# taxon_id=485 #Neisseria gonorrhoeae
ref_accession=GCF_000006845.1
mkdir -p $species_snvs_dir 


echo downloading genomes for $species
# need ncbi datasets tool available in PATH: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/
datasets download genome accession $genome_accession --dehydrated --annotated --assembly-source 'RefSeq'
# another example 
# datasets download genome taxon $taxon_id --dehydrated --annotated --assembly-source 'RefSeq'
mv ncbi_dataset.zip  $species_snvs_dir/ 
unzip $species_snvs_dir/ncbi_dataset.zip -d $species_snvs_dir/accession_genomes
datasets rehydrate --directory $species_snvs_dir/accession_genomes 

echo download ref strain for $species
#download ref strain (for Neisseria gonorrhoeae: FA1090)
datasets download genome accession $ref_accession --annotated --assembly-source 'RefSeq' --include genome,gff3,protein
mv ncbi_dataset.zip $species_snvs_dir/ref_ncbi_dataset.zip 
unzip $species_snvs_dir/ref_ncbi_dataset.zip -d $species_snvs_dir/ref_genomes


### IF USING OWN DATA, STILL RUN BELOW AND POPULATE DIRECTORIES ACCORDINGLY

mkdir -p $species_snvs_dir/${species}_fastas
mkdir -p $species_snvs_dir/${species}_ref_files
mkdir -p $species_snvs_dir/${species}_catalog
cp $species_snvs_dir/*genomes/ncbi_dataset/data/*/*.fna $species_snvs_dir/${species}_fastas/
cp $species_snvs_dir/ref_genomes/ncbi_dataset/data/*/* $species_snvs_dir/${species}_ref_files/

echo generating $species pairs
ref_genome=`ls $species_snvs_dir/${species}_ref_files/*.fna | rev | cut -f 1 -d / | rev | sed 's/.fna//'`
len_isolates=`ls -1 $species_snvs_dir/${species}_fastas | wc -l`
paste <( yes "$ref_genome" | head -n $len_isolates) <(ls -1 $species_snvs_dir/${species}_fastas | sed 's/.fna//') > ${species_snvs_dir}/input_pairs.txt


