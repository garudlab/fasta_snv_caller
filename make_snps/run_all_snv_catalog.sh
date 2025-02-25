# Peter Laurin
# Finalized Feb 19, 2025
# bash script to generate snp catalog from isolate fasta files
# requires mummer (I use version 4.0.0rc1) and python2

SPECIES=Neisseria_gonorrhoeae
ref_name=fa_1090 # whatever you want to call it, in case you try multiple refs
ref_fasta=GCF_000006845.1_ASM684v1_genomic
species_snvs_dir=/u/scratch/p/peterlau/m12/${SPECIES}_aligned/snps/


GENOMES_DIR=$species_snvs_dir/${SPECIES}_fastas
PAIRS_IN=${species_snvs_dir}/input_pairs.txt

OUTPUT_DIR=$species_snvs_dir/${SPECIES}_catalog

while IFS= read -r pair
do
    genome1=$(echo $pair | cut -d' ' -f1)
    genome2=$(echo $pair | cut -d' ' -f2)

    GB_PATH_1=$GENOMES_DIR/$genome1.fna
    GB_PATH_2=$GENOMES_DIR/$genome2.fna


    nucmer $GENOMES_DIR/$genome1.fna $GENOMES_DIR/$genome2.fna --prefix $OUTPUT_DIR/$genome1-$genome2 --threads 1
    delta-filter -q -r $OUTPUT_DIR/$genome1-$genome2.delta > $OUTPUT_DIR/$genome1-$genome2.filter.delta
    show-coords $OUTPUT_DIR/$genome1-$genome2.filter.delta > $OUTPUT_DIR/$genome1-$genome2.coords
    show-snps $OUTPUT_DIR/$genome1-$genome2.filter.delta > $OUTPUT_DIR/$genome1-$genome2.snps
    show-diff $OUTPUT_DIR/$genome1-$genome2.filter.delta > $OUTPUT_DIR/$genome1-$genome2.diff

	echo "done whole genome alignment for $genome1 and $genome2"
done < $PAIRS_IN

# get all snps among all alignments
# uses 15G memory -- might want to change!
ls $OUTPUT_DIR/ | grep snps | sed 's/.snps//' | xargs -I[] bash -c 'sed "1,5d" '$OUTPUT_DIR'/[].snps | awk "$0"' '$2 != "." && $3 != "." {printf "%s\t%s||%s||%s||%s\n", "[]", $14, $1, $2, $3}' | cut -f2 | LC_ALL=C sort -k2,2 -S15G --parallel=4 | uniq -c | awk '$1 > 1 {print $2}' > ${species_snvs_dir}/$SPECIES.snps.list

python3 generate_catalog.py --shared ${species_snvs_dir}/$SPECIES.snps.list --in-list <(ls $OUTPUT_DIR/*.snps) --out ${species_snvs_dir}/${ref_name}.catalog.tsv

sed -i "1s/${ref_fasta}-//g" ${species_snvs_dir}/${ref_name}.catalog.tsv

# don't think I need this since I fixed generate catalog error?
#echo $ref_name | xargs -I[] bash -c 'head -n 1 [].catalog.tsv | awk "$0"' '{for (i=1; i <= NF; ++i) {if ($i == "[]") {print "[]",i}}}' | awk '{printf "cut -f1-%s,%s- %s.catalog.tsv > ${SPECIES}.catalog.noAuto.tsv\n", $2-1, $2+1, $1, $1}' | xargs -I[] bash -c "[]"

#python bin_snps_by_source.py --catalog $SPECIES.catalog.noAuto.tsv --name $SPECIES --genome-info <(cut -f1,22 ./genomes_metadata.tsv) > $SPECIES.continentBin.tsv

#python2 identify_ref_allele.py --catalog $ref_name.catalog.noAuto.tsv --name $SPECIES --coords-dir $OUTPUT_DIR --out $ref_name.catalog.noAuto.wtRef.tsv
python3 identify_ref_allele.py --catalog ${species_snvs_dir}/${ref_name}.catalog.tsv --coords-dir $OUTPUT_DIR --out ${species_snvs_dir}/${ref_name}.catalog.final.tsv

rm ${species_snvs_dir}/$SPECIES.snps.list
rm ${species_snvs_dir}/${ref_name}.catalog.tsv
#rm $SPECIES.catalog.noAuto.tsv
