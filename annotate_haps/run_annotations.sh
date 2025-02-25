#!/bin/bash
# runs annotation

SPECIES=Neisseria_gonorrhoeae
species_snvs_dir=/u/scratch/p/peterlau/m12/${SPECIES}_aligned/snps/
hap_dir=/u/scratch/p/peterlau/m12/${SPECIES}_aligned/

python3 annotation_utils.py --species $SPECIES --snps_file_path $species_snvs_dir --raw_dir $hap_dir

