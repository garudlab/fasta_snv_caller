#!/bin/bash
# runs annotation

SPECIES=Neisseria_gonorrhoeae
species_snvs_dir=../${SPECIES}_aligned/snps/
hap_dir=../${SPECIES}_aligned/

python3 annotation_utils.py --species $SPECIES --snps_file_path $species_snvs_dir --raw_dir $hap_dir

