# Peter Laurin
# Feb 24, 2025
# code adapted from Wolff and Laurin "UHGG" codebase
# converts 'Almeida style' snps file into pandas-df style used in Garud lab, 
# with 'contig', 'gene', 'site_pos' and 'site_type' fields, etc. 

import pandas as pd
import sys
import numpy as np
import os
import glob
from collections import defaultdict
import fastparquet as fp

def nuc_map(nuc):
    
    if nuc == "A":
        return "T"
    
    elif nuc == "C":
        return "G"
        
    elif nuc == "T":
        return "A"
    
    elif nuc == "G":
        return "C"
    
    elif nuc == "N":
        return "N"

def fetch_codon_table():
    
    codon_table = pd.read_csv("codons_table.txt",index_col=0)
    
    return(codon_table)
    
def return_gff(gff_path):
    

    df_gff = pd.read_csv(gff_path,comment='#',sep="\t",header=None)

    df_gff = df_gff.dropna()
    df_gff = df_gff[[0,2,3,4,5,6,7,8]]
    df_gff.columns = ["scaffold_id","gene_type","start","end","score","strand","phase","gene_description"]
    #need to account for phase!
    df_gff.loc[:,"start"] = df_gff["start"].astype(int)
    df_gff.loc[:,"end"] = df_gff["end"].astype(int)
    df_gff = df_gff[["scaffold_id","gene_type","start","end","strand","gene_description"]]

    df_gff_desc_lists = df_gff["gene_description"].str.split(";")
    all_gene_dics = [{g.split("=")[0]:g.split("=")[1] for g in df_gff_desc_lists[i]} for i in df_gff.index]
    df_descs = pd.DataFrame(all_gene_dics)
    df_gff.loc[:,df_descs.columns] = df_descs
    df_gff.index = [x.strip('cds-') for x in df_gff["ID"]]
    df_gff = df_gff.loc[(df_gff["gene_type"] == "CDS") & (df_gff["pseudo"] != "true")]
    
    return(df_gff)

def fetch_contig_seq(fasta_path, contig_id):
    #changed to fetching individual contigs (rather than dict) to save memory
        
    contig_seq = []
    add_to_seq = False
    with open(fasta_path) as file:    

        for line in file:

            if line[0] == ">":
                new_contig = line.strip("\n").strip(">")
                if contig_id in new_contig:
                    add_to_seq = True
                else:
                    add_to_seq = False

            elif add_to_seq:

                seq = line.strip("\n")
                contig_seq.extend(seq)

    
    return(contig_seq)

def fetch_aa_sequence(aa_path):
    
    aa_seq = {}
    with open(aa_path) as file:    

        for line in file:

            if line[0] == ">":
                new_cds = line.strip("\n").strip(">").split(" ",maxsplit=2)
                wp_ref = new_cds[0]
                description = new_cds[2]
                aa_seq[wp_ref] = [description,[]]
            else:
                seq = line.strip("\n")
                aa_seq[wp_ref][1].extend(seq)

    
    return(aa_seq)



## takes raw SNV catalogue and splits it by contig
## then saves the split contigs as a pickle for efficient retrieval
def split_contigs(snv_init_path, raw_dir, chunksize=10**5):
    
    
    sys.stdout.write("\tReading in SNV catalogue\n\n")

    os.makedirs(f"{raw_dir}/contigs",exist_ok=True)
    
    #remove files on rerun so doesn't keep appending
    already_written = glob.glob(f"{raw_dir}/contigs/*[0-9].tmp.pq")
    for file in already_written:
        os.remove(file)

    i = 0

    types = defaultdict(pd.Int8Dtype, snp_id=str)
    
    #slow reading csvs - sep by contigs and write to parquets
    current_contig = ""
    for chunk in pd.read_csv(snv_init_path, index_col=0,chunksize=chunksize,sep="\t",na_values=255, dtype=types, engine="c"):

        #assign new index from SNV output from pipeline
        chunk = chunk.reset_index()
        new_ind = chunk["snp_id"].str.split("||",regex=False,expand=True)
        chunk[["Contig", "Pos", "Ref", "Alt"]] = new_ind
        chunk.set_index(["Contig", "Pos", "Ref", "Alt"], inplace=True)
        chunk = chunk.drop("snp_id",axis=1)

        chunk=chunk.loc[chunk.index.dropna()]
        dup_pos = chunk.index.get_level_values("Pos").duplicated()
        chunk = chunk.loc[~dup_pos]
        new_levels = [level.astype(desired_type) for level, desired_type in zip(chunk.index.levels, [str, pd.Int64Dtype(), str, str])]
        chunk.index = chunk.index.set_levels(new_levels)

        chunk_snvs_by_contig = chunk.groupby("Contig")
        for group in chunk_snvs_by_contig:
            contig,df_contig = group
            filename = f"{raw_dir}/contigs/{contig}.tmp.pq"
            if not os.path.isfile(filename):
                df_contig.to_parquet(filename, engine="fastparquet",index=True)
            else:
                df_contig.to_parquet(filename, engine="fastparquet", append=True,index=True) 
            if contig != current_contig:
                sys.stdout.write(f"{current_contig} split and written\n")
                current_contig = contig
        i+=1
        sys.stdout.write(f"\t\tChunk {i} processed\n")



def rewrite_contig_with_genes(contig_id, raw_dir, ref_paths):

    sys.stdout.write(f" \n writing contig {contig_id} with genes! \n")
    
    #remove files on rerun so doesn't keep appending
    already_written = glob.glob(f"{raw_dir}/contigs/*w_genes.tmp.pq")
    for file in already_written:
        os.remove(file)

    df_gff = return_gff(ref_paths["gff"])
    gene_bin_edges_df = df_gff.loc[df_gff.scaffold_id == contig_id][["start","end"]]
    filename = f"{raw_dir}/contigs/{contig_id}.tmp.pq"
    contig_seq = fetch_contig_seq(ref_paths["fna"], contig_id)
    i = 0

    # read in contig

    p_file = fp.ParquetFile(filename)
    #chunk = p_file.to_pandas()
    for chunk in p_file.iter_row_groups(): #filename, chunksize=chunksize, index_col=["Contig", "Pos", "Ref", "Alt"], dtypes=types, engine="c"):
        contig_snv_pos_redux_idx = []
        contig_snv_pos_redux_rows = []
        
        for idx,row in chunk.iterrows():
            contig,pos,ref,alt = idx

            #quick check
            check_pos = int(pos) - 1
            if contig_seq[check_pos] != ref:
                raise Exception(f"for {contig} {pos} snv catalogue state does not match .fna-inferred reference state for contig")

            #get genes in pos 
            pos = int(pos)
            rng_true = (gene_bin_edges_df["start"] < pos)&(gene_bin_edges_df["end"] > pos)
            gffs_true = rng_true.loc[rng_true]
        
            ## if site falls in a gene
            if gffs_true.sum() == 1:
                site_gene_id = gene_bin_edges_df.loc[gffs_true.index[0]].name.strip("cds-")
                idx_new = [contig,str(pos),ref,alt,site_gene_id]
                contig_snv_pos_redux_idx.append(idx_new)
                contig_snv_pos_redux_rows.append(row)     
            
            ## if site does not fall in a gene
            elif gffs_true.sum() == 0:
                idx_new = [contig,str(pos),ref,alt,"non_coding"]
                contig_snv_pos_redux_idx.append(idx_new)
                contig_snv_pos_redux_rows.append(row)     
            
            ## if site falls in multiple genes
            else:
                continue
                #skip for now
                count = gffs_true.sum()    
                print("mult gene")
                for i in range(count):

                    site_gene_id = gffs_true.index[i] 

                    idx_new = [contig,str(pos),ref,alt,site_gene_id]

                    # add index multiple times -- might want to filter to one site
                    contig_snv_pos_redux_idx.append(idx_new)
                    contig_snv_pos_redux_rows.append(row)     

        contig_snv_pos_redux_idx = np.array(contig_snv_pos_redux_idx).T
        idx_redux = pd.MultiIndex.from_arrays(contig_snv_pos_redux_idx,names=["contig","site_pos","ref","alt","gene_id"])

        df_snvs_redux = pd.DataFrame(np.array(contig_snv_pos_redux_rows),index=idx_redux,columns=chunk.columns, dtype=pd.Int8Dtype())
        
        write_filename = f"{raw_dir}/contigs/{contig_id}_w_genes.tmp.pq"

        if not os.path.isfile(write_filename):
            df_snvs_redux.to_parquet(write_filename, engine = "fastparquet")
        else:
            df_snvs_redux.to_parquet(write_filename, engine = "fastparquet", append=True)
        i += 1
        sys.stdout.write(f"Chunk {i} (with genes) processed! \n")

    sys.stdout.write(f" contig {contig_id} with genes processed! \n")

def get_aa_and_substitution(site_pos, ref, alt, start_pos, gene_seq, codon_table,aa_seq, gene, strand):

    pos_index = int(site_pos)
    aa_num = (pos_index - start_pos) // 3
    mut_pos = (pos_index - start_pos) % 3

    codon = gene_seq[aa_num*3:aa_num*3+3]
    
    if strand == "-":
        codon = codon[::-1]
        codon = [nuc_map(n) for n in codon]
        ref = nuc_map(ref)
        alt = nuc_map(alt)
        if mut_pos == 2:
            mut_pos = 0
        elif mut_pos == 0:
            mut_pos = 2
        aa_num = len(gene_seq)//3 - aa_num - 1 
    
    assert codon[mut_pos] == ref, f"\t reference state {ref} nor alt state {alt} match reference {codon[mut_pos]}, in codon {codon} at site {site_pos}. Strand = '{strand}' \n"
     
        #if codon[mut_pos] != nuc_map(ref) and codon[mut_pos] != nuc_map(alt):
        #print(f"\t reference state {ref} nor alt state {alt} match reference {codon[mut_pos]}, in codon {codon} at site {site_pos}. Strand = '{strand}' \n")
        #print(f"\t reference state {ref} nor alt state {alt} match reference {codon[mut_pos]}, in codon {codon} at site {site_pos}. Strand = '{strand}' \n")
        #with open("substitution_error_A_muciniphilia_C_contig_1.err", "a") as err_file:
        #    err_file.write(f"\t reference state {ref} nor alt state {alt} match reference {codon[mut_pos]}, in codon {codon} at site {site_pos}. Strand = '{strand}' \n")
    codon_mut = codon.copy()
    
    codon_mut[mut_pos] = alt
    

    codon_str = "".join(codon)
    codon_mut_str = "".join(codon_mut)

    if "N" not in codon_str:
        ref_aa = codon_table.loc[codon_str].AA
    else:
        ref_aa = "X"
    if "N" not in codon_mut_str:
        alt_aa = codon_table.loc[codon_mut_str].AA
    else:
        alt_aa = "X"
    
    if gene in aa_seq and ref_aa != "*":
        if ref_aa != aa_seq[gene][1][aa_num]:
            if (aa_num == 0) & (aa_seq[gene][1][aa_num] == "M") & (ref_aa in ["V", "L", "I"]): #alternative start codon
                if alt_aa in ["M", "V", "L", "I"]:
                    return(ref_aa, alt_aa, "syn") 
                else:
                    return(ref_aa, alt_aa, "nonsyn")
            else:  
                print(f"\t reference amino acid {aa_seq[gene][1][aa_num]} does not match calculated amino acid {ref_aa}")
                return("NC","NC","NC")


    if ref_aa != alt_aa:
        sub_stat = "nonsyn"
    else:
        sub_stat = "syn"

    return(ref_aa, alt_aa, sub_stat) 


def write_contig_haplotype(contig_id, raw_dir, ref_paths):
    
    annotation_dir = f"{raw_dir}/snp_annotations/"
    hap_dir = f"{raw_dir}/haplotypes/"
    
    os.makedirs(hap_dir,exist_ok=True)
    os.makedirs(annotation_dir,exist_ok=True)

    #check to see if files already exist
    files_in_hap_dir = glob.glob(f"{hap_dir}/*{contig_id}*")
    files_in_ann_dir = glob.glob(f"{annotation_dir}/*{contig_id}*")
    for file in (files_in_hap_dir + files_in_ann_dir):
        os.remove(file)


    df_genes_filename = f"{raw_dir}/contigs/{contig_id}_w_genes.tmp.pq"

    codon_table = fetch_codon_table()
    df_gff = return_gff(ref_paths["gff"])
    contig_nuc_seq = fetch_contig_seq(ref_paths["fna"], contig_id)
    aa_seq = fetch_aa_sequence(ref_paths["faa"])

    p_file = fp.ParquetFile(df_genes_filename)
    # chunk = p_file.to_pandas()
    for chunk in p_file.iter_row_groups():
        chunk = chunk.astype(pd.Int8Dtype())
        chunk_gb = chunk.groupby("gene_id")

        ann_dic = {}
        change_dic = {}

        ## assign syn/non syn status to coding sites, record change 
        for gene_id in chunk_gb.groups.keys():
            
            #gene_id = "GUT_GENOME143197_00002"
            #sys.stdout.write(f"\t{gene_id}\n")

            if gene_id == "non_coding" or "ncRNA" in gene_id:
                
                for idx in chunk_gb.get_group(gene_id).index:

                    contig,pos,ref,alt,gene_id = idx

                    idx = (contig,pos,ref,alt,gene_id)

                    ann_dic[idx] = "NC"

                    change_dic[idx] = ["NC","NC"]

            else:   #coding

                strand = df_gff.loc[gene_id].strand

                start_pos = df_gff.loc[gene_id].start   # PL 0 based indexing?
                end_pos = df_gff.loc[gene_id].end
                gene_seq = contig_nuc_seq[start_pos-1:end_pos] 

                #if strand == "-":
                #    #reverse seq
                #    gene_seq = gene_seq[::-1]
                #    gene_seq = [nuc_map(n) for n in gene_seq]

                for idx in chunk_gb.get_group(gene_id).index:
                        
                    contig,pos,ref,alt,gene = idx

                    idx = (contig,pos,ref,alt,gene)

                    ref_aa, alt_aa, ann = get_aa_and_substitution(pos, ref, alt, start_pos, gene_seq, codon_table, aa_seq, gene, strand)
                    ann_dic[idx] = ann
                    change_dic[idx] = [ref_aa,alt_aa]

                ## utility to check if amino acids are correct against .faa file
                ## need to have written gene_AA_dic first
                #sys.stderr.write(f"Reference AA: {codon_table.loc[aa_str].AA} — True Reference AA: {gene_AA_dic[gene_id][aa_num]} — Alternate AA: {codon_table.loc[aa_mut_str].AA}\n")

                #sys.stderr.write(f"Reference AA: {codon_table.loc[aa_str].AA} — True Reference AA: {gene_AA_dic[gene_id][aa_num]} — Alternate AA: {codon_table.loc[aa_mut_str].AA}\n")


        sys.stdout.write("\nChunk annotated\n\n")            

        #format files for writing

        ann = pd.Series(ann_dic)            

        idxs = change_dic.keys()
        idxs = np.array(list(idxs)).T
        idxs_arr = pd.MultiIndex.from_arrays(idxs,names=["contig","site_pos","ref","alt","gene_id"])
        df_change = pd.DataFrame.from_dict(change_dic,orient="index",columns=["ref_aa","alt_aa"])
        df_change.index = idxs_arr

        chunk["site_type"] = ann
        chunk.set_index('site_type', append=True, inplace=True)

        df_change["site_type"] = ann
        df_change.set_index('site_type', append=True, inplace=True)

        level_to_change = 1
        chunk.index = chunk.index.set_levels(chunk.index.levels[level_to_change].astype(int), level=level_to_change)
        df_change.index = df_change.index.set_levels(df_change.index.levels[level_to_change].astype(int), level=level_to_change)
    
        chunk = chunk.sort_index(level="site_pos")
    
        df_change = df_change.loc[chunk.index]
    
        chunk = chunk.droplevel(["ref","alt"])

        df_change = df_change.reorder_levels(["contig","gene_id","site_pos","site_type","ref","alt"])
        chunk = chunk.reorder_levels(["contig","gene_id","site_pos","site_type"])

        hap_filename = f"{hap_dir}/{contig_id}_haplotypes.csv"
        ann_filename =  f"{annotation_dir}/{contig_id}_annotations.csv"
        if not os.path.isfile(hap_filename):
                chunk.to_csv(hap_filename)
                df_change.to_csv(ann_filename)
        else:
            chunk.to_csv(hap_filename, mode="a",header=False)
            df_change.to_csv(ann_filename, mode="a",header=False)
    
    sys.stdout.write(f"{contig} written\n\n")


    ####################################
    ## Deprecated: polarize by consensus
    ## Now, polarize downstream
    ####################################
    
##    repolarize by consensus
#     to_repolarize_idx = df_snvs_redux.loc[df_snvs_redux.mean(axis=1) > 0.5].index
#     df_snvs_redux.loc[to_repolarize_idx] = 1 - df_snvs_redux.loc[to_repolarize_idx]
#     to_repolarize_AA = df_change.loc[to_repolarize_idx].values
#     repolarized_AA = to_repolarize_AA[:,[1,0]]
#     df_change.loc[to_repolarize_idx] = repolarized_AA 
#    sys.stderr.write("Contig polarized\n\n")  

    
if __name__ == "__main__":
    
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--species',
                        help="",
                        type=str, default="Neisseria_gonorrhoeae")
    parser.add_argument("--snps_file_path",
                        type=str, required=True)
    parser.add_argument("--raw_dir",
                        type=str, required=True)

    args = parser.parse_args()
    species = args.species
    snps_file_path = args.snps_file_path
    raw_dir = args.raw_dir

    snps_files =  glob.glob(f"{snps_file_path}/*")
    al_contig = next((i for i in snps_files if i.endswith("final.tsv")),None) 

    ref_files = glob.glob(f"{snps_file_path}/{species}_ref_files/*")
    ref_paths = {"gff":"", "fna":"", "faa":""}
    for extension in ref_paths.keys():
        ref_paths[extension] += next((i for i in ref_files if i.endswith(extension)),None)


    
    sys.stdout.write(f"\nProcessing {species}\n\n")
    
    sys.stdout.write("Splitting contigs\n\n")
    split_contigs(al_contig, raw_dir)
    sys.stdout.write("\n\n\nContigs split, now performing annotation\n\n")
    
    contig_list = glob.glob(f"{raw_dir}/contigs/*[0-9].tmp.pq")
    contig_list = [os.path.basename(x) for x in contig_list]
    contig_list = [c[:-7] for c in contig_list] 
    
    for contig in contig_list:
        
        sys.stdout.write(f"Annotating {contig} contig")
        rewrite_contig_with_genes(contig, raw_dir, ref_paths)
        write_contig_haplotype(contig, raw_dir, ref_paths)
    
    sys.stdout.write("\nAll contigs processed")