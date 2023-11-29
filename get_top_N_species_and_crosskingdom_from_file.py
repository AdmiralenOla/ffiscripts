#!/usr/bin/env python

'''
Script for extracting the top X species from an AGGREGATED DECONTAMMED Bracken OTU matrix, with avg. abundance and prevalence.
'''

import argparse, os, sys
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description='Get top species from OTU matrix')
parser.add_argument("-i", "--input_file", help="Path to count matrix TSV file.", required=True)
#parser.add_argument("-o", "--output_file", help='Output file', required=True,type=str)
parser.add_argument("-n", "--num_species", help='Number of species in results',default=20,type=int)
parser.add_argument("-t", "--taxids", help="Taxid to name table", required=True)
parser.add_argument("--print_individuals", help="Whether to output individual top files", default=False, action='store_true')
parser.add_argument("--print_crosskingdom", help="Whether to output cross-kingdom files", default=False, action='store_true')
parser.add_argument("--print_topspecies", help="Whether to output top N files", default=False, action='store_true')
parser.add_argument("-l", "--library", help="Path to NCBIs taxidlineage.dmp file. (For cross-kingdom only).", required=False)
parser.add_argument("-m", "--merged", help="Path to NCBIs merged.dmp file (For cross-kingdom only).", required=False)
args = parser.parse_args()


def get_top_species_per_sample(col, num_species, colname):
    '''
    Takes in a pd.Series and the num_species number of top species to return. Then writes that file to output dir
    '''
    sorted_col = col.sort_values(ascending=False)
    sorted_perc = sorted_col/sorted_col.sum()*100
    top_n = sorted_perc[:num_species]
    # Take percentages
    top_df = pd.DataFrame(top_n)
    #print(top_df)
    this_file_name = "INDIVIDUAL/" + str(colname) + "_top_" + str(num_species) + ".tsv"
    #print(this_file_name)
    top_df.to_csv(this_file_name,sep="\t") # This function should support writing file names from string and not just handle
    
def process_NCBI_files(taxidlineage, merged):
    taxidmap = {}
    with open(taxidlineage,'r') as libraryfile:
        my_library = libraryfile.readlines()
        for row in my_library:
            taxid = row.split("\t|\t")[0]
            lineage = row.split("\t|\t")[1].split(" ")
            # i = taxid, j = the row
            if "4751" in lineage:
                taxidmap[taxid] = "Fungi"
            if "2" in lineage:
                taxidmap[taxid] = "Bacteria"
            if "2157" in lineage:
                taxidmap[taxid] = "Archaea"
            if "10239" in lineage:
                taxidmap[taxid] = "Viruses"
            if "33090" in lineage:
                taxidmap[taxid] = "Viridiplantae"
            if "33208" in lineage:
                taxidmap[taxid] = "Metazoa"
    mergedmap = {}
    with open(merged, 'r') as mergedfile:
        my_merged = mergedfile.readlines()
        for row in my_merged:
            taxid_old = row.split("\t|\t")[0]
            taxid_new = row.split("\t|\t")[1].split("\t")[0]
            mergedmap[taxid_old] = taxid_new
    return taxidmap, mergedmap

def main():
    if not any([args.print_crosskingdom, args.print_individuals, args.print_topspecies]):
        sys.exit("Error: Must output either crosskingdom, individuals or topspecies")
    if args.print_crosskingdom and (args.library is None or args.merged is None):
        sys.exit("Error: print_crosskingdom required arguments library and merged")
    with open (args.input_file,"r") as my_infile:
        my_infile_pd = pd.read_csv(my_infile, sep="\t",index_col=0, skip_blank_lines=True,keep_default_na=False) # squeeze = True
        # Normalize to 10.000.000 OTU counts
        my_df = (my_infile_pd/my_infile_pd.sum())*10000000 # Considering np.ceil, but this will make the counts not be exactly 10M
        #colsums = my_infile_pd.sum(axis=0)
        #rowsums = my_infile_pd.sum(axis=1)
        # TODO:
        # Convert taxid labels to strings
        translation = pd.read_csv(open(args.taxids,'r'),sep="\t",index_col=0)
        scientific_names = []
        for taxa in my_df.index:
            scientific_names.append(translation.loc[taxa]["name"])
        my_df_taxid = my_df.copy(deep=True)
        my_df.index = scientific_names
        # Write top N species per column (= per sample)
        if args.print_individuals:
            for col in my_df:
                get_top_species_per_sample(my_df[col],args.num_species, col)
       
        # Get MATRIX-wide top N
        taxasums = my_df.sum(axis=1)
        taxasums_id = my_df_taxid.sum(axis=1)
        taxasort = taxasums.sort_values(ascending=False)
        taxasort_id = taxasums_id.sort_values(ascending=False)
        taxaperc = taxasort/taxasort.sum()*100
        top_n_matrix = taxaperc[:args.num_species]
        
        # Get MATRIX-wide prev. per row in top_n_matrix
        prevalences = []
        for taxa in top_n_matrix.index:
            prevalences.append(sum(my_df.loc[taxa] > 0)/len(my_df.columns)*100)
        
        # Get taxonomic info from NCBI files
        output_filename_base = os.path.splitext(os.path.basename(args.input_file))[0]
        if args.print_crosskingdom:
            taxidmap, mergedmap = process_NCBI_files(args.library, args.merged)
            counts = {"bacteria" : 0.0, "virus" : 0.0, "archaea" : 0.0, "eukaryota" : 0.0, "fungi" : 0.0, "viridiplantae" : 0.0, "metazoa" : 0.0}
            denominator = taxasort_id.sum()
            for taxa in taxasort_id.index:
                taxa_use = str(taxa)
                if taxa_use in mergedmap:
                    taxa_use = mergedmap[taxa_use]
                if taxa_use in taxidmap:
                    if taxidmap[taxa_use] == "Bacteria":
                        counts["bacteria"] += taxasort_id[taxa]
                    elif taxidmap[taxa_use] == "Fungi":
                        counts["fungi"] += taxasort_id[taxa]
                        counts["eukaryota"] += taxasort_id[taxa]
                    elif taxidmap[taxa_use] == "Archaea":
                        counts["archaea"] += taxasort_id[taxa]
                    elif taxidmap[taxa_use] == "Viruses":
                        counts["virus"] += taxasort_id[taxa]
                    elif taxidmap[taxa_use] == "Metazoa":
                        counts["metazoa"] += taxasort_id[taxa]
                        counts["eukaryota"] += taxasort_id[taxa]
                    elif taxidmap[taxa_use] == "Viridiplantae":
                        counts["viridiplantae"] += taxasort_id[taxa]
                        counts["eukaryota"] += taxasort_id[taxa]
                    else:
                        print("Not assigned: %s" % taxa_use)
                else:
                    print("Not found: %s" % taxa_use)
                    sys.exit()
            counts_perc = {a: [(counts[a]/denominator*100)] for a in counts}
            crosskingdom_results = pd.DataFrame(counts_perc)
            crosskingdom_results.to_csv("CROSSKINGDOM/" + output_filename_base + "_CROSSKINGDOM.tsv", sep="\t", index=False)
        
        if args.print_topspecies:
            overall_results = pd.DataFrame({"Percentage_of_total": top_n_matrix, "Prevalence": prevalences}) # "Mean abundance": mean_abundances
            overall_results.to_csv("OVERALL/" + os.path.splitext(os.path.basename(args.input_file))[0] + "_RESULTS.tsv",sep="\t")


if __name__ == "__main__":
    main()
