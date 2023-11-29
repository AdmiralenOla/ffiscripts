#!/usr/bin/env python

'''Script for separating OTU counts into Bacteria (taxid: 2) and Fungi (taxid: 4751)'''

import sys, argparse, os
import pandas as pd

parser = argparse.ArgumentParser(description='Split OTU counts tsv file into Bacteria and Fungi')
parser.add_argument("-i", "--input_file",help="Path to OTU counts file (TSV format)", required=True)
parser.add_argument("-l", "--library", help="Path to NCBIs taxidlineage.dmp file", required=True)
parser.add_argument("-m", "--merged", help="Path to NCBIs merged.dmp file", required=True)
#parser.add_argument("-o", "--outfile", required=False,help="Path to output file.", default="Aggregated_Bracken_output.tsv")
#parser.add_argument("-l", "--level", help="Taxonomic level.",required=True,type=str,choices=['Domain','Phylum','Class','Order','Family','Genus','Species'])
args = parser.parse_args()


def main():
    taxidmap = {}
    output_filename_base = os.path.splitext(os.path.basename(args.input_file))[0]
    with open(args.library,'r') as libraryfile:
        #my_library_pd = pd.read_csv(libraryfile,sep="\s+",header=None,index_col=0)
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
                taxidmap[taxid] = "Archea"
            if "10239" in lineage:
                taxidmap[taxid] = "Viruses"
            if "33090" in lineage:
                taxidmap[taxid] = "Viridiplantae"
            if "33208" in lineage:
                taxidmap[taxid] = "Metazoa"
    mergedmap = {}
    with open(args.merged, 'r') as mergedfile:
        my_merged = mergedfile.readlines()
        for row in my_merged:
            taxid_old = row.split("\t|\t")[0]
            taxid_new = row.split("\t|\t")[1].split("\t")[0]
            mergedmap[taxid_old] = taxid_new
            
    with open(args.input_file,'r') as infile:
        my_infile_pd = pd.read_csv(infile, sep="\t",header=0,index_col=0) # Setting tax_id numeric col as index
        data_frame_fungi = my_infile_pd
        #fungidic = {}
        #bacteriadic = {}
        data_frame_bacteria = my_infile_pd
        for i,j in my_infile_pd.iterrows():
            k = i # Set up k to be the same as i, except in missing cases
            if str(i) not in taxidmap:
                print("Not found: %s" % i)
                if str(i) in mergedmap:
                    k = mergedmap[str(i)]
#                elif i == 206324:
#                    k = 2773346
#                elif i == 478866:
#                    k = 2932460
                # elif i == 2656914:
                    # k = 2219224
                # elif i == 2761535:
                    # k = 37331
                # elif i == 2823898:
                    # k = 585529
                # elif i == 2823899:
                    # k = 38303
                # elif i == 2831617:
                    # k = 2960088
                # elif i == 2840473:
                    # k = 2840469
                # elif i == 2841037:
                    # k = 3028070
            if taxidmap[str(k)] != "Fungi":
                data_frame_fungi = data_frame_fungi.drop(i, axis=0)
            if taxidmap[str(k)] != "Bacteria":
                data_frame_bacteria = data_frame_bacteria.drop(i,axis=0)
    with open(output_filename_base + "_fungi.tsv",'w') as outfile:
        data_frame_fungi.to_csv(outfile, sep="\t")
    with open(output_filename_base + "_bacteria.tsv",'w') as outfile:
        data_frame_bacteria.to_csv(outfile,sep="\t")





if __name__ == '__main__':
    main()
