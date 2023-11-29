#!/usr/bin/env python

'''
Script for reading in Bracken classifications from Kraken2 and generating aggregated sample outputs.
The preferred file type for this script is *_KRAKEN2_FBAV2/XX-SEA-YY_{G,S}_OUT.bracken
Usage: python generate_aggregated_krona_output.py <input1.mpa> <input2.mpa> ... <inputN.mpa>
'''

import csv, argparse, os
import pandas as pd

parser = argparse.ArgumentParser(description='Generate aggregated Bracken output')
parser.add_argument("-i", "--input_files",nargs='+',help="Path to Bracken file(s). (Accepts multiple)", required=True)
parser.add_argument("-o", "--outfile", required=False,help="Path to output file.", default="Aggregated_Bracken_output.tsv")
#parser.add_argument("-l", "--level", help="Taxonomic level.",required=True,type=str,choices=['Domain','Phylum','Class','Order','Family','Genus','Species'])
args = parser.parse_args()

def find_level(txstring,level):
    txsplit = txstring.split("|")
    l = txsplit[-1]
    if level == "Domain" and l.startswith("d__"):
        return l[3:]
    elif level == "Phylum" and l.startswith("p__"):
        return l[3:]
    elif level == "Class" and l.startswith("c__"):
        return l[3:]
    elif level == "Order" and l.startswith("o__"):
        return l[3:]
    elif level == "Family" and l.startswith("f__"):
        return l[3:]
    elif level == "Genus" and l.startswith("g__"):
        return l[3:]
    elif level == "Species" and l.startswith("s__"):
        return l[3:]
    else:
        return None

def main():
    results = {}
    allnames = pd.Series()
    for infile in args.input_files:
        # Get samplename from first part of filename
        #samplename = infile.split("/")[0
        samplename = os.path.basename(infile).split("_")[0]
        with open (infile,"r") as my_infile:
            my_infile_pd = pd.read_csv(my_infile, sep="\t",header=0, index_col=1, squeeze = True) # Setting tax_id numeric col as index
            my_infile_dict = my_infile_pd.to_dict() # NOTE - NO real reason to convert via dictionary. Let pandas handle everything
            #text = csv.reader(my_infile,delimiter="\t")
            #sampleresults = {}
            #header = next(text)
            #for line in text:
            #    taxa_string = line[0]
            #    for elem in line:
            #        sampleresults[taxa_string][
             #    taxa = find_level(taxa_string,args.level)
            #    count = line[-1]
            #    if taxa is not None:
            #        sampleresults[taxa] = count
            print("Number of counts in file %s: %s" % (samplename, str(my_infile_pd.shape[0])))
            # NOTE: This will truncate so that only index keys are used
            results[samplename] = my_infile_pd["new_est_reads"] # This is already a pd.Series, with taxid as key
            if allnames.empty:
                allnames = my_infile_pd["name"]
            else:
                allnames = pd.concat([allnames,my_infile_pd["name"]]).drop_duplicates()
            #results[samplename] = pd.Series(sampleresults) 
    resultsdf = pd.DataFrame(results).fillna(0)
    print(resultsdf)
    with open(args.outfile,'w') as outfile:
        resultsdf.to_csv(outfile, sep="\t")
    with open("taxid_to_name.tsv",'w') as translationfile:
        allnames.to_csv(translationfile,sep="\t")

if __name__ == '__main__':
    main()
