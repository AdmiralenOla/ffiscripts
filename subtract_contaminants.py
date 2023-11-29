#!/usr/bin/env python

'''
Script for subtracting contaminants from a read counts matrix
'''

import argparse, os, sys, subprocess, yaml
import pandas as pd

parser = argparse.ArgumentParser(description='Set contaminant taxa to 0')
parser.add_argument("-i", "--input_file", help="Path to count matrix TSV file.", required=True)
parser.add_argument("-c", "--contaminants", help="YAML file with contaminants.", required=True)
parser.add_argument("-t", "--type", help="Type of contaminants to exclude. Default=All", default="All", type=str)
parser.add_argument("-o", "--output_file", help='Output file', required=True,type=str)
args = parser.parse_args()

def main():
    with open (args.input_file,"r") as my_infile:
        my_infile_pd = pd.read_csv(my_infile, sep="\t",index_col=0) # squeeze = True
        with open(args.contaminants,'r') as my_yaml_input:
            my_yaml = yaml.safe_load(my_yaml_input)["Karis Pankitome organisms"]
            if args.type not in my_yaml:
                sys.exit("Invalid type specified. Options with this file are %s." % ",".join(k for k in my_yaml.keys()))
            else:
                use_yaml = my_yaml[args.type]
            print(my_infile_pd)
            for taxa in use_yaml["ids"]:
                # Blank input matrix for this taxa
                # NOTE! Only blank taxa IF they exist in the input!!
                if taxa in my_infile_pd.index:
                    my_infile_pd.loc[taxa] = 0
            with open(args.output_file,'w') as my_outfile:
                #pd.to_csv(my_outfile,sep="\t")
                my_infile_pd.to_csv(my_outfile,sep="\t")

if __name__ == "__main__":
    main()
