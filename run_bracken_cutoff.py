#!/usr/bin/env python
'''
Script that accepts Kraken2 file(s) and runs Bracken with an x% threshold. Default x is 0.005%. NOTE: You must have Bracken installed in a conda environment named "bracken".
'''

import argparse, os, sys, subprocess
import pandas as pd

def range_type(astr, min=0.0, max=100.0):
    value = float(astr)
    if min<= value <= max:
        return value
    else:
        raise argparse.ArgumentTypeError('value not in range %s-%s'%(min,max))


parser = argparse.ArgumentParser(description='Generate Bracken output with x% cutoff')
parser.add_argument("-i", "--input_files",nargs='+',help="Path to Kraken2 file(s). (Accepts multiple)", required=True)
parser.add_argument("-c", "--cutoff", help='Cutoff in Bracken, by percent. Default 0.005.', metavar="[0-100]",type=range_type, default=0.005)
parser.add_argument("--db", help='Kraken2 database. Default=/media/ubuntu/Elements/NEWPIPELINE_MetaAIR/reference_genomes/Kraken2/NCBInr_7_22', default='/media/ubuntu/Elements/NEWPIPELINE_MetaAIR/reference_genomes/Kraken2/NCBInr_7_22',type=str)
args = parser.parse_args()


def main():
    for infile in args.input_files:
        with open (infile,"r") as my_infile:
            my_infile_pd = pd.read_csv(my_infile, sep="\t",header=None, names=['Percent','TotalReads','SpecificReads','Label','taxid','name'],index_col=4) # squeeze = True
            #print(my_infile_pd["TotalReads"][0]) # Unclassified reads
            #print(my_infile_pd["TotalReads"][1]) # Root reads
            #print(my_infile_pd["TotalReads"][2])
            total_reads = int(my_infile_pd["TotalReads"][0]) + int(my_infile_pd["TotalReads"][1])
            outfile = os.path.splitext(infile)[0] + "_" + str(args.cutoff) + ".bracken"
            outreport = os.path.splitext(infile)[0] + "_" + str(args.cutoff) + ".bracken_report"
            threshold = int(total_reads * args.cutoff / 100)
            print("Output file: " + outfile)
            print("Minimum number of reads: " + str(threshold))
            activate_command = f"conda activate bracken"
            deactivate_command =f"conda deactivate"
            #subprocess.run(activate_command, shell=True, check=True)
            
            cmd = "mamba run -n bracken bracken -d %s -i %s -o %s -w %s -r 150 -l S -t %s" % (args.db, infile, outfile, outreport, threshold)
            try:
                subprocess.run(cmd, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error: {e}")
                
            #subprocess.run(deactivate_command, shell=True, check=True)
            
            #input_file1 = subprocess.run(["mamba", "activate" "bracken"], stdout=subprocess.PIPE, text=True, input="Hello How are you?")

if __name__ == '__main__':
    main()