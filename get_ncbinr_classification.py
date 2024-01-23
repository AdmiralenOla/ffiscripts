#!/usr/bin/env python


import os, sys, argparse
import pandas as pd


parser = argparse.ArgumentParser(description='Get crosskingdom info from Kraken2 reports')
parser.add_argument("-i", "--input_file", nargs='+', help="Path to Kraken2 report file. (Accepts multiple)", required=True)
parser.add_argument("-o", "--output_file", help='Output file', default='crosskingdom_output.tsv',type=str)
args = parser.parse_args()

samples_to_include = ["SL310935", "SL310936", "SL310937", "SL310938", "SL310939", "SL310940", "SL310941", "SL310942", "SL310943", "SL310944", "SL310945", "SL310946", 
"SL310947", "SL310948", "SL310949", "SL310950", "SL310951", "SL310952", "SL310953", "SL310954", "SL310955", "SL310956", "SL310957", "SL310958", "SL310959", "SL310960", 
"SL310961", "SL310962", "SL310963", "SL310964", "SL310965", "SL310966", "SL310967", "SL310968", "SL310969", "SL310970", "SL310971", "SL310972", "SL310973", "SL310974", 
"SL310975", "SL310976", "SL310977", "SL310978", "SL310979", "SL310980", "SL310981", "SL310982", "SL310983", "SL310984", "SL310985", "SL310986", "SL310987", "SL310988", 
"SL310989", "SL310990", "SL310991", "SL310992", "SL310993", "SL310994", "SL310995", "SL310996", "SL310997", "SL310998", "SL310999", "SL311000", "SL311001", "SL311002", 
"SL311003", "SL311004", "SL311005", "SL311006", "SL311007", "SL311008", "SL311009", "SL311010", "SL311011", "SL311012", "SL311013", "SL311014", "SL311015", "SL311016", 
"SL311017", "SL311018", "SL311019", "SL311020", "SL311021", "SL311022", "SL311023", "SL311024", "SL311025", "SL311026", "SL311027", "SL311028", "SL311029", "SL311030", 
"SL335676", "SL335677", "SL335678", "SL335679", "SL335680", "SL335681", "SL335682", "SL335683", "SL335684", "SL335685", "SL335686", "SL335687", "SL335688", "SL335689", 
"SL335690", "SL335691", "SL335692", "SL335693", "SL335694", "SL335695", "SL335696", "SL335697", "SL335698", "SL335699", "SL335700", "SL335701", "SL335702", "SL335703", 
"SL335704", "SL335705", "SL335706", "SL335707", "SL335708", "SL335709", "SL335710", "SL335711", "SL335712", "SL335713", "SL335714", "SL335715", "SL335716", "SL335717", 
"SL335718", "SL335719", "SL335720", "SL335721", "SL335722", "SL335723", "SL335739", "SL335742", "SL335743", "SL335744", "SL335745", "SL335746", "SL335747", "SL335748", 
"SL335751", "SL335752", "SL335753", "SL335755", "SL335756", "SL335757", "SL335758", "SL335759", "SL335760", "SL335761", "SL335762", "SL335763", "SL335764", "SL335765", 
"SL335766", "SL335767", "SL335768", "SL335769", "SL335770", "SL335771", "SL342363", "SL342364", "SL342365", "SL342366", "SL342367", "SL342368", "SL342369", "SL342370", 
"SL342371", "SL342372", "SL342373", "SL342374", "SL342375", "SL342376", "SL342377", "SL342378", "SL342379", "SL342380", "SL342381", "SL342382", "SL342383", "SL342384", 
"SL342385", "SL342386", "SL342387", "SL342388", "SL342389", "SL342390", "SL342391", "SL342392", "SL342393", "SL342394", "SL342395", "SL342396", "SL342397", "SL342398", 
"SL342399", "SL342400", "SL342401", "SL342402", "SL342403", "SL342404", "SL342405", "SL342406", "SL342407", "SL342408", "SL342409", "SL342410", "SL342411", "SL342412", 
"SL342413", "SL342414", "SL342415", "SL342416", "SL342417", "SL342418", "SL342419", "SL342420", "SL342421", "SL342422", "SL342423", "SL342424", "SL342425", "SL342426", 
"SL342427", "SL342428", "SL342429", "SL342430", "SL342431", "SL342432", "SL342433", "SL342434", "SL342435", "SL342436", "SL342437", "SL342438", "SL342439", "SL342440", 
"SL467134", "SL467135", "SL467136", "SL467137", "SL467138", "SL467139", "SL467140", "SL467141", "SL467142", "SL467143", "SL467144", "SL467145", "SL467146", "SL467147", 
"SL467148", "SL467149", "SL467150", "SL467151", "SL467152", "SL467153", "SL467154", "SL467155", "SL467156", "SL467157", "SL467158", "SL467159", "SL467160", "SL467161", 
"SL467162", "SL467163", "SL467164", "SL467165", "SL467166", "SL467167", "SL467168", "SL467169", "SL467170", "SL467171", "SL467172", "SL467173", "SL467174", "SL467175", 
"SL467176", "SL467177", "SL467178", "SL467179", "SL467180", "SL467182", "SL467183", "SL467184", "SL467185", "SL467186", "SL467187", "SL467188", "SL467189", "SL467190", 
"SL467191", "SL467192", "SL467193", "SL467194", "SL467195", "SL467196", "SL467197", "SL467198", "SL467199", "SL467200", "SL467201", "SL467202", "SL467203", "SL467204", 
"SL467205", "SL467206", "SL467207", "SL467208", "SL467209", "SL467210", "SL467211", "SL467212", "SL467213", "SL467214", "SL467215", "SL467216", "SL467217", "SL467218", 
"SL467219", "SL467220", "SL467221", "SL467222", "SL467223", "SL467224", "SL467225", "SL467226", "SL467227", "SL469677", "SL469678", "SL469679", "SL469680", "SL469681", 
"SL469682", "SL469683", "SL469684", "SL469685", "SL469686", "SL469687", "SL469688", "SL469689", "SL469690", "SL469691", "SL469692", "SL469693", "SL469694", "SL469695", 
"SL469696", "SL469697", "SL469698", "SL469699", "SL469700", "SL469701", "SL469702", "SL469703", "SL469704", "SL469705", "SL469706", "SL469707", "SL469708", "SL469709", 
"SL469710", "SL469711", "SL469712", "SL469713", "SL469714", "SL469715", "SL469716", "SL469717", "SL469718", "SL469719", "SL469720", "SL469721", "SL469722", "SL469723", 
"SL469724", "SL469725", "SL469726", "SL469727", "SL469728", "SL469729", "SL469730", "SL469731", "SL469732", "SL469733", "SL469734", "SL469735", "SL469736", "SL469737", 
"SL469738", "SL469739", "SL469740", "SL469741", "SL469742", "SL469743", "SL469744", "SL469745", "SL469746", "SL469747", "SL469748", "SL469749", "SL469750", "SL469751", 
"SL469752", "SL469753", "SL469754", "SL469755", "SL469756", "SL469757", "SL469758", "SL469759", "SL469760", "SL469761", "SL469762", "SL469763", "SL469764", "SL469765", 
"SL469766", "SL469767", "SL469768", "SL469769", "SL469770", "SL469771", "SL470203", "SL470205", "SL470206", "SL470209", "SL470210", "SL470211", "SL470212", "SL470213", 
"SL470214", "SL470215", "SL470216", "SL470217", "SL470218", "SL470219", "SL470220", "SL470221", "SL470222", "SL470223", "SL470224", "SL470225", "SL470226", "SL470227", 
"SL470228", "SL470229", "SL470230", "SL470231", "SL470232", "SL470233", "SL470234", "SL470235", "SL470236", "SL470237", "SL470238", "SL470239", "SL470240", "SL470241", 
"SL470242", "SL470243", "SL470244", "SL470245", "SL470246", "SL470247", "SL470248", "SL470249", "SL470250", "SL470251", "SL470252", "SL470253", "SL470254", "SL470255", 
"SL470256", "SL470257", "SL470258", "SL470259", "SL470260", "SL470261", "SL470262", "SL470263", "SL470264", "SL470265", "SL470266", "SL470267", "SL470268", "SL470269", 
"SL470270", "SL470271", "SL470272", "SL470273", "SL470274", "SL470275", "SL470296", "SL469867", "SL469868", "SL469869", "SL469870", "SL469871", "SL469872", "SL469873", 
"SL469874", "SL469875", "SL469876", "SL469877", "SL469878", "SL469879", "SL469880", "SL469881", "SL469882", "SL469883", "SL469884", "SL469885", "SL469886", "SL469887", 
"SL469888", "SL469889", "SL469890", "SL469891", "SL469892", "SL469893", "SL469894", "SL469895", "SL469896", "SL469897", "SL469898", "SL469899", "SL469900", "SL469901", 
"SL469902", "SL469903", "SL469904", "SL469905", "SL469906", "SL469907", "SL469908", "SL469909", "SL469910", "SL469911", "SL469912", "SL469913", "SL469914", "SL469915", 
"SL469916", "SL469917", "SL469918", "SL469919", "SL469920", "SL469921", "SL469922", "SL469923", "SL469924", "SL469925", "SL469926", "SL469927", "SL469928", "SL469929", 
"SL469930", "SL469931", "SL469932", "SL469933", "SL469934", "SL469935", "SL469936", "SL469937", "SL469938", "SL469939", "SL469940", "SL469941", "SL469942", "SL469943", 
"SL469944", "SL469945", "SL469946", "SL469947", "SL469948", "SL469949", "SL469950", "SL469951", "SL469952", "SL469953", "SL469954", "SL469955", "SL469956", "SL469957", 
"SL469958", "SL469959", "SL469960", "SL469961", "SL469962", "SL470282", "SL470283", "SL470284", "SL470285", "SL470286", "SL470287", "SL470288", "SL470289", "SL470290", 
"SL470291", "SL470292", "SL470293", "SL470299", "SL470300", "SL470301", "SL470302", "SL470303", "SL470304", "SL470305", "SL470306", "SL470307", "SL470308", "SL470309", 
"SL470310", "SL470311", "SL470312", "SL470313", "SL470314", "SL470315", "SL470317", "SL470318", "SL470319", "SL470320", "SL470321", "SL470322", "SL470323", "SL470324", 
"SL470325", "SL470326", "SL470327", "SL470328", "SL470329", "SL470330", "SL470332", "SL470333", "SL470336", "SL470412", "SL470413", "SL470414", "SL470415", "SL470416", 
"SL470417", "SL470418", "SL470419", "SL470420", "SL470421", "SL470422", "SL470423", "SL470424", "SL470425", "SL470426", "SL470427", "SL470428", "SL470429", "SL470430", 
"SL470431", "SL470432", "SL470433", "SL470434", "SL470435", "SL470436", "SL470437", "SL470438", "SL470439", "SL470440", "SL470441", "SL470442", "SL470443", "SL470444", 
"SL470445", "SL470446", "SL470447", "SL470448", "SL470449", "SL470450", "SL470451", "SL470452", "SL470453", "SL470454", "SL470455", "SL470456", "SL470457", "SL470458", 
"SL470459", "SL470460", "SL470461", "SL470462", "SL470463", "SL470464", "SL470465", "SL470466", "SL470467", "SL470468", "SL470469", "SL470470", "SL470471", "SL470472", 
"SL470473", "SL470474", "SL470475", "SL470476", "SL470477", "SL470478", "SL470479", "SL470480", "SL470481", "SL470482", "SL470483", "SL470484", "SL470485", "SL470486", 
"SL470487", "SL470488", "SL470489", "SL470490", "SL470491", "SL470492", "SL470493", "SL470494", "SL470495", "SL470496", "SL470497", "SL470498", "SL470499", "SL470500", 
"SL470501", "SL470502", "SL470503", "SL470504", "SL470505", "SL470506", "SL470507"]

def main():
    results = {}
    for report in args.input_file: # Repeat for each input file
        samplename = os.path.basename(report).split("_")[0]
        if samplename not in samples_to_include:
            continue
        with open(report,'r') as inhandle:
            myf = pd.read_csv(inhandle,sep="\t",index_col=4,names=['Percent','ReadsTotal','ReadsExact','Classification','taxid','Name'])
            bacteria = 0 if 2 not in myf.index else myf.loc[2]["ReadsTotal"]
            archaea = 0 if 2157 not in myf.index else myf.loc[2157]["ReadsTotal"]
            unclassified = 0 if 0 not in myf.index else myf.loc[0]["ReadsTotal"]
            fungi = 0 if 4751 not in myf.index else myf.loc[4751]["ReadsTotal"]
            metazoa = 0 if 33208 not in myf.index else myf.loc[33208]["ReadsTotal"] 
            virus = 0 if 10239 not in myf.index else myf.loc[10239]["ReadsTotal"]
            viridiplantae = 0 if 33090 not in myf.index else myf.loc[33090]["ReadsTotal"]
            resultsdic = {"Unclassified": unclassified, "Bacteria":bacteria,"Archaea":archaea, "Fungi":fungi, "Metazoa":metazoa, "Viruses": virus, "Viridiplantae":viridiplantae}
            results[samplename] = pd.Series(resultsdic)
    
    resultsdf = pd.DataFrame(results)

    with open(args.output_file,'w') as outfile:
        resultsdf.to_csv(outfile, sep="\t")
    
    
if __name__ == '__main__':
    main()
