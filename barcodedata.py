import sys
import os
import glob
from Bio import SeqIO
import pandas


def constitantwithconsistantregion(consistantregion, expseq):
    if expseq==consistantregion:
        return 0
    else:
        errors = 0
        expseql = list(expseq)
        consistantregionl = list(consistantregion)
        for i, n in enumerate(expseql):
            n2 = consistantregionl[i]
            if n!=n2:
                errors = errors+1
        return errors

#"exp6":"TGAGCGTCTGAACTCCAGTCACGTCCTGCATCTGATCATCTCGTATGCCGTCTTCTGCTTG"

def parsetofile(fastqfile, outdirpluspre):
    reads = [rec for rec in SeqIO.parse(fastqfile, "fastq")]
    totallines = len(reads)
    readswithconstrantdf = outdirpluspre+".perfectdata.txt"
    readswithoutconstrantdf = outdirpluspre+".notperfectdata.txt"
    wf = open(readswithconstrantdf, "w")
    wf2 = open(readswithoutconstrantdf, "w")
    locexpl = [["exp1", 5, "ATC"], ["exp2", 10, "GAT"], ["exp3", 15, "AAA"],["exp4", 20, "GGT"],["exp5", 25,"ACCCAGCTTTCTTGTACAAAGTGGTTGATCGATGCGATGTACGGGCCAGATATACGCGTATCTGAGGGGACTAGGGTGTGTTTAGGCGAAAAGCGG"],["exp6", 129,"TGAGCGTCTGAACTCCAGTCAC"]]
    for readen in reads:
        read = readen.seq
        errors = False
        for exp in locexpl:
            name, pos, seq = exp
            stop = pos+len(seq)
            thisreadseq = read[pos:stop]
            if constitantwithconsistantregion(thisreadseq, seq)!=0:
                errors= True
        if errors==False:
            pcrbarcode = str(read[121:129])
            var1 = str(read[0:5])
            var2 = str(read[8:10])
            var3 = str(read[13:15])
            var4 = str(read[18:20])
            var5 = str(read[23:25])
            barcodeparts = [var1,var2,var3,var4,var5]
            biobarcode = "".join(barcodeparts)
            information = [str(readen.name), str(read), biobarcode, pcrbarcode]
            wf.write("\t".join(information)+"\n")
        else:
            information = [str(readen.name), str(read)]
            wf2.write("\t".join(information)+"\n")
    wf.close()
    wf2.close()    
    return readswithconstrantdf,readswithoutconstrantdf, totallines


def countsonperfectfastqreads(filename,outdirpluspre):
    pcrbarcodesfile, cellcountsfile = outdirpluspre+"pcrdups.txt", outdirpluspre+"cellcounts.txt"
    wf = open(pcrbarcodesfile, "w")
    df = pandas.read_csv(filename, sep="\t", names=["name", "seq", "biobarcode", "pcrbarcode"])
    totalperfectcontraintlines, colcount = df.shape
    g = df.groupby(['biobarcode', 'pcrbarcode'])
    line = [[name[0], name[1],len(group)] for name, group in g]
    print (line)
    for subline in line:
        wf.write("\t".join(map(str,subline))+"\n")
    wf.close()
    df2 = pandas.read_csv(pcrbarcodesfile, sep="\t", names=["biobarcode", "pcrbarcode", "countsofeachpcrdup"])
    totalpcrdupslines, colcount = df2.shape
    c = df2["biobarcode"].value_counts()
    totalcelllines = len(c)
    print (totalcelllines)
    c.to_csv(cellcountsfile, sep="\t")
    return totalperfectcontraintlines, totalperfectcontraintlines, totalpcrdupslines, pcrbarcodesfile, cellcountsfile, totalcelllines 
    


def main(indir, outdir):
	wfinto = open(outdir+"summary.csv", "w")
	labelline = []
	wfinto.write(",".join(labelline)+"\n")
	fastqfiles = sorted([filename for filename in glob.glob(os.path.join(indir, '*.fastq'))])
	for fastqfile in fastqfiles:
		print (fastqfile)
		fastqroot = fastqfile.split("/")[-1]
		outdirpluspre = outdir+fastqroot
		readswithconstrantdf, readswithoutconstrantdf, totallines = parsetofile(fastqfile, outdirpluspre)
		totalperfectcontraintlines, totalperfectcontraintlines, totalpcrdupslines, pcrbarcodesfile, cellcountsfile, totalcelllines = countsonperfectfastqreads(readswithconstrantdf, outdirpluspre)		
		line = [fastqfile, readswithconstrantdf, readswithoutconstrantdf, pcrbarcodesfile, cellcountsfile, totallines, totalperfectcontraintlines, totalpcrdupslines, totalcelllines]	
		wfinto.write(",".join(map(str,line))+"\n")
	wfinto.close()	

#def makegraphs(summaryfile):
#    df2 = pandas.read_csv(pcrbarcodesfile, sep="\t", names=["biobarcode", "pcrbarcode", "countsofeachpcrdup"])
#    countsforpcrdups= df2["countsofeachpcrdup"].value_counts()
#    df3 = pandas.read_csv(cellcountsfile, sep="\t", names=["biobarcode", "countafterremovepcrdups"])	
#    countsforcellbarcodes = df3[countafterremovepcrdups].value_counts()


if __name__=="__main__":
    #directorywithfastq = "/Users/allenma/olwin/" 
    #outdir = "/Users/allenma/olwin/" 
    directorywithfastq = sys.argv[1]
    outdir = sys.argv[2]
    main(directorywithfastq, outdir)

