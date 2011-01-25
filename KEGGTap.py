##################################################################
# @Program: KEGGTap.py                                           #
# @Version: 2                                                    #
# @Author: Chris Plaisier                                        #
# @License: Licensed under the Academic Free License version 2.1 #
# @Sponsored by:                                                 #
# Paivi Pajukanta Lab, UCLA                                      #
# University of California Department of Human Genetics          #
# David Geffen School of Medicine at UCLA                        #
# Gonda Neuroscience and Genetics Research Center                #
# 695 Charles E. Young Drive South, Box 708822                   #
# Los Angeles, CA 90095-7088, USA                                #
# (310) 794-9802                                                 #
# @Also Sponsored by:                                            #
# Genomics Analysis and Interpretation Training Grant            #
#                                                                #
# If this program is used in your analysis please mention who    #
# built it. Thanks. :-)                                          #
#                                                                #
# Copyrighted by Chris Plaisier  6/7/2005                        #
##################################################################

import sys
import os
import math
import cmath
import cPickle
import time
from optparse import OptionParser
from ftplib import FTP
from xml.dom import minidom, Node

print "###############################################################"
print "#                                                             #"
print "# KEGGTap 2.0 developed by Chris Plaisier (plaisier@ucla.edu) #"
print "# Licensed under the Academic Free License version 2.1        #"
print "#                                                             #"
print "###############################################################"

#### Option Parsing ####
usage = "%prog [options]"
parser = OptionParser(usage)
parser.add_option("-o", "--org", action="store", type="string", dest="org", default="NA", help="Three letter organism name from KEGG, http://www.genome.jp/kegg/catalog/org_list.html")
parser.add_option("-d", "--data", action="store", type="int", dest="data", default=2, help="How to obtain the KEGG data; 0 = Cached, 1 = Download XML Files, 2 = Get Directly From Kegg Database")
parser.add_option("-i", "--input", action="store", type="string", dest="input", default="input.csv", help="Name of input file, default is input.csv")
parser.add_option("-s", "--output", action="store", type="string", dest="output", default="KEGGTapOutput.csv", help="Name of output file, default is KEGGTapOutput.csv")
parser.add_option("-p", "--pvalue", action="store", type="float", dest="pvalue", default=0.05, help="p-value to use as threshold for calling a pathway over-represented")
(options, args) = parser.parse_args()

if options.org=="NA":
    parser.print_help()
    sys.exit(0)

#### Organism Name ####
# Need to get this from a file or from the command prompt sys.argv
orgName = options.org

#### Download Clean ####
# If this flag is set then go ahead and download the most recent data, otherwise if the files are downloaded then don't bother
# 2 = Yes, WSDL; 1 = Yes, XML; 0 = no
downClean = options.data

#### Input File ####
# Need to get this from a file or from the command prompt sys.argv
inputFile = options.input

#### Output File ####
# Need to get this from a file or from the command prompt sys.argv
outputFile = options.output

#### p-value Threshold ####
# Need to get this from a file or from the command prompt sys.argv
threshold = options.pvalue

def getWSDLConn():
    ## May not need SOAPpy after all for this program
    from SOAPpy import WSDL
    # Open a connection with the KEGG server
    wsdl = 'http://soap.genome.jp/KEGG.wsdl'
    serv = WSDL.Proxy(wsdl)
    return serv

def getHTMLPathway(pathway, objects, serv):
    queryMe = pathway
    doItAgain = True
    while doItAgain:
        try:
            doItAgain = False
            results = serv.get_html_of_marked_pathway_by_objects(pathway,objects)
            return results
        except:
            doItAgain = True
            time.sleep(1)

def getFTPText(ftp, filename, outfile=None):
    # fetch a text file
    if outfile is None:
        outfile = sys.stdout
    # use a lambda to add newlines to the lines read from the server
    ftp.retrlines("RETR " + filename, lambda s, w=outfile.write: w(s+"\n"))

def dlOrgFiles(orgName):
    ftp = FTP('ftp.genome.jp')
    ftp.login()
    ftp.cwd('pub/kegg/xml/current/'+orgName)
    files = ftp.nlst()

    # Create directory if it doesn't exist
    if os.path.exists(os.path.abspath(orgName)):
        dels=os.listdir(os.path.abspath(orgName))
        for i in dels:
            os.remove(os.path.abspath(orgName)+"\\"+str(i))
        os.removedirs(os.path.abspath(orgName))
    os.mkdir(orgName)
    
    # Selects only those files with xml in the name
    for x in files:
        if not str(x).find('xml')==-1:
            xmlFile = file(str(orgName)+'/'+str(x),'w')
            getFTPText(ftp, x, xmlFile)
            xmlFile.close()

# Good old fashioned factorial, uses a natural logarithm to allow for larger values to be used
def logFact(n):
    if n==0:
        return 0
    tot = 0
    for i in range(n):
        i += 1
        tot += cmath.log(i)
    return abs(tot)

# Binomial uses logarithm factorials
def lognCk(n,k):
    return logFact(n)-logFact(k)-logFact(n-k)

# N = number of objects in Urn, n = number of possible sucesses, N-n = number of other objects in Urn, m = sample size
# taken out at random, k = number of success in sample size
def hypergeo(k,m,n,N):
    a = lognCk(n,k)
    b = lognCk(N-n,m-k)
    c = lognCk(N,m)
    prob = abs(cmath.e**(a+b-c))
    return prob

# Take the current plus all those more extreme
def hypergeoPValue(k,m,n,N):
    pvalue = 0
    pvalue = hypergeo(k,m,n,N)
    if not k==n and not k==m:
        k = k+1
        while not k>n and not k>m:
            pvalue = pvalue + hypergeo(k,m,n,N)
            k = k+1
    return pvalue

def computeZScore(a1,a2,b1,b2):
    n1 = float(a1+a2)
    n2 = float(b1+b2)
    p1 = float(a1/n1)
    p2 = float(b1/n2)
    num = p1-p2
    denomA = ((p1*(1-p1))/n1)
    denomB = ((p2*(1-p2))/n2)
    denom = abs(cmath.sqrt(denomA+denomB))
    zscore = num/denom
    return zscore

def qsort(L):
    if L == []: return []
    return qsort([x for x in L[1:] if x< L[0]]) + L[0:1] + qsort([x for x in L[1:] if x>=L[0]])

def cutoff(perc, data):
    lenData = len(data)
    data = qsort(data)
    # multiply the length of the data by cutoff percentage
    co = math.floor(float(lenData) * (float(perc) / float(100)))
    # get connectivity for gene at number gotten above
    returnIndex = ((lenData - 1) - (co - 1))
    return data[int(returnIndex)]

# Query the KEGG database
def soapQuery(query):
    while doItAgain:
        try:
            doItAgain = False
            results = serv.get_linked_pathways(queryMe)
        except:
            doItAgain = True
            time.sleep(1)
        return results

# Get the number of genes in a pathway
def numGenes(pathway):
    result = soapQuery(pathway)
    return len(results)

## First thing to do is build a loci pathway dictionary
## So first need to download the XML files
lociPathDict = {}
if downClean==1:
    print "Downloading XML files, this is dependent on your network connection:"
    dlOrgFiles(orgName)
    print "downloading complete."
    ## Now compile the lociPathDict
    # Build a dictionary of the loci and the pathways
    # <key = locus> = <value = [path1, path2, ..., pathn]>
    print ""
    print "Parsing XML files:"
    pathway = ""
    listDir = os.listdir(orgName)
    for x in listDir:
        inFile = open(str(orgName)+"/"+str(x),"r")
        doc = minidom.parse(inFile)
        parent = doc.documentElement
        level = 0
        if parent.nodeName=="pathway":
            pathway = parent.attributes.get("name").value.encode()
        for node in parent.childNodes:
            if node.nodeType == Node.ELEMENT_NODE:
                # Write out the attributes.
                if node.nodeName == "entry":
                    attrs = node.attributes
                    attrType = attrs.get("type")
                    if attrType.value=="gene":
                        attrSep = (str(attrs.get("name").value)).split(" ")
                        for x in range(len(attrSep)):
                            gene = attrSep[x].encode()
                            if lociPathDict.has_key(gene):
                                tempMe = lociPathDict[gene]
                                tempMe.append(pathway)
                                lociPathDict[gene] = tempMe
                            else:
                                lociPathDict[gene] = [pathway]
    print "finished parsing XML files."
    ## Write out the data so don't have to do it everytime
    cPickle.dump(lociPathDict,open(str(orgName)+"LociPathDict.pkl","w"))
elif downClean==2:
    serv = getWSDLConn()
    # Load the Gene Names File. Example of gene name = "<org>:<geneName>"
    geneNamesFile = file(inputFile,"r")
    lines = geneNamesFile.readlines()
    numLines = len(lines)
    print "Starting run for all the pathways in dataset:"
    for i in range(numLines):
	queryMe = str(orgName)+":"+str((lines[i].split(","))[0])
        queryThis = [queryMe]
        doItAgain = True
        while doItAgain:
            try:
                doItAgain = False
                results = serv.get_pathways_by_genes(queryThis)
            except:
                doItAgain = True
                time.sleep(1)
	if i%10==0:
	    sys.stdout.write(str(i))
	elif ((str(results)).count("path"))!=0:
            locus = queryMe
            pathways = []
            for j in range(len(results)):
                pathways.append(str(results[j]))
            lociPathDict[locus] = pathways
	    sys.stdout.write("|")
	else:
	    sys.stdout.write(".")  	
    ## Write out the data so don't have to do it everytime
    cPickle.dump(lociPathDict,open(str(orgName)+"LociPathDict.pkl","w"))
    sys.stdout.write("\n")
    print "finished gathering information from database."
else:
    print "Using cached data source."
    lociPathFile = file(str(orgName)+"LociPathDict.pkl","r")
#    lociPathFile = file("lociPathDict.pkl","r")
    lociPathDict = cPickle.load(lociPathFile)

## Hypergeometric distribution for this application:
# N = the N is the number of genes in the network
# n = the number of genes in a specific pathway that are in the module
# m = the number of genes in a module
# k = the number of genes from the pathway in the module
# This probability from this will be the probability of getting a specific number of genes from a pathway given the
# total number of gene in the network and the size of the module along with how many genes from the pathway
# there are in the network
#
# Should use hypergeo(N, n, m, k)

# Therefore will need a tally for each dataset the nubmer of genes that are found for each pathway, and how many genes
# comprise the pathway

## Dictionaries for the data
# <key = moduleId> = <value = pathwayDict>
moduleDict = {}

# <key = locus> = <value = connectivity>
locusDict = {}

# <key = pathwayId> = <value = [locusDict, p-value, URL]>
pathwayDict = {}

# Input file
# <locus>, <module #>[, <connectivity>]
lociModFile = open(inputFile,"r")

## Now read in the inputFile and build the lociModDict
# can be modified to add the connectivity later
inputLine = lociModFile.readline()
N = 0 # Number of genes in the dataset
print "Annotating genes..."
while inputLine:
    inputSep = (inputLine.split("\n"))[0].split(",")
    module = inputSep[1]
    locus = str(orgName)+":"+str(inputSep[0])
    if lociPathDict.has_key(locus):	
        N += 1
	# Make Module portion
        if moduleDict.has_key(module):
           loci = moduleDict[module]
           loci[locus] = lociPathDict[locus]
        else:
	    loci = {}
	pathways = {}
	paths = lociPathDict[locus]
	for i in range(len(paths)):
	    pathways[paths[i]] = [666.666,"NA"]
   	loci[locus] = pathways
        # Add the module if it doesn't exist
	moduleDict[module] = loci
    inputLine = lociModFile.readline()
lociModFile.close()
print "KEGG annotates", N, "genes."

## Sum the number of genes for each pathway found in the dataset
# Calculates all the n's
sumPathwaysDict = {}
for i in moduleDict:
    lociDict = moduleDict[i]
    for j in lociDict:
	pathDict = lociDict[j]
	for k in pathDict:
	    if sumPathwaysDict.has_key(k):
		sumPathwaysDict[k] = sumPathwaysDict[k]+1
	    else:
		sumPathwaysDict[k] = 1

## Now for each module get a p-value for each pathway
print "Calculating p-Values..."
writeMe = "Module,Pathway,p-Value,k,m,n,N,KEGG Pathway URL,Genes\n"
serv = getWSDLConn()
for i in moduleDict:
    pathDone = {}
    lociDict = moduleDict[i]
    print "Starting Module #:", i
    for j in lociDict:
	pathDict = lociDict[j]
	for p in pathDict:
	    loci2HTML = []
	    if not pathDone.has_key(p):
		loci2HTML = [j]
	        k = 0 # Number of genes from a pathway in a module
		for l in lociDict:
		    if not l==j:
                        path2Dict = lociDict[l]
			if path2Dict.has_key(p):
			    k += 1
			    loci2HTML.append(l)
	        pValue = hypergeoPValue(k,len(lociDict),sumPathwaysDict[p],N)
		if pValue<threshold:
		    url = getHTMLPathway(p,loci2HTML,serv)
		else:
		    url = "NA"
		## Now need to set the pathway value really quick
		pathDict[p] = [pValue,url]
		pathDone[p] = [pValue,url]
		if pValue<threshold:
		    writeMe = writeMe + str(i) + "," + str(p) + "," + str(pValue) + "," + str(k) + "," + str(len(lociDict)) + "," + str(sumPathwaysDict[p]) + "," + str(N) + "," + url + ","
		    for bitty in range(len(loci2HTML)):
		        if not bitty==0:
			    writeMe = writeMe + "|"
			writeMe = writeMe + loci2HTML[bitty]
		    writeMe = writeMe + "\n"
	    else:
		pathDict[p] = pathDone[p]
	lociDict[j] = pathDict
    moduleDict[i] = lociDict
print "finished calculating p-Values."

## Write out a csv file
outFile = open(outputFile,"w")
outFile.write(writeMe)
outFile.close()
