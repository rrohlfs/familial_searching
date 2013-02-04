#this script should read in the Y_STR_ABI file and print out files for each pop group with alleles as plain ints

#set constants
nloci = 16;
datafile = "Y_STR_ABI.txt";
popnames = ["African_American", "Asian", "Asian_Indian", "Caucasian", "Chinese", "Filipino", "Hispanic", "Japanese", "Malay", "Native_American", "Sub-sharan_African", "Thai", "Vietnamese"];

#set up refHaps data struct
refHaps = {}
for i in range(len(popnames)):
    refHaps[popnames[i]] = []

#read in ref haps, saving them appropriately
infile = open(datafile, 'r');
headerline = infile.readline();
headerline = infile.readline();

alleles = [[] for row in range(nloci)]
for rawline in infile:
    rawhap = str.split(rawline)
    curpop = rawhap[0]
    curhap = rawhap[1:(len(rawhap)+1)]
    
    #track what alleles have been observed
    for i in range(nloci):
        newallele = 1
        for j in range(len(alleles[i])):
            if alleles[i][j] == curhap[i]:
                newallele = 0
	if newallele == 1:
            alleles[i].append(curhap[i])
	
    #save hap in right place
    refHaps[curpop].append(curhap)
    
    

infile.close()


#recode alleles as ints
cRefHaps = {}
for i in range(len(popnames)):
    cRefHaps[popnames[i]] = refHaps[popnames[i]]

for i in range(nloci):
    curalleles = sorted(alleles[i])
    curcode = {}
    for j in range(len(curalleles)):
        if curalleles[j] == '0':
            curcode[curalleles[j]] = -1
        else:
            curcode[curalleles[j]] = j
    
    for j in range(len(popnames)):
        for k in range(len(refHaps[popnames[j]])):
            cRefHaps[popnames[j]][k][i] = curcode[refHaps[popnames[j]][k][i]]


#remove profiles with missing data
for i in range(len(popnames)):
    indivs = range(len(refHaps[popnames[i]]))
    indivs.reverse()
    for j in indivs:
        nomissing = 1
        for k in range(nloci):
            if cRefHaps[popnames[i]][j][k] == -1:
                nomissing = 0
        
        if nomissing == 0:
            cRefHaps[popnames[i]].pop(j)


#write refHaps to files
for i in range(len(popnames)):
    outfile = open(''.join(["Y_STR_ABI_clean/",popnames[i],".dat"]), 'w')
    for j in range(len(cRefHaps[popnames[i]])):
        for k in range(nloci):
            outfile.write("%d " % cRefHaps[popnames[i]][j][k])
        outfile.write("\n")
    outfile.close()

