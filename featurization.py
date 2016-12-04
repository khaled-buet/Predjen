import csv
import os
import itertools
from itertools import product

def findIndexes(s, ch):
    return [i for i, letter in enumerate(s) if letter == ch]

def findAllIndexes(string, substr):
    ind = []
    index = 0;
    while(index<len(string)):
        index = string.find(substr, index)
        if index == -1:
            break
        ind.append(index)
        index = index + len(substr)
    return ind

def isAllExist(a, b):
    return all(x in a for x in b)
        
def loadViennaData(filename):
    filep = open(filename)
    csv_reader = csv.reader(filep)
    sgrna = []
    alldata = []
    for row in csv_reader:
        #print row[2]
        sgrna.append(row[1])
        alldata.append(row)
    #sgrna = sgrna[-1]
    sgrna = sgrna[1:]
    #print(sgrna)
    txtfile = open('rnaseq.txt', 'w')
    #writer = csv.writer(txtfile);
    #print(sgrna)
    for row in sgrna:
        txtfile.write("%s\n" % row)
    txtfile.close()
    os.system("RNAfold < rnaseq.txt > output_mfe.txt")
    os.system("RNAheat --Tmin=50 --Tmax=50 < rnaseq.txt > output_heat.txt")
    os.system("rm rnaseq.txt")
    os.system("rm rna.ps")
    
    outfile = open('output_mfe.txt')
    mfestr = []
    for line in outfile:
        mfestr.append(line)
    mfestr = ''.join(mfestr[1:len(mfestr):2])
    outfile.close()
    os.system('rm output_mfe.txt')
    #print(mfestr)
    mfesplit = mfestr.splitlines()
    #print(mfesplit)
    mfenumlist = []
    #mfenumlist.append('MFE')
    for i in range(0, len(mfesplit)):
        strnum = (''.join(mfesplit[i])).split(' ')
        mfenumlist.append(''.join(strnum[len(strnum)-1]).strip('()'))
        #print(mfenumlist)
    mfenumlist = list(map(float, mfenumlist))
    #print(mfenumlist)
    #print(len(alldata[0]))
    alldata[0].insert(len(alldata[0]), 'MFE')
    for i in range(1, len(mfenumlist) + 1):
        alldata[i].insert(len(alldata[i]), mfenumlist[i - 1])
    
    outfile2 = open('output_heat.txt')
    heat = []
    for line in outfile2:
        heat.append(line)
    #heat = ''.join(heat[1:len(heat):2])
    outfile2.close()
    #os.system('rm output_heat.txt')
    #print(heat)
    heatlist = []
    for i in range(0, len(heat)):
        tmp = heat[i].splitlines()
        heatlist.append(tmp[0].rstrip('\n').split(' ')[3])
    heatlist = list(map(float, heatlist))
    #print(heatlist)
        
    alldata[0].insert(len(alldata[0]), 'Heat')
    for i in range(1, len(heatlist) + 1):
        alldata[i].insert(len(alldata[i]), heatlist[i - 1])
    
    file1 = open(filename, 'wb')
    csv_writer = csv.writer(file1)
    csv_writer.writerows(alldata)
    

def checkextrafeatures(sgrnas, feature):
    dlist = []
    if(feature == 'A'):
        for sgrna in sgrnas:
            dlist.append(len(findIndexes(sgrna,'A')))
    elif(feature == 'B'):
        for sgrna in sgrnas:
            dlist.append(len(findIndexes(sgrna,'B')))
    elif(feature == 'C'):
        for sgrna in sgrnas:
            dlist.append(len(findIndexes(sgrna,'C')))
    elif(feature == 'D'):
        for sgrna in sgrnas:
            dlist.append(len(findIndexes(sgrna,'D')))
    elif(feature == 'E'):
        for sgrna in sgrnas:
            dlist.append(len(findIndexes(sgrna,'E')))
    elif(feature == 'F'):
        for sgrna in sgrnas:
            dlist.append(len(findIndexes(sgrna,'F')))
    elif(feature == 'G'):
        for sgrna in sgrnas:
            dlist.append(len(findIndexes(sgrna,'G')))
    elif(feature == 'H'):
        for sgrna in sgrnas:
            dlist.append(len(findIndexes(sgrna,'H')))
    elif(feature == 'I'):
        for sgrna in sgrnas:
            dlist.append(len(findIndexes(sgrna,'I')))
    elif(feature == 'J'):
        for sgrna in sgrnas:
            dlist.append(len(findIndexes(sgrna,'J')))
    elif(feature == 'K'):
        for sgrna in sgrnas:
            dlist.append(len(findIndexes(sgrna,'K')))
    elif(feature == 'L'):
        for sgrna in sgrnas:
            dlist.append(len(findIndexes(sgrna,'L')))
    elif(feature == 'M'):
        for sgrna in sgrnas:
            dlist.append(len(findIndexes(sgrna,'M')))
    elif(feature == 'N'):
        for sgrna in sgrnas:
            dlist.append(len(findIndexes(sgrna,'N')))
    elif(feature == 'O'):
        for sgrna in sgrnas:
            dlist.append(len(findIndexes(sgrna,'O')))
    elif(feature == 'P'):
        for sgrna in sgrnas:
            dlist.append(len(findIndexes(sgrna,'P')))
    elif(feature == 'Q'):
        for sgrna in sgrnas:
            dlist.append(len(findIndexes(sgrna,'Q')))
    elif(feature == 'R'):
        for sgrna in sgrnas:
            dlist.append(len(findIndexes(sgrna,'R')))
    elif(feature == 'S'):
        for sgrna in sgrnas:
            dlist.append(len(findIndexes(sgrna,'S')))
    elif(feature == 'T'):
        for sgrna in sgrnas:
            dlist.append(len(findIndexes(sgrna,'T')))
    elif(feature == 'U'):
        for sgrna in sgrnas:
            dlist.append(len(findIndexes(sgrna,'U')))
    elif(feature == 'V'):
        for sgrna in sgrnas:
            dlist.append(len(findIndexes(sgrna,'V')))
    elif(feature == 'W'):
        for sgrna in sgrnas:
            dlist.append(len(findIndexes(sgrna,'W')))
    elif(feature == 'Y'):
        for sgrna in sgrnas:
            dlist.append(len(findIndexes(sgrna,'Y')))
    elif(feature == 'Z'):
        for sgrna in sgrnas:
            dlist.append(len(findIndexes(sgrna,'Z')))
    
    return dlist


def appendtoCSV(dname, datalist, filename):
    filep = open(filename)
    csv_reader = csv.reader(filep)
    alldata = []
    for row in csv_reader:
        alldata.append(row)
    print(alldata)
    alldata[0].insert(len(alldata[0]), dname)
    for i in range(1, len(datalist)):
        alldata[i].insert(len(alldata[i]), datalist[i - 1])
    file1 = open(filename, 'w')
    csv_writer = csv.writer(file1)
    csv_writer.writerows(alldata)
    filep.close()
    file1.close()

    
def n_order(sgrnas, substr, pos):
    dlist = []
    for sgrna in sgrnas:
        #print(sgrnas)
        if(pos in findAllIndexes(sgrna, substr)):
            dlist.append(1)
        else:
            dlist.append(0)
    return dlist

if __name__ == "__main__":
    filename = "vaxijensequence.csv"
    #extrafeaturelist = ["GC_Count","AT_Count","A_Count","C_Count","G_Count","T_Count"]
    totalfeature = 0;
    filep = open(filename)
    csv_reader = csv.reader(filep)
    proteinseq = []
    for row in csv_reader:
        proteinseq.append(row[3])
    aminoacid = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "Y", "Z"]
    firstorderfeature = 0;
    #print(proteinseq)
    #for f in proteinseq:
        #print(f)
    for feature in aminoacid:
        dat = checkextrafeatures(proteinseq, feature)
        dat[0] = feature
        print(dat)
        #appendtoCSV(feature, dat, filename)
        break
##    ########first order##########

##    for acid in aminoacid:
##        for i in range(0,8):
##            feature = acid + "_" + str(i+1)
##            dat = n_order(protein, acid, i)
##            appendtoCSV(feature, dat, filename)
##            firstorderfeature = firstorderfeature + 1
##    print("First Order Features:",firstorderfeature)
##    for acid in aminoacid:
##        for i in range(0,30):
##            feature = nucleotide + "_" + str(i+1)
##            dat = n_order(sgrnas, nucleotide, i)
##            appendtoCSV(feature, dat, filename)
##            firstorderfeature = firstorderfeature + 1
##    print("First Order Features:",firstorderfeature)
##    
##    ########second order#########
##    secondorderfeature = 0;
##    secondorder = keywords = [''.join(i) for i in product(nucleotides, repeat = 2)]
##    for nucleotide in secondorder:
##        for i in range(0,29):
##            feature = nucleotide + "_" + str(i+1)
##            dat = n_order(sgrnas, nucleotide, i)
##            appendtoCSV(feature, dat, filename)
##            secondorderfeature = secondorderfeature + 1
##    print("Second Order Features:",secondorderfeature)
##            
##################third order######################
##    thirdorderfeature = 0;
##    thirdorder = keywords = [''.join(i) for i in product(nucleotides, repeat = 3)]
##    for nucleotide in thirdorder:
##        for i in range(0,28):
##            feature = nucleotide + "_" + str(i+1)
##            dat = n_order(sgrnas, nucleotide, i)
##            appendtoCSV(feature, dat, filename)
##            thirdorderfeature = thirdorderfeature + 1
##    print("Third Order Features:",thirdorderfeature)
##            
#################fourth order######################
##    fourthorderfeature = 0;        
##    fourthorder = keywords = [''.join(i) for i in product(nucleotides, repeat = 4)]
##    for nucleotide in fourthorder:
##        for i in range(0,27):
##            feature = nucleotide + "_" + str(i+1)
##            dat = n_order(sgrnas, nucleotide, i)
##            appendtoCSV(feature, dat, filename)
##            fourthorderfeature = fourthorderfeature + 1
##    print("Fourth Order Features:",fourthorderfeature)
##
##

##
##    loadViennaData(filename)
##
##    #totalfeatures = firstorderfeature + secondorderfeature + thirdorderfeature + fourthorderfeature + 8
##    #print('Total Features:',totalfeatures)
##
