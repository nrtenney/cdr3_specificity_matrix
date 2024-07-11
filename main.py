import sys

peptideListFile = "input/CDR3/" + sys.argv[1]
proteaseFile = "input/protease/" + sys.argv[2]
outputFile = "output/" + sys.argv[3]

#correlates the 3 letter amino acid abreviations from specificity matrix
#to one letter abbreviations from CDR3 file
def getAminoAcidAcronyms(input):
    aminoAcids = {}
    with open(input) as inFile:
        lines = inFile.readlines()
        for line in lines:
            line = line.strip()
            aminoAcid = line.split('\t')
            aminoAcids[aminoAcid[1]] = aminoAcid[0]
    return aminoAcids

#gets specificity matrix for certain protease from input file
def getProteaseTable(input, aminoAcids):
    proteaseTable = {}
    with open(input) as inFile:
        lines = inFile.readlines()
        for line in lines:
            line = line.strip()
            aminoAcidScore = line.split(' \t')
            if len(aminoAcidScore[0]) == 3:
                for i in range(1, 9):
                    aminoAcidScore[i] = int(aminoAcidScore[i])
                proteaseTable[aminoAcids[aminoAcidScore[0]]] = aminoAcidScore[1:]
    return proteaseTable

#calculate max possible score of specificity matrix
def getMaxScore(proteaseTable):
    keys = proteaseTable.keys()
    score = 0
    for i in range(8):
        maxPosition = 0
        for key in keys:
            value = int(proteaseTable[key][i])
            if value > maxPosition:
                maxPosition = value
        score += maxPosition
    return score

#input of cdr3 list from input file
def getPeptideList(peptideListFile):
    peptideList = []
    with open(peptideListFile) as inFile:
        lines = inFile.readlines()
        for line in lines[1:]:
            line = line.strip().split('\t')
            id = line[0]
            peptide = line[1]
            if len(peptide) > 7:
                peptideList.append('\t'.join([id, peptide]))
    peptideList = list(set(peptideList))
    noDuplicatesList = []
    for peptide in peptideList:
        noDuplicatesList.append(peptide.split('\t'))
    return noDuplicatesList

#calculates score of each length of 8 peptides from cdr3
def getScore(proteaseTable, peptide):
    score = 0
    for i in range(len(peptide)):
        score += proteaseTable[peptide[i]][i]
    return score

#calcultes the ratio of maxium score from cdr3 to the max possible score from the protease
def getMaxRatio(protein, proteaseTable, maxProteaseScore):
    maxProteinScore = 0
    for i in range(len(protein) - 7):
        seq = protein[i:i+8]
        score = getScore(proteaseTable, seq)
        if score > maxProteinScore:
            maxProteinScore = score
    ratio = maxProteinScore / maxProteaseScore
    return ratio

def printPeptideList(peptideList, outputFile):
    with open(outputFile, 'w') as outFile:
        outFile.write("id\tcdr3\tRatio\n")
        for peptide in peptideList:
            outFile.write(peptide[0] + '\t' + peptide[1] + '\t' + str(peptide[2]) + '\n')

aminoAcids = getAminoAcidAcronyms("input/aminoAcids.txt")
proteaseTable = getProteaseTable(proteaseFile, aminoAcids)
maxProteaseScore = getMaxScore(proteaseTable)
peptideList = getPeptideList(peptideListFile)
for id in peptideList:
    id.append(getMaxRatio(id[1], proteaseTable, maxProteaseScore))
peptideList = sorted(peptideList, key= lambda x:x[2])
printPeptideList(peptideList, outputFile)