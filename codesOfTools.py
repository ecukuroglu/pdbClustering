#################
## Created by Engin Cukuroglu
#################

import os

def interfaceNameSorter(interfaceName):
    splittedInterface = interfaceName.split('_')
    pdbName = splittedInterface[0].upper()
    chain_1_sortedList = sorted(splittedInterface[1])
    chain_1 = ''
    for chain_1_item in chain_1_sortedList:
        chain_1 = '%s%s' %(chain_1, chain_1_item)
    chain_2_sortedList = sorted(splittedInterface[2])
    chain_2 = ''
    for chain_2_item in chain_2_sortedList:
        chain_2 = '%s%s' %(chain_2, chain_2_item)
    allChains_sortedList = sorted('%s%s' %(chain_1, chain_2))
    allChains = ''
    for allChains_item in allChains_sortedList:
        allChains = '%s%s' %(allChains, allChains_item)
    commonChainCounter = 0
    for i in list(range(len(chain_1))):
        for j in list(range(len(chain_2))):
            if chain_1[i] == chain_2[j]:
                commonChainCounter = commonChainCounter + 1
    if chain_1[0] > chain_2[0]:
        tempChain = chain_1
        chain_1 = chain_2
        chain_2 = tempChain
    interfaceResidueStatusOfChains = splittedInterface[3]
    interface = '%s_%s_%s_%s' %(pdbName, chain_1, chain_2, interfaceResidueStatusOfChains)
    return pdbName, chain_1, chain_2, allChains, interfaceResidueStatusOfChains, interface, commonChainCounter

def getPDB(pdbID, pdbDir):
    pdb = pdbID.lower()
    try:
        os.system("rsync -rlpt -v -z --delete --port=33444 rsync.wwpdb.org::ftp_data/structures/divided/pdb/%s/pdb%s.ent.gz %s/ >/dev/null 2>&1" %(pdb[1:3], pdb, pdbDir))
        tempDownloadedPDBDirectory = '%s/pdb%s.ent.gz' %(pdbDir, pdb)
        if os.path.exists(tempDownloadedPDBDirectory):
            os.system("gunzip %s/pdb%s.ent.gz" %(pdbDir, pdb))
            os.system("mv %s/pdb%s.ent %s/%s.pdb" %(pdbDir, pdb, pdbDir, pdbID))
        else:
            print('%s does not in the ftp server or an error occurred during download process.\n' %(pdbID))
    except:
        print('Exception handler: %s could not be downloaded.\n' %(pdbID))

def pdbChainListExtractor(pdbFileDirectory):
    chainDict = {}
    chainList = []
    pdbFile = open(pdbFileDirectory, 'r')
    while 1:
        pdbLine = pdbFile.readline()
        if pdbLine == '' or pdbLine[0:3] == 'END':
            break
        if pdbLine[0:4] == 'ATOM':
            chainID = pdbLine[21]
            if not chainID in chainDict:
                chainDict[chainID] = chainID
                chainList.append(chainID)
    pdbFile.close()
    return chainList

def three2One(residue):
    if residue == 'ALA': return 'A'
    elif residue == 'CYS': return 'C'
    elif residue == 'ASP': return 'D'
    elif residue == 'GLU': return 'E'
    elif residue == 'PHE': return 'F'
    elif residue == 'GLY': return 'G'
    elif residue == 'HIS': return 'H'
    elif residue == 'ILE': return 'I'
    elif residue == 'LYS': return 'K'
    elif residue == 'LEU': return 'L'
    elif residue == 'MET': return 'M'
    elif residue == 'ASN': return 'N'
    elif residue == 'PRO': return 'P'
    elif residue == 'GLN': return 'Q'
    elif residue == 'ARG': return 'R'
    elif residue == 'SER': return 'S'
    elif residue == 'THR': return 'T'
    elif residue == 'VAL': return 'V'
    elif residue == 'TRP': return 'W'
    elif residue == 'TYR': return 'Y'
    else: return 'X'

def naccessRSAFileReader(naccessResultDirectory, fileName):
    rsaFile = open('%s/%s.rsa' %(naccessResultDirectory, fileName), 'r')
    residueData = []
    while 1:
        line = rsaFile.readline()
        if line == '':
            break
        if line[0:3] == 'RES':
            tempList = []
            tempList.append(line[0:3].strip())
            tempList.append(line[4:7].strip())
            tempList.append(line[8].strip())
            tempList.append(line[9:13].strip())
            tempList.append(line[15:22].strip())
            tempList.append(line[22:28].strip())
            tempList.append(line[28:35].strip())
            tempList.append(line[35:41].strip())
            tempList.append(line[41:48].strip())
            tempList.append(line[48:54].strip())
            tempList.append(line[54:61].strip())
            tempList.append(line[61:67].strip())
            tempList.append(line[67:74].strip())
            tempList.append(line[74:80].strip())

            residueData.append(tempList)
            if len(tempList) != 14:
                print('Error in %s' %(fileName))
    rsaFile.close()
    return residueData

def naccessRSAFileReaderOnlyAbsASADictionary(naccessResultDirectory, fileName):
    rsaFile = open('%s/%s.rsa' %(naccessResultDirectory, fileName), 'r')
    asaDictionary = {}
    while 1:
        line = rsaFile.readline()
        if line == '':
            break
        if line[0:3] == 'RES':
            dictionaryKey = '%s%s%s' %(line[4:7].strip(), line[9:13].strip(), line[8])
            asaDictionary[dictionaryKey] = line[15:22].strip()
    rsaFile.close()
    return asaDictionary

def naccessRSAFileReaderReturnOnlyAbsASADictionary(naccessResultDirectory):
    asaDictionary = {}
    try:
        rsaFile = open('%s' %(naccessResultDirectory), 'r')
        while 1:
            line = rsaFile.readline()
            if line == '':
                break
            if line[0:3] == 'RES':
                dictionaryKey = '%s_%s_%s' %(line[4:7].strip(), line[9:13].strip(), line[8])
                asaDictionary[dictionaryKey] = line[15:22].strip()
        rsaFile.close()
    except:
        pass
    return asaDictionary

def distanceCalculator(coordinates1, coordinates2):
    x1 = coordinates1[0]
    y1 = coordinates1[1]
    z1 = coordinates1[2]
    x2 = coordinates2[0]
    y2 = coordinates2[1]
    z2 = coordinates2[2]
    return (((x2-x1)**2)+((y2-y1)**2)+((z2-z1)**2))**0.5

def dictionaryVDWRadii():
    VDWRadii = {'C':1.76, 'N':1.65, 'O':1.40, 'CA':1.87, 'H':1.20, 'S':1.85, 'CB':1.87, 'CZ':1.76, 'NZ':1.50, 'CD':1.81, 'CE':1.81, 'CG':1.81, 'C1':1.80, 'P':1.90}
    return VDWRadii

def dictionaryVDWRadiiExtended():
    VDWRadiiExtended = { 'ALA': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'OXT': 1.40 },\
                         'ARG': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.87, 'CD': 1.87, 'NE': 1.65, 'CZ': 1.76, 'NH1': 1.65, 'NH2': 1.65, 'OXT': 1.40 },\
                         'ASP': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.76, 'OD1': 1.40, 'OD2': 1.40, 'OXT': 1.40  },\
                         'ASN': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.76, 'OD1': 1.40, 'ND2': 1.65, 'OXT': 1.40 },\
                         'CYS': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'SG': 1.85, 'OXT': 1.40 },\
                         'GLU': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.87, 'CD': 1.76, 'OE1': 1.40, 'OE2': 1.40 , 'OXT': 1.40 },\
                         'GLN': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.87, 'CD': 1.76, 'OE1': 1.40, 'NE2': 1.65, 'OXT': 1.40  },\
                         'GLY': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'OXT': 1.40 },\
                         'HIS': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.76, 'ND1': 1.65, 'CD2': 1.76, 'CE1': 1.76, 'NE2': 1.65, 'OXT': 1.40  },\
                         'ILE': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG1': 1.87, 'CG2': 1.87, 'CD1': 1.87 , 'OXT': 1.40 },\
                         'LEU': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.87, 'CD1': 1.87, 'CD2': 1.87, 'OXT': 1.40 },\
                         'LYS': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.87, 'CD': 1.87, 'CE': 1.87, 'NZ': 1.50 , 'OXT': 1.40 },\
                         'MET': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.87, 'SD': 1.85, 'CE': 1.87, 'OXT': 1.40 },\
                         'PHE': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.76, 'CD1': 1.76, 'CD2': 1.76, 'CE1': 1.76, 'CE2': 1.76, 'CZ': 1.76, 'OXT': 1.40  },\
                         'PRO': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.87, 'CD': 1.87, 'OXT': 1.40  },\
                         'SER': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'OG': 1.40, 'OXT': 1.40 },\
                         'THR': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'OG1': 1.40, 'CG2': 1.87, 'OXT': 1.40 },\
                         'TRP': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.76, 'CD1': 1.76, 'CD2': 1.76, 'NE1': 1.65, 'CE2': 1.76, 'CE3': 1.76, 'CZ2': 1.76, 'CZ3': 1.76, 'CH2': 1.76, 'OXT': 1.40  },\
                         'TYR': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.76, 'CD1': 1.76, 'CD2': 1.76, 'CE1': 1.76, 'CE2': 1.76, 'CZ': 1.76, 'OH': 1.40, 'OXT': 1.40 },\
                         'VAL': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG1': 1.87, 'CG2': 1.87, 'OXT': 1.40  },\
                         'ASX': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.76, 'AD1': 1.50, 'AD2': 1.50 },\
                         'GLX': { 'N': 1.65, 'CA': 1.87, 'C': 1.76, 'O': 1.40, 'CB': 1.87, 'CG': 1.76, 'CD': 1.87, 'AE1': 1.50, 'AE2': 1.50, 'OXT': 1.40  },\
                         'ACE': { 'C': 1.76, 'O': 1.40, 'CA': 1.87 },\
                         'PCA': dictionaryVDWRadii() }
    return VDWRadiiExtended

def readPDBReturnAtomProperties(PDBDirectory, chainList):
    pdbFile = open(PDBDirectory,'r')
    atomPropertyList = []
    chainType = 'Protein'
    coordinateDictionaryCA = {}
    VDWRadiiExtended = dictionaryVDWRadiiExtended()
    dummyResTypeDict = {}
    while 1:
        pdbLine = pdbFile.readline()
        if pdbLine == '' or pdbLine[0:3] == 'END':
            break
        if pdbLine[0:4] == 'ATOM':
            if pdbLine[21] in chainList:
                resType = pdbLine[17:20].strip()
                if not chainType == 'Unknown':
                    if resType == 'DG' or resType == 'DT' or resType == 'DA' or resType == 'DC':
                        chainType = 'DNA'
                    elif resType == 'G' or resType == 'U' or resType == 'A' or resType == 'C':
                        chainType = 'RNA'
                    elif three2One(resType) == 'X':
                        chainType = 'Unknown'
                        break
                    else:
                        atomType = pdbLine[12:16].strip()
                        if atomType in VDWRadiiExtended[resType]:
                            VDWDistance = VDWRadiiExtended[resType][atomType]
                            resNo = pdbLine[22:26].strip()
                            coordinates = [float(pdbLine[30:38].strip()), float(pdbLine[38:46].strip()), float(pdbLine[46:54].strip())]
                            atomPropertyList.append([pdbLine[21], resNo, resType, atomType, coordinates, VDWDistance])
                            if pdbLine[13:15] == 'CA':
                                resKey = '%s%s%s' %(resType, resNo, pdbLine[21])
                                dummyResKey = '%s_%s' %(resNo, pdbLine[21])
                                if not dummyResKey in dummyResTypeDict:
                                    coordinateDictionaryCA[resKey] = coordinates
                                    dummyResTypeDict[dummyResKey] = 1
    pdbFile.close()
    return atomPropertyList, chainType, coordinateDictionaryCA
                
def readPDBReturnAtomProperties_withMultipleChains(PDBDirectory, chainDict):
    chainTypeDict = {}
    for chainID in chainDict:
        chainTypeDict[chainID] = {}
        chainTypeDict[chainID]['type'] = 'Protein'
        chainTypeDict[chainID]['residueNumber'] = 0
    atomPropertyList = []
    coordinateDictionaryCA = {}
    VDWRadiiExtended = dictionaryVDWRadiiExtended()
    pdbResidueDict = {}
    pdbAtomDict = {}
    if os.path.exists(PDBDirectory):
        pdbFile = open(PDBDirectory,'r')
        while 1:
            pdbLine = pdbFile.readline()
            if pdbLine == '' or pdbLine[0:3] == 'END':
                break
            if pdbLine[0:4] == 'ATOM':
                chainID = pdbLine[21]
                if chainID in chainDict:
                    resType = pdbLine[17:20].strip()
                    if resType == 'DG' or resType == 'DT' or resType == 'DA' or resType == 'DC':
                        chainType = 'DNA'
                        if chainTypeDict[chainID]['type'] == 'Protein':
                            chainTypeDict[chainID]['type'] = chainType
                        elif not chainTypeDict[chainID]['type'] == chainType:
                            chainTypeDict[chainID]['type'] = 'Conflict'
                    elif resType == 'G' or resType == 'U' or resType == 'A' or resType == 'C':
                        chainType = 'RNA'
                        if chainTypeDict[chainID]['type'] == 'Protein':
                            chainTypeDict[chainID]['type'] = chainType
                        elif not chainTypeDict[chainID]['type'] == chainType:
                            chainTypeDict[chainID]['type'] = 'Conflict'
                    elif three2One(resType) == 'X':
                        chainType = 'Unknown'
                        if chainTypeDict[chainID]['type'] == 'Protein':
                            chainTypeDict[chainID]['type'] = chainType
                        elif not chainTypeDict[chainID]['type'] == chainType:
                            chainTypeDict[chainID]['type'] = 'Conflict'
                    else:
                        atomType = pdbLine[12:16].strip()
                        if atomType in VDWRadiiExtended[resType]:
                            VDWDistance = VDWRadiiExtended[resType][atomType]
                            resNo = pdbLine[22:26].strip()
                            resICode = pdbLine[26]
                            atomAlternativeLocationIndicator = pdbLine[16]
                            resDictKey = '%s_%s' %(resNo, chainID)
                            resDictValue = '%s_%s_%s' %(resType, resNo, chainID)
                            atomDictKey = '%s_%s_%s_%s' %(resType, resNo, chainID, atomType)
                            if not resDictKey in pdbResidueDict:
                                pdbResidueDict[resDictKey] = [resDictValue, resICode, atomAlternativeLocationIndicator]
                                pdbAtomDict[atomDictKey] = 1
                                coordinates = [float(pdbLine[30:38].strip()), float(pdbLine[38:46].strip()), float(pdbLine[46:54].strip())]
                                tempFactor = pdbLine[60:66].strip()
                                atomPropertyList.append([pdbLine[21], resNo, resType, atomType, coordinates, VDWDistance, tempFactor])
                                if atomType == 'CA':
                                    coordinateDictionaryCA[resDictValue] = coordinates
                                    chainTypeDict[chainID]['residueNumber'] = chainTypeDict[chainID]['residueNumber'] + 1 
                            else:
                                if pdbResidueDict[resDictKey][0] == resDictValue:
                                    if pdbResidueDict[resDictKey][1] == resICode:
                                        if not atomDictKey in pdbAtomDict:
                                            pdbAtomDict[atomDictKey] = 1
                                            coordinates = [float(pdbLine[30:38].strip()), float(pdbLine[38:46].strip()), float(pdbLine[46:54].strip())]
                                            tempFactor = pdbLine[60:66].strip()
                                            atomPropertyList.append([pdbLine[21], resNo, resType, atomType, coordinates, VDWDistance, tempFactor])
                                            if atomType == 'CA':
                                                coordinateDictionaryCA[resDictValue] = coordinates
                                                chainTypeDict[chainID]['residueNumber'] = chainTypeDict[chainID]['residueNumber'] + 1 
                                            
        pdbFile.close()
    else:
        print('%s does not present' %(PDBDirectory))
    return atomPropertyList, chainTypeDict, coordinateDictionaryCA
                
def vdwInterfaceDictionaryCreator(interfaceFileDirectory):
    interfaceDict = {}
    try:
        interfaceFile = open(interfaceFileDirectory, 'r')
        for line in interfaceFile:
            line = line.strip()
            if line[0:4] == 'ATOM':
                chainID = line[21]
                resType = line[17:20].strip()
                resNo = line[22:26].strip()
                resKey = '%s_%s_%s' %(resType, resNo, chainID)
                if not resKey in interfaceDict:
                    interfaceDict[resKey] = 1
    except:
        pass
    return interfaceDict 

def popsResultFileReader(popsResultFileDirectory):
    sasaDict = {}
    if os.path.exists(popsResultFileDirectory):
        popsResultFile = open(popsResultFileDirectory, 'r')
        while 1:
            line = popsResultFile.readline()
            if line == '':
                break
            splittedLine = line.strip().split()
            if len(splittedLine) > 5:
                resKey = '%s_%s_%s' %(splittedLine[0], splittedLine[2], splittedLine[1])
                sasaDict[resKey] = splittedLine
        popsResultFile.close()
    else:
        print('%s does not exist' %(popsResultFileDirectory))
    return sasaDict

def one2Three(residue):
        if residue == 'A': return 'ALA'
        elif residue == 'C': return 'CYS'
        elif residue == 'D': return 'ASP'
        elif residue == 'E': return 'GLU'
        elif residue == 'F': return 'PHE'
        elif residue == 'G': return 'GLY'
        elif residue == 'H': return 'HIS'
        elif residue == 'I': return 'ILE'
        elif residue == 'K': return 'LYS'
        elif residue == 'L': return 'LEU'
        elif residue == 'M': return 'MET'
        elif residue == 'N': return 'ASN'
        elif residue == 'P': return 'PRO'
        elif residue == 'Q': return 'GLN'
        elif residue == 'R': return 'ARG'
        elif residue == 'S': return 'SER'
        elif residue == 'T': return 'THR'
        elif residue == 'V': return 'VAL'
        elif residue == 'W': return 'TRP'
        elif residue == 'Y': return 'TYR'
        else: return 'XXX'

def multiprotSolutionReader(multiprotSolutionFileDirectory, solutionNo):
    startParsingSolutionFileIdentifier = 'Solution Num : %d' %(solutionNo)
    endParsingSolutionFileIdentifier = 'Solution Num : %d' %(solutionNo + 1)
    readingStatus = 2
    alignmentDict = {}
    referenceMonomer = ''
    alignedMonomer = ''
    rmsd = -1.0
    transMatrix = ''
    multiprotSolutionFile = open(multiprotSolutionFileDirectory, 'r')
    moleculeNameDict = {}
    while readingStatus > 0:
        solutionLine = multiprotSolutionFile.readline()
        if solutionLine == '':
            break
        solutionLine = solutionLine.strip()
        if solutionLine == startParsingSolutionFileIdentifier:
            readingStatus = 1
        elif solutionLine == endParsingSolutionFileIdentifier:
            readingStatus = 0
        elif solutionLine[0:4] == 'Mol-':
            splittedSolutionLine = solutionLine.split(':')
            moleculeNo = int(splittedSolutionLine[0].replace('Mol-','').strip())
            moleculePDBName = splittedSolutionLine[1].split('/')[-1]
            moleculeName = moleculePDBName.replace('.pdb', '')
            moleculeNameDict[moleculeNo] = moleculeName
        if readingStatus == 1:
            if solutionLine[0:9] == 'Reference':
                referenceMonomerIndex = int(solutionLine.strip().split()[-1])
                referenceMonomer = moleculeNameDict[referenceMonomerIndex]
            elif solutionLine[0:10] == 'Molecule :':
                alignedMonomerIndex = int(solutionLine.strip().split()[-1])
                alignedMonomer = moleculeNameDict[alignedMonomerIndex]
            elif solutionLine[0:7] == 'Trans :':
                transMatrix = solutionLine[7:].strip()
                    elif solutionLine[0:6] == 'RMSD :':
                            rmsd = float(solutionLine.strip().split()[-1])
            elif solutionLine[0:10] == 'Match List':
                while 1:
                    solutionLine = multiprotSolutionFile.readline()
                    if solutionLine[0:3] == 'End' or solutionLine == '':
                        break
                    tempResKeyList = solutionLine.strip().split()
                    referenceMoleculeItems = tempResKeyList[0].split('.')
                    alignedMoleculeItems = tempResKeyList[1].split('.')
                    referenceDictKey = '%s_%s_%s' %(one2Three(referenceMoleculeItems[1]), referenceMoleculeItems[2], referenceMoleculeItems[0])
                    alignedDictKey = '%s_%s_%s' %(one2Three(alignedMoleculeItems[1]), alignedMoleculeItems[2], alignedMoleculeItems[0])
                    alignmentDict[referenceDictKey] = alignedDictKey
                if solutionLine[0:3] == 'End' or solutionLine == '':
                    readingStatus = 0
    multiprotSolutionFile.close()
    return alignmentDict, referenceMonomer, alignedMonomer, rmsd, transMatrix

