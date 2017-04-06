#################
## Created by Engin Cukuroglu
#################
##
##usage: nohup python pythonCodes/proteinMonomerInfoExtractor.py & 
##
#################

import os, sys, time

t1 = time.time()
pdbChainListFileDirectory = 'fullLists/fullPDBChainList.txt'
#pdbChainListFileDirectory = 'fullLists/denemefullPDBChainList.txt'
pdbFileDirectory = 'pdbFiles'
pdbChainSpecificationFileDirectory = 'fullLists/fullPDBChainSpecifications.txt'

pdbChainListFile = open(pdbChainListFileDirectory, 'r')
pdbChainSpecificationFile = open(pdbChainSpecificationFileDirectory, 'w')

for infoLine in pdbChainListFile:
    splittedInfoLine = infoLine.strip().split('\t')
    pdbID = splittedInfoLine[0]
    moleculeSpecificationDict = {}
    moleculeSpecificationDict['resolution'] = ''
    moleculeSpecificationDict['experimentalType'] = ''
    moleculeSpecificationDict['temperature'] = ''
    moleculeSpecificationDict['PH'] = ''
    chainMoleculeIDMappingDict = {}
    pdbFile = open('%s/%s.pdb' %(pdbFileDirectory, pdbID.upper()), 'r')
    
    while 1:
        pdbLine = pdbFile.readline()
        if pdbLine == '' or pdbLine[0:3] == 'END':
            break
        if pdbLine[0:4] == 'ATOM':
            break
        if pdbLine[0:6] == 'COMPND':
            while not pdbLine.strip()[-1] == ';':
                lastPositionInFile = pdbFile.tell()
                newPDBLine = pdbFile.readline()
                if not newPDBLine[0:6] == 'COMPND':
                    pdbFile.seek(lastPositionInFile)
                    break
                else:
                    pdbLine = '%s %s' %(pdbLine.strip(), newPDBLine[10:])
            compoundLine = pdbLine[10:].strip()
            splittedCompoundLine = compoundLine.replace(';','').split(':')
            compoundSpecification = splittedCompoundLine[0].strip()
            if compoundSpecification == 'MOL_ID':
                compoundDefinition = splittedCompoundLine[1].strip()
                moleculeID = compoundDefinition
                if not moleculeID in moleculeSpecificationDict:
                    moleculeSpecificationDict[moleculeID] = {}
                    moleculeSpecificationDict[moleculeID]['definition'] = ''
                    moleculeSpecificationDict[moleculeID]['synonym'] = ''
                    moleculeSpecificationDict[moleculeID]['organismScientific'] = ''
                    moleculeSpecificationDict[moleculeID]['organismCommon'] = ''
                    moleculeSpecificationDict[moleculeID]['organismTaxID'] = ''
            elif compoundSpecification == 'MOLECULE':
                compoundDefinition = splittedCompoundLine[1].strip()
                moleculeSpecificationDict[moleculeID]['definition'] = compoundDefinition
            elif compoundSpecification == 'CHAIN':
                compoundDefinition = splittedCompoundLine[1].strip()
                splittedCompoundDefinition = compoundDefinition.strip().split(',')
                for tempChain in splittedCompoundDefinition:
                    tempChain = tempChain.strip()
                    chainMoleculeIDMappingDict[tempChain] = moleculeID
            elif compoundSpecification == 'SYNONYM':
                compoundDefinition = splittedCompoundLine[1].strip()
                moleculeSpecificationDict[moleculeID]['synonym'] = compoundDefinition
        elif pdbLine[0:6] == 'SOURCE':
            while not pdbLine.strip()[-1] == ';':
                lastPositionInFile = pdbFile.tell()
                newPDBLine = pdbFile.readline()
                if not newPDBLine[0:6] == 'SOURCE':
                    pdbFile.seek(lastPositionInFile)
                    break
                else:
                    pdbLine = '%s %s' %(pdbLine.strip(), newPDBLine[10:])
            sourceLine = pdbLine[10:].strip()
            splittedSourceLine = sourceLine.replace(';','').split(':')
            sourceSpecification = splittedSourceLine[0].strip()
            if sourceSpecification == 'MOL_ID':
                sourceDefinition = splittedSourceLine[1].strip()
                moleculeID = sourceDefinition
                if not moleculeID in moleculeSpecificationDict:
                    moleculeSpecificationDict[moleculeID] = {}
                    moleculeSpecificationDict[moleculeID]['definition'] = ''
                    moleculeSpecificationDict[moleculeID]['synonym'] = ''
                    moleculeSpecificationDict[moleculeID]['organismScientific'] = ''
                    moleculeSpecificationDict[moleculeID]['organismCommon'] = ''
                    moleculeSpecificationDict[moleculeID]['organismTaxID'] = ''
            elif sourceSpecification == 'ORGANISM_SCIENTIFIC':
                sourceDefinition = splittedSourceLine[1].strip()
                moleculeSpecificationDict[moleculeID]['organismScientific'] = sourceDefinition
            elif sourceSpecification == 'ORGANISM_COMMON':
                sourceDefinition = splittedSourceLine[1].strip()
                moleculeSpecificationDict[moleculeID]['organismCommon'] = sourceDefinition
            elif sourceSpecification == 'ORGANISM_TAXID':
                sourceDefinition = splittedSourceLine[1].strip()
                moleculeSpecificationDict[moleculeID]['organismTaxID'] = sourceDefinition
        elif pdbLine[0:6] == 'REMARK':
            if pdbLine[6:10].strip() == '2':
                if pdbLine[11:22] == 'RESOLUTION.':
                    if pdbLine[22:].count('NOT'):
                        moleculeSpecificationDict['resolution'] = ''
                    else:
                        if pdbLine[23:30].strip() == 'NULL':
                            moleculeSpecificationDict['resolution'] = ''
                        else:
                            moleculeSpecificationDict['resolution'] = pdbLine[23:30].strip()
            elif int(pdbLine[6:10].strip()) < 300 and int(pdbLine[6:10].strip()) >= 200:
                if pdbLine[12:27] == 'EXPERIMENT TYPE':
                    splittedRemarkLine = pdbLine[12:].strip().split(':')
                    if len(splittedRemarkLine) < 2:
                        continue
                    splittedExperimentalType = splittedRemarkLine[1].split(';')
                    experimentalType = splittedExperimentalType[0].strip()
                    if experimentalType == 'NULL':
                        experimentalType = ''
                    moleculeSpecificationDict['experimentalType'] = experimentalType
                if pdbLine[13:26] == 'SPECIMEN TYPE':
                    remarkNumber = pdbLine[6:10].strip()
                    if remarkNumber == '240':
                        experimentalType = 'ELECTRON CRYSTALLOGRAPHY'
                    elif remarkNumber == '245':
                        experimentalType = 'ELECTRON MICROSCOPY'
                    else:
                        experimentalType = ''
                    moleculeSpecificationDict['experimentalType'] = experimentalType
                if pdbLine[12:23] == 'TEMPERATURE':
                    splittedRemarkLine = pdbLine[12:].strip().split(':')
                    if len(splittedRemarkLine) < 2:
                        continue
                    splittedTemperature = splittedRemarkLine[1].split(';')
                    temperature = splittedTemperature[0].strip()
                    if temperature == 'NULL':
                        temperature = ''
                    moleculeSpecificationDict['temperature'] = temperature
                elif pdbLine[12:14] == 'PH':
                    splittedRemarkLine = pdbLine[12:].strip().split(':')
                    if len(splittedRemarkLine) < 2:
                        continue
                    splittedPH = splittedRemarkLine[1].split(';')
                    PH = splittedPH[0].strip()
                    if PH == 'NULL':
                        PH = ''
                    moleculeSpecificationDict['PH'] = PH

    pdbFile.close()

    resolution = moleculeSpecificationDict['resolution']
    PH = moleculeSpecificationDict['PH']
    temperature = moleculeSpecificationDict['temperature']
    experimentalType = moleculeSpecificationDict['experimentalType']
    for chainID in splittedInfoLine[1:]:
        monomerID = '%s_%s' %(pdbID, chainID)
        moleculeID = chainMoleculeIDMappingDict[chainID]
        moleculeDefinition = moleculeSpecificationDict[moleculeID]['definition']
        moleculeSynonym = moleculeSpecificationDict[moleculeID]['synonym']
        moleculeOrganismScientific = moleculeSpecificationDict[moleculeID]['organismScientific']
        moleculeOrganismCommon = moleculeSpecificationDict[moleculeID]['organismCommon']
        moleculeOrganismTaxID = moleculeSpecificationDict[moleculeID]['organismTaxID']
        pdbChainSpecificationFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(monomerID, moleculeDefinition, moleculeSynonym, moleculeOrganismScientific, moleculeOrganismCommon, moleculeOrganismTaxID, experimentalType, resolution, PH, temperature))
                
pdbChainListFile.close()
pdbChainSpecificationFile.close()
t2 = time.time()
print('Elapsed time = %f' %(t2-t1))
