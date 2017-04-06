#################
## Created by Engin Cukuroglu
#################

import os, sys, glob, time
from codesOfTools import pdbChainListExtractor


def mainPossibleInterfaceExtractor(pdbFileDirectory, possibleInterfaceListFileDirectory, monomerListFileDirectory, complexListFileDirectory, pdbChainListFileDirectory, fullPDBIDListFileDirectory, fullPossibleInterfaceListFileDirectory, fullMonomerListFileDirectory, fullComplexListFileDirectory, fullPDBChainListFileDirectory, interfaceResidueStatusOfChains, possibleInterfaceLogFileDirectory):
    t1 = time.time()
    print('\n* POSSIBLE INTERFACE FINDER STARTED *\n')
    print('Time stamp : %s' %(time.asctime()))

    possibleInterfaceLogFile = open(possibleInterfaceLogFileDirectory, 'w')
    exPossibleInterfaceListDict = {}
    if os.path.exists(fullPossibleInterfaceListFileDirectory):
        fullPossibleInterfaceListFile = open(fullPossibleInterfaceListFileDirectory, 'r')
        for possibleInterfaceListLine in fullPossibleInterfaceListFile:
            possibleInterface = possibleInterfaceListLine.strip()
            splittedPossibleInterface = possibleInterface.split('_')
            pdbID = splittedPossibleInterface[0]
            if pdbID in exPossibleInterfaceListDict:
                exPossibleInterfaceListDict[pdbID].append(possibleInterface)
            else:
                exPossibleInterfaceListDict[pdbID] = []
                exPossibleInterfaceListDict[pdbID].append(possibleInterface)
        fullPossibleInterfaceListFile.close()
    else:
        possibleInterfaceLogFile.write('%s does not present' %(fullPossibleInterfaceListFileDirectory))

    exPDBChainListDict = {}
    if os.path.exists(fullPDBChainListFileDirectory):
        fullPDBChainListFile = open(fullPDBChainListFileDirectory, 'r')
        for pdbChainListLine in fullPDBChainListFile:
            splittedPDBChainListLine = pdbChainListLine.strip().split('\t')
            if len(splittedPDBChainListLine) > 1:
                exPDBChainListDict[splittedPDBChainListLine[0]] = splittedPDBChainListLine[1:]
            else:
                exPDBChainListDict[splittedPDBChainListLine[0]] = []
        fullPDBChainListFile.close()
    else:
        possibleInterfaceLogFile.write('%s does not present' %(fullPDBChainListFileDirectory))
        
    tempFullPDBChainListFileDirectory = '%s_temp' %(fullPDBChainListFileDirectory)
    tempFullPDBChainListFile = open(tempFullPDBChainListFileDirectory, 'w')
    tempFullPossibleInterfaceListFileDirectory = '%s_temp' %(fullPossibleInterfaceListFileDirectory)
    tempFullPossibleInterfaceListFile = open(tempFullPossibleInterfaceListFileDirectory, 'w')
    tempFullMonomerListFileDirectory = '%s_temp' %(fullMonomerListFileDirectory)
    tempFullMonomerListFile = open(tempFullMonomerListFileDirectory, 'w')
    tempFullComplexListFileDirectory = '%s_temp' %(fullComplexListFileDirectory)
    tempFullComplexListFile = open(tempFullComplexListFileDirectory, 'w')
    possibleInterfaceListFile = open(possibleInterfaceListFileDirectory, 'w')
    monomerListFile = open(monomerListFileDirectory, 'w')
    complexListFile = open(complexListFileDirectory, 'w')
    pdbChainListFile = open(pdbChainListFileDirectory, 'w')

    fullPDBIDListFile = open(fullPDBIDListFileDirectory, 'r')
    for pdbID in fullPDBIDListFile:
        pdbID = pdbID.strip()
        if pdbID in exPDBChainListDict:
            tempFullPDBChainListFile.write('%s' %(pdbID))
            pdbChainListFile.write('%s' %(pdbID))
            chainList = exPDBChainListDict[pdbID]
            for chainID in chainList:
                tempFullPDBChainListFile.write('\t%s' %(chainID))
                pdbChainListFile.write('\t%s' %(chainID))
            tempFullPDBChainListFile.write('\n')
            pdbChainListFile.write('\n')
            if len(chainList) > 1:
                tempFullComplexListFile.write('%s\n' %(pdbID))
                complexListFile.write('%s\n' %(pdbID))
                if pdbID in exPossibleInterfaceListDict:
                    for exPossibleInterface in exPossibleInterfaceListDict[pdbID]:
                        tempFullPossibleInterfaceListFile.write('%s\n' %(exPossibleInterface))
                        possibleInterfaceListFile.write('%s\n' %(exPossibleInterface))
                else:
                    for i in list(range(len(chainList))):
                        for j in list(range(i+1, len(chainList),1)):
                            tempFullPossibleInterfaceListFile.write('%s_%s_%s_%d\n' %(pdbID,chainList[i],chainList[j], interfaceResidueStatusOfChains))
                            possibleInterfaceListFile.write('%s_%s_%s_%d\n' %(pdbID,chainList[i],chainList[j], interfaceResidueStatusOfChains))
            elif len(chainList) == 1: 
                tempFullMonomerListFile.write('%s\n' %(pdbID))
                monomerListFile.write('%s\n' %(pdbID))
            else:
                possibleInterfaceLogFile.write('There is no chain in %s\n' %(pdbID))

        else:
            currentPDBFileDirectory = '%s/%s.pdb' %(pdbFileDirectory, pdbID)
            if os.path.exists(currentPDBFileDirectory):
                chainList = pdbChainListExtractor(currentPDBFileDirectory)
                tempFullPDBChainListFile.write('%s' %(pdbID))
                pdbChainListFile.write('%s' %(pdbID))
                for chainID in chainList:
                    tempFullPDBChainListFile.write('\t%s' %(chainID))
                    pdbChainListFile.write('\t%s' %(chainID))
                tempFullPDBChainListFile.write('\n')
                pdbChainListFile.write('\n')
                if len(chainList) > 1:
                    tempFullComplexListFile.write('%s\n' %(pdbID))
                    complexListFile.write('%s\n' %(pdbID))
                    if pdbID in exPossibleInterfaceListDict:
                        for exPossibleInterface in exPossibleInterfaceListDict[pdbID]:
                            tempFullPossibleInterfaceListFile.write('%s\n' %(exPossibleInterface))
                            possibleInterfaceListFile.write('%s\n' %(exPossibleInterface))
                    else:
                        for i in list(range(len(chainList))):
                            for j in list(range(i+1, len(chainList),1)):
                                tempFullPossibleInterfaceListFile.write('%s_%s_%s_%d\n' %(pdbID,chainList[i],chainList[j], interfaceResidueStatusOfChains))
                                possibleInterfaceListFile.write('%s_%s_%s_%d\n' %(pdbID,chainList[i],chainList[j], interfaceResidueStatusOfChains))
                elif len(chainList) == 1: 
                    tempFullMonomerListFile.write('%s\n' %(pdbID))
                    monomerListFile.write('%s\n' %(pdbID))
                else:
                    possibleInterfaceLogFile.write('There is no chain in %s\n' %(pdbID))
            else:
                possibleInterfaceLogFile.write('%s does not present in %s\n' %(pdbID, pdbFileDirectory))
    
    tempFullPDBChainListFile.close()
    tempFullPossibleInterfaceListFile.close()
    tempFullMonomerListFile.close()
    tempFullComplexListFile.close()
    possibleInterfaceListFile.close()
    monomerListFile.close()
    complexListFile.close()
    pdbChainListFile.close()
    fullPDBIDListFile.close()
    
    os.system('mv %s %s' %(tempFullPDBChainListFileDirectory, fullPDBChainListFileDirectory))    
    os.system('mv %s %s' %(tempFullPossibleInterfaceListFileDirectory, fullPossibleInterfaceListFileDirectory))    
    os.system('mv %s %s' %(tempFullMonomerListFileDirectory, fullMonomerListFileDirectory))    
    os.system('mv %s %s' %(tempFullComplexListFileDirectory, fullComplexListFileDirectory))    
    possibleInterfaceLogFile.close()
    t2 = time.time()
    print('\nElapsed time = %f seconds\n' %(t2-t1))
    print('Time stamp : %s' %(time.asctime()))
    print('\n* END POSSIBLE INTERFACE COMPLETED *\n')
    
