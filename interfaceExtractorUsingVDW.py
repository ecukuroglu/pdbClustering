#################
## Created by Engin Cukuroglu
#################

import multiprocessing
from codesOfTools import interfaceNameSorter
from multiprocessing import Queue, Process
import os, sys, time
from codesOfTools import *

def interfaceExtractorWithInterfaceResidueStatusOfChains_On(fullChainAtomPropertyList, fullChainCoordinatesCA, fullChainDict, VDWCriteria):
    interfaceDict = {}
    skeletonDict = {}
    contactDict = {}
    tempFactorDict = {}
    for i in list(range(0,len(fullChainAtomPropertyList))):
        atom_1 = fullChainAtomPropertyList[i]
        resDictKey_1 = '%s_%s_%s' %(atom_1[2], atom_1[1], atom_1[0])
        if atom_1[3] == 'CA':
            tempFactorDict[resDictKey_1] = atom_1[6]
        elif not resDictKey_1 in tempFactorDict:
            tempFactorDict[resDictKey_1] = ''
        for j in list(range(i+1,len(fullChainAtomPropertyList))):
            atom_2 = fullChainAtomPropertyList[j]
            resDictKey_2 = '%s_%s_%s' %(atom_2[2], atom_2[1], atom_2[0])
            if not resDictKey_1 == resDictKey_2:
                distanceThreshold = atom_1[5] + atom_2[5] + VDWCriteria
                distanceBetweenAtoms = distanceCalculator(atom_1[4], atom_2[4])
                if distanceBetweenAtoms < distanceThreshold:
                    if not resDictKey_1 in contactDict:
                        contactDict[resDictKey_1] = {}
                    if not resDictKey_2 in contactDict:
                        contactDict[resDictKey_2] = {}
                    contactDict[resDictKey_1][resDictKey_2] = 1
                    contactDict[resDictKey_2][resDictKey_1] = 1
                    if not atom_1[0] == atom_2[0]:
                        if resDictKey_1 in fullChainCoordinatesCA:
                            if resDictKey_2 in fullChainCoordinatesCA:
                                interfaceDict[resDictKey_1] = 1
                                interfaceDict[resDictKey_2] = 1
    for interfaceResidue in interfaceDict:
        splittedInterfaceResidue = interfaceResidue.split('_')
        resChain = splittedInterfaceResidue[2]
        fullChainDict[resChain]['interfaceResidues'] = fullChainDict[resChain]['interfaceResidues'] + 1 
        if not interfaceResidue in skeletonDict:
            skeletonDict[interfaceResidue] = 1
            fullChainDict[resChain]['skeletonResidues'] = fullChainDict[resChain]['skeletonResidues'] + 1 
        for contactResidue in contactDict[interfaceResidue]:
            splittedContactResidue = contactResidue.split('_')
            if splittedContactResidue[2] == resChain:
                if not contactResidue in skeletonDict:
                    skeletonDict[contactResidue] = 1
                    fullChainDict[resChain]['skeletonResidues'] = fullChainDict[resChain]['skeletonResidues'] + 1 
    return interfaceDict, skeletonDict, contactDict, tempFactorDict, fullChainDict

def interfaceExtractorWithInterfaceResidueStatusOfChains_Off(fullChainAtomPropertyList, fullChainCoordinatesCA, fullChainDict, VDWCriteria):
    interfaceDict = {}
    skeletonDict = {}
    contactDict = {}
    tempFactorDict = {}
    for i in list(range(0,len(fullChainAtomPropertyList))):
        atom_1 = fullChainAtomPropertyList[i]
        resDictKey_1 = '%s_%s_%s' %(atom_1[2], atom_1[1], atom_1[0])
        tempFactorDict[resDictKey_1] = atom_1[6]
        if atom_1[3] == 'CA':
            tempFactorDict[resDictKey_1] = atom_1[6]
        elif not resDictKey_1 in tempFactorDict:
            tempFactorDict[resDictKey_1] = ''
        for j in list(range(i+1,len(fullChainAtomPropertyList))):
            atom_2 = fullChainAtomPropertyList[j]
            resDictKey_2 = '%s_%s_%s' %(atom_2[2], atom_2[1], atom_2[0])
            if not resDictKey_1 == resDictKey_2:
                distanceThreshold = atom_1[5] + atom_2[5] + VDWCriteria
                distanceBetweenAtoms = distanceCalculator(atom_1[4], atom_2[4])
                if distanceBetweenAtoms < distanceThreshold:
                    if not resDictKey_1 in contactDict:
                        contactDict[resDictKey_1] = {}
                    if not resDictKey_2 in contactDict:
                        contactDict[resDictKey_2] = {}
                    contactDict[resDictKey_1][resDictKey_2] = 1
                    contactDict[resDictKey_2][resDictKey_1] = 1
                    if not fullChainDict[atom_1[0]]['chain'] == fullChainDict[atom_2[0]]['chain']:
                        if resDictKey_1 in fullChainCoordinatesCA:
                            if resDictKey_2 in fullChainCoordinatesCA:
                                interfaceDict[resDictKey_1] = 1
                                interfaceDict[resDictKey_2] = 1
    for interfaceResidue in interfaceDict:
        splittedInterfaceResidue = interfaceResidue.split('_')
        resChain = splittedInterfaceResidue[2]
        fullChainDict[resChain]['interfaceResidues'] = fullChainDict[resChain]['interfaceResidues'] + 1 
        if not interfaceResidue in skeletonDict:
            skeletonDict[interfaceResidue] = 1
            fullChainDict[resChain]['skeletonResidues'] = fullChainDict[resChain]['skeletonResidues'] + 1 
        for contactResidue in contactDict[interfaceResidue]:
            splittedContactResidue = contactResidue.split('_')
            if splittedContactResidue[2] == resChain:
                if not contactResidue in skeletonDict:
                    skeletonDict[contactResidue] = 1
                    fullChainDict[resChain]['skeletonResidues'] = fullChainDict[resChain]['skeletonResidues'] + 1 
    return interfaceDict, skeletonDict, contactDict, tempFactorDict, fullChainDict

def interfaceExtractionWork(taskQueue_interface, interfaceResultsQueue, allPDBFilesDirectory, VDWCriteria):
    while True:
        interface, interfaceFileDirectories = taskQueue_interface.get()
        taskQueue_interface.task_done()
        if interface == None:
            break
        pdbName, chain_1, chain_2, allChains, interfaceResidueStatusOfChains, interface, commonChainCounter = interfaceNameSorter(interface)
            
        interfacePDBFileDirectory = interfaceFileDirectories[0]
        interfaceSkeletonPDBFileDirectory = interfaceFileDirectories[1]
        interfaceResidueContactInfoFileDirectory = interfaceFileDirectories[2]
        fullChainDict = {}
        chain_1_dict = {}
        for i in list(range(len(chain_1))):
            chainID = chain_1[i]
            chain_1_dict[chainID] = 1
            fullChainDict[chainID] = {}
            fullChainDict[chainID]['chain'] = 1
            fullChainDict[chainID]['interfaceResidues'] = 0
            fullChainDict[chainID]['skeletonResidues'] = 0
        chain_2_dict = {}
        for i in list(range(len(chain_2))):
            chainID = chain_2[i]
            chain_2_dict[chainID] = 2
            fullChainDict[chainID] = {}
            fullChainDict[chainID]['chain'] = 2
            fullChainDict[chainID]['interfaceResidues'] = 0
            fullChainDict[chainID]['skeletonResidues'] = 0
        
        pdbFileDirectory = '%s/%s.pdb' %(allPDBFilesDirectory, pdbName)
        fullChainAtomPropertyList, fullChainTypeDict, fullChainCoordinatesCA = readPDBReturnAtomProperties_withMultipleChains(pdbFileDirectory, fullChainDict)

        proteinInterfaceCalculationStatus = 1
        
        for chainID in fullChainTypeDict:
            if not fullChainTypeDict[chainID]['type'] == 'Protein':
                proteinInterfaceCalculationStatus = 0
                chainLengthString = '%d' %(fullChainTypeDict[chain_1[0]]['residueNumber'])
                chainTypeString = '%s' %(fullChainTypeDict[chain_1[0]]['type'])
                for i in list(range(1,len(chain_1))):
                    chainLengthString = '%s,%d' %(chainLengthString,fullChainTypeDict[chain_1[i]]['residueNumber'])
                    chainTypeString = '%s,%s' %(chainTypeString,fullChainTypeDict[chain_1[i]]['type'])
                chainLengthString = '%s\t%d' %(chainLengthString,fullChainTypeDict[chain_2[0]]['residueNumber'])
                chainTypeString = '%s\t%s' %(chainTypeString,fullChainTypeDict[chain_2[0]]['type'])
                for i in list(range(1,len(chain_2))):
                    chainLengthString = '%s,%d' %(chainLengthString,fullChainTypeDict[chain_2[i]]['residueNumber'])
                    chainTypeString = '%s,%s' %(chainTypeString,fullChainTypeDict[chain_2[i]]['type'])
                interfaceExtractionResultText = '%s\t1\t%s\t%s\t%s\t%s\tAt least one of the chain has a residue which is not in common 20 amino acids.\n' %(interface, chain_1, chain_2, chainLengthString, chainTypeString)
                break
            if fullChainTypeDict[chainID]['residueNumber'] == 0:
                proteinInterfaceCalculationStatus = 0
                chainTypeString = '%s' %(fullChainTypeDict[chain_1[0]]['type'])
                chainLengthString = '%d' %(fullChainTypeDict[chain_1[0]]['residueNumber'])
                for i in list(range(1,len(chain_1))):
                    chainTypeString = '%s,%s' %(chainTypeString,fullChainTypeDict[chain_1[i]]['type'])
                    chainLengthString = '%s,%d' %(chainLengthString,fullChainTypeDict[chain_1[i]]['residueNumber'])
                chainTypeString = '%s\t%s' %(chainTypeString,fullChainTypeDict[chain_2[0]]['type'])
                chainLengthString = '%s\t%d' %(chainLengthString,fullChainTypeDict[chain_2[0]]['residueNumber'])
                for i in list(range(1,len(chain_2))):
                    chainTypeString = '%s,%s' %(chainTypeString,fullChainTypeDict[chain_2[i]]['type'])
                    chainLengthString = '%s,%d' %(chainLengthString,fullChainTypeDict[chain_2[i]]['residueNumber'])
                interfaceExtractionResultText = '%s\t2\t%s\t%s\t%s\t%s\tAt least one of the chain does not present in pdb file.\n' %(interface, chain_1, chain_2, chainLengthString, chainTypeString)
                break
        if proteinInterfaceCalculationStatus == 1:
            if interfaceResidueStatusOfChains == '0':
                interfaceDict, interfaceSkeletonDict, contactDict, tempFactorDict, fullChainDict = interfaceExtractorWithInterfaceResidueStatusOfChains_On(fullChainAtomPropertyList, fullChainCoordinatesCA, fullChainDict, VDWCriteria)
            elif interfaceResidueStatusOfChains == '1':
                interfaceDict, interfaceSkeletonDict, contactDict, tempFactorDict, fullChainDict = interfaceExtractorWithInterfaceResidueStatusOfChains_Off(fullChainAtomPropertyList, fullChainCoordinatesCA, fullChainDict, VDWCriteria)
            else:
                print('%s does not have interfaceResidueStatusOfChains options' %(interface))
                continue
            pdbFile = open(pdbFileDirectory, 'r')
            interfacePDBFile = open(interfacePDBFileDirectory, 'w')
            interfaceSkeletonPDBFile = open(interfaceSkeletonPDBFileDirectory, 'w')
            writtenCALineCheckerDictionary = {}
            while 1:
                pdbLine = pdbFile.readline()
                if pdbLine == '' or pdbLine[0:3] == 'END':
                    break
                if pdbLine[0:4] == 'ATOM':
                    atomType = pdbLine[12:16].strip()
                    if atomType == 'CA':
                        chainID = pdbLine[21]
                        if chainID in fullChainDict:
                            resType = pdbLine[17:20].strip()
                            resNo = pdbLine[22:26].strip()
                            resDictKey = '%s_%s_%s' %(resType, resNo, chainID)
                            if not resDictKey in writtenCALineCheckerDictionary:
                                writtenCALineCheckerDictionary[resDictKey] = 1
                                if resDictKey in interfaceDict:
                                    interfacePDBFile.write('%s' %(pdbLine))
                                if resDictKey in interfaceSkeletonDict:
                                    interfaceSkeletonPDBFile.write('%s' %(pdbLine))
            pdbFile.close()
            interfacePDBFile.close()
            interfaceSkeletonPDBFile.close()
            
            interfaceResidueContactInfoFile = open(interfaceResidueContactInfoFileDirectory, 'w')
            for residue in contactDict:
                tempFactorNeighborTotal = 0.0
                tempFactorNeighborNumber = 0
                tempFactorNeighborAveString = ''
                contactResidueString = ''
                contactResidueList = contactDict[residue].keys()
                if len(contactResidueList) > 0:
                                        contactResidue = contactResidueList[0]
                                        contactResidueString = '%s' %(contactResidue)
                                        if not tempFactorDict[contactResidue] == '':
                                                tempFactorNeighborTotal = tempFactorNeighborTotal + float(tempFactorDict[contactResidue])
                                                tempFactorNeighborNumber = tempFactorNeighborNumber + 1
                                        for i in list(range(1,len(contactResidueList))):
                                                contactResidue = contactResidueList[i]
                                                contactResidueString = '%s,%s' %(contactResidueString, contactResidue)
                                                if not tempFactorDict[contactResidue] == '':
                                                        tempFactorNeighborTotal = tempFactorNeighborTotal + float(tempFactorDict[contactResidue])
                                                        tempFactorNeighborNumber = tempFactorNeighborNumber + 1
                if not tempFactorNeighborNumber == 0:
                    tempFactorNeighborAveString = '%.3f' %(tempFactorNeighborTotal/tempFactorNeighborNumber)
                interfaceResidueContactInfoFile.write('%s\t%s\t%s\t%s\t%s\n' %(residue, tempFactorDict[residue], tempFactorNeighborAveString, len(contactDict[residue]), contactResidueString))
            interfaceResidueContactInfoFile.close()

            chain_1_lengthString = '%d' %(fullChainTypeDict[chain_1[0]]['residueNumber'])
            chain_1_typeString = '%s' %(fullChainTypeDict[chain_1[0]]['type'])
            chain_1_interfaceLengthString = '%d' %(fullChainDict[chain_1[0]]['interfaceResidues'])
            chain_1_skeletonLengthString = '%d' %(fullChainDict[chain_1[0]]['skeletonResidues'])
            overallChain_1_interfaceLength = fullChainDict[chain_1[0]]['interfaceResidues']
            overallChain_1_skeletonLength = fullChainDict[chain_1[0]]['skeletonResidues']
            for i in list(range(1,len(chain_1))):
                chain_1_lengthString = '%s,%d' %(chain_1_lengthString,fullChainTypeDict[chain_1[i]]['residueNumber'])
                chain_1_typeString = '%s,%s' %(chain_1_typeString, fullChainTypeDict[chain_1[i]]['type'])
                chain_1_interfaceLengthString = '%s,%d' %(chain_1_interfaceLengthString, fullChainDict[chain_1[i]]['interfaceResidues'])
                chain_1_skeletonLengthString = '%s,%d' %(chain_1_skeletonLengthString, fullChainDict[chain_1[i]]['skeletonResidues'])
                overallChain_1_interfaceLength = overallChain_1_interfaceLength + fullChainDict[chain_1[i]]['interfaceResidues']
                overallChain_1_skeletonLength = overallChain_1_skeletonLength + fullChainDict[chain_1[i]]['skeletonResidues']
            chain_2_lengthString = '%d' %(fullChainTypeDict[chain_2[0]]['residueNumber'])
            chain_2_typeString = '%s' %(fullChainTypeDict[chain_2[0]]['type'])
            chain_2_interfaceLengthString = '%d' %(fullChainDict[chain_2[0]]['interfaceResidues'])
            chain_2_skeletonLengthString = '%d' %(fullChainDict[chain_2[0]]['skeletonResidues'])
            overallChain_2_interfaceLength = fullChainDict[chain_2[0]]['interfaceResidues']
            overallChain_2_skeletonLength = fullChainDict[chain_2[0]]['skeletonResidues']
            for i in list(range(1,len(chain_2))):
                chain_2_lengthString = '%s,%d' %(chain_2_lengthString,fullChainTypeDict[chain_2[i]]['residueNumber'])
                chain_2_typeString = '%s,%s' %(chain_2_typeString, fullChainTypeDict[chain_2[i]]['type'])
                chain_2_interfaceLengthString = '%s,%d' %(chain_2_interfaceLengthString, fullChainDict[chain_2[i]]['interfaceResidues'])
                chain_2_skeletonLengthString = '%s,%d' %(chain_2_skeletonLengthString, fullChainDict[chain_2[i]]['skeletonResidues'])
                overallChain_2_interfaceLength = overallChain_2_interfaceLength + fullChainDict[chain_2[i]]['interfaceResidues']
                overallChain_2_skeletonLength = overallChain_2_skeletonLength + fullChainDict[chain_2[i]]['skeletonResidues']

            interfaceExtractionResultText = '%s\t0\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\n' %(interface, chain_1, chain_2, chain_1_lengthString, chain_2_lengthString, chain_1_typeString, chain_2_typeString, chain_1_interfaceLengthString, chain_2_interfaceLengthString, chain_1_skeletonLengthString, chain_2_skeletonLengthString, overallChain_1_interfaceLength, overallChain_2_interfaceLength, overallChain_1_skeletonLength, overallChain_2_skeletonLength)
            interfaceResultsQueue.put(interfaceExtractionResultText)
        else:
            interfaceResultsQueue.put(interfaceExtractionResultText)
    return

def mainInterfaceExtractorUsingVDW(interfaceListFileDirectory, allPDBFilesDirectory, VDWCriteria, VDWInterfaceFileDirectory, VDWInterfaceSkeletonFileDirectory, residueContactInfoFileDirectory, VDWInterfaceTypeListFileDirectory, fullVDWInterfaceTypeListFileDirectory, fullPDBIDListFileDirectory, numberOfProcesses):
    
    print('\n* INTERFACE EXTRACTION USING VDW STARTED *\n')
    print('Time stamp : %s' %(time.asctime()))
    t1 = time.time()

    if not os.path.exists(interfaceListFileDirectory):
        sys.exit('\nThe interfaceListFileDirectory does not exist\n')

    if not os.path.exists(allPDBFilesDirectory):
        sys.exit('\nThe pdbFile path is not correct\n')

    if not os.path.exists(VDWInterfaceFileDirectory):
        os.system('mkdir %s' %(VDWInterfaceFileDirectory))

    if not os.path.exists(VDWInterfaceSkeletonFileDirectory):
        os.system('mkdir %s' %(VDWInterfaceSkeletonFileDirectory))

    if not os.path.exists(residueContactInfoFileDirectory):
        os.system('mkdir %s' %(residueContactInfoFileDirectory))

    fullVDWInterfaceResultDict = {}
    if os.path.exists(fullVDWInterfaceTypeListFileDirectory):
        fullVDWInterfaceTypeListFile = open(fullVDWInterfaceTypeListFileDirectory, 'r')
        for vdwInterfaceResultEntry in fullVDWInterfaceTypeListFile:
            splittedVDWInterfaceResultEntry = vdwInterfaceResultEntry.strip().split()
            fullVDWInterfaceResultDict[splittedVDWInterfaceResultEntry[0]] = [splittedVDWInterfaceResultEntry[1], vdwInterfaceResultEntry]
        fullVDWInterfaceTypeListFile.close()

    totalNumberOfInterfaces = 0
    taskDict = {}
    interfaceResultsQueue = multiprocessing.Queue()
    interfaceListFile = open(interfaceListFileDirectory, 'r')
    for interface in interfaceListFile:
        interface = interface.strip()
        if interface == '':
            continue
        totalNumberOfInterfaces = totalNumberOfInterfaces + 1
        pdbName, chain_1, chain_2, allChains, interfaceResidueStatusOfChains, interface, commonChainCounter = interfaceNameSorter(interface)
        if len(chain_1) == 0:
            interfaceExtractionResultText = '%s\t2\tChain_1 does not have any monomer.\n' %(interface)
            interfaceResultsQueue.put(interfaceExtractionResultText)
            continue
        if len(chain_2) == 0:
            interfaceExtractionResultText = '%s\t2\tChain_2 does not have any monomer.\n' %(interface)
            interfaceResultsQueue.put(interfaceExtractionResultText)
            continue
        if commonChainCounter > 0:
            interfaceExtractionResultText = '%s\t2\tChain_1 and Chain_2 have %d common monomers.\n' %(interface, commonChainCounter)
            interfaceResultsQueue.put(interfaceExtractionResultText)
            continue
            
        interfacePDBFileDirectory = '%s/%s.pdb' %(VDWInterfaceFileDirectory, interface)
        interfaceSkeletonPDBFileDirectory = '%s/%s.pdb' %(VDWInterfaceSkeletonFileDirectory, interface)
        interfaceResidueContactInfoFileDirectory = '%s/%s_contacts.txt' %(residueContactInfoFileDirectory, interface)

        if interface in fullVDWInterfaceResultDict:
            if int(fullVDWInterfaceResultDict[interface][0]) == 0:
                if os.path.exists(interfacePDBFileDirectory):
                    if os.path.exists(interfaceSkeletonPDBFileDirectory):
                        if os.path.exists(interfaceResidueContactInfoFileDirectory):
                            interfaceExtractionResultText = fullVDWInterfaceResultDict[interface][1]
                            interfaceResultsQueue.put(interfaceExtractionResultText)
                            continue
        taskDict[interface] = [interfacePDBFileDirectory, interfaceSkeletonPDBFileDirectory, interfaceResidueContactInfoFileDirectory]
    interfaceListFile.close()

    taskQueue_interface = multiprocessing.JoinableQueue()
    interfaceWorkers = [Process(target=interfaceExtractionWork, args=(taskQueue_interface, interfaceResultsQueue, allPDBFilesDirectory, VDWCriteria)) for i in range(numberOfProcesses)]
    for tempWorkers in interfaceWorkers:
        tempWorkers.start()

    interfaceListFile = open(interfaceListFileDirectory, 'r')
    for interface in taskDict:
        taskQueue_interface.put([interface, taskDict[interface]])
    interfaceListFile.close()

    for i in range(numberOfProcesses):
        taskQueue_interface.put([None, None])

    taskQueue_interface.join()
    VDWInterfaceTypeListFile = open(VDWInterfaceTypeListFileDirectory, 'w')
    interfaceTypeDict = {}
    processedTotalNumberOfInterfaces = totalNumberOfInterfaces
    while totalNumberOfInterfaces:
        interfaceResultString = interfaceResultsQueue.get()
        splittedInterfaceResultString = interfaceResultString.strip().split('\t')
        interfaceTypeDict[splittedInterfaceResultString[0]] = [0, splittedInterfaceResultString[1], interfaceResultString]
        VDWInterfaceTypeListFile.write(interfaceResultString)
        totalNumberOfInterfaces = totalNumberOfInterfaces - 1
    VDWInterfaceTypeListFile.close()

    fullPDBIDDict = {}
    fullPDBIDListFile = open(fullPDBIDListFileDirectory, 'r')
    for pdbID in fullPDBIDListFile:
        pdbID = pdbID.strip()
        fullPDBIDDict[pdbID] = 1
    fullPDBIDListFile.close()

    tempFullVDWInterfaceTypeListFileDirectory = '%s_temp' %(fullVDWInterfaceTypeListFileDirectory)
    tempFullVDWInterfaceTypeListFile = open(tempFullVDWInterfaceTypeListFileDirectory, 'w')
    if os.path.exists(fullVDWInterfaceTypeListFileDirectory):
        fullVDWInterfaceTypeListFile = open(fullVDWInterfaceTypeListFileDirectory, 'r')
        for interfaceTypeLine in fullVDWInterfaceTypeListFile:
            splittedInterfaceTypeLine = interfaceTypeLine.strip().split('\t')
            interfaceName = splittedInterfaceTypeLine[0]
            splittedInterfaceName = interfaceName.split('_')
            pdbID = splittedInterfaceName[0]
            if pdbID in fullPDBIDDict:
                if interfaceName in interfaceTypeDict:
                    interfaceTypeDict[interfaceName][0] = 1
                    if splittedInterfaceTypeLine[1] >= interfaceTypeDict[interfaceName][1]:
                        tempFullVDWInterfaceTypeListFile.write(interfaceTypeDict[interfaceName][2])
                    else:
                        tempFullVDWInterfaceTypeListFile.write(interfaceTypeLine)
                else:
                    tempFullVDWInterfaceTypeListFile.write(interfaceTypeLine)
        fullVDWInterfaceTypeListFile.close()
        for interfaceName in interfaceTypeDict:
            if interfaceTypeDict[interfaceName][0] == 0:
                tempFullVDWInterfaceTypeListFile.write(interfaceTypeDict[interfaceName][2])
    else:
        for interfaceName in interfaceTypeDict:
            tempFullVDWInterfaceTypeListFile.write(interfaceTypeDict[interfaceName][2])
    tempFullVDWInterfaceTypeListFile.close()
    os.system('mv %s %s' %(tempFullVDWInterfaceTypeListFileDirectory, fullVDWInterfaceTypeListFileDirectory))
        
    t2 = time.time()
    print('\nTotal number of processed interfaces = %d\n' %(processedTotalNumberOfInterfaces))
    print('\nElapsed time = %f seconds\n' %(t2-t1))
    print('Time stamp : %s' %(time.asctime()))
    print('\n* INTERFACE EXTRACTION USING VDW COMPLETED *\n')
