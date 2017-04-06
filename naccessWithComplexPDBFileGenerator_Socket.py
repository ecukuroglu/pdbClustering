#################
## Created by Engin Cukuroglu
#################

import multiprocessing
from codesOfTools import interfaceNameSorter, naccessRSAFileReaderReturnOnlyAbsASADictionary
from multiprocessing import Queue, Process
import os, sys, time
import subprocess
import tempfile
import socket

def generateComplexPDBFilesAndRunNaccessWork(taskQueue_pdbDict, naccessResultsQueue, allPDBFilesDirectory, minInterfaceResidueCriteria, differenceBetweenComplexAndMonomerASACriteria, naccessRunFileDirectory, rsaFileKeeperDict, naccessErrorQueue):
    while True:
        pdbName, pdbDict = taskQueue_pdbDict.get()
        taskQueue_pdbDict.task_done()
        if pdbName is None:
            break
        pdbFileDirectory = '%s/%s.pdb' %(allPDBFilesDirectory, pdbName)
        if os.path.exists(pdbFileDirectory):
            createPDBFileList = []
            for chainInfo in pdbDict['chains']:
                createPDBFileList.append([chainInfo, '', pdbDict['chains'][chainInfo][0], pdbDict['chains'][chainInfo][1], pdbDict['chains'][chainInfo][2], pdbDict['chains'][chainInfo][3]])
            pdbFile = open(pdbFileDirectory, 'r')
            pdbResidueDict = {}
            pdbAtomDict = {}
            pdbChainIDDict = {}
            while 1:
                pdbLine = pdbFile.readline()
                if pdbLine == '' or pdbLine[0:3] == 'END':
                    break
                if pdbLine[0:4] == 'ATOM':
                    chainID = pdbLine[21]
                    pdbChainIDDict[chainID] = 1
                    resNo = pdbLine[22:26].strip()
                    resICode = pdbLine[26]
                    atomAlternativeLocationIndicator = pdbLine[16]
                    resType = pdbLine[17:20].strip()
                    atomType = pdbLine[12:16].strip()
                    resDictKey = '%s_%s' %(resNo, chainID)
                    resDictValue = '%s_%s_%s' %(resType, resNo, chainID)
                    atomDictKey = '%s_%s_%s_%s' %(resType, resNo, chainID, atomType)
                    if not resDictKey in pdbResidueDict:
                        pdbResidueDict[resDictKey] = [resDictValue, resICode, atomAlternativeLocationIndicator]
                        pdbAtomDict[atomDictKey] = 1
                        for tempPDBFileProperty in createPDBFileList:
                            if chainID in tempPDBFileProperty[0]:
                                tempPDBFileProperty[1] = '%s%s' %(tempPDBFileProperty[1], pdbLine)
                    else:
                        if pdbResidueDict[resDictKey][0] == resDictValue:
                            if pdbResidueDict[resDictKey][1] == resICode:
                                if not atomDictKey in pdbAtomDict:
                                    pdbAtomDict[atomDictKey] = 1
                                    for tempPDBFileProperty in createPDBFileList:
                                        if chainID in tempPDBFileProperty[0]:
                                            tempPDBFileProperty[1] = '%s%s' %(tempPDBFileProperty[1], pdbLine)
            pdbFile.close()
            for tempPDBFileProperty in createPDBFileList:
                chainInfo = tempPDBFileProperty[0]
                chainExistenceStatus = 1
                for chainID in chainInfo:
                    if not chainID in pdbChainIDDict:
                        chainExistenceStatus = 0
                        break
                if chainExistenceStatus == 1:
                    tempPDBFile = open(tempPDBFileProperty[3], 'w')
                    tempPDBFile.write(tempPDBFileProperty[1])
                    tempPDBFile.close()
#                    tempNaccessErrorFile = tempfile.NamedTemporaryFile(mode='w+b')
#                    naccessRunString = '%s %s' %(naccessRunFileDirectory, tempPDBFileProperty[3])
#                    proc = subprocess.Popen([naccessRunFileDirectory, tempPDBFileProperty[3]], stdout=None, stderr=tempNaccessErrorFile, close_fds=True)
#                    ret_code = proc.wait()
#                    naccessErr = False
#                    if not ret_code == None:
#                        tempNaccessErrorFile.flush()
#                        tempNaccessErrorFile.seek(0)
#                        naccessErr = tempNaccessErrorFile.readlines()
#                        tempNaccessErrorFile.close()
#                    if naccessErr:
#                        naccessErrorQueue.put('Error: %s' %(tempPDBFileProperty[3]))
#                    try:
#                        ret_code =  subprocess.check_call([naccessRunFileDirectory, tempPDBFileProperty[3]])
#                    except:
#                        naccessErrorQueue.put('Error: %s' %(tempPDBFileProperty[3]))
                    stderr = socket.socketpair()
                    stderr[0].settimeout(0.01)
                    errMessage = ''
                    proc = subprocess.Popen([naccessRunFileDirectory, tempPDBFileProperty[3]], stdout=None, stderr=stderr[1], close_fds=True)
                    err = u''
                    while True:
                        proc.poll()
                        try:
                            errtmp = stderr[0].recv(4096)
                        except socket.timeout as exc:
                            errtmp = ''
                        err = err + errtmp
                        if len(err) > 4096:
                            proc.kill()
                            proc.wait()
                        if proc.returncode != None:
                            returnCode = proc.returncode
                            break
                    if err:
                        naccessErrorQueue.put('Error: %s, %s' %(tempPDBFileProperty[3], err))

                    rsaFileDirectory = '%s.rsa' %(tempPDBFileProperty[2])
                    asaFileDirectory = '%s.asa' %(tempPDBFileProperty[2])
                    logFileDirectory = '%s.log' %(tempPDBFileProperty[2])
                    
#                    os.system('%s %s > %s' %(naccessRunFileDirectory, tempPDBFileProperty[3], naccessRunOutputFileDirectory))
                    if os.path.exists(rsaFileDirectory):
                        os.system('mv %s %s' %(rsaFileDirectory, tempPDBFileProperty[4]))
                        if os.path.exists(asaFileDirectory):
                            os.system('mv %s %s' %(asaFileDirectory, tempPDBFileProperty[5]))
                    else:
                        if os.path.exists(asaFileDirectory):
                            os.system('rm %s' %(asaFileDirectory))
                    if os.path.exists(logFileDirectory):
                        os.system('rm %s' %(logFileDirectory))
                #    if os.path.exists(naccessRunOutputFileDirectory):
                #        os.system('rm %s' %(naccessRunOutputFileDirectory))
                    if os.path.exists(tempPDBFileProperty[3]):
                        os.system('rm %s' %(tempPDBFileProperty[3]))
            naccessSolutionUsageDictionary = {}
            for interfaceList in pdbDict['interfaces']:
                chain_1_rsaFileDirectory = interfaceList[1]
                chain_2_rsaFileDirectory = interfaceList[2]
                allChains_rsaFileDirectory = interfaceList[3]
                chain_1_rsaDict = {}
                chain_2_rsaDict = {}
                allChains_rsaDict = {}
                if not chain_1_rsaFileDirectory in naccessSolutionUsageDictionary:
                    naccessSolutionUsageDictionary[chain_1_rsaFileDirectory] = [0, interfaceList[4]]
                if not chain_2_rsaFileDirectory in naccessSolutionUsageDictionary:
                    naccessSolutionUsageDictionary[chain_2_rsaFileDirectory] = [0, interfaceList[5]]
                if not allChains_rsaFileDirectory in naccessSolutionUsageDictionary:
                    naccessSolutionUsageDictionary[allChains_rsaFileDirectory] = [0, interfaceList[6]]
                if os.path.exists(chain_1_rsaFileDirectory):
                    chain_1_rsaDict = naccessRSAFileReaderReturnOnlyAbsASADictionary(chain_1_rsaFileDirectory)
                    if len(chain_1_rsaDict) == 0:
                        interfaceNaccessResultText = '%s\t2\t%s does not exist\n' %(interfaceList[0], chain_1_rsaFileDirectory)
                        naccessResultsQueue.put(interfaceNaccessResultText)    
                        continue
                else:
                    interfaceNaccessResultText = '%s\t2\t%s does not exist\n' %(interfaceList[0], chain_1_rsaFileDirectory)
                    naccessResultsQueue.put(interfaceNaccessResultText)    
                    continue
                if os.path.exists(chain_2_rsaFileDirectory):
                    chain_2_rsaDict = naccessRSAFileReaderReturnOnlyAbsASADictionary(chain_2_rsaFileDirectory)
                    if len(chain_2_rsaDict) == 0:
                        interfaceNaccessResultText = '%s\t2\t%s does not exist\n' %(interfaceList[0], chain_2_rsaFileDirectory)
                        naccessResultsQueue.put(interfaceNaccessResultText)    
                        continue
                else:
                    interfaceNaccessResultText = '%s\t2\t%s does not exist\n' %(interfaceList[0], chain_2_rsaFileDirectory)
                    naccessResultsQueue.put(interfaceNaccessResultText)    
                    continue
                if os.path.exists(allChains_rsaFileDirectory):
                    allChains_rsaDict = naccessRSAFileReaderReturnOnlyAbsASADictionary(allChains_rsaFileDirectory)
                    if len(allChains_rsaDict) == 0:
                        interfaceNaccessResultText = '%s\t2\t%s does not exist\n' %(interfaceList[0], allChains_rsaFileDirectory)
                        naccessResultsQueue.put(interfaceNaccessResultText)    
                        continue
                else:
                    interfaceNaccessResultText = '%s\t2\t%s does not exist\n' %(interfaceList[0], allChains_rsaFileDirectory)
                    naccessResultsQueue.put(interfaceNaccessResultText)    
                    continue
                residueASAChangesCounter_chain_1 = 0
                residueASAChangesCounter_chain_2 = 0
                for resKey in allChains_rsaDict:
                    if resKey in chain_1_rsaDict:
                        if (float(chain_1_rsaDict[resKey]) - float(allChains_rsaDict[resKey])) > differenceBetweenComplexAndMonomerASACriteria:
                            residueASAChangesCounter_chain_1 = residueASAChangesCounter_chain_1 + 1
                    elif resKey in chain_2_rsaDict:
                        if (float(chain_2_rsaDict[resKey]) - float(allChains_rsaDict[resKey])) > differenceBetweenComplexAndMonomerASACriteria:
                            residueASAChangesCounter_chain_2 = residueASAChangesCounter_chain_2 + 1
                if residueASAChangesCounter_chain_1 >= minInterfaceResidueCriteria and residueASAChangesCounter_chain_2 >= minInterfaceResidueCriteria:
                    naccessSolutionUsageDictionary[chain_1_rsaFileDirectory][0] = 1
                    naccessSolutionUsageDictionary[chain_2_rsaFileDirectory][0] = 1
                    naccessSolutionUsageDictionary[allChains_rsaFileDirectory][0] = 1
                interfaceNaccessResultText = '%s\t0\t%d\t%d\n' %(interfaceList[0], residueASAChangesCounter_chain_1, residueASAChangesCounter_chain_2)
                naccessResultsQueue.put(interfaceNaccessResultText)    
            for naccessRSAFileDirectory in naccessSolutionUsageDictionary:
                if naccessSolutionUsageDictionary[naccessRSAFileDirectory][0] == 0:
                    if not naccessRSAFileDirectory in rsaFileKeeperDict:
                        if os.path.exists(naccessRSAFileDirectory):
                            os.system('rm %s' %(naccessRSAFileDirectory))
                        if os.path.exists(naccessSolutionUsageDictionary[naccessRSAFileDirectory][1]):
                            os.system('rm %s' %(naccessSolutionUsageDictionary[naccessRSAFileDirectory][1]))
                
        else:
            for interfaceList in pdbDict['interfaces']:
                interfaceNaccessResultText = '%s\t2\tPDB file does not present in %s.\n' %(interfaceList[0], pdbFileDirectory)
                naccessResultsQueue.put(interfaceNaccessResultText)    

def mainNaccessWithComplexPDBFileGenerator(interfaceListFile, allPDBFilesDirectory, generatedPDBFilesDirectory, fullPDBIDListFileDirectory, fullInterfaceNaccessResultsFileDirectory, interfaceNaccessResultsFileDirectory, naccessResultsFileDirectory, naccessRunFileDirectory, naccessLogFileDirectory, minInterfaceResidueCriteria, differenceBetweenComplexAndMonomerASACriteria, fullInterfaceListAfterNaccessFileDirectory, numberOfProcesses):
    
    print('\n* GENERATE COMPLEX PDB FILES STARTED *\n')
    print('Time stamp : %s' %(time.asctime()))
    t1 = time.time()

    if not os.path.exists(interfaceListFile):
        sys.exit('\nThe %s does not exist\n' %(interfaceListFile))

    if not os.path.exists(allPDBFilesDirectory):
        sys.exit('\nThe %s path is not correct.\n' %(allPDBFilesDirectory))

    if not os.path.exists(generatedPDBFilesDirectory):
        os.system('mkdir %s' %(generatedPDBFilesDirectory))

    if not os.path.exists(fullPDBIDListFileDirectory):
        sys.exit('\nThe %s does not exist\n' %(fullPDBIDListFileDirectory))

    if not os.path.exists(naccessResultsFileDirectory):
        os.system('mkdir %s' %(naccessResultsFileDirectory))

    if not os.path.exists(naccessRunFileDirectory):
        sys.exit('\nThe %s does not exist\n' %(naccessRunFileDirectory))



    fullNaccessResultDict = {}
    if os.path.exists(fullInterfaceNaccessResultsFileDirectory):
        fullInterfaceNaccessResultFile = open(fullInterfaceNaccessResultsFileDirectory, 'r')
        for naccessResultEntry in fullInterfaceNaccessResultFile:
            splittedNaccessResultEntry = naccessResultEntry.strip().split('\t')
            fullNaccessResultDict[splittedNaccessResultEntry[0]] = splittedNaccessResultEntry[1:]
        fullInterfaceNaccessResultFile.close()

    interfaceList = open(interfaceListFile,'r')
    numberOfCreatedPDBFile = 0
    existedPDBFiles = 0
    totalNumberOfInterface = 0
    taskDict = {}
    naccessResultsQueue = multiprocessing.Queue()
    rsaFileKeeperDict = {}
    for interface in interfaceList:
        interface = interface.strip()
        if interface == '':
            continue
        totalNumberOfInterface = totalNumberOfInterface + 1
        pdbName, chain_1, chain_2, allChains, interfaceResidueStatusOfChains, interface, commonChainCounter = interfaceNameSorter(interface)
        
        if len(chain_1) == 0:
            interfaceNaccessResultText = '%s\t2\tChain 1 does not have any monomer.\n' %(interface)
            naccessResultsQueue.put(interfaceNaccessResultText)    
            continue
        if len(chain_2) == 0:
            interfaceNaccessResultText = '%s\t2\tChain 2 does not have any monomer.\n' %(interface)
            naccessResultsQueue.put(interfaceNaccessResultText)    
            continue
        if commonChainCounter > 0:
            interfaceNaccessResultText = '%s\t2\tChain 1 and Chain 2 have %d common monomer.\n' %(interface, commonChainCounter)
            naccessResultsQueue.put(interfaceNaccessResultText)    
            continue
    
        chain_1_name = '%s_%s' %(pdbName, chain_1)
        chain_2_name = '%s_%s' %(pdbName, chain_2)
        allChainsName = '%s_%s' %(pdbName, allChains)

        pdbFileDirectoryWith_chain_1 = '%s/%s.pdb' %(generatedPDBFilesDirectory, chain_1_name)
        pdbFileDirectoryWith_chain_2 = '%s/%s.pdb' %(generatedPDBFilesDirectory, chain_2_name)
        pdbFileDirectoryWith_allChains = '%s/%s.pdb' %(generatedPDBFilesDirectory, allChainsName)
    
        naccessRSAFileDirectoryWith_chain_1 = '%s/%s.rsa' %(naccessResultsFileDirectory, chain_1_name)
        naccessRSAFileDirectoryWith_chain_2 = '%s/%s.rsa' %(naccessResultsFileDirectory, chain_2_name)
        naccessRSAFileDirectoryWith_allChains = '%s/%s.rsa' %(naccessResultsFileDirectory, allChainsName)

        naccessASAFileDirectoryWith_chain_1 = '%s/%s.asa' %(naccessResultsFileDirectory, chain_1_name)
        naccessASAFileDirectoryWith_chain_2 = '%s/%s.asa' %(naccessResultsFileDirectory, chain_2_name)
        naccessASAFileDirectoryWith_allChains = '%s/%s.asa' %(naccessResultsFileDirectory, allChainsName)

        if interface in fullNaccessResultDict:
            if int(fullNaccessResultDict[interface][0]) == 0:
                if int(fullNaccessResultDict[interface][1]) < minInterfaceResidueCriteria or int(fullNaccessResultDict[interface][2]) < minInterfaceResidueCriteria:
                    interfaceNaccessResultText = '%s\t%s\t%s\t%s\n' %(interface, fullNaccessResultDict[interface][0], fullNaccessResultDict[interface][1], fullNaccessResultDict[interface][2])
                    naccessResultsQueue.put(interfaceNaccessResultText)    
                    continue
                else:
                    if (os.path.exists(naccessRSAFileDirectoryWith_chain_1) and os.path.exists(naccessASAFileDirectoryWith_chain_1)):
                        if (os.path.exists(naccessRSAFileDirectoryWith_chain_2) and os.path.exists(naccessASAFileDirectoryWith_chain_2)):
                            if (os.path.exists(naccessRSAFileDirectoryWith_allChains) and os.path.exists(naccessASAFileDirectoryWith_allChains)):
                                interfaceNaccessResultText = '%s\t%s\t%s\t%s\n' %(interface, fullNaccessResultDict[interface][0], fullNaccessResultDict[interface][1], fullNaccessResultDict[interface][2])
                                naccessResultsQueue.put(interfaceNaccessResultText)
                                rsaFileKeeperDict[naccessRSAFileDirectoryWith_chain_1] = 1
                                rsaFileKeeperDict[naccessRSAFileDirectoryWith_chain_2] = 1
                                rsaFileKeeperDict[naccessRSAFileDirectoryWith_allChains] = 1
                                continue
                    
        if pdbName in taskDict:
            taskDict[pdbName]['interfaces'].append([interface, naccessRSAFileDirectoryWith_chain_1, naccessRSAFileDirectoryWith_chain_2, naccessRSAFileDirectoryWith_allChains, naccessASAFileDirectoryWith_chain_1, naccessASAFileDirectoryWith_chain_2, naccessASAFileDirectoryWith_allChains])
        else:
            taskDict[pdbName] = {}
            taskDict[pdbName]['chains'] = {}
            taskDict[pdbName]['interfaces'] = []
            taskDict[pdbName]['interfaces'].append([interface, naccessRSAFileDirectoryWith_chain_1, naccessRSAFileDirectoryWith_chain_2, naccessRSAFileDirectoryWith_allChains, naccessASAFileDirectoryWith_chain_1, naccessASAFileDirectoryWith_chain_2, naccessASAFileDirectoryWith_allChains])

        if not (os.path.exists(naccessRSAFileDirectoryWith_chain_1) and os.path.exists(naccessASAFileDirectoryWith_chain_1)):
            if chain_1 in taskDict[pdbName]['chains']:
                existedPDBFiles = existedPDBFiles + 1
            else:
                taskDict[pdbName]['chains'][chain_1] = [chain_1_name, pdbFileDirectoryWith_chain_1, naccessRSAFileDirectoryWith_chain_1, naccessASAFileDirectoryWith_chain_1]
                numberOfCreatedPDBFile = numberOfCreatedPDBFile + 1
        else:
            existedPDBFiles = existedPDBFiles + 1
        if not (os.path.exists(naccessRSAFileDirectoryWith_chain_2) and os.path.exists(naccessASAFileDirectoryWith_chain_2)):
            if chain_2 in taskDict[pdbName]['chains']:
                existedPDBFiles = existedPDBFiles + 1
            else:
                taskDict[pdbName]['chains'][chain_2] = [chain_2_name, pdbFileDirectoryWith_chain_2, naccessRSAFileDirectoryWith_chain_2, naccessASAFileDirectoryWith_chain_2]
                numberOfCreatedPDBFile = numberOfCreatedPDBFile + 1
        else:
            existedPDBFiles = existedPDBFiles + 1
        if not (os.path.exists(naccessRSAFileDirectoryWith_allChains) and os.path.exists(naccessASAFileDirectoryWith_allChains)):
            if allChains in taskDict[pdbName]['chains']:
                existedPDBFiles = existedPDBFiles + 1
            else:
                taskDict[pdbName]['chains'][allChains] = [allChainsName, pdbFileDirectoryWith_allChains, naccessRSAFileDirectoryWith_allChains, naccessASAFileDirectoryWith_allChains]
                numberOfCreatedPDBFile = numberOfCreatedPDBFile + 1
        else:
            existedPDBFiles = existedPDBFiles + 1
        
    interfaceList.close()
    
    taskQueue_pdbDict = multiprocessing.JoinableQueue()
    naccessErrorQueue = multiprocessing.Manager().Queue()
    generateComplexWorkers = [Process(target=generateComplexPDBFilesAndRunNaccessWork, args=(taskQueue_pdbDict, naccessResultsQueue, allPDBFilesDirectory, minInterfaceResidueCriteria, differenceBetweenComplexAndMonomerASACriteria, naccessRunFileDirectory, rsaFileKeeperDict, naccessErrorQueue)) for i in range(numberOfProcesses)]
    for tempWorkers in generateComplexWorkers:
        tempWorkers.start()

    interfaceTotal = 0
    for pdbName in taskDict:
        print('%s task added' %(pdbName))
        interfaceTotal = interfaceTotal + len(taskDict[pdbName]['interfaces'])
        taskQueue_pdbDict.put([pdbName, taskDict[pdbName]])
    print('total interfaces in taskDict = %d' %(interfaceTotal))

    for i in range(numberOfProcesses):
        taskQueue_pdbDict.put([None,None])

    taskQueue_pdbDict.join()
    print('taskQueue join oldu')
    interfaceNaccessResultsFile = open(interfaceNaccessResultsFileDirectory, 'w')
    interfaceNaccessResultDict = {}
    tempTotalNumberOfInterface = totalNumberOfInterface
    while tempTotalNumberOfInterface:
        interfaceNaccessResultString = naccessResultsQueue.get()
        splittedInterfaceNaccessResultString = interfaceNaccessResultString.strip().split('\t')
        interfaceNaccessResultDict[splittedInterfaceNaccessResultString[0]] = [0, splittedInterfaceNaccessResultString[1], interfaceNaccessResultString] 
        interfaceNaccessResultsFile.write(interfaceNaccessResultString)
        tempTotalNumberOfInterface = tempTotalNumberOfInterface - 1
    interfaceNaccessResultsFile.close()

    naccessErrorQueue.put('Done')
    naccessLogFile = open(naccessLogFileDirectory, 'w')
    while True:
        naccessErrorEntry = naccessErrorQueue.get()
        if naccessErrorEntry == 'Done':
            break
        naccessLogFile.write('%s\n' %(naccessErrorEntry))
    naccessLogFile.close()

    fullPDBIDDict = {}
    fullPDBIDListFile = open(fullPDBIDListFileDirectory, 'r')
    for pdbID in fullPDBIDListFile:
        pdbID = pdbID.strip()
        fullPDBIDDict[pdbID] = 1
    fullPDBIDListFile.close()
    tempFullInterfaceNaccessResultsFileDirectory = '%s_temp' %(fullInterfaceNaccessResultsFileDirectory)
    tempFullInterfaceNaccessResultsFile = open(tempFullInterfaceNaccessResultsFileDirectory,'w')
    tempFullInterfaceListAfterNaccessFileDirectory = '%s_temp' %(fullInterfaceListAfterNaccessFileDirectory)
    tempFullInterfaceListAfterNaccessFile = open(tempFullInterfaceListAfterNaccessFileDirectory, 'w')
    if os.path.exists(fullInterfaceNaccessResultsFileDirectory):
        fullInterfaceNaccessResultsFile = open(fullInterfaceNaccessResultsFileDirectory,'r')
        for naccessResultLine in fullInterfaceNaccessResultsFile:
            splittedNaccessResultLine = naccessResultLine.strip().split('\t')
            interfaceName = splittedNaccessResultLine[0]
            splittedInterfaceName = interfaceName.split('_')
            pdbID = splittedInterfaceName[0]
            if pdbID in fullPDBIDDict:
                if interfaceName in interfaceNaccessResultDict:
                    interfaceNaccessResultDict[interfaceName][0] = 1
                    if splittedNaccessResultLine[1] > interfaceNaccessResultDict[interfaceName][1]:
                        tempFullInterfaceNaccessResultsFile.write(interfaceNaccessResultDict[interfaceName][2])
                        if int(interfaceNaccessResultDict[interfaceName][1]) == 0:
                            splittedInterfaceNaccessResultString = interfaceNaccessResultDict[interfacename][2].strip().split('\t')
                            if int(splittedInterfaceNaccessResultString[2]) >= minInterfaceResidueCriteria and int(splittedInterfaceNaccessResultString[3]) >= minInterfaceResidueCriteria:
                                tempFullInterfaceListAfterNaccessFile.write('%s\n' %(interfaceName))
                    else:
                        tempFullInterfaceNaccessResultsFile.write(naccessResultLine)
                        if int(splittedNaccessResultLine[1]) == 0:
                            if int(splittedNaccessResultLine[2]) >= minInterfaceResidueCriteria and int(splittedNaccessResultLine[3]) >= minInterfaceResidueCriteria:
                                tempFullInterfaceListAfterNaccessFile.write('%s\n' %(interfaceName))
                else:
                    tempFullInterfaceNaccessResultsFile.write(naccessResultLine)
                    if int(splittedNaccessResultLine[1]) == 0:
                        if int(splittedNaccessResultLine[2]) >= minInterfaceResidueCriteria and int(splittedNaccessResultLine[3]) >= minInterfaceResidueCriteria:
                            tempFullInterfaceListAfterNaccessFile.write('%s\n' %(interfaceName))
        fullInterfaceNaccessResultsFile.close()
        for interfaceName in interfaceNaccessResultDict:
            if interfaceNaccessResultDict[interfaceName][0] == 0:
                tempFullInterfaceNaccessResultsFile.write(interfaceNaccessResultDict[interfaceName][2])
                if int(interfaceNaccessResultDict[interfaceName][1]) == 0:
                    splittedInterfaceNaccessResultString = interfaceNaccessResultDict[interfaceName][2].strip().split('\t')
                    if int(splittedInterfaceNaccessResultString[2]) >= minInterfaceResidueCriteria and int(splittedInterfaceNaccessResultString[3]) >= minInterfaceResidueCriteria:
                        tempFullInterfaceListAfterNaccessFile.write('%s\n' %(interfaceName))
    else:
        for interfaceName in interfaceNaccessResultDict:
            tempFullInterfaceNaccessResultsFile.write(interfaceNaccessResultDict[interfaceName][2])
            if int(interfaceNaccessResultDict[interfaceName][1]) == 0:
                splittedInterfaceNaccessResultString = interfaceNaccessResultDict[interfacename][2].strip().split('\t')
                if int(splittedInterfaceNaccessResultString[2]) >= minInterfaceResidueCriteria and int(splittedInterfaceNaccessResultString[3]) >= minInterfaceResidueCriteria:
                    tempFullInterfaceListAfterNaccessFile.write('%s\n' %(interfaceName))
    tempFullInterfaceNaccessResultsFile.close()
    tempFullInterfaceListAfterNaccessFile.close()
    os.system('mv %s %s' %(tempFullInterfaceNaccessResultsFileDirectory, fullInterfaceNaccessResultsFileDirectory))            
    os.system('mv %s %s' %(tempFullInterfaceListAfterNaccessFileDirectory, fullInterfaceListAfterNaccessFileDirectory))            


    t2 = time.time()
    print('\nTotal number of interface = %d\n' %(totalNumberOfInterface))
    print('\nExisted PDB Files = %d\n' %(existedPDBFiles))
    print('\nNumber of created PDB File = %d\n' %(numberOfCreatedPDBFile))
    print('\nElapsed time = %f seconds\n' %(t2-t1))    
    print('Time stamp : %s' %(time.asctime()))
    print('\n* GENERATE COMPLEX PDB FILES COMPLETED *\n')
        
