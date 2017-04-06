#################
## Created by Engin Cukuroglu
#################

import os,sys,time
import multiprocessing
from operator import itemgetter
from multiprocessing import Queue, Process

def multiprotInterfaceListWork(taskQueue_multiprotInterfaceList, multiprotInterfaceResultQueue, interfaceList, interfaceListLength, multiprotInterfaceListFileDirectory, multiprotListSuffix, multiprotResultsForInterfaceFilesDirectory, multiprotResultsSuffix, maximumToMinimumInterfaceExcessRatio, maximumToMinimumSkeletonExcessRatio):
    while True:
        interface, interfaceNumber, interfaceSize, skeletonSize = taskQueue_multiprotInterfaceList.get()
        taskQueue_multiprotInterfaceList.task_done()
        if interface == None:
            break
        interfaceMultiprotResultsDict= {}
        interfaceMultiprotResultsFileDirectory = '%s/%s_%s' %(multiprotResultsForInterfaceFilesDirectory, interface, multiprotResultsSuffix)
        if os.path.exists(interfaceMultiprotResultsFileDirectory):
            interfaceMultiprotResultsFile = open(interfaceMultiprotResultsFileDirectory, 'r')
            for interfaceMultiprotResults in interfaceMultiprotResultsFile:
                splittedInterfaceMultiprotResults = interfaceMultiprotResults.strip().split('\t')
                interfaceMultiprotResultsDict[splittedInterfaceMultiprotResults[1]] = 1
            interfaceMultiprotResultsFile.close()
            
        interfaceMultiprotListFileDirectory = '%s/%s_%s' %(multiprotInterfaceListFileDirectory, interface, multiprotListSuffix)
        interfaceMultiprotListFile = open(interfaceMultiprotListFileDirectory, 'w')
        numberOfMultiprotComparison = 0
        numberOfStoredMultiprotComparison = 0
        for i in xrange(interfaceNumber+1, interfaceListLength):
            tempInterfaceInfo = interfaceList[i]
            tempInterfaceName = tempInterfaceInfo[0]
            if tempInterfaceName in interfaceMultiprotResultsDict:
                numberOfStoredMultiprotComparison = numberOfStoredMultiprotComparison + 1
                continue
            tempInterfaceSize = float(tempInterfaceInfo[1])
            tempSkeletonSize = float(tempInterfaceInfo[2])
            if interfaceSize < tempInterfaceSize:
                minInterfaceSize = interfaceSize
                maxInterfaceSize = tempInterfaceSize
            else:    
                minInterfaceSize = tempInterfaceSize
                maxInterfaceSize = interfaceSize
            if skeletonSize < tempSkeletonSize:
                minSkeletonSize = skeletonSize
                maxSkeletonSize = tempSkeletonSize
            else:
                minSkeletonSize = tempSkeletonSize
                maxSkeletonSize = skeletonSize
            interfaceRatio = (maxInterfaceSize/minInterfaceSize) - 1
            skeletonRatio = (maxSkeletonSize/minSkeletonSize) - 1
            if interfaceRatio <= maximumToMinimumInterfaceExcessRatio:
                if skeletonRatio <= maximumToMinimumSkeletonExcessRatio:
                    interfaceMultiprotListFile.write('%s\t%s\t%.0f\t%.0f\t%.0f\t%.0f\n' %(interface, tempInterfaceName, interfaceSize, tempInterfaceSize, skeletonSize, tempSkeletonSize))
                    numberOfMultiprotComparison = numberOfMultiprotComparison + 1
        interfaceMultiprotListFile.close()
        interfaceMultiprotResultText = '%s\t%d\t%.0f\t%.0f\t%d\t%d\n' %(interface, interfaceNumber, interfaceSize, skeletonSize, numberOfMultiprotComparison, numberOfStoredMultiprotComparison)
        multiprotInterfaceResultQueue.put([interfaceMultiprotResultText, numberOfMultiprotComparison])

def mainMultiprotInterfaceListExtractor(interfaceTypeListFileDirectory, minInterfaceResidueCriteria, maximumToMinimumInterfaceExcessRatio, maximumToMinimumSkeletonExcessRatio, fullInterfaceTypeListForMultiprotFileDirectory, multiprotInterfaceListFileDirectory, multiprotListSuffix, multiprotResultsForInterfaceFilesDirectory, multiprotResultsSuffix, fullMultiprotStatisticsFileDirectory,  numberOfProcesses):

    print('\n* MULTIPROT INTERFACE LIST EXTRACTOR STARTED *\n')
    print('Time stamp: %s' %(time.asctime()))
    t1 = time.time()

    if not os.path.exists(interfaceTypeListFileDirectory):
        sys.exit('\nThe %s does not exist\n' %(interfaceTypeListFileDirectory))
    if os.path.exists(multiprotInterfaceListFileDirectory):
        os.system('rm -r %s' %(multiprotInterfaceListFileDirectory))
        os.system('mkdir %s' %(multiprotInterfaceListFileDirectory))
    else:
        os.system('mkdir %s' %(multiprotInterfaceListFileDirectory))

    interfaceTypeListFile = open(interfaceTypeListFileDirectory, 'r')
    fullInterfaceTypeListForMultiprotFile = open(fullInterfaceTypeListForMultiprotFileDirectory, 'w')
    interfaceList = []
    interfaceDict = {}
    for interfaceProperty in interfaceTypeListFile:
        splittedInterfaceProperty = interfaceProperty.strip().split('\t')
        interfaceName = splittedInterfaceProperty[0]
        if int(splittedInterfaceProperty[1]) == 0:
            if not interfaceName in interfaceDict:
                interfaceStatus = 1
                chain_1_interfaceResidueList = splittedInterfaceProperty[8].split(',')
                for interfaceResidueNumber in chain_1_interfaceResidueList:
                    if int(interfaceResidueNumber) < minInterfaceResidueCriteria:
                        interfaceStatus = 0
                        continue
                if interfaceStatus == 1:
                    chain_2_interfaceResidueList = splittedInterfaceProperty[9].split(',')
                    for interfaceResidueNumber in chain_2_interfaceResidueList:
                        if int(interfaceResidueNumber) < minInterfaceResidueCriteria:
                            interfaceStatus = 0
                            continue
                if interfaceStatus == 1:
                    interfaceList.append([interfaceName, int(splittedInterfaceProperty[12])+int(splittedInterfaceProperty[13]), int(splittedInterfaceProperty[14])+int(splittedInterfaceProperty[15])])
                    fullInterfaceTypeListForMultiprotFile.write(interfaceProperty)
                    interfaceDict[interfaceName] = 1
    interfaceTypeListFile.close()
    fullInterfaceTypeListForMultiprotFile.close()
    interfaceList = sorted(interfaceList, key=itemgetter(0))
    interfaceListLength = len(interfaceList)
    multiprotInterfaceResultQueue = multiprocessing.Queue()
    taskQueue_multiprotInterfaceList = multiprocessing.JoinableQueue()
    
    generateMultiprotInterfaceListWorkers = [Process(target=multiprotInterfaceListWork, args=(taskQueue_multiprotInterfaceList, multiprotInterfaceResultQueue, interfaceList, interfaceListLength, multiprotInterfaceListFileDirectory, multiprotListSuffix, multiprotResultsForInterfaceFilesDirectory, multiprotResultsSuffix, maximumToMinimumInterfaceExcessRatio, maximumToMinimumSkeletonExcessRatio)) for i in range(numberOfProcesses)]

    for tempWorkers in generateMultiprotInterfaceListWorkers:
        tempWorkers.start()
    
    interfaceCounter = 0
    for interfaceInfo in interfaceList:
        taskQueue_multiprotInterfaceList.put([interfaceInfo[0], interfaceCounter, float(interfaceInfo[1]), float(interfaceInfo[2])])
        interfaceCounter = interfaceCounter + 1

    for i in range(numberOfProcesses):
        taskQueue_multiprotInterfaceList.put([None, None, None, None])
        
    taskQueue_multiprotInterfaceList.join()
    totalMultiprotRunNumber = 0
    fullMultiprotStatisticsFile = open(fullMultiprotStatisticsFileDirectory, 'w')
    while interfaceCounter:
        multiprotInterfaceList = multiprotInterfaceResultQueue.get()
        fullMultiprotStatisticsFile.write(multiprotInterfaceList[0])
        totalMultiprotRunNumber = totalMultiprotRunNumber + multiprotInterfaceList[1]
        interfaceCounter = interfaceCounter - 1
    fullMultiprotStatisticsFile.close()
            
    t2 = time.time()
    print('Total number of multiprot runs = %d\n' %(totalMultiprotRunNumber))
    print('\nElapsed time = %f seconds\n' %(t2-t1))
    print('Time stamp: %s' %(time.asctime()))
    print('\n* MULTIPROT INTERFACE LIST EXTRACTOR COMPLETED *\n')
