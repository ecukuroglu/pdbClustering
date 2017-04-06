#################
## Created by Engin Cukuroglu
#################

import os,sys,time
import multiprocessing
from multiprocessing import Queue, Process
import subprocess
import socket
from codesOfTools import multiprotSolutionReader, vdwInterfaceDictionaryCreator

def multiprotRunListWork(taskQueue_multiprotRunList, multiprotRunResultQueue, multiprotRunFileDirectory, multiprotRunForProcessesFileDirectory, multiprotInterfaceListFileDirectory, multiprotListSuffix, multiprotResultsForInterfaceFilesDirectory, multiprotResultsSuffix, interfaceFilesDirectory, interfaceSkeletonFilesDirectory, multiprotLogFileDirectory, i):
    processWorkingDirectory = '%s/process_%s' %(multiprotRunForProcessesFileDirectory, i)
    if not os.path.exists(processWorkingDirectory):
        os.system('mkdir %s' %(processWorkingDirectory))
    os.chdir(processWorkingDirectory)
    while True:
        interface, interfaceIndex, interfaceSize, skeletonSize, numberOfMultiprotRuns = taskQueue_multiprotRunList.get()
        taskQueue_multiprotRunList.task_done()
        if interface == None:
            break
        if numberOfMultiprotRuns > 0:
            interfaceMultiprotResultsFileDirectory = '../../%s/%s_%s' %(multiprotResultsForInterfaceFilesDirectory, interface, multiprotResultsSuffix)
            tempInterfaceMultiprotResultsFileDirectory = '../../%s/%s_%s_temp' %(multiprotResultsForInterfaceFilesDirectory, interface, multiprotResultsSuffix)
            interfaceMultiprotListFileDirectory = '../../%s/%s_%s' %(multiprotInterfaceListFileDirectory, interface, multiprotListSuffix)
            ##### start changes for yunus cluster ##### 
            interfaceMultiprotListFileDirectory = '%s/%s_%s' %(multiprotInterfaceListFileDirectory, interface, multiprotListSuffix)
            ##### end changes for yunus cluster ##### 
            if os.path.exists(tempInterfaceMultiprotResultsFileDirectory):
                os.system('rm %s' %(tempInterfaceMultiprotResultsFileDirectory))
            interfaceMultiprotResultsDict = {}
            if os.path.exists(interfaceMultiprotResultsFileDirectory):
                interfaceMultiprotResultsFile = open(interfaceMultiprotResultsFileDirectory, 'r')
                tempInterfaceMultiprotResultsFile = open(tempInterfaceMultiprotResultsFileDirectory, 'w')
                for interfaceMultiprotResults in interfaceMultiprotResultsFile:
                    tempInterfaceMultiprotResultsFile.write(interfaceMultiprotResults)    
                    splittedInterfaceMultiprotResults = interfaceMultiprotResults.strip().split('\t')
                    interfaceMultiprotResultsDict[splittedInterfaceMultiprotResults[1]] = 1
                interfaceMultiprotResultsFile.close()
                tempInterfaceMultiprotResultsFile.close()
            if os.path.exists(interfaceMultiprotListFileDirectory):
                interfaceMultiprotListFile = open(interfaceMultiprotListFileDirectory, 'r')
                tempInterfaceMultiprotResultsFile = open(tempInterfaceMultiprotResultsFileDirectory, 'a')
                multiprotErrorString = ''
                multiprotStoredRunsCounter = 0
                multiprotSuccessfulRunsCounter = 0
                multiprotUnsuccessfulRunsCounter = 0
                for interfaceMultiprotEntry in interfaceMultiprotListFile:
                    splittedInterfaceMultiprotEntry = interfaceMultiprotEntry.strip().split('\t')
                    tempInterface = splittedInterfaceMultiprotEntry[1]
                    tempInterfaceSize = splittedInterfaceMultiprotEntry[3]
                    tempSkeletonSize = splittedInterfaceMultiprotEntry[5]
                    if not tempInterface in interfaceMultiprotResultsDict:
                        if os.path.exists('2_sol.res'):
                            os.system('rm 2_sol.res')
                        if os.path.exists('2_sets.res'):
                            os.system('rm 2_sets.res')
                        if os.path.exists('log_multiprot.txt'):
                            os.system('rm log_multiprot.txt')

                        stderr = socket.socketpair()
                        stderr[0].settimeout(0.01)
                        stdout = socket.socketpair()
                        stdout[0].settimeout(0.01)
                        errMessage = ''
                        multiprotRunDirectory_forProcess = '../../%s' %(multiprotRunFileDirectory)
                        skeleton_1_pdbDirectory = '../../%s/%s.pdb' %(interfaceSkeletonFilesDirectory, interface)
                        skeleton_2_pdbDirectory = '../../%s/%s.pdb' %(interfaceSkeletonFilesDirectory, tempInterface)
                        proc = subprocess.Popen([multiprotRunDirectory_forProcess, skeleton_1_pdbDirectory, skeleton_2_pdbDirectory], stdout=stdout[1], stderr=stderr[1], close_fds=True)
                        err = u''
                        while True:
                            proc.poll()
                            try:
                                errtmp = stderr[0].recv(4096)
                            except socket.timeout as exc:
                                errtmp = ''
                            if len(err) > 4096:
                                proc.kill()
                                proc.wait()
                            if proc.returncode != None:
                                returnCode = proc.returncode
                                break
                        if err:
                            multiprotErrorString = '%sError: %s and %s.\t%s\n' %(multiprotErrorString, interface, tempInterface, err.strip())
                            multiprotUnsuccessfulRunsCounter = multiprotUnsuccessfulRunsCounter + 1
                        else:
                            if os.path.exists('2_sol.res'):
                                alignmentDict, referenceMonomer, alignedMonomer, rmsd, transMatrix = multiprotSolutionReader('2_sol.res', 0)
                                interface_1_pdbDirectory = '../../%s/%s.pdb' %(interfaceFilesDirectory, interface)
                                interface_2_pdbDirectory = '../../%s/%s.pdb' %(interfaceFilesDirectory, tempInterface)
                                interface_1_dict = vdwInterfaceDictionaryCreator(interface_1_pdbDirectory)
                                interface_2_dict = vdwInterfaceDictionaryCreator(interface_2_pdbDirectory)
                                interfaceMatchingSize = 0
                                skeletonMatchingSize = 0
                                if referenceMonomer == interface:
                                    for referenceDictKey in alignmentDict:
                                        skeletonMatchingSize = skeletonMatchingSize + 1
                                        if referenceDictKey in interface_1_dict:
                                            if alignmentDict[referenceDictKey] in interface_2_dict:
                                                interfaceMatchingSize = interfaceMatchingSize + 1
                                else:
                                    for referenceDictKey in alignmentDict:
                                        skeletonMatchingSize = skeletonMatchingSize + 1
                                        if referenceDictKey in interface_2_dict:
                                            if alignmentDict[referenceDictKey] in interface_1_dict:
                                                interfaceMatchingSize = interfaceMatchingSize + 1
                                if interfaceSize < tempInterfaceSize:
                                    minInterfaceSize = interfaceSize
                                else:
                                    minInterfaceSize = tempInterfaceSize
                                if skeletonSize < tempSkeletonSize:
                                    minSkeletonSize = skeletonSize
                                else:
                                    minSkeletonSize = tempSkeletonSize
                                interfaceMatchingRatio = float(interfaceMatchingSize) / minInterfaceSize
                                skeletonMatchingRatio = float(skeletonMatchingSize) / minSkeletonSize
                                tempInterfaceMultiprotResultsFile.write('%s\t%d\t%d\t%.2f\t%.2f\t%.2f\t%s\n' %(interfaceMultiprotEntry.strip(), interfaceMatchingSize, skeletonMatchingSize, interfaceMatchingRatio, skeletonMatchingRatio, rmsd, transMatrix))
                                multiprotSuccessfulRunsCounter = multiprotSuccessfulRunsCounter + 1
                            else:
                                multiprotUnsuccessfulRunsCounter = multiprotUnsuccessfulRunsCounter + 1
                        if os.path.exists('2_sol.res'):
                            os.system('rm 2_sol.res')
                        if os.path.exists('2_sets.res'):
                            os.system('rm 2_sets.res')
                        if os.path.exists('log_multiprot.txt'):
                            os.system('rm log_multiprot.txt')
                    else:
                        multiprotStoredRunsCounter = multiprotStoredRunsCounter + 1
    
                interfaceMultiprotListFile.close()
                tempInterfaceMultiprotResultsFile.close()
                os.system('mv %s %s' %(tempInterfaceMultiprotResultsFileDirectory, interfaceMultiprotResultsFileDirectory))
                if not multiprotErrorString == '':
                    interfaceMultiprotLogFileDirectory = '../../%s/%s_multiprotErrorLogs.txt' %(multiprotLogFileDirectory, interface)
                    interfaceMultiprotLogFile = open(interfaceMultiprotLogFileDirectory, 'w')
                    interfaceMultiprotLogFile.write(multiprotErrorString)
                multiprotRunResultString = '%s\t0\t%d\t%d\t%d\t%d\t%d\n' %(interface, interfaceIndex, numberOfMultiprotRuns, multiprotStoredRunsCounter, multiprotSuccessfulRunsCounter, multiprotUnsuccessfulRunsCounter)
                multiprotRunResultQueue.put(multiprotRunResultString)
            else:
                multiprotRunResultString = '%s\t2\t%d\t%d\tError: The %s does not exist\n' %(interface, interfaceIndex, numberOfMultiprotRuns, interfaceMultiprotListFileDirectory)
                multiprotRunResultQueue.put(multiprotRunResultString)
            
        else:
            multiprotRunResultString = '%s\t0\t%d\t%d\t0\t0\t0\n' %(interface, interfaceIndex, numberOfMultiprotRuns)
            multiprotRunResultQueue.put(multiprotRunResultString)
    

def mainMultiprotDriver(fullMultiprotStatisticsFileDirectory, interfaceStartIndex, interfaceEndIndex, multiprotRunFileDirectory, multiprotRunForProcessesFileDirectory, multiprotInterfaceListFileDirectory, multiprotListSuffix, multiprotResultsForInterfaceFilesDirectory, multiprotResultsSuffix, multiprotRunListLogFileDirectory, multiprotRunListResultFileDirectory, interfaceFilesDirectory, interfaceSkeletonFilesDirectory, multiprotLogFileDirectory, numberOfProcesses):
    
    print('\n* MULTIPROT DRIVER STARTED *\n')
    print('Time stamp: %s' %(time.asctime()))
    t1 = time.time()

    if not os.path.exists(fullMultiprotStatisticsFileDirectory):
        sys.exit('\nThe %s does not exist.\n' %(fullMultiprotStatisticsFileDirectory))

    if not os.path.exists(multiprotResultsForInterfaceFilesDirectory):
        os.system('mkdir %s' %(multiprotResultsForInterfaceFilesDirectory))

    if os.path.exists(multiprotRunForProcessesFileDirectory):
        os.system('rm -r %s' %(multiprotRunForProcessesFileDirectory))
        os.system('mkdir %s' %(multiprotRunForProcessesFileDirectory))
    else:    
        os.system('mkdir %s' %(multiprotRunForProcessesFileDirectory))
    
    if not os.path.exists(multiprotLogFileDirectory):
        os.system('mkdir %s' %(multiprotLogFileDirectory))    
    taskList = []
    taskDict = {}
    totalNumberOfMultiprotRuns = 0
    totalNumberOfInterfaces = 0
    fullMultiprotStatisticsFile = open(fullMultiprotStatisticsFileDirectory, 'r')
    multiprotRunListLogFile = open(multiprotRunListLogFileDirectory, 'w')
    for interfaceMultiprotInfo in fullMultiprotStatisticsFile:
        splittedInterfaceMultiprotInfo = interfaceMultiprotInfo.strip().split('\t')
        interfaceIndex = int(splittedInterfaceMultiprotInfo[1])
        interfaceStatus = 0
        if interfaceStartIndex < 0:
            if interfaceEndIndex < 0:
                interfaceStatus = 1
            elif interfaceIndex <= interfaceEndIndex:
                interfaceStatus = 1
        elif interfaceIndex >= interfaceStartIndex:
            if interfaceEndIndex < 0:
                interfaceStatus = 1
            elif interfaceIndex <= interfaceEndIndex:
                interfaceStatus = 1
        if interfaceStatus == 1:
            interface = splittedInterfaceMultiprotInfo[0]
            if not interface in taskDict:
                taskDict[interface] = 1
                totalNumberOfMultiprotRuns = totalNumberOfMultiprotRuns + int(splittedInterfaceMultiprotInfo[4])
                totalNumberOfInterfaces = totalNumberOfInterfaces + 1
                taskList.append([interface, interfaceIndex, int(splittedInterfaceMultiprotInfo[2]), int(splittedInterfaceMultiprotInfo[3]), int(splittedInterfaceMultiprotInfo[4])])
                multiprotRunListLogFile.write(interfaceMultiprotInfo)
    
    fullMultiprotStatisticsFile.close()
    multiprotRunListLogFile.write('---\nTotal Number of Interfaces: %d\n' %(totalNumberOfInterfaces))
    multiprotRunListLogFile.write('Total Number of Multiprot Runs: %d\n---\n' %(totalNumberOfMultiprotRuns))
    multiprotRunListLogFile.close()

    taskQueue_multiprotRunList = multiprocessing.JoinableQueue()
    multiprotRunResultQueue = multiprocessing.Queue()
    generateMultiprotRunListWorkers = [Process(target=multiprotRunListWork, args=(taskQueue_multiprotRunList, multiprotRunResultQueue, multiprotRunFileDirectory, multiprotRunForProcessesFileDirectory, multiprotInterfaceListFileDirectory, multiprotListSuffix, multiprotResultsForInterfaceFilesDirectory, multiprotResultsSuffix, interfaceFilesDirectory, interfaceSkeletonFilesDirectory, multiprotLogFileDirectory, i)) for i in range(numberOfProcesses)]

    for tempWorkers in generateMultiprotRunListWorkers:
        tempWorkers.start()

    for interfaceInfoList in taskList:
        taskQueue_multiprotRunList.put(interfaceInfoList)    

    for i in range(numberOfProcesses):
        taskQueue_multiprotRunList.put([None, None, None, None, None])    

    taskQueue_multiprotRunList.join()
    totalInterfaceResults = totalNumberOfInterfaces
    multiprotRunListResultFile = open(multiprotRunListResultFileDirectory, 'w')
    while totalInterfaceResults:
        interfaceMultiprotResultString = multiprotRunResultQueue.get()
        multiprotRunListResultFile.write(interfaceMultiprotResultString)
        totalInterfaceResults = totalInterfaceResults - 1
    
    t2 = time.time()
    print('\nElapsed time = %f seconds\n' %(t2-t1))
    print('Time stamp: %s' %(time.asctime()))
    print('\n* MULTIPROT DRIVER COMPLETED *\n')
