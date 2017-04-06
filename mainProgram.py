#################
## Created by Engin Cukuroglu
#################
##
##usage: nohup python pythonCodes/mainProgram.py & 
##
#################

import os, sys, time
from pdbSync import mainPDBSync
from pdbUpdater import mainUpdater
from possibleInterfaceExtractor import mainPossibleInterfaceExtractor
from naccessWithComplexPDBFileGenerator import mainNaccessWithComplexPDBFileGenerator
from interfaceExtractorUsingVDW import mainInterfaceExtractorUsingVDW
from interfaceASAValueFinder import mainInterfaceASAValueFinder
from multiprotInterfaceListExtractor import mainMultiprotInterfaceListExtractor
from multiprotDriver import mainMultiprotDriver

if 1:
    #### start yunus cluster run strategy ####
    multiprotInterfaceListFileDirectory = '/mnt/kufs/scratch/ecukuroglu/project/multiprotInterfaceLists'
    pdbSyncTarFileDirectory = 'PDBSync.tar.gz'
    currentDirectory = os.getcwd()
    ownerDirectory = '/state/partition1/ecukuroglu'
    if not os.path.exists(ownerDirectory):
        os.system('mkdir %s' %(ownerDirectory))
    os.system('cp %s %s' %(pdbSyncTarFileDirectory, ownerDirectory))
    os.chdir(ownerDirectory)
    os.system('tar -xzf %s' %(pdbSyncTarFileDirectory))
    os.system('rm %s' %(pdbSyncTarFileDirectory))
    workingDirectory = '/state/partition1/ecukuroglu/PDBSync'
    if not os.path.exists(workingDirectory):
        sys.exit('\nThe %s does not exist.\n' %(workingDirectory))
    os.chdir(workingDirectory)
    #### end yunus cluster run strategy ####

    print('\n*** MAIN PROGRAM STARTED ***\n')
    print('Time stamp : %s' %(time.asctime()))
    t1 = time.time()
    timeStamp = '%s' %(time.strftime('%Y_%B_%d'))
    numberOfProcesses = 24
    archivedFilesDirectory = 'archivedFiles'
    fullListFileDirectory = 'fullLists'
    allPDBFilesDirectory = 'pdbFiles'
    mainLogFileDirectory = 'logFiles'
    logFileDirectory = '%s/log_%s' %(mainLogFileDirectory, timeStamp)
    archivedPDBFileDirectory = 'archivedPDBFiles'

    if not os.path.exists(allPDBFilesDirectory):
        os.system('mkdir %s' %(allPDBFilesDirectory))

    if not os.path.exists(archivedPDBFileDirectory):
        os.system('mkdir %s' %(archivedPDBFileDirectory))

    if not os.path.exists(archivedFilesDirectory):
        os.system('mkdir %s' %(archivedFilesDirectory))

    if not os.path.exists(fullListFileDirectory):
        os.system('mkdir %s' %(fullListFileDirectory))

    if not os.path.exists(mainLogFileDirectory):
        os.system('mkdir %s' %(mainLogFileDirectory))

    if not os.path.exists(logFileDirectory):
        os.system('mkdir %s' %(logFileDirectory))

    ### PDB check new pdb files ###
    archivedPDBIDListFileDirectory = '%s/pdbIDList_%s.txt' %(archivedFilesDirectory, timeStamp)
    archivedPDBIDEntriesFileDirectory = '%s/pdbIDEntries_%s.txt' %(archivedFilesDirectory, timeStamp)
#    mainPDBSync(archivedPDBIDListFileDirectory, archivedPDBIDEntriesFileDirectory)

    ### download missing pdb files from PDB database ###
    fullPDBIDListFileDirectory = '%s/fullPDBIDList.txt' %(fullListFileDirectory)
    pdbDownloadLogFileDirectory = '%s/pdbDownloadLog_%s.txt' %(logFileDirectory, timeStamp)
#    mainUpdater(allPDBFilesDirectory, archivedPDBIDListFileDirectory, fullPDBIDListFileDirectory, archivedPDBFileDirectory, pdbDownloadLogFileDirectory)

    ### extract possible interfaces ###
    possibleInterfaceListFile = '%s/possibleInterfaceList_%s.txt' %(archivedFilesDirectory, timeStamp)
    fullPossibleInterfaceListFile = '%s/fullPossibleInterfaceList.txt' %(fullListFileDirectory)
    monomerListFile = '%s/monomerList_%s.txt' %(archivedFilesDirectory, timeStamp)
    fullMonomerListFile = '%s/fullMonomerList.txt' %(fullListFileDirectory)
    complexListFile = '%s/complexList_%s.txt' %(archivedFilesDirectory, timeStamp)
    fullComplexListFile = '%s/fullComplexList.txt' %(fullListFileDirectory)
    pdbChainListFile = '%s/pdbChainList_%s.txt' %(archivedFilesDirectory, timeStamp)
    fullPDBChainListFile = '%s/fullPDBChainList.txt' %(fullListFileDirectory)
    interfaceResidueStatusOfChains = 0
    possibleInterfaceLogFileDirectory = '%s/possibleInterfaceLog_%s.txt' %(logFileDirectory, timeStamp)
#    mainPossibleInterfaceExtractor(allPDBFilesDirectory, possibleInterfaceListFile, monomerListFile, complexListFile, pdbChainListFile, fullPDBIDListFileDirectory, fullPossibleInterfaceListFile, fullMonomerListFile, fullComplexListFile, fullPDBChainListFile, interfaceResidueStatusOfChains, possibleInterfaceLogFileDirectory)

    ### run naccess for accessible surface area calculations ###    
    generatedPDBFilesDirectory = 'generatedPDBFiles'
    interfaceNaccessResultsFileDirectory = '%s/interfaceNaccessResults_%s.txt' %(archivedFilesDirectory, timeStamp)
    fullInterfaceNaccessResultsFileDirectory = '%s/fullInterfaceNaccessResults.txt' %(fullListFileDirectory)
    naccessResultsFileDirectory = 'naccessResults'
    naccessRunFileDirectory = 'externalTools/naccess/naccess'
    naccessLogFileDirectory = '%s/naccessLog_%s.txt' %(logFileDirectory, timeStamp)
    minInterfaceResidueCriteria = 5
    differenceBetweenComplexAndMonomerASACriteria = 1.0 
    fullInterfaceListAfterNaccessFileDirectory = '%s/fullInterfaceListAfterNaccess.txt' %(fullListFileDirectory)
#    mainNaccessWithComplexPDBFileGenerator(fullPossibleInterfaceListFile, allPDBFilesDirectory, generatedPDBFilesDirectory, fullPDBIDListFileDirectory, fullInterfaceNaccessResultsFileDirectory, interfaceNaccessResultsFileDirectory, naccessResultsFileDirectory, naccessRunFileDirectory, naccessLogFileDirectory, minInterfaceResidueCriteria, differenceBetweenComplexAndMonomerASACriteria, fullInterfaceListAfterNaccessFileDirectory, numberOfProcesses)

    ### interface extraction by using VDW ###
    VDWInterfaceFileDirectory = 'interfaceFilesGeneratedByVDW'
    VDWInterfaceSkeletonDirectory = 'interfaceSkeletonFilesGeneratedByVDW'
    VDWResidueContactInfoFileDirectory = 'residueContactInfoGeneratedByVDW'
    VDWInterfaceList = '%s/VDWInterfaceList_%s.txt' %(archivedFilesDirectory, timeStamp)
    VDWInterfaceTypeListFileDirectory = '%s/VDWInterfaceTypeList_%s.txt' %(archivedFilesDirectory, timeStamp)
    fullVDWInterfaceTypeListFileDirectory = '%s/fullVDWInterfaceTypeList.txt' %(fullListFileDirectory)
    VDWCriteria = 0.5
    fullNaccessInterfaceList = possibleInterfaceListFile
#    mainInterfaceExtractorUsingVDW(fullInterfaceListAfterNaccessFileDirectory, allPDBFilesDirectory, VDWCriteria, VDWInterfaceFileDirectory, VDWInterfaceSkeletonDirectory, VDWResidueContactInfoFileDirectory, VDWInterfaceTypeListFileDirectory, fullVDWInterfaceTypeListFileDirectory, fullPDBIDListFileDirectory, numberOfProcesses)

    ### interface ASA extraction ###
    vdwInterfaceASAValuesFileDirectory = '%s/VDWInterfaceASAValues_%s.txt' %(archivedFilesDirectory, timeStamp)
    fullVDWInterfaceASAValuesFileDirectory = '%s/fullVDWInterfaceASAValues.txt' %(fullListFileDirectory)
    interfaceASAExtractionLogFileDirectory = '%s/interfaceASAExtractionLog_%s.txt' %(logFileDirectory, timeStamp)
#    mainInterfaceASAValueFinder(fullInterfaceListAfterNaccessFileDirectory, VDWInterfaceFileDirectory, naccessResultsFileDirectory, vdwInterfaceASAValuesFileDirectory, fullVDWInterfaceASAValuesFileDirectory, interfaceASAExtractionLogFileDirectory, numberOfProcesses)

    ### multiprot interface list extractor ###
    maximumToMinimumInterfaceExcessRatio = 0.25
    maximumToMinimumSkeletonExcessRatio = 0.50
    fullInterfaceTypeListForMultiprotFileDirectory = '%s/fullInterfaceTypeListForMultiprot.txt' %(fullListFileDirectory)
    multiprotInterfaceListFileDirectory = 'multiprotInterfaceLists'
    multiprotListSuffix = 'multiprotList.txt'
    multiprotResultsForInterfaceFilesDirectory = 'multiprotResultsForInterfaces'
    multiprotResultsSuffix = 'multiprotResults.txt'
    fullMultiprotStatisticsFileDirectory = '%s/fullMultiprotStatistics.txt' %(fullListFileDirectory)
#    mainMultiprotInterfaceListExtractor(fullVDWInterfaceTypeListFileDirectory, minInterfaceResidueCriteria, maximumToMinimumInterfaceExcessRatio, maximumToMinimumSkeletonExcessRatio, fullInterfaceTypeListForMultiprotFileDirectory, multiprotInterfaceListFileDirectory, multiprotListSuffix, multiprotResultsForInterfaceFilesDirectory, multiprotResultsSuffix, fullMultiprotStatisticsFileDirectory,  numberOfProcesses)

    ### multiprot driver ###
    interfaceStartIndex = 173001
    interfaceEndIndex = 173024
    multiprotRunFileDirectory = 'externalTools/multiprot/multiprot.Linux'
    multiprotRunForProcessesFileDirectory = 'multiprotRuns'
    multiprotRunListLogFileDirectory = '%s/multiprotRunListLogs_%s.txt' %(logFileDirectory, timeStamp)
    multiprotRunListResultFileDirectory = '%s/multiprotRunListResults_%s.txt' %(archivedFilesDirectory, timeStamp)
    multiprotLogFileDirectory = '%s/multiprotLogs' %(logFileDirectory)

    #### start yunus cluster run strategy ####
    multiprotInterfaceListFileDirectory = '/mnt/kufs/scratch/ecukuroglu/project/multiprotInterfaceLists'
    #### end yunus cluster run strategy ####

    mainMultiprotDriver(fullMultiprotStatisticsFileDirectory, interfaceStartIndex, interfaceEndIndex, multiprotRunFileDirectory, multiprotRunForProcessesFileDirectory, multiprotInterfaceListFileDirectory, multiprotListSuffix, multiprotResultsForInterfaceFilesDirectory, multiprotResultsSuffix, multiprotRunListLogFileDirectory, multiprotRunListResultFileDirectory, VDWInterfaceFileDirectory, VDWInterfaceSkeletonDirectory, multiprotLogFileDirectory, numberOfProcesses)

    #### start yunus cluster run strategy ####
    os.system('tar czf %s.tar.gz %s' %(multiprotResultsForInterfaceFilesDirectory, multiprotResultsForInterfaceFilesDirectory))    
    os.system('mv %s.tar.gz %s' %(multiprotResultsForInterfaceFilesDirectory, currentDirectory))
    os.system('rm -r %s' %(ownerDirectory))
    #### end yunus cluster run strategy ####
    t2 = time.time()
    print('\nTotal elapsed time = %f seconds\n' %(t2-t1))
    print('Time stamp : %s' %(time.asctime()))
    print('\n*** MAIN PROGRAM COMPLETED ***\n')
else:
    print('Enter the name of the PDB depositor file and pdbIdListFile')

