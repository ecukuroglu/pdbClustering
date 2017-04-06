#################
## Created by Engin Cukuroglu
#################

import os, sys, glob, time

from codesOfTools import getPDB


def mainUpdater(pdbFileDirectory, pdbIDListFileDirectory, fullPDBIDListFileDirectory, archivedPDBFileDirectory, pdbDownloadLogFileDirectory):
    t1 = time.time()
    print('\n* PDB UPDATER STARTED *\n')
    print('Time stamp : %s' %(time.asctime()))
    pdbIDListFile = open(pdbIDListFileDirectory, 'r')
    pdbDownloadLogFile = open(pdbDownloadLogFileDirectory, 'w')
    downloadCounter = 0
    pdbCounter = 0
    prePDBIDListDict = {}
    while 1:
        line = pdbIDListFile.readline()
        pdbID = line.strip()
        if pdbID == '':
            break
        pdbCounter = pdbCounter + 1
        prePDBIDListDict[pdbID] = 1
        if not os.path.exists('%s/%s.pdb' %(pdbFileDirectory, pdbID)):
            getPDB(pdbID, pdbFileDirectory)
            if os.path.exists('%s/%s.pdb' %(pdbFileDirectory, pdbID)):    
                downloadCounter = downloadCounter + 1
            else:
                pdbDownloadLogFile.write('%s is not in the ftp\n' %(pdbID))
    
    pdbIDListFile.close()

    pdbsInPDBFile = glob.glob('%s/*.pdb' %(pdbFileDirectory))
    tempFullPDBIDListFileDirectory = '%s_temp' %(fullPDBIDListFileDirectory)
    tempFullPDBIDListFile = open(tempFullPDBIDListFileDirectory, 'w')
    fullPDBIDListDict = {}
    for pdbInDirectory in pdbsInPDBFile:
        splittedDir = pdbInDirectory.split('/')
        pdbInDirectory = splittedDir[-1].replace('.pdb', '')
        if pdbInDirectory in prePDBIDListDict:
            tempFullPDBIDListFile.write('%s\n' %(pdbInDirectory))
            fullPDBIDListDict[pdbInDirectory] = 1
        else:
            os.system('mv %s/%s.pdb %s/' %(pdbFileDirectory, pdbInDirectory, archivedPDBFileDirectory))
            pdbDownloadLogFile.write('Now, %s is an archived pdb.\n' %(pdbInDirectory))
            
    tempFullPDBIDListFile.close()    
    pdbDownloadLogFile.close()
    
    os.system('mv %s %s' %(tempFullPDBIDListFileDirectory, fullPDBIDListFileDirectory))

    t2 = time.time()
    print('Number of PDB in the list = %d\n' %(pdbCounter))
    print('Number of PDB Downloaded = %d\n' %(downloadCounter))
    print('Number of PDB in %s directory = %d\n' %(pdbFileDirectory, len(fullPDBIDListDict)))
    print('\nElapsed time = %f seconds\n' %(t2-t1))
    print('Time stamp : %s' %(time.asctime()))
    print('\n* PDB UPDATER COMPLETED *\n')

