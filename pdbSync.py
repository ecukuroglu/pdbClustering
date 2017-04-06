#################
## Created by Engin Cukuroglu
#################

import os, sys, glob, time

def mainPDBSync(pdbIDListFileDirectory, pdbIDEntriesFileDirectory):
    print('\n\n* PDB SYNC STARTED *\n\n')
    print('Time stamp : %s' %(time.asctime()))

    t1 = time.time()    

    entriesFileDirectory = 'entries.idx'
    os.system("rsync -rlpt -v -z --delete --port=33444 rsync.wwpdb.org::ftp_derived/index/%s ./ > rsync_ftp.log 2>/dev/null" %(entriesFileDirectory))
    os.system('mv %s %s' %(entriesFileDirectory, pdbIDEntriesFileDirectory))

    pdbIDEntriesFile = open(pdbIDEntriesFileDirectory, 'r')
    pdbIDListFile = open(pdbIDListFileDirectory, 'w')
    for line in pdbIDEntriesFile:
        splittedLine = line.strip().split('\t')
        if len(splittedLine[0]) == 4:
            pdbName = splittedLine[0].upper()
            pdbIDListFile.write('%s\n' %(pdbName))
    pdbIDEntriesFile.close()
    pdbIDListFile.close()

    t2 = time.time()
    print('Elapsed time = %f seconds' %(t2-t1))
    print('Time stamp : %s' %(time.asctime()))
    print('\n\n* PDB SYNC COMPLETED *\n\n')

