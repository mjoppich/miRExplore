import os, sys
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")

from ftplib import FTP
from urllib import request

from utils.parallel import MapReduce

if __name__ == '__main__':

    ftpBase = 'ftp.ncbi.nlm.nih.gov'
    pdfBase = 'pub/pmc/oa_pdf'
    saveDirectory = "/mnt/raidtmpbio2/joppich/pmc_apr2020/"

    ftp = FTP(ftpBase)
    ftp.login()
    ftp.cwd(pdfBase)


    allFolders = []
    ftp.dir('-d','*/',lambda L:allFolders.append(L.split()[-1]))

    ftp.quit()

    def downloadFiles(allPaths, env):

        for folder in allPaths:

            cmd = "bash -c \"mkdir -p {localFolder}; cd {localFolder}; lftp {ftpServer} -e 'cd {ftpBaseFolder}; mirror --only-newer . {localFolder}; exit'\"".format(localFolder=saveDirectory+"/"+folder,
                                                                                                                                            ftpServer = ftpBase,
                                                                                                                                            ftpBaseFolder = pdfBase + "/" + folder)

            print(cmd)
            os.system(cmd)

    ll = MapReduce(procs=16)
    result = ll.exec(allFolders, downloadFiles, None, 1, None)


