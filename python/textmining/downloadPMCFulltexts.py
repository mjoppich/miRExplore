import os
from ftplib import FTP
from urllib import request

from utils.parallel import MapReduce

if __name__ == '__main__':

    ftpBase = 'ftp.ncbi.nlm.nih.gov'
    pdfBase = 'pub/pmc/oa_pdf'
    saveDirectory = "/media/disk/pmc/"

    ftp = FTP(ftpBase)
    ftp.login()
    ftp.cwd(pdfBase)


    allFolders = []
    ftp.dir('-d','*/',lambda L:allFolders.append(L.split()[-1]))

    print(allFolders)

    def downloadFiles(allURLFilePaths, env):


        for (fileURL, filePATH) in allURLFilePaths:

            if os.path.exists(filePATH) and os.path.isfile(filePATH):
                print("EXISTS:", fileURL, filePATH)
                continue

            print(fileURL, filePATH)

            request.urlretrieve(fileURL, filePATH)

    ftp.close()

    for primFolder in allFolders:

        if primFolder in ['00/', '01/', '02/', '03/', '04/']:
            print("Skipping Primary Folder", primFolder)
            continue

        ftp = FTP(ftpBase)
        ftp.login()
        ftp.cwd(pdfBase + "/" + primFolder)

        print("from base folder",ftp.pwd())
        #ftp.cwd(primFolder)
        print("to prim folder", ftp.pwd())

        allSecFolders = []
        ftp.dir('-d', '*/', lambda L: allSecFolders.append(L.split()[-1]))

        if not os.path.exists(saveDirectory + "/" + primFolder):
            os.makedirs(saveDirectory + "/" + primFolder)

        procFolders = [x for x in allSecFolders]
        toDownloadFiles = []

        for folderPath in procFolders:

            commonFilePath = "/"+primFolder + "/" + folderPath
            saveFolderPath = saveDirectory + "/" + commonFilePath

            if not os.path.exists(saveFolderPath):
                os.makedirs(saveFolderPath)

            #print(ftp.pwd())
            ftp.cwd(folderPath)

            allFiles = ftp.nlst()

            print(folderPath, "Files", len(allFiles))


            for downloadFile in allFiles:
                fileURL = "ftp://" + ftpBase + "/" + pdfBase + "/" + commonFilePath + "/" + downloadFile
                filePATH = saveFolderPath + "/" + downloadFile
                toDownloadFiles.append((fileURL, filePATH))

                print("Adding", fileURL, filePATH)

            ftp.cwd("..")

            if len(toDownloadFiles) > 50:
                break
            #print(ftp.pwd())



        print("Downloading", len(toDownloadFiles), "Files")
        ftp.close()

        ll = MapReduce(2)
        result = ll.exec(toDownloadFiles, downloadFiles, None, 1, None)


    ftp.quit()