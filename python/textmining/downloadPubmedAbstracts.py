import os
from urllib import request
from urllib.error import URLError

downloadBase = True
downloadUpdates = False
updateStart = 928
medlineBase="pubmed18n"

downloadLocation = "/local/storage/pubmed18/"

directory = os.path.dirname(downloadLocation)
if not os.path.exists(directory):
    os.makedirs(directory)

if downloadBase:
    for i in range(1, 928):

        fileID = "{:>4}".format(i).replace(" ", "0")
        downloadFile = medlineBase+fileID+".xml.gz"
        print(downloadFile)

        request.urlretrieve("ftp://ftp.ncbi.nlm.nih.gov/pubmed/baseline/"+downloadFile, downloadLocation + "/" + downloadFile)

if downloadUpdates:

    for i in range(updateStart, 10000):
        fileID = "{:>4}".format(i).replace(" ", "0")
        downloadFile = medlineBase + fileID + ".xml.gz"
        print(downloadFile)

        try:
            (url, resp) = request.urlretrieve("ftp://ftp.ncbi.nlm.nih.gov/pubmed/updatefiles/" + downloadFile,
                                downloadLocation + "/" + downloadFile)
        except URLError as e:
            print(downloadFile, "not existing")
            break