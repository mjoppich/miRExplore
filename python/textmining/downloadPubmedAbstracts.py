from urllib import request
from urllib.error import URLError

downloadBase = False
downloadUpdates = True
updateStart = 1289

if downloadBase:
    for i in range(1, 893):

        fileID = "{:>4}".format(i).replace(" ", "0")
        downloadFile = "medline17n"+fileID+".xml.gz"
        downloadLocation = "/local/storage/pubmed/"
        print(downloadFile)

        request.urlretrieve("ftp://ftp.ncbi.nlm.nih.gov/pubmed/baseline/"+downloadFile, downloadLocation + "/" + downloadFile)

if downloadUpdates:

    for i in range(updateStart, 10000):
        fileID = "{:>4}".format(i).replace(" ", "0")
        downloadFile = "medline17n" + fileID + ".xml.gz"
        downloadLocation = "/local/storage/pubmed/"
        print(downloadFile)

        try:
            (url, resp) = request.urlretrieve("ftp://ftp.ncbi.nlm.nih.gov/pubmed/updatefiles/" + downloadFile,
                                downloadLocation + "/" + downloadFile)
        except URLError as e:
            print(downloadFile, "not existing")
            break