import os, sys

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")

from urllib import request
from urllib.error import URLError

from utils.parallel import MapReduce

downloadBase = True
downloadUpdates = False
downloadUpdatesParallel = True
updateStart = 1016
updateEnd = 1284
medlineBase="pubmed20n"

downloadLocation = "/mnt/d/dev/data/pmid_jul2020/"

directory = os.path.dirname(downloadLocation)
if not os.path.exists(directory):
    os.makedirs(directory)

onlyNew = True

def downloadDataBase(data, env):

    for i in data:

        fileID = "{:>4}".format(i).replace(" ", "0")
        downloadFile = medlineBase+fileID+".xml.gz"
        print(downloadFile)

        if onlyNew:
            if os.path.exists(downloadLocation + "/" + downloadFile):
                continue

        request.urlretrieve("ftp://ftp.ncbi.nlm.nih.gov/pubmed/baseline/"+downloadFile, downloadLocation + "/" + downloadFile)

def downloadDataUpdate(data, env):

    for i in data:

        fileID = "{:>4}".format(i).replace(" ", "0")
        downloadFile = medlineBase+fileID+".xml.gz"
        print(downloadFile)

        if onlyNew:
            if os.path.exists(downloadLocation + "/" + downloadFile):
                continue

        request.urlretrieve("ftp://ftp.ncbi.nlm.nih.gov/pubmed/updatefiles/"+downloadFile, downloadLocation + "/" + downloadFile)

if downloadBase:

    ll = MapReduce(8)
    ll.exec([i for i in range(1,updateStart)], downloadDataBase, None)

if downloadUpdatesParallel:

    ll = MapReduce(8)
    ll.exec([i for i in range(updateStart, updateEnd)], downloadDataUpdate, None)


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
