import glob
import os

sentDir = "/mnt/c/dev/data/pmidp/"
nSentDir = "/mnt/c/dev/data/pmid/"

for file in glob.glob(sentDir + "pubmed*.sent"):

    fileName = os.path.basename(file)
    print(fileName)

    with open(file, 'r') as file:

        fileLines = []

        for line in file:

            line = line.strip().rstrip(".")
            fileLines.append(line)


        with open(nSentDir + "/" + fileName, 'w') as fout:

            for line in fileLines:
                fout.write(line + "\n")




