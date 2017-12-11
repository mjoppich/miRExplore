from collections import defaultdict

import editdistance
from nertoolkit.geneontology.GeneOntology import GeneOntology
from porestat.utils.DataFrame import DataFrame

from utils.idutils import miRExploreDir

dbData = DataFrame.parseFromFile(miRExploreDir + "/miR2Disease/mirna_disease.tsv", bConvertTextToNumber=False)

allDiseases = set()
for row in dbData:

    disease = row['disease']

    if disease == 'None':
        continue

    allDiseases.add(disease.upper())

print(len(allDiseases))

diseaseObo = GeneOntology(miRExploreDir + "/doid.obo")

disease2obo = defaultdict(set)

"""

find perfect matches

"""
for oboID in diseaseObo.dTerms:
    oboNode = diseaseObo.dTerms[oboID]
    doidName = oboNode.name

    upperDiseaseName = doidName.upper()

    if upperDiseaseName in allDiseases:
        disease2obo[upperDiseaseName].add(oboID)
    else:
        pass

"""

Find duplicates

"""

for disease in allDiseases:

    if disease[-1] == ')':
        actDisease = disease.split('(')[0].strip()
        if actDisease in disease2obo:
            disease2obo[disease] = disease2obo[actDisease]


for x in disease2obo:
    print(x, disease2obo[x])

    if x in allDiseases:
        allDiseases.remove(x)

print(len(allDiseases))
print(allDiseases)

"""

find approximative matches

"""

oboname2node = {}

for oboID in diseaseObo.dTerms:
    oboNode = diseaseObo.dTerms[oboID]
    doidName = oboNode.name.upper()

    oboname2node[doidName] = oboNode

for disease in allDiseases:

    if disease[-1] == ')':
        actDisease = disease.split('(')[0].strip().upper()

        if actDisease.startswith('ORAL SQUAMOUS'):
            print(actDisease)

        if actDisease in oboname2node:
            disease2obo[disease].add(oboname2node[actDisease].id)


subname2node = defaultdict(set)
for disease in allDiseases:

    if disease[-1] == ')':
        actDisease = disease.split('(')[0].strip().upper()
    else:
        actDisease = disease

    if actDisease == 'ACUTE LYMPHOBLASTIC LEUKEMIA':
        print(actDisease)

    obo2dist = {}

    allDiseaseWords = actDisease.split(' ')

    for oboname in oboname2node:

        foundOne = True
        for word in allDiseaseWords:
            if len(word) > 4 and word in oboname:
                foundOne = True
                break

        if not foundOne:
            continue

        dist = editdistance.eval(oboname, actDisease)
        obo2dist[oboname] = dist

    if len(obo2dist) == 0:
        continue

    bestDistances = sorted(obo2dist.items(), key=lambda x: x[1])

    print(disease, bestDistances[0])

    if len(disease) < 8 and actDisease not in bestDistances[0][0]:
        continue

    # add best distance
    disease2obo[disease].add(oboname2node[bestDistances[0][0]].id)

    if bestDistances[0][1] == 0:
        continue

    # add other good distances
    for distancePairs in bestDistances:

        if distancePairs[1] < 5:
            disease2obo[disease].add( oboname2node[ distancePairs[0] ].id )

for x in disease2obo:
    print(x, disease2obo[x])

    if x in allDiseases:
        allDiseases.remove(x)

print(len(allDiseases))
print(allDiseases)


for x in disease2obo:
    print(x, disease2obo[x], [diseaseObo.getID(y).name for y in disease2obo[x]])

doidColIdx = dbData.addColumn("doid", None)
diseaseColIdx = dbData.getColumnIndex('disease')

def addDoidToRow(oldrow):

    ddata = oldrow[diseaseColIdx]

    if ddata == None:
        return oldrow

    if not ddata.upper() in disease2obo:
        return oldrow

    doids = disease2obo[ddata.upper()]

    oldrow[doidColIdx] = ", ".join([str(x) for x in doids])

    return oldrow

dbData.applyToRow( addDoidToRow )


dbData.export(miRExploreDir + "/miR2Disease/mirna_disease.doid.tsv")
