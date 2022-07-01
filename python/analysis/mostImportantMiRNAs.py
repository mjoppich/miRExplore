import matplotlib
from collections import defaultdict, OrderedDict

from plots.DotSetPlot import DotSetPlot

processToTitle = {
    "targetMirsECA": "EC activation and\n inflammation",
    "targetMirsMonocyte": "Monocyte diff. &\nMacrophage act.",
    "targetMirsFCF": "Foam cell formation",
    "targetMirsAngio": "Angiogenesis",
    "targetMirsVasRemod": "Vascular remodeling",
    "targetMirsTCell": "T cell differentiation &\n activation",
    "targetMirsCholEfflux": "Cholesterol efflux",
    "targetMirsSMCProlif": "SMC proliferation &\n SMC migration"
}
network2nicename = {
"CV-IPN-Plaque_destabilization_1": "(VI) Plaque destabilization",
"CV-IPN-Platelet_activation_1": "(V) Platelet activation",
"CV-IPN-Smooth_muscle_cell_activation_1": "(IV) SMC activation",
"CV-IPN-Foam_cell_formation_1": "(III) Foam cell formation",
"CV-IPN-Endothelial_cell-monocyte_interaction_1": "(II) EC/MC interaction",
"CV-IPN-Endothelial_cell_activation_1": "(I) EC activation",
}

celltype2nicename = {
    'SMC': "Smooth muscle cell",
    'EC': "Endothelial cell",
    "MC": "Macrophage/Monocyte",
    "FC": "Foam cell"
}

def source2index( sname ):

    if sname != None and sname.startswith("CV-IPN"):
        return 0

    return 1


mirna2evidenceCellT = defaultdict(lambda: defaultdict(set))
mirna2evidenceCBN = defaultdict(lambda: defaultdict(set))
mirna2evidenceProcess = defaultdict(lambda: defaultdict(set))

pubmed2tuples = defaultdict(set)
mirna2evflows = defaultdict(set)

dataLabels = defaultdict(set)
#"miR-98",  "miR-125a"
manuMirnas = ["miR-98",  "miR-125a","miR-21", "miR-34a", "miR-93", "miR-125b", "miR-126", "miR-146a", "miR-155", "miR-370"]
#manuMirnas = ['miR-126', 'miR-21', 'miR-155', 'miR-146a', 'miR-125b', 'miR-34a', 'miR-499', 'miR-221', 'miR-370', 'miR-504']
#manuMirnas = ['miR-181c', 'miR-222', 'miR-126', 'miR-155', 'miR-125b', 'miR-34a', 'miR-370', 'miR-146a', 'miR-21', 'miR-93']

manuMirnas = list({'miR-155', 'miR-93', 'miR-181c', 'miR-370', 'miR-222', 'miR-125b', 'miR-34a', 'miR-146a', 'miR-126', 'miR-21'})
manuMirnas = ["miR-98",  "miR-125a","miR-21", "miR-34a", "miR-93", "miR-125b", "miR-126", "miR-146a", "miR-155", "miR-370"]

miRNA2InteractionPartner = defaultdict(set)
miRNA2Evidences = defaultdict(set)

with open("/mnt/d/yanc_network/disease_pw_important_cbn.txt", 'r') as fin:

    for line in fin:

        line = line.strip().split("\t")

        #CV-IPN-Endothelial_cell-monocyte_interaction_1	VEGFA	miR-140	EC	27035554

        cbn = network2nicename.get(line[0], line[0])
        gene = line[1]
        miRNA = line[2]
        cellT = celltype2nicename.get(line[3], line[3])
        evidence = line[4]

        if "US" in miRNA:
            continue

        miRNA2InteractionPartner[miRNA].add(gene)
        miRNA2Evidences[miRNA].add(evidence)

        dataLabels["Cell-Type"].add(cellT)
        dataLabels["CBN"].add(cbn)

        mirna2evidenceCellT[miRNA][evidence].add(cellT)
        mirna2evidenceCBN[miRNA][evidence].add(cbn)

#important_process
with open("/mnt/d/yanc_network/pathway_important_process.txt", 'r') as fin:

    for line in fin:

        line = line.strip().split("\t")

        #CV-IPN-Endothelial_cell-monocyte_interaction_1	VEGFA	miR-140	EC	27035554

        process = processToTitle.get(line[0], line[0])
        gene = line[1]
        miRNA = line[2]
        cellT = celltype2nicename.get(line[3], line[3])
        evidence = line[4]

        if "US" in miRNA:
            continue

        miRNA2InteractionPartner[miRNA].add(gene)
        miRNA2Evidences[miRNA].add(evidence)

        dataLabels["Cell-Type"].add(cellT)
        dataLabels["Process"].add(process)

        mirna2evidenceCellT[miRNA][evidence].add(cellT)
        mirna2evidenceProcess[miRNA][evidence].add(process)

for x in manuMirnas:
    print(x, miRNA2InteractionPartner[x], miRNA2Evidences[x])

allMiRNA = set()
for x in mirna2evidenceCellT:
    allMiRNA.add(x)

for x in mirna2evidenceProcess:
    allMiRNA.add(x)

for x in mirna2evidenceCBN:
    allMiRNA.add(x)


dataUpPlot = {}

for miRNA in allMiRNA:

    miRNAEvs = set()
    for x in mirna2evidenceCBN.get(miRNA, []):
        miRNAEvs.add(x)

    for x in mirna2evidenceProcess.get(miRNA, []):
        miRNAEvs.add(x)

    for x in mirna2evidenceCellT.get(miRNA, []):
        miRNAEvs.add(x)

    miRNAData = {
        "CBN": set(),
        "Process": set(),
        "Cell-Type": set()
    }

    for ev in miRNAEvs:

        cellT = mirna2evidenceCellT[miRNA].get(ev, None)
        cbns = mirna2evidenceCBN[miRNA].get(ev, None)
        process = mirna2evidenceProcess[miRNA].get(ev, None)

        if cellT != None:
            miRNAData['Cell-Type'] = miRNAData['Cell-Type'].union(cellT)

        if cbns != None:
            miRNAData['CBN'] = miRNAData['CBN'].union(cbns)

        if process != None:
            miRNAData['Process'] = miRNAData['Process'].union(process)

    dataUpPlot[miRNA] = miRNAData

orderDict = OrderedDict()

for type in ["CBN", "Process", "Cell-Type"]:

    orderDict[type] = sorted(dataLabels[type])

def makeMIRNAName(miRNA):
    return miRNA
    return miRNA + " (" + str(len(miRNA2InteractionPartner[miRNA])) + ","+ str(len(miRNA2Evidences[miRNA]))+")"

filteredData = OrderedDict()

for miRNA in manuMirnas:
    if miRNA in dataUpPlot:
        filteredData[makeMIRNAName(miRNA)] = dataUpPlot[miRNA]
    else:
        print("Missing manu", miRNA)

stages2 = 0
stages0 = 0

from natsort import natsorted
for miRNA in natsorted(dataUpPlot, key=lambda x: x.split("-")[1]):

    stages = dataUpPlot[miRNA]['CBN']

    if len(miRNA2Evidences[miRNA]) <= 0:
        continue

    if len(dataUpPlot[miRNA]['Process']) == 0:
        pass#continue
    if len(dataUpPlot[miRNA]['CBN']) == 0:
        continue

    filteredData[makeMIRNAName(miRNA)] = dataUpPlot[miRNA]


print(len(dataUpPlot))
print(len(filteredData))
print(stages2)
print(stages0)

fout = open("/mnt/c/Users/mjopp/Desktop/d3-parsets-d3v5/titanic.csv", "w")
print("miRNA", "CBN", "PROCESS", "CELLTYPE", sep=",", file=fout)


mirna2printTuple = defaultdict(list)

for miRNA in allMiRNA:

    miRNAEvs = set()
    for x in mirna2evidenceCBN.get(miRNA, []):
        miRNAEvs.add(x)

    for x in mirna2evidenceProcess.get(miRNA, []):
        miRNAEvs.add(x)

    for x in mirna2evidenceCellT.get(miRNA, []):
        miRNAEvs.add(x)

    miRNAData = {
        "CBN": set(),
        "Process": set(),
        "Cell-Type": set()
    }

    for ev in miRNAEvs:

        cellT = mirna2evidenceCellT[miRNA].get(ev, ["None"])
        cbns = mirna2evidenceCBN[miRNA].get(ev, ["None"])
        processes = mirna2evidenceProcess[miRNA].get(ev, ["None"])

        if miRNA == "miR-98":
            print(ev, cbns, celltype, process)

        if "None" in cbns:# or "None" in processes:
            continue

        for celltype in cellT:
            for cbn in cbns:
                for process in processes:
                    mirna2printTuple[miRNA].append( (cbn, process, celltype) )


selMirnas = sorted([x for x in mirna2printTuple], reverse=True, key=lambda x: len(mirna2printTuple[x]))
print(selMirnas[0:10])


for miRNA in manuMirnas:
    for (cbn, process, celltype) in mirna2printTuple[miRNA]:
        print(miRNA, cbn.replace("\n", " ").replace("  ", " "), process.replace("\n", " ").replace("  ", " "), celltype, sep=",", file=fout)


interactorCounts = [len(miRNA2InteractionPartner[miRNA]) for miRNA in filteredData]
pubmedCounts = [len(miRNA2Evidences[miRNA]) for miRNA in filteredData]

DotSetPlot().plot(dataLabels, filteredData, numbers={"Interactor Count":interactorCounts , "PubMed Evidence Count": pubmedCounts },sortData=False,order=orderDict)#, max=30)
matplotlib.pyplot.savefig("/mnt/d/owncloud/markus/uni/publications/miReview/dotset_important.pdf")
matplotlib.pyplot.show()
