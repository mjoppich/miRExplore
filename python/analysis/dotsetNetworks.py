import matplotlib
from collections import defaultdict, OrderedDict

from plots.DotSetPlot import DotSetPlot

network2title = {

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

pubmed2tuples = defaultdict(set)
mirna2evflows = defaultdict(set)

dataLabels = defaultdict(set)
#"miR-98",  "miR-125a"
#manuMirnas = ["miR-98",  "miR-125a","miR-21", "miR-34a", "miR-93", "miR-125b", "miR-126", "miR-146a", "miR-155", "miR-370"]
#manuMirnas = ['miR-126', 'miR-21', 'miR-155', 'miR-146a', 'miR-125b', 'miR-34a', 'miR-499', 'miR-221', 'miR-370', 'miR-504']
#manuMirnas = ['miR-181c', 'miR-222', 'miR-126', 'miR-155', 'miR-125b', 'miR-34a', 'miR-370', 'miR-146a', 'miR-21', 'miR-93']

miRNA2InteractionPartner = defaultdict(set)
miRNA2Evidences = defaultdict(set)

network2nicename = {
    'andreou_fig2_athero': "miRNA-Gene interactions in\natherosclerotic plaque destabilization\n(Fig.3, Andreou et al.)",
    'andreou_fig1_athero':"Endothelial miRNA-Gene\ninteractions in Atherosclerosis\n(Fig.2, Andreou et al.)",
    'andreou_table1_athero': "miRNA-Gene interactions involved in\ninitiation, progression and thrombotic\ncomplications of atherosclerosis\n(Tab.1, Andreou et al.)",
    'inflammatory_ec_athero': "miRNA-mediated inflammatory\nresponse in endothelial cells\n(Hartmann et. al.)",
    'macrophages_athero': "miRNA-mediated CCL2\nexpression in macrophages\n(Hartmann et. al.)"
}

dataLabels['network'] = [
    network2nicename[x] for x in ['andreou_table1_athero', 'andreou_fig2_athero', 'andreou_fig1_athero', 'inflammatory_ec_athero', 'macrophages_athero']
]
with open("/mnt/d/yanc_network/important_networks.txt", 'r') as fin:

    for line in fin:

        line = line.strip().split("\t")

        #CV-IPN-Endothelial_cell-monocyte_interaction_1	VEGFA	miR-140	EC	27035554
        cbn = line[0]
        if not cbn in ['andreou_fig1_athero', 'andreou_fig2_athero', 'andreou_table1_athero',
                       'inflammatory_ec_athero', 'macrophages_athero']:
            continue

        cbn = network2nicename.get(cbn, cbn)

        gene = line[1]
        miRNA = line[2]
        cellT = celltype2nicename.get(line[3], line[3])
        evidence = line[4]

        if "US" in miRNA:
            continue

        miRNA2InteractionPartner[miRNA].add(gene)
        miRNA2Evidences[miRNA].add(evidence)

        dataLabels["Cell-Type"].add(cellT)
        #dataLabels["network"].add(cbn)

        mirna2evidenceCellT[miRNA][evidence].add(cellT)
        mirna2evidenceCBN[miRNA][evidence].add(cbn)


#for x in manuMirnas:
#    print(x, miRNA2InteractionPartner[x], miRNA2Evidences[x])

allMiRNA = set()
for x in mirna2evidenceCellT:
    allMiRNA.add(x)


for x in mirna2evidenceCBN:
    allMiRNA.add(x)


dataUpPlot = {}

for miRNA in allMiRNA:

    miRNAEvs = set()
    for x in mirna2evidenceCBN.get(miRNA, []):
        miRNAEvs.add(x)

    for x in mirna2evidenceCellT.get(miRNA, []):
        miRNAEvs.add(x)

    miRNAData = {
        "network": set(),
        "Cell-Type": set()
    }

    for ev in miRNAEvs:

        cellT = mirna2evidenceCellT[miRNA].get(ev, None)
        cbns = mirna2evidenceCBN[miRNA].get(ev, None)

        if cellT != None:
            miRNAData['Cell-Type'] = miRNAData['Cell-Type'].union(cellT)

        if cbns != None:
            miRNAData['network'] = miRNAData['network'].union(cbns)


    dataUpPlot[miRNA] = miRNAData

orderDict = OrderedDict()

dataLabels['Cell-Type'] = sorted(dataLabels['Cell-Type'])

for type in ["network", "Cell-Type"]:
    orderDict[type] = dataLabels[type]

def makeMIRNAName(miRNA):
    return miRNA
    return miRNA + " (" + str(len(miRNA2InteractionPartner[miRNA])) + ","+ str(len(miRNA2Evidences[miRNA]))+")"

filteredData = OrderedDict()

#for miRNA in manuMirnas:
#    if miRNA in dataUpPlot:
#        filteredData[makeMIRNAName(miRNA)] = dataUpPlot[miRNA]
#    else:
#        print("Missing manu", miRNA)

stages2 = 0
stages0 = 0

from natsort import natsorted
for miRNA in natsorted(dataUpPlot, key=lambda x: x.split("-")[1]):

    stages = dataUpPlot[miRNA]['network']

    if len(miRNA2Evidences[miRNA]) <= 0:
        continue

    if len(dataUpPlot[miRNA]['network']) == 0:
        pass

    filteredData[makeMIRNAName(miRNA)] = dataUpPlot[miRNA]


print(len(dataUpPlot))
print(len(filteredData))
print(stages2)
print(stages0)

interactorCounts = [len(miRNA2InteractionPartner[miRNA]) for miRNA in filteredData]
pubmedCounts = [len(miRNA2Evidences[miRNA]) for miRNA in filteredData]

DotSetPlot().plot(dataLabels, filteredData,
                  numbers={"Interactor Count":interactorCounts , "PubMed Evidence Count": pubmedCounts },
                  sortData=False,order=orderDict,
label2orderspan = {
    'network': 3,
    'Cell-Type': 2
}
)#, max=30)
matplotlib.pyplot.savefig("/mnt/d/owncloud/markus/uni/publications/miReview/dotset_networks.pdf")
matplotlib.pyplot.show()
