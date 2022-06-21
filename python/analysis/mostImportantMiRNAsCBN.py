import matplotlib
from collections import defaultdict, OrderedDict, Counter

from plots.DotSetPlot import DotSetPlot

processToTitle = {
    "targetMirsECA": "EC activation and inflammation",
    "targetMirsMonocyte": "Monocyte diff. & Macrophage act.",
    "targetMirsFCF": "Foam cell formation",
    "targetMirsAngio": "Angiogenesis",
    "targetMirsVasRemod": "Vascular remodeling",
    "targetMirsTCell": "T-cell differentiation & activation",
    "targetMirsCholEfflux": "Cholesterol efflux",
    "targetMirsSMCProlif": "SMC proliferation & SMC migration"
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

manuMirnas = ["miR-21", "miR-34a", "miR-93", "miR-98", "miR-125a", "miR-125b", "miR-126", "miR-146a", "miR-155", "miR-370"]
manuMirnas = ['miR-126', 'miR-21', 'miR-155', 'miR-146a', 'miR-125b', 'miR-34a', 'miR-499', 'miR-221', 'miR-370', 'miR-504']
manuMirnas = ['miR-181c', 'miR-222', 'miR-126', 'miR-155', 'miR-125b', 'miR-34a', 'miR-370', 'miR-146a', 'miR-21', 'miR-93']
# old version
manuMirnas = ["miR-98",  "miR-125a","miR-21", "miR-34a", "miR-93", "miR-125b", "miR-126", "miR-146a", "miR-155", "miR-370"]
#manuMirnas = list({'miR-155', 'miR-93', 'miR-181c', 'miR-370', 'miR-222', 'miR-125b', 'miR-34a', 'miR-146a', 'miR-126', 'miR-21'})
miRNA2InteractionPartner = defaultdict(set)
miRNA2Evidences = defaultdict(set)

with open("/mnt/d/yanc_network/disease_important_cbn.txt", 'r') as fin:

    for line in fin:

        line = line.strip().split("\t")

        #CV-IPN-Endothelial_cell-monocyte_interaction_1	VEGFA	miR-140	EC	27035554

        cbn = network2nicename.get(line[0], line[0])
        gene = line[1]
        miRNA = line[2]
        cellT = celltype2nicename.get(line[3], line[3])
        evidence = line[4]

        miRNA2InteractionPartner[miRNA].add(gene)
        miRNA2Evidences[miRNA].add(evidence)

        dataLabels["Cell-Type"].add(cellT)
        dataLabels["CBN"].add(cbn)

        mirna2evidenceCellT[miRNA][evidence].add(cellT)
        mirna2evidenceCBN[miRNA][evidence].add(cbn)


allMiRNA = set()
for x in mirna2evidenceCellT:
    allMiRNA.add(x)

for x in mirna2evidenceProcess:
    allMiRNA.add(x)

for x in mirna2evidenceCBN:
    allMiRNA.add(x)



fout = open("/mnt/c/Users/mjopp/Desktop/d3-parsets-d3v5/cbn_plot.csv", "w")
print("miRNA", "CBN", "CELLTYPE", sep=",", file=fout)


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
        "Cell-Type": set()
    }

    for ev in miRNAEvs:

        cellT = mirna2evidenceCellT[miRNA].get(ev, ["None"])
        cbns = mirna2evidenceCBN[miRNA].get(ev, ["None"])

        if "None" in cbns:
            continue

        print(len(cellT), len(cbns))

        #mirna2printTuple[miRNA].append((cbns[0], cellT[0]))
        for celltype in cellT:
            for cbn in cbns:
                mirna2printTuple[miRNA].append( (cbn, celltype) )




for miRNA in manuMirnas:
    for (cbn, celltype) in mirna2printTuple[miRNA]:
        print(miRNA, cbn, celltype, sep=",", file=fout)

