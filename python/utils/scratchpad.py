from synonymes.GeneOntology import GeneOntology

cellsObo = GeneOntology("/home/mjoppich/dev/data/tm_soehnlein/obos/cl.obo")


ecTerm = cellsObo.dTerms['CL:0000255']

ac = ecTerm.getAllChildren(withLevel=True)

for child, l in ac:
    print(child.term.id, child.term.name, l, sep="\t")