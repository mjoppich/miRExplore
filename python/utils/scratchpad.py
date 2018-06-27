from synonymes.GeneOntology import GeneOntology

cellsObo = GeneOntology("/home/mjoppich/dev/data/tm_soehnlein/obos/cl.obo")


ecTerm = cellsObo.dTerms['CL:0000084']

ac = ecTerm.getAllChildren()

for child in ac:
    print(child.term.id, child.term.name)