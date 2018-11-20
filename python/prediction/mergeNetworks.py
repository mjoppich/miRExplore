from collections import defaultdict

if __name__ == "__main__":

    base = "/home/mjoppich/ownCloud/data/mirpredict/"

    allRels = set()

    gene2rels = defaultdict(set)

    with open(base + "transmir_rels.tsv", 'r') as fin:

        for line in fin:

            line = line.strip().split("\t")

            miRNA = line[0]
            gene = line[1]

            action = line[2].upper()

            if not "ACTIV" in action:
                continue

            allRels.add((gene, miRNA))

            gene2rels[gene].add( (miRNA, action, "transmir") )

    with open(base + "mirexplore_rels.tsv", 'r') as fin:

        for line in fin:

            line = line.strip().split("\t")

            miRNA = line[1]
            gene = line[0]

            action = line[2].upper()


            if (gene, miRNA) in allRels:
                #print(gene, miRNA)
                continue

            gene2rels[gene].add((miRNA, action, "mirexplore"))

    for gene in gene2rels:

        allrels = gene2rels[gene]

        for rel in allrels:

            pass
            print(gene, rel[0], rel[1], rel[2], sep="\t")