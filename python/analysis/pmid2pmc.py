


pmid2pmc = {}
with open('/home/mjoppich/evi_ocr/pmid_pmc', 'r') as fin:

    for line in fin:

        line = line.strip().split('\t')

        if len(line) != 2:
            continue

        pmid2pmc[line[1]] = line[0]




allPMIDs2PMC = set()

with open('/home/mjoppich/evi_ocr/pmids_sorted', 'r') as fin:

    for line in fin:

        line = line.strip()
        aline = line.split('\t')

        pmid = aline[0]

        if pmid in pmid2pmc:
            allPMIDs2PMC.add((pmid, pmid2pmc[pmid]))



for elem in allPMIDs2PMC:

    print("\t".join(elem))
