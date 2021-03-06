import re

interactions = {
    'CCL9': ['miR-30d-3p', 'miR-3473c', 'let-7g-5p'],
    'CXCL5': ['miR-204-5p', 'let-7g-5p', 'miR-362-3p', 'miR-155-5p'],
    'CXCL1': ['miR-194-2-3p', 'miR-128-3p', 'miR-194-5p', 'miR-199b-5p', 'miR-467g', 'miR-122-5p'],
    'CXCL13': ['miR-122-5p'],
    'CXCL14': ['miR-301b-3p'],
    'CXCR2': ['let-7g-5p', 'let-7b-5p', 'let-7f-5p', 'let-7c-5p', 'let-7a-5p', 'let-7i-5p', 'miR-98-5p'],
    'CXCL7': ['let-7g-5p'],
    'CCL2': ['let-7a-5p', 'let-7b-5p', 'let-7f-5p', 'let-7c-5p', 'let-7g-5p', 'let-7i-5p', 'miR-181a-5p'],
    'CXCL9': ['miR-1935'],
    'CCL3': ['miR-30a-5p', 'miR-30b-5p', 'miR-30c-5p', 'miR-30d-5p', 'miR-30e-5p'],
    'CCL7': ['miR-181a-5p', 'miR-322-5p', 'miR-29a-5p', 'miR-29b-1-5p'],
    'CCL22': ['miR-34a-5p'],
    'CXCL10': ['miR-503-3p', 'miR-186-5p'],
    'CCR5': ['miR-186-5p', 'miR-669j', 'miR-21-5p', 'miR-146a-5p', 'miR-150-5p', 'miR-146b-5p', 'miR-669k-3p',
             'miR-142-3p', 'miR-34a-5p'],
    'CCL4': ['miR-27b-3p', 'miR-27a-3p', 'miR-21-3p', 'miR-467f'],
    'CX3CL1': ['miR-15a-5p', 'miR-322-5p', 'miR-706', 'miR-762', 'miR-665-3p', 'miR-758-3p', 'miR-381-3p'],
    'CXCR4': ['miR-381-3p', 'miR-21-3p', 'miR-467a-5p', 'miR-467h', 'miR-218-5p', 'miR-1a-3p', 'miR-181d-5p',
              'miR-206-3p', 'miR-181b-5p', 'miR-9-5p', 'miR-132-3p', 'miR-25-3p', 'miR-467d-5p', 'miR-669k-3p',
              'miR-146b-5p', 'miR-467b-5p', 'miR-467e-5p', 'miR-467f', 'miR-146a-5p'],
    'CCR7': ['let-7g-5p', 'miR-23b-3p', 'miR-669p-5p', 'miR-23a-5p', 'let-7e-5p', 'miR-669l-5p', 'miR-15a-5p',
             'miR-467e-5p', 'miR-21-5p', 'miR-16-5p', 'let-7d-5p', 'miR-669n', 'miR-98-5p', 'let-7b-5p', 'let-7a-5p',
             'let-7i-5p', 'let-7c-5p', 'miR-15b-5p', 'miR-467h'],
    'CXCL12': [
        'miR-532-5p', 'miR-130b-3p', 'miR-222-3p', 'miR144-3p', 'miR-542-3p', 'miR-149-5p', 'miR-330-3p', 'miR-532-3p',
        'miR-3470b', 'miR-125b-5p', 'miR-221-3p', 'miR-19b-3p', 'miR-301b-3p',
        'miR-34b-5p', 'miR-125a-3p', 'miR-126-3p', 'miR-16-1-3p', 'miR-882', 'miR-497-5p', 'miR-26a-5p', 'miR-124-3p',
        'miR-26b-5p', 'miR-5620-3p', 'mIR-19a-3p', 'miR-130a-3p', 'miR-690',
        'miR-185-5p', 'miR-31-5p', 'miR-340-5p', 'miR-1843-5p', 'miR-466f-3p', 'miR-301a-3p', 'miR-101a-3p',
        'miR-210-3p', 'miR-107-3p', 'miR-706', 'miR-23b-3p', 'miR-146a-5p', 'miR-467f',
        'miR-322-5p', 'miR-15a-5p', 'miR-29b-1-5p', 'let-7e-5p', 'miR-23a-3p', 'miR-338-3p', 'miR-103-3p', 'miR-362-3p',
        'let-7g-5p', 'miR-155-5p', 'miR-140-5p', 'miR-122-5p', 'miR-22-3p', 'miR-3470a', 'let-7d-5p'
    ]

}



fetchedInteractions = set()


with open('/tmp/mirtex/mirel', 'r') as fin:

    for iline in fin:

        if 'something' in iline:
            continue

        alllines = []

        spl = re.search(r"\]([A-Za-z])", iline)

        if spl != None:



            alllines.append(iline[0:spl.regs[0][0]+1])
            alllines.append(iline[spl.regs[0][0]+1:])
        else:
            alllines.append(iline)

        for line in alllines:

            aline = line.strip().split('\t')

            if len(aline) < 3:
                continue

            gene = aline[0]
            intertype = "chemokine"

            if not gene in interactions:
                intertype = "other"

            mirna = aline[1]

            try:
                evidences = eval(aline[3])

            except:

                spl = re.match(r"\][A-Za-z]", iline)

                exit(-1)

            for x in evidences:

                pmid = x[0]

                fetchedInteractions.add(

                    (
                        pmid,
                        intertype,
                        gene,
                        mirna
                    )

                )


for elem in sorted(fetchedInteractions, key=lambda x: x[0]):

    print("\t".join(elem))