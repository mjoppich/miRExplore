from collections import Counter

import spacy
from spacy import displacy
import textacy


#alldeps = [(t.idx, t.text, t.dep_, t.pos_, t.head.text) for t in doc]


def getStack(t):
    h = t
    stack = []
    while h != None:
        stack.append( (h, h.dep_, h.pos_) )

        h = h.head if h != h.head else None

    return stack

def getAllChildren(gen):

    allBase = [x for x in gen]
    allElems = []

    for x in allBase:

        allElems += getAllChildren(x.children)

    return allBase + allElems


def analyseStacks(stackL, stackR):

    tokensL = [x[0] for x in stackL]
    tokensR = [x[0] for x in stackR]

    intersection = set(tokensL).intersection(set(tokensR))

    ret = []

    for elem in intersection:

        if elem.pos_ == "VERB":

            lBase = stackL[0][0]
            rBase = stackR[0][0]

            allLefts = getAllChildren(elem.lefts)
            allRights = getAllChildren(elem.rights)

            lPart = lBase in allLefts
            rPart = rBase in allRights

            if lPart == True and rPart == True:
                ret.append((lPart, rPart, elem))

    return ret

def analyseVerbs(doc, lWord, rWord):

    ret = []
    for verb in [x for x in doc if x.pos_ == "VERB" and x.dep_ != "amod"]:

        lelems = getAllChildren(verb.lefts)
        relems = getAllChildren(verb.rights)

        if lWord in lelems and rWord in relems:
            ret.append(verb)

    return ret

def analyseConjunction(stackL, stackR):

    tokensL = [x[0] for x in stackL]
    tokensR = [x[0] for x in stackR]

    retL = [stackL[0]]
    retR = [stackR[0]]

    for elem in stackL:

        if elem[0].dep_ == "conj":
            if elem[0].head in tokensL:
                if not elem in retL:
                    retL.append(elem)

        if elem[0] in tokensR:
            break

    for elem in stackR:

        if elem[0].dep_ == "conj":
            if elem[0].head in tokensR:

                if not elem in retR:
                    retR.append(elem)

        if elem[0] in tokensL:
            break

    return set(retL).intersection(set(retR))




#displacy.serve(doc, style='dep', port=5005)

testCases = [
    (u'UBI may enhance the phagocytic capacity of various phagocytic cells (neutrophils and dendritic cells), inhibit lymphocytes, and oxidize blood lipids', 11, 14),
    (u'Only rapamycin proved effective against acute Th17-dependent airway inflammation, accompanied by increased plasmacytoid dendritic cells (pDCs) and reduced neutrophils as well as reduced CXCL-1 levels in BAL', 21, 6),
    (u'Only rapamycin proved effective against acute Th17-dependent airway inflammation, accompanied by increased plasmacytoid dendritic cells (pDCs) and reduced neutrophils as well as reduced CXCL-1 levels in BAL', 21, 15),
    (u'Neutrophils are often regarded as cells causing tissue damage, but in recent years it has become clear that a subset of human neutrophils is capable of suppressing T-cells, which is dependent on Mac-1 (CD11b/CD18)', 0, 30),
    (u'Neutrophils are often regarded as cells causing tissue damage, but in recent years it has become clear that a subset of human neutrophils is capable of suppressing T-cells, which is dependent on Mac-1 (CD11b/CD18)', 23, 30),
    (u'Our data show that neutrophil CD11b/CD18 limits pathology in influenza-induced, T-cell mediated disease', 4, 17),
    (u'Ectosomes from neutrophil-like cells down-regulate nickel-induced dendritic cell maturation and promote Th2 polarization.', 2, 13),
    (u'Microvesicles released by apoptotic human neutrophils suppress proliferation and IL-2/IL-2 receptor expression of resting T helper cells.', 5,17),
    (u'Both events converged to trigger the recruitment of 2 waves of immune cells: a swift, massive recruitment of neutrophils, followed by a delayed rise in monocytes and CD8 T cells in the tumor mass', 18,25),
    (u'The authors of this study aimed to elucidate the efficiency of preoperative inflammatory markers, including neutrophil/lymphocyte ratio (NLR), derived NLR (dNLR), platelet/lymphocyte ratio (PLR), lymphocyte/monocyte ratio (LMR), and prognostic nutritional index (PNI), and their paired combinations as tools for the preoperative diagnosis of glioma, with particular interest in its most aggressive form, glioblastoma (GBM)', 0,0),
    (u'Multivariable-adjusted linear regression models showed that neutrophil concentration was associated with FEV1 decline rate (1.14 ml/year decline per 1000 neutrophils/µl, 95% CI: 0.69-1.60 ml/year, p<0.001), while eosinophil concentration was associated with FEV1 decline rate in ever-smokers (1.46 ml/year decline per 100 eosinophils/µl, 95% CI: 0.65-2.26 ml/year, p<0.001) but not in never-smokers (p for interaction=0.004)', 0,0),
    (u'Gene expression levels were strongly associated with cell differentials, explaining 71% of variation in eosinophil counts and 64% of variation in neutrophil counts',21 , 14),
    (u'TLR7/8-activated neutrophils promoted cleavage of FcgRIIA on plasmacytoid dendritic cells and monocytes, resulting in impaired overall clearance of ICs and increased complement C5a generation', 1, 11),
    (u'Lupus neutrophils produce elevated levels of factors known to fuel autoantibody production, including IL-6 and B cell survival factors, but also reactive oxygen intermediates, which can suppress lymphocyte proliferation', 1, 16),
    (u'Human neutrophils mediate trogocytosis rather than phagocytosis of CLL B cells opsonized with anti-CD20 antibodies', 1,9),
    (u'Mycobacterium tuberculosis-induced neutrophil ectosomes decrease macrophage activation.', 4,7),
    (u'Furthermore, KLF12, HMGB1 and CIT mRNAs were confirmed as direct targets of the p53-induced miR-34a, miR-205 and miR-486-5p, respectively', 2, 16),
    (u'It has been shown that macrophage can communicate with endothelial cells via ICAM1 and miR-98', 12, 14),
(u'It has been shown that macrophage can communicate with endothelial cells via ICAM1 and miR-98', 5, 14),
    (u'MicroRNA-495 regulates the proliferation and apoptosis of human umbilical vein endothelial cells by targeting chemokine CCL2', 0, 15),
    (u'miR-100 suppresses the migration of phages by targeting FZD-8', 0, 8),
    (u'Furthermore, luciferase reporter assay analysis identified ZEB2 as a direct target of miR-215', [(7, 13)]),
    (u'Furthermore, some chemical influences both ZEB2 and miR-215', [(6, 8)]),
    (u'Dual-luciferase reporter assays showed that miR-195-5p binds the 3\'-untranslated region (UTR) of CDK8, suggesting that CDK8 should be a direct target of miR-195-5p', 7,16),
    (u'Antitumor miR-150-5p and miR-150-3p inhibit cancer cell aggressiveness by targeting SPOCK1 in head and neck squamous cell carcinoma', 1, 10),
    (u'These findings demonstrated that miR‑186 acted as a tumor suppressor by targeting IGF‑1R in glioma, suggesting miR‑186 may be a potential therapeutic target for the treatment of this disease.', 4,12),
    (u'One potential microRNA that regulates Bcan is miR-9 and overexpression of miR-9 can partly rescue the effects of Dicer1 deletion on the MG phenotype', [(5, 7), (12, 19)]),
    (u'The levels of miR-330-3p were positively correlated with the status of TNM stage (p = 0.011) and lymph node metastasis', [(3, 12)]),
    (u'In conclusion, our findings show that miR-140 acts as a tumor suppressor in OS by targeting HDAC4', [(17, 7)]),
    (u'Moreover, a bioinformatics prediction indicated that the histone deacetylase 4 (HDAC4) is a target gene of miR-140 and is involved in miR-140-mediated suppressive effects', [(9, 19)])
    #(u'One potential microRNA that regulates Bcan is miR-9', [(5, 7)])
]

"We provide evidence that exposure of monocyte-derived dendritic cells (MDDCs) to recombinant HIV-1 R5 gp120, but not to CCR5 natural ligand CCL4, influences the expression of a panel of miRs"


#29099614.2.11
userMirs = ['miR-150-5p', 'miR-150-3p', 'miR-195-5p', 'miR-486-5p', 'miR‑186', 'miR-9', 'miR-330-3p', 'miR-140']
nlp = spacy.load('/mnt/d/spacy/models/en_core_web_lg-2.2.0/en_core_web_lg/en_core_web_lg-2.2.0/')  # create blank Language class #en_core_web_lg
#nlp = spacy.load('/mnt/d/spacy/models/en_core_web_sm-2.2.0/en_core_web_sm/en_core_web_sm-2.2.0/')

for mir in userMirs:

    umir = mir

    nlp.vocab[umir]
    nlp.tokenizer.add_special_case(umir, [{'ORTH': mir, 'POS': 'NOUN', 'TAG': 'NN', "NER": "MIRNA"}])


testCases = [testCases[-2]]


def get_sdp_path(doc, subj, obj, lca_matrix):
    lca = lca_matrix[subj, obj]

    current_node = doc[subj]
    subj_path = [current_node]

    seenNodes = Counter()

    if lca != -1:
        if lca != subj:
            while current_node.head.i != lca:
                current_node = current_node.head
                subj_path.append(current_node)
                seenNodes[current_node] += 1

                if seenNodes[current_node] > 15:
                    return subj_path + [obj]

            subj_path.append(current_node.head)
            current_node = doc[obj]

    seenNodes = Counter()

    obj_path = [current_node]
    if lca != -1:
        if lca != obj:
            while current_node.head.i != lca:
                current_node = current_node.head
                obj_path.append(current_node)
                seenNodes[current_node] += 1

                if seenNodes[current_node] > 15:
                    return subj_path + [obj]

            obj_path.append(current_node.head)

    return subj_path + obj_path[::-1][1:]

for testCase in testCases:

    for lwIdx, rwIdx in testCase[1]:

        doc = nlp(testCase[0])
        lWord = doc[lwIdx]
        rWord = doc[rwIdx]

        for t in doc:
            x = (t.idx, t.text, t.dep_, t.pos_, t.head.text)
            print(x)

        print(testCase[0])
        print(lWord, rWord, lwIdx, rwIdx)

        for chunk in doc.noun_chunks:
            print(chunk.text, chunk.root.text, chunk.root.dep_,
                  chunk.root.head.text)

        sdp = get_sdp_path(doc, lwIdx, rwIdx, doc.get_lca_matrix())
        print(sdp)

        tuples = textacy.extract.subject_verb_object_triples(doc)
        ntuples = []
        if tuples:
            ntuples = list(tuples)

        print("svo", ntuples)

        #for x in [(t.idx, t.text, t.dep_, t.pos_, t.head.text) for t in doc]:
        #    print(x)

        stacks = []
        for token in [lWord, rWord]:
            stack = getStack(token)
            stacks.append(stack)
            print("Stack", token, ":", stack)

        stackRes = analyseStacks(stacks[0], stacks[1])
        print("analyseStacks", stackRes, [(t.idx, t.text, t.dep_, t.pos_, t.head.text) for (l, r, t) in stackRes])
        # analyseConjugation(doc[22], doc[27])
        verbRes = analyseVerbs(doc, lWord, rWord)
        print("analyseVerbs", verbRes, [(t.idx, t.text, t.dep_, t.pos_, t.head.text) for t in verbRes])

        conjRes = analyseConjunction(stacks[0], stacks[1])
        print("analyseConjunction", conjRes, [(t.idx, t.text, t.dep_, t.pos_, t.head.text) for (t, s, r) in conjRes])

    if len(testCases) == 1:
        displacy.serve(doc, style="dep", port=5005)
    print()
    print()
    print()