from collections import defaultdict
import os, sys
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")


import scispacy
import spacy
from spacy import displacy

from textmining.MirGeneRelCheck import MirGeneRelCheck, nlp

#nlp = spacy.load('/mnt/d/spacy/models/en_core_web_lg-2.2.0/en_core_web_lg/en_core_web_lg-2.2.0/')
#nlp = spacy.load('/mnt/d/spacy/models/en_core_sci_lg-0.2.4/en_core_sci_lg/en_core_sci_lg-0.2.4/')
#nlp = spacy.load("/mnt/d/spacy/models/en_ner_bionlp13cg_md-0.2.4/en_ner_bionlp13cg_md/en_ner_bionlp13cg_md-0.2.4")


# create blank Language class #en_core_web_lg
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

    ret = []

    for elem in stackL:

        if elem[0].dep_ == "conj":
            if elem[0].head in tokensL:
                if not elem in ret:
                    ret.append(elem)

        if elem[0] in tokensR:
            break

    for elem in stackR:

        if elem[0].dep_ == "conj":
            if elem[0].head in tokensR:

                if not elem in ret:
                    ret.append(elem)

        if elem[0] in tokensL:
            break

    return ret

doc = nlp(u'Mycobacterium tuberculosis neutrophil ectosomes decrease macrophage activation.')
doc = nlp(u'Neutrophils are often regarded as cells causing tissue damage, but in recent years it has become clear that a subset of human neutrophils is capable of suppressing T-cells, which is dependent on Mac-1 (CD11b/CD18)')
doc = nlp(u'UBI may enhance the phagocytic capacity of various phagocytic cells (neutrophils and dendritic cells), inhibit lymphocytes, and oxidize blood lipids')



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
    (u'Furthermore, luciferase reporter assay analysis identified ZEB2 as a direct target of miR-215', 7, 13),
    (u'One potential microRNA that regulates Bcan is miR-9, and overexpression of miR-9 can partly rescue the effects of Dicer1 deletion on the MG phenotype', 7, 19),
    (u'Furthermore, ectopic expression of LINC00152 partially halted the decrease in CCND1 expression and cell proliferation capacity induced by miR-193a/b-3p overexpression', 11, 19),
    #(u'Thus, LINC00152 acts as a competing endogenous RNA (ceRNA) by sponging miR-193a/b-3p to modulate its target gene, CCND1', 2, 14),
    #(u'MS2-RIP analysis indicated that LINC00152 binds directly to miR-193a/b-3p, as confirmed by luciferase reporter assays', 4,8),
    (u'Furthermore, ectopic expression of LINC00152 partially halted the decrease in CCND1 expression and cell proliferation capacity induced by miR-193a/b-3p overexpression', 5, 19), # should be false???,
    (u'We also found MVs downregulated RUNX3 expression, and inhibition of miR-210 upregulated RUNX3 expression', 5, 11), #check!
    (u'The data provide the solid evidence that miR-375 plays a tumor-suppressive role in ccRCC progression, partially through regulating YWHAZ', 7, 21), #should pass
    (u'Overall we can see that miR-200, miR-201 and miR-202 all target the genes KLF12, KFL13 and CXCR4', 5, 18), #should pass
    (u'It has been shown that macrophage can communicate with endothelial cells via ICAM1 and miR-98.',12, 14),
    (u'In addition, miR-146a and miR-125a-3p/5p were not regulated at transcriptional levels upon infection, and miR-125a-3p/5p were found to be TLR2 responsive', 3, 21),
    (u'miR-98 is actively inhibiting ICAM1 activity.',0, 4),
    (u'Computational analysis predicted ICAM1 as the strongest target gene for miR-98.', 3, 10 ),
    (u'We found that hematopoietic stem/progenitor cells, stromal cells and endothelial cells readily incorporate these miR-98-loaded microvescicles, and that miR-98 represses ICAM1 expression on bone marrow hematopoietic stem/progenitor cells, stromal cells and endothelial cells.', 17, 24),
    (u'Among these miRNAs, miR-139, miR-101, miR-29A, and miR-181 could target FOS; miR-200 and miR-29A could target KLF4; miR-101 and miR-200 could target DUSP1; miR-139 and miR-181 could target JUN and CD163, respectively', 14, 19),
    (u'miR-200 induction was inversely correlated with expression of known targets, transcription factors ZEB1/2 and TGF-β2.', 0, 13),
    #(u'Western blotting indicated that protein levels of pERK/ERK mediated by SPRED1 were significantly higher in miR-126-treated mice (control vs miR-126)', 12, 21),
    (u'KLF13 messenger RNA expression levels were significantly lower in miR-126-treated mice.', 0, 9),
    (u' miR-200 and miR-29A could target KLF4', 0, 5 ),
    (u'miR-126 has been found to promote angiogenesis and inhibit vascular inflammation in endothelial cells by repressing three target genes Sprouty-related EVH1 domain-containing protein (SPRED1), phosphoinositol-3 kinase regulatory subunit (PIK3R2), and vascular cell adhesion molecule (VCAM1)', 0, 42),
    #(u'Moreover, overexpression of miR-126 significantly decreased (p < 0.05) the protein levels of Crk, which has previously been identified as a direct target of miR-126.', 4 ,16),
    (u'Luciferase activity assays confirmed that miR-200b regulated adipocyte differentiation by directly targeting KLF4, and miR-200 suppressed the expression of KLF4 in both mRNA level and protein level.', 15, 20),
    #(u'These results suggest that AKT2 can regulate miR-200a in a histology- or stage-specific manner and that this regulation is independent of subsequent involvement of miR-200a in epithelial-mesenchymal transition.', 4, 7),
    #(u'These results suggest that AKT2 can regulate miR-200a in a histology- or stage-specific manner and that this regulation is independent of subsequent involvement of miR-200a in epithelial-mesenchymal transition.', 4, 24),
    #(u'We demonstrated that AKT2 is a direct target of miR-200, Spearman\'s rank correlation analysis showed that the expression levels of AKT2 and miR-200c in 35 pairs of osteosarcoma specimens were inversely correlated.', 3, 9),
    #(u'Zinc-finger enhancer binding (ZEB) transcription factors (ZEB1 and ZEB2) are crucial EMT activators, whereas members of the miR-200 family induce epithelial differentiation.', 9,22),
    (u'Thus, miR-200a appears to act as a multifunctional tumor suppressor miRNA in meningiomas through effects on the E-cadherin and Wnt/beta-catenin signaling pathways.', 2, 21),
    (u'Downregulated miR-200 in meningiomas promotes tumor growth by reducing E-cadherin and activating the Wnt/beta-catenin signaling pathway.', 1, 9),
    (u'Strikingly, cells with high levels of PKM2 expressed lower levels of miR-326, suggestive of endogenous regulation of PKM2 by miR-326.', 7, 21),
    #(u'Lack of miRNAs in post-mitotic neurons in vivo is associated with APP exons 7 and 8 inclusion, while ectopic expression of miR-124, an abundant neuronal-specific miRNA, reversed these effects in cultured neurons.', 11, 22),
    #(u'By incorporating miR-126 target sequences into a GALC-expressing vector, we suppressed GALC expression in HSCs while maintaining robust expression in mature hematopoietic cells.', 2, 13),
    #(u"Association studies using all common variants detected in the 3' UTR of BACE1 and the miR-29 gene cluster did not identify an association with AD risk.", 12, 15),
    #(u"Strikingly, cells with high levels of PKM2 expressed lower levels of miR-326, suggestive of endogenous regulation of PKM2 by miR-326.", 12, 19),
    #(u"Interaction of microRNA-338 and its potential targeting protein eiF4E3", 2, 8),
    #(u"Our results showed that the expression of miR-106b was inversely correlated with TβR II protein levels and miR-106b can directly inhibit the TβR II translation in vitro", 17, 22),
    #can the below be fixed? somehow?
    #(u"Genetic deletion of let-7, antagomir-mediated blockage of let-7 and miR-184* action, transgenic expression of dp target protector, or replacement of endogenous dp with a dp transgene non-responsive to let-7 each had toxic effects similar to those of pathogenic LRRK2.", 32, 42),
    #(u"These findings suggested that bcl2 is an important functional target for miR-34a, and the abnormal expression of miR-34a may contribute to the pathogenesis of Alzheimer's disease, at least in part by affecting the expression of bcl2.", 18, 38),
    #(u"Incubation of a protected antisense miRNA-146a was found to inhibit miRNA-146a and restore IRAK-1, whereas IRAK-2 remained unaffected.", 10, 16),
    #(u"Up-regulation of miRNA-146a coupled to down-regulation of CFH was observed in AD brain and in interleukin-1beta, Abeta42, and/or oxidatively stressed human neural (HN) cells in primary culture. ", 2, 17),
    #(u"MicroRNA-338 can recognize the 3'UTR of eiF4E3 while it has no significant effect on the expression of eiF4E3.", 0, 17),
    #(u"Pharmacologic and genetic inhibition of p53 corroborated the hypothesis that pre-miR-630 (but not pre-miR-181a) blocks the upstream signaling pathways that are ignited by DNA damage and converge on p53 activation.", 10, 30),
    #(u"miR-15a and miR-16-1 function by targeting multiple oncogenes, including BCL2, MCL1, CCND1, and WNT3A.", 2, 10),
    (u"miR-29a, decreased in Alzheimer disease brains, targets neurone navigator 3.  AIMS.", 0, 9),
    #(u"Using miRNA-146a-, IRAK-1-, or IRAK-2 promoter-luciferase reporter constructs, we observe decreases in IRAK-1 and increases in miRNA-146a and IRAK-2 expression in interleukin-1β (IL-1β) and amyloid-β-42 (Aβ42) peptide-stressed HAG cells.", 19,29),
    #(u"In human cell lines, we found that miR-205 down-regulates the expression of LRP1 by targeting sequences in the 3'UTR of LRP1 mRNA.", 8, 21),
    #(u"Using transiently transfected murine neuronal N2a cells in culture, in parallel to a mouse model of AD, we were able to demonstrate a role for two miRNAs (miR-298 and miR-328) in the regulation of beta-amyloid (Abeta) precursor protein (APP)-converting enzyme (BACE) messenger RNA (mRNA) translation, thereby providing key insights into the molecular basis underlying BACE deregulation in AD.", 32, 65),
    #(u"Association studies using all common variants detected in the 3' UTR of BACE1 and the miR-29 gene cluster did not identify an association with AD risk.", 13, 16),
    #(u"Our results showed that the expression of miR-106b was inversely correlated with TβR II protein levels and miR-106b can directly inhibit the TβR II translation in vitro.", 17, 23),
    #(u"A direct correlation was found between the downregulation of miR-200a and the upregulation of beta-catenin in human meningioma samples.", 9, 14),
    #(u"miR-9 targets REST and miR-9* targets CoREST.", 0,7),
    #(u"Adjusted miR-107 and BACE1 mRNA levels tended to correlate negatively (trend with regression P< 0.07).", 1, 3),
    (u"The database search on TargetScan, PicTar and miRBase Target identified neurone navigator 3 (NAV3), a regulator of axon guidance, as a principal target of miR-29a, and actually NAV3 mRNA levels were elevated in AD brains.", 29, 33),
    (u"BACE1 mRNA levels tended to increase as miR-107 levels decreased in the progression of AD. Cell culture reporter assays performed with a subset of the predicted miR-107 binding sites indicate the presence of at least one physiological miR-107 miRNA recognition sequence in the 3'-UTR of BACE1 mRNA.", 27, 46),
    (u"The expression of microRNA miR-107 decreases early in Alzheimer's disease and may accelerate disease progression through regulation of beta-site amyloid precursor protein-cleaving enzyme 1.", 4, 24),
    (u"Common variation in the miR-659 binding-site of GRN is a major risk factor for TDP43-positive frontotemporal dementia.", 4, 7),
    (u"We further demonstrate that miR-659 can regulate GRN expression in vitro, with miR-659 binding more efficiently to the high risk T-allele of rs5848 resulting in augmented translational inhibition of GRN.", 7, 13),
    (u"These findings support a causal link between let-7 and high-mobility group A2 whereby loss of let-7 expression induces high-mobility group A2 upregulation that represents an important mechanism in pituitary tumorigenesis and progression.", 7, 20),
    (u"Here we found that miR-106b and TGF-β type II receptor (TβR II) were aberrantly expressed in APPswe/PS∆E9 mice (a double transgenic mouse model for AD).", 4, 12),
    (u"Interestingly, both miRNAs are capable of binding directly to TDP-43 in different positions: within the miRNA sequence itself (let-7b) or in the hairpin precursor (miR-663).", 10, 21),
    (u"We further demonstrate that miRNAs belonging to the miR-15 family are potent regulators of ERK1 expression in mouse neuronal cells and co-expressed with ERK1/2 in vivo.", 8, 23),
    (u"A direct correlation was found between the downregulation of miR-200a and the upregulation of beta-catenin in human meningioma samples.", 9, 14),
    (u"Three bioinformatics-verified miRNA-128 targets, angiopoietin-related growth factor protein 5 (ARP5; ANGPTL6), a transcription suppressor that promotes stem cell renewal and inhibits the expression of known tumor suppressor genes involved in senescence and differentiation, Bmi-1, and a transcription factor critical for the control of cell-cycle progression, E2F-3a, were found to be up-regulated.", 2, 39),
(u"Reduction of miR-21 induces glioma cell apoptosis via activating caspase 9 and 3.  Extensive data indicate that miR-21 plays a critical role in gliomagenesis, however, knowledge is limited on the mechanism of action of miR-21, including cell proliferation, apoptosis, and migration.", 2, 10),


]

#29099614.2.11
mir = 'miR-486-5p'
lex = nlp.vocab[mir]
print(nlp.vocab[mir])

genes = ["neurone navigator 3", "Sprouty-related EVH1 domain-containing protein 1", "phosphoinositol-3 kinase regulatory subunit 2", "vascular cell adhesion molecule 1", "beta-catenin"]

for entity in genes+["miR-29a", "miR-200", "miR-201", "miR-202", "KLF12", "KLF13", "CXCR4", "miR-98", "ICAM1", "miR-200", "KLF4", "miR-101", "DUSP1", "miR-139", "miR-126", "miR-181", "JUN", "CD163", "vascular cell adhesion molecule 1"]:
    nlp.tokenizer.add_special_case(entity, [{'ORTH': entity, 'TAG': 'NNP'}])

nlp.tokenizer.add_special_case(u'miR-486-5p', [{'ORTH': 'miR-486-5p', 'TAG': 'NNP'}])
nlp.tokenizer.add_special_case(u'miR-100', [{'ORTH': 'miR-100', 'TAG': 'NNP'}])

testCases = [testCases[-1]]


relCheck = MirGeneRelCheck()


for testCase in testCases:

    if testCase[1] == 0 and testCase[2] == 0:
        continue

    doc = nlp(testCase[0])
    lWord = doc[testCase[1]]
    rWord = doc[testCase[2]]

    for t in doc:
        print(t.idx, t, t.pos_, t.dep_, t.head, t.head.idx)

    print(lWord, rWord)

    checkResults = relCheck.checkRelation(doc, lWord, rWord, verbose=True)
    print(checkResults)

    if len(testCases) == 1:
        displacy.serve(doc, style="dep", port=5005)
    print()
    print()
    print()