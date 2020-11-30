import os,sys

sys.path.insert(0, "/mnt/f/dev/git/miRExplore/python/")

import time

from textdb.MiGenRelDB import MiGenRelDB
from textdb.SentenceDB import SentenceDB

from collections import defaultdict
from natsort import natsorted

sentDB, _ = SentenceDB.loadFromFile("./test/", "./development/pmid2sent", returnAll=True, redoPmid2Sent=True)
mmuDB = MiGenRelDB.loadFromFile("./aggregated_test/mirna_gene.mmu.pmid", ltype="mirna", rtype="gene")
hsaDB = MiGenRelDB.loadFromFile("./aggregated_test/mirna_gene.hsa.pmid", ltype="mirna", rtype="gene")


referenceSolution = [

('16831872','miR-9','ONECUT2','MIR_GENE'),
('16831872','miR-9','SYTL4','MIR_GENE'),

#('17438130','miR-17-92','MYC','GENE_MIR'), # not a miRNA (miR-17-92 cluster)
('17438130','let-7c','MYC','MIR_GENE'),
('17438130','let-7c','MIR17HG','MIR_GENE'), #The PPARalpha-mediated induction of c-myc via let-7C subsequently increased expression of the oncogenic mir-17-92 cluster; these events did not occur in Pparalpha-null mice.


('18185580','miR-335','SOX4','MIR_GENE'),
('18185580','miR-335','TNC','MIR_GENE'),
('18755897','miR-34','TP53','GENE_MIR'), # wrong: acetylated TP53 // correct: Finally, miR-34a itself is a transcriptional target of p53, suggesting a positive feedback loop between p53 and miR-34a.
('18755897','miR-34a','TP53','GENE_MIR'),
('18755897','miR-34a','SIRT1','MIR_GENE'),
('18755897','miR-34','SIRT1','MIR_GENE'),
('18755897','miR-34','TP53','MIR_GENE'), # wrong: acetylated TP53 // correct: Finally, miR-34a itself is a transcriptional target of p53, suggesting a positive feedback loop between p53 and miR-34a.
('18755897','miR-34','CDKN1A','MIR_GENE'), #p21 #miR-34a
('18755897','miR-34','TPT1','MIR_GENE'), # p21
('18755897','miR-34','NSG1','MIR_GENE'), # p21
('18755897','miR-34','H3F3AP6','MIR_GENE'), # p21
('18755897','miR-34','TCEAL1','MIR_GENE'), # p21
('18755897','miR-34','BBC3','MIR_GENE'), #34a, has syn PUMA


('19059913','miR-223','SPI1','GENE_MIR'),
('19059913','miR-223','NFIA','MIR_GENE'),
('19059913','miR-223','NFIC','MIR_GENE'), # TM error, means NFIA
('19059913','miR-223','CSF1R','MIR_GENE'),

('19073597','miR-133a','MYOD1','GENE_MIR'),
('19073597','miR-133a','UCP2','MIR_GENE'),
('19073597','miR-133a','BMIQ4','MIR_GENE'), # has UCP-2 as syn

('19073597','miR-133a-mediated','UCP2','MIR_GENE'),
('19073597','miR-133a-mediated','BMIQ4','MIR_GENE'), # has UCP-2 as syn


('22066022', 'miR-21', 'GPT', 'MIR_GENE'), # missing mirtex: Serum miR-21 levels correlated with histological activity index (HAI) in the liver, alanine aminotransferase (ALT), aspartate aminotransferase , bilirubin, international normalized ratio and gamma-glutamyltransferase.

('19158092','miR-21','PDCD4','MIR_GENE'),
('19158092','miR-21-mediated','PDCD4','MIR_GENE'),

#('19378336','miR-145','KRT7','MIR_GENE'), #no interaction
('19378336','miR-30','KRT7','MIR_GENE'),
('19378336','miR-133a','KRT7','MIR_GENE'),
#('19378336','miR-133b','KRT7','MIR_GENE'), #no interaction in text
#('19378336','miR-195','KRT7','MIR_GENE'), #no interaction
#('19378336','miR-125b','KRT7','MIR_GENE'), #no interaction in text
('19378336','miR-199a','KRT7','MIR_GENE'),

('19524507','miR-31','RHOA','MIR_GENE'),
('19544458','miR-92b','CORO1A','MIR_GENE'), # was p57
('19625769','miR-101','EZH2','MIR_GENE'),

('19723773','miR-290','CDKN2A','MIR_GENE'), # p16

('19839716','miR-205','ERBB3','MIR_GENE'),
('19839716','miR-205','ZEB1','MIR_GENE'),



('19956414','miR-29b','COL1A1','MIR_GENE'),
('19956414','miR-29b','OI4','MIR_GENE'), #COL1A1
('19956414','miR-29b','COL1A2','MIR_GENE'),
('19956414','miR-29b','COL4A1','MIR_GENE'),
('19956414','miR-29b','COL5A1','MIR_GENE'),
('19956414','miR-29b','COL5A2','MIR_GENE'),
('19956414','miR-29b','COL3A1','MIR_GENE'),
('19956414','miR-29b','EDS4A','MIR_GENE'), #COL3A1
('19956414','miR-29b','LAMC1','MIR_GENE'),
('19956414','miR-29b','FBN1','MIR_GENE'),
('19956414','miR-29b','SPARC','MIR_GENE'),
('19956414','miR-29b','ON','MIR_GENE'), #osteonectin
('19956414','miR-29b','BMP1','MIR_GENE'),
('19956414','miR-29b','PCOLC','MIR_GENE'), #BMP1
('19956414','miR-29b','ADAM12','MIR_GENE'),
('19956414','miR-29b','NKIRAS2','MIR_GENE'),

('20012062','miR-221','PSMD9','MIR_GENE'), #p27
('20012062','miR-222','PSMD9','MIR_GENE'),
('20012062','miR-221','SSSCA1','MIR_GENE'), #p27
('20012062','miR-222','SSSCA1','MIR_GENE'),

('20017139','miR-146a','CNTN2','GENE_MIR'),
('20017139','miR-146a','NFKB1','GENE_MIR'),

#('20046097', 'miR-449', 'CDK', 'MIR_GENE'), # CDK not a gene symbol, not in mirtex
('20046097', 'miR-449', 'E2F1', 'MIR_GENE'), #not in mirtex :miR-449 regulates CDK-Rb-E2F1 through an auto-regulatory feedback circuit.
('20046097', 'miR-449', 'RB1', 'MIR_GENE'), #not in mirtex

('20103675','miR-222','PPP2R2A','MIR_GENE'),
('20143188','miR-21','PDCD4','MIR_GENE'),

('20299489','miR-34a','ERK','GENE_MIR'),
('20299489','miR-34a','EPHB2','GENE_MIR'),# ERK syn
('20299489','miR-34a','MAPK1','GENE_MIR'),# ERK syn
('20299489','miR-34a','MAP2K1','MIR_GENE'),
('20299489','miR-221','FOS','MIR_GENE'),
('20299489','miR-222','FOS','MIR_GENE'),
('20299489','miR-34a','FOSB','GENE_MIR'), #mirtex missing: induced miR-34a expression by transactivation via the activator protein-1 binding site in the upstream region of the miR-34a gene.
('20299489','miR-34a','JUND','GENE_MIR'), # activator protein 1 syn
('20299489','miR-34a','JUN','GENE_MIR'), #  induced miR-34a expression by transactivation via the activator protein-1 binding site


('20462046','miR-21','PDCD4','MIR_GENE'),
('20478254','miR-183','SLC1A1','MIR_GENE'),
('20478254','miR-96','SLC1A1','MIR_GENE'),
('20478254','miR-182','SLC1A1','MIR_GENE'),
('20498046','miR-200b','ATP2A2','MIR_GENE'),
('20498046','miR-214','ATP2A2','MIR_GENE'),

('20603081','miR-150','MYB','MIR_GENE'),

('20606648', 'miR-34a', 'BIRC5', 'MIR_GENE'), # missing in mirtex, miRNA-34a (miR-34a) induced apoptosis, inhibited survivin expression, and downregulated MAPK pathway in B16F10 cells.

('20620960','miR-200c','FAP','MIR_GENE'),
('20620960','miR-200','FAP','MIR_GENE'),
('20620960','miR-200c','GLMN','MIR_GENE'), # has FAP as syn
('20620960','miR-200','GLMN','MIR_GENE'), # has FAP as syn
('20620960','miR-200','FAS','MIR_GENE'), # CD95; quite indirect though. miR-200c regulates induction of apoptosis through CD95 by targeting FAP-1.
('20620960','miR-200c','FAS','MIR_GENE'), # CD95; quite indirect though. miR-200c regulates induction of apoptosis through CD95 by targeting FAP-1.

('20620960','miR-200','ZEB1','MIR_GENE'), # 200c
('20620960','miR-200','ZEB2','MIR_GENE'), # 200c
('20620960','miR-200','PPCD3','MIR_GENE'), # ZEB1

('20676061','miR-29c','WNT5A','MIR_GENE'),
('20676061','miR-130b','WNT5A','MIR_GENE'),
('20676061','miR-101','WNT5A','MIR_GENE'),
('20676061','miR-30b','WNT5A','MIR_GENE'),
('20676061','miR-140','WNT5A','MIR_GENE'),
('20676061','miR-29c','ZIC1','MIR_GENE'),
('20676061','miR-130b','ZIC1','MIR_GENE'),
('20676061','miR-101','ZIC1','MIR_GENE'),
('20676061','miR-30b','ZIC1','MIR_GENE'),
('20676061','miR-140','ZIC1','MIR_GENE'),
('20676061','miR-29c','TGFB1','MIR_GENE'),
('20676061','miR-130b','TGFB1','MIR_GENE'),
('20676061','miR-101','TGFB1','MIR_GENE'),
('20676061','miR-30b','TGFB1','MIR_GENE'),
('20676061','miR-140','TGFB1','MIR_GENE'),

('20676061','miR-29c','DPD1','MIR_GENE'), # has TGFB1 as syn
('20676061','miR-130b','DPD1','MIR_GENE'),
('20676061','miR-101','DPD1','MIR_GENE'),
('20676061','miR-30b','DPD1','MIR_GENE'),
('20676061','miR-140','DPD1','MIR_GENE'),

('20736365','miR-196','HOXC8','MIR_GENE'),
('20736365','miR-196','HOX3A','MIR_GENE'), # has syn HOXC8

('20859756', 'miR-126', 'TMEM8B', 'GENE_MIR'), # missing mirtex: In particular, miR-126, miR-142-3p, miR-155, miR-552, and miR-630 were all upregulated, whereas miR-146a, miR-152, miR-205, miR-365, miR-449, miR-518c, miR-584, miR-615, and miR-622 were downregulated after NGX6 transfection.
('20859756', 'miR-142', 'TMEM8B', 'GENE_MIR'),
('20859756', 'miR-155', 'TMEM8B', 'GENE_MIR'),
('20859756', 'miR-552', 'TMEM8B', 'GENE_MIR'),
('20859756', 'miR-630', 'TMEM8B', 'GENE_MIR'),
('20859756', 'miR-146a', 'TMEM8B', 'GENE_MIR'),
('20859756', 'miR-152', 'TMEM8B', 'GENE_MIR'),
('20859756', 'miR-205', 'TMEM8B', 'GENE_MIR'),
('20859756', 'miR-365', 'TMEM8B', 'GENE_MIR'),
('20859756', 'miR-449', 'TMEM8B', 'GENE_MIR'),
('20859756', 'miR-518c', 'TMEM8B', 'GENE_MIR'),
('20859756', 'miR-584', 'TMEM8B', 'GENE_MIR'),
('20859756', 'miR-615', 'TMEM8B', 'GENE_MIR'),
('20859756', 'miR-622', 'TMEM8B', 'GENE_MIR'),




('20945501', 'miR-141', 'AR', 'MIR_GENE'), # missing mirtex,  inhibition of miR-141 by anti-miR-141 suppressed the growth of the LNCaP subline overexpressing AR.
('20945501', 'miR-141', 'SBMA', 'MIR_GENE'), # has AR as syn
('20945501', 'miR-141', 'DHTR', 'MIR_GENE'), # has AR as syn
('20945501', 'miR-141', 'AKR1B3', 'MIR_GENE'), # has AR as syn
('20945501', 'miR-141', 'AKR1B7', 'MIR_GENE'), # has AR as syn
('20945501', 'miR-141', 'AKR1B8', 'MIR_GENE'), # has AR as syn
('20945501', 'miR-141', 'SBMA', 'MIR_GENE'), # has AR as syn
('20945501', 'miR-141', 'AREG', 'MIR_GENE'), # has AR as syn
('20945501', 'miR-141', 'FDXR', 'MIR_GENE'), # has AR as syn

('20947507','miR-155','NFKB1','GENE_MIR'),
('20947507','miR-155','CARD11','GENE_MIR'),
('20947507','miR-155','SPI1','MIR_GENE'), # PU.1
#('20947507','miR-155','CD10','MIR_GENE'),
('20947507','miR-155','MME','MIR_GENE'), # CD10


('21088996','miR-21','PDCD4','MIR_GENE'),
('21276775','miR-145','ROBO2','MIR_GENE'),
('21276775','miR-145','SRGAP2','MIR_GENE'),
('21276775','miR-145','SRGAP3','MIR_GENE'), # has SRGAP2 syn
('21276775','miR-214','ROBO2','MIR_GENE'),
('21276775','miR-214','SRGAP2','MIR_GENE'),
('21276775','miR-214','SRGAP3','MIR_GENE'),# has SRGAP2 syn

('21285947','miR-24','INS','MIR_GENE'),
('21285947','miR-26','INS','MIR_GENE'),
('21285947','miR-182','INS','MIR_GENE'),
('21285947','miR-148','INS','MIR_GENE'),

#('21347332','miR-21','serum','GENE_MIR'), # not a gene
('21347332','miR-21','FGF2','GENE_MIR'),
('21347332','miR-21','RHOB','MIR_GENE'),
('21415212','miR-486','OLFM4','MIR_GENE'),
('21454627','mmu-miR-183','mSEL-1L','MIR_GENE'),



('21609717','miR-98-mediated','IL10','MIR_GENE'),
('21609717','miR-98','IL10','MIR_GENE'),
('21609717','miR-98','PTGS2','MIR_GENE'), #COX-2
('21609717','miR-98', 'LPS', 'GENE_MIR'),
('21609717','miR-98', 'IRF6', 'GENE_MIR'), #missing mirtex, MicroRNA-98 negatively regulates IL-10 production and endotoxin tolerance in macrophages after LPS stimulation.

('21666774','miR-21','LH (luteinizing hormone)','GENE_MIR'),
('21666774','miR-132','LH (luteinizing hormone)','GENE_MIR'),
('21666774','miR-212','LH (luteinizing hormone)','GENE_MIR'),

('21685392','miR-143','NOTCH1','MIR_GENE'), #should be N1ICD, TM issue
('21685392','miR-145','NOTCH1','MIR_GENE'),
('21685392','miR-143','TAN1','MIR_GENE'), #should be N1ICD, TM issue
('21685392','miR-145','TAN1','MIR_GENE'),

('21685392','miR-143','RBPJ','GENE_MIR'), #We also identified N1ICD complex binding to CBF1 sites within the endogenous human miR-143/145 promoter.
('21685392','miR-145','RBPJ','GENE_MIR'),

('21685392','miR-143','JAG1','GENE_MIR'), #Using SRF knockdown, we found that Jag-1/Notch induction of miR-143/145 is SRF independent, although full acquisition of contractile markers requires SRF.
('21685392','miR-145','JAG1','GENE_MIR'),

('21685392','miR-143','SRF','GENE_MIR'), #Using SRF knockdown, we found that Jag-1/Notch induction of miR-143/145 is SRF independent, although full acquisition of contractile markers requires SRF.
('21685392','miR-145','SRF','GENE_MIR'),
('21685392','miR-145','MYOCD','GENE_MIR'), #The serum response factor (SRF)/myocardin complex binds to CArG sequences to activate miR-143/145 transcription
('21685392','miR-143','MYOCD','GENE_MIR'),

('21693621','miR-21','MYC','GENE_MIR'),
('21693621','miR-29a','MYC','GENE_MIR'),

('21898400','miR-520c','MTOR','MIR_GENE'),
('21898400','miR-373','MTOR','MIR_GENE'),
('21898400','miR-520c','SIRT1','MIR_GENE'),
('21898400','miR-373','SIRT1','MIR_GENE'),
('21898400','miR-520c','MMP9','MIR_GENE'),
('21898400','miR-373','MMP9','MIR_GENE'),
('21898400','miR-520c','CLG4B','MIR_GENE'), # MMP-9
('21898400','miR-373','CLG4B','MIR_GENE'),

('22123611','miR-195','BCL2','MIR_GENE'),
('22123611','miR-195','CASP3','MIR_GENE'),
('22123611','miR-195','WT1','MIR_GENE'), #missing mirtex: miR-195-treated podocytes underwent actin rearrangement and failed to synthesize sufficient levels of WT-1 and synaptopodin proteins, which suggests that the cells had suffered injuries similar to those observed in diabetic nephropathy in both humans and animal models.
('22123611','miR-195','SYNPO','MIR_GENE'),
('22123611','miR-195','GUD','MIR_GENE'),

('22139444','miR-30c','MTA1','MIR_GENE'),

('22249219','miR-214','ADORA2A','MIR_GENE'),
('22249219','miR-15','ADORA2A','MIR_GENE'),
('22249219','miR-16','ADORA2A','MIR_GENE'),

('22269326','miR-29b','COL1A1','MIR_GENE'),
('22269326','miR-29b','COL3A1','MIR_GENE'), 
('22269326','miR-29b','EDS4A','MIR_GENE'), #COL3A1 syn 
('22269326','miR-29b','COL5A1','MIR_GENE'),
('22269326','miR-29b','ELN','MIR_GENE'),


('22286762','miR-21','NF-kappaB','GENE_MIR'),
('22286762','miR-10b','NF-kappaB','GENE_MIR'),
('22286762','miR-17','NF-kappaB','GENE_MIR'),
('22286762','miR-9','NF-kappaB','GENE_MIR'),

('22569260','miR-223','FOXO1','MIR_GENE'),
('22569260','miR-223','FKHR','MIR_GENE'), # FOXO1 syn

('22634495','miR-10a','CHL1','MIR_GENE'),
('22634495','miR-10a','DDX11','MIR_GENE'), # CHL1 syn

('22698995','let-7','BACH1','MIR_GENE'),
('22698995','let-7b','BACH1','MIR_GENE'),
('22698995','let-7c','BACH1','MIR_GENE'),
('22698995','miR-98','BACH1','MIR_GENE'),
#same gene symbol
('22698995','let-7','BRIP1','MIR_GENE'),
('22698995','let-7b','BRIP1','MIR_GENE'),
('22698995','let-7c','BRIP1','MIR_GENE'),
('22698995','miR-98','BRIP1','MIR_GENE'),

('22698995','let-7','HMOX1','MIR_GENE'),

('22761336','miR-96','REV1','MIR_GENE'),
('22761336','miR-96','RAD51','MIR_GENE'),
('22761336','miR-96','RECA','MIR_GENE'), #RAD51
('22761336','miR-96','RAD51A','MIR_GENE'),

('22847613','miR-130b','TP53','GENE_MIR'),
('22847613','miR-130b','ZEB1','MIR_GENE'),
('22847613','miR-130b','PPCD3','MIR_GENE'), # zeb1

('22891274','miR-146a','NFKB1','GENE_MIR'),
('22891274','miR-146a','NFKB1','MIR_GENE'),
('22891274','miR-146a','TRAF6','MIR_GENE'),
('22891274','miR-146a','IRAK1','MIR_GENE'),

('22925189','miR-30c','ERBB2','MIR_GENE'), #Her-2
('22925189','miR-30d','ERBB2','MIR_GENE'),
('22925189','miR-30e','ERBB2','MIR_GENE'),
('22925189','miR-532','ERBB2','MIR_GENE'),

('22955854','miR-144','ZFX','MIR_GENE'),
('22956424','miR-21','PTEN','MIR_GENE'),
('22956424','miR-21','MHAM','MIR_GENE'), # has PTEN syn
('22956424','miR-21','BZS','MIR_GENE'), # has PTEN syn


('22982443','miR-200c','BMI1','MIR_GENE'),
('22982443','miR-200c','ABCG2','MIR_GENE'),
('22982443','miR-200c','ABCG5','MIR_GENE'),
('22982443','miR-200c','MDR1','MIR_GENE'),
('22982443','miR-200c','TBC1D9','MIR_GENE'), # has syn MDR1
('22982443','miR-200c','ABCB1','MIR_GENE'), # has syn MDR1
('22982443','miR-200c','CDH1','MIR_GENE'),


('23010597','miR-134','FOXM1','MIR_GENE'),
('23010597','miR-134','FKHL16','MIR_GENE'), # has FOXM1 as syn
('23010597','miR-134','ITK','MIR_GENE'), # has EMT as syn; mirtex missing: Functional assays demonstrated that miR-134 inhibited EMT in NSCLC cells.
('23010597','miR-134','SLC22A3','MIR_GENE'), # has EMT as syn

('23041385','miR-21','CRP','MIR_GENE'),
#('23041385','miR-21','fibrinogen','MIR_GENE'), # not a gene
('23041385','miR-21','TGFB2','MIR_GENE'),
('23097316','miR-34c','RARg','MIR_GENE'),

('23113351','miR-29','TP53','MIR_GENE'), # missing in mirtex: While miRNA-29 members induced apoptosis through p53 gene activation, the effect of miRNA-29a on osteoblastic cells was independent on p53 expression level.
('23113351','miR-29a','TP53','MIR_GENE'),
('23113351','miR-29','BCL2','MIR_GENE'),
('23113351','miR-29','MCL1','MIR_GENE'),
('23113351','miR-29','CLEC4D','MIR_GENE'), # CLEC4D has syn mcl
('23113351','miR-29a','CLEC4D','MIR_GENE'), # CLEC4D has syn mcl
#('23113351','miR-29','E2F1','MIR_GENE'),
#('23113351','miR-29','E2F3','MIR_GENE'),

('23113351','miR-29a','BCL2','MIR_GENE'),
('23113351','miR-29a','MCL1','MIR_GENE'),
('23113351','miR-29a','E2F1','MIR_GENE'),
('23113351','miR-29a','E2F3','MIR_GENE'),

('23113351','miR-29a','E2F1','MIR_GENE'), # possibly too
('23113351','miR-29a','E2F3','MIR_GENE'), # possibly too

('23148210','miR-210','HIF1A','GENE_MIR'), # was actived in ...-dependant

('23169590','miR-451','IL6','GENE_MIR'),
('23169590','miR-451','IFNA1','GENE_MIR'), #type I IFN
('23169590','miR-451','YWHAZ','MIR_GENE'),
('23169590','miR-451','YWHAD','MIR_GENE'), # has syn YWHAZ
#('23169590','miR-451','14-3-3zeta','MIR_GENE'), # is YWHAZ
('23169590','miR-451','ZFP36','MIR_GENE'),

#Three types of primary DCs treated with antisense RNA antagomirs directed against miR-451 secreted elevated levels of IL-6, TNF, CCL5/RANTES, and CCL3/MIP1alpha, and these results were confirmed using miR-451(null) cells.
#this suggests that miR-451 suppresses these genes normalls
('23169590','miR-451','IL6','MIR_GENE'),
('23169590','miR-451','CCL3','MIR_GENE'),
('23169590','miR-451','CCL5','MIR_GENE'),
('23169590','miR-451','IL6','MIR_GENE'),
('23169590','miR-451','IFNB2','MIR_GENE'), #IL6
('23169590','miR-451','TNF','MIR_GENE'),
('23169590','miR-451','TNFA','MIR_GENE'), #TNF

#miR-451 levels are themselves increased by IL-6 and type I IFN, potentially forming a regulatory loop.
('23169590','miR-451','IL6','GENE_MIR'), #IL6
('23169590','miR-451','IFNA1','GENE_MIR'), #type I IFN
('23169590','miR-451','IFNB2','GENE_MIR'), #IL6




('23190607','miR-203','RAN','MIR_GENE'),
('23190607','miR-203','RAPH1','MIR_GENE'),
('23190607','miR-203','ALS2CR18','MIR_GENE'), # synonym

('23190608','miR-29b','SP1','GENE_MIR'),
('23190608','miR-29b','SP1','MIR_GENE'), # mirtex missing: miR-29b sensitizes multiple myeloma cells to bortezomib-induced apoptosis through the activation of a feedback loop with the transcription factor Sp1.
('23190608','miR-29b','DAND5','GENE_MIR'), # DAND5 has SP1 as syn
('23190608','miR-29b','DAND5','MIR_GENE'),


('23206698','miR-7','IRS2','MIR_GENE'),
('23396109','miR-17','PTEN','MIR_GENE'), #miR-17~92
('23396109','miR-17','MHAM','MIR_GENE'), # MHAM syn is PTEN too
('23396109','miR-17','BZS','MIR_GENE'), # BZS syn is PTEN too
#('23396109','miR-17','BIM','MIR_GENE'), # miR-17~92 is pten
('23396109','miR-17','BCL2L11','MIR_GENE'), # BCL2L11 syn is BIM

('23472202','miR-183','TAOK1','MIR_GENE'),
('23516615','miR-143','ERK5','MIR_GENE'),
('23516615','miR-143','PPARg','MIR_GENE'),
('23516615','miR-204','ERK5','MIR_GENE'),
('23516615','miR-204','PPARg','MIR_GENE'),
('23519125','miR-125a','erbB2','MIR_GENE'),
('23519125','miR-125a','erbB3','MIR_GENE'),
('23519125','miR-125b','erbB2','MIR_GENE'),
('23519125','miR-125b','erbB3','MIR_GENE'),
('23519125','miR-205','erbB2','MIR_GENE'),
('23519125','miR-205','erbB3','MIR_GENE'),

# these are no genes!
#('23527070','miR-21','collagen I','MIR_GENE'),
#('23527070','miR-21','collagen III','MIR_GENE'),

('23527070','miR-21','ELN','MIR_GENE'),
('23527070','miR-21','SMAD7','MIR_GENE'),
('23527070','miR-21','SMAD5','MIR_GENE'),
('23527070','miR-21','SMAD2','MIR_GENE'),

#('23534973','miR-152','HLA-G','MIR_GENE'), not a specific miRNA: miR-152 family

('23579289','miR-214','SP7','MIR_GENE'),
('23583389','miR-96','IRS1','MIR_GENE'),

('23592910','miR-146a','IL1B','GENE_MIR'),
('23592910','miR-146a','IFNG','GENE_MIR'),
#('23592910','miR-146a','TNFA','GENE_MIR'), # is TNF
('23592910','miR-146b','IL1B','GENE_MIR'),
('23592910','miR-146b','IFNG','GENE_MIR'),
('23592910','miR-146b','TNF','GENE_MIR'),
('23592910','miR-146a','IRAK','MIR_GENE'),
('23592910','miR-146b','IRAK','MIR_GENE'),

('23611780','miR-106b','FBXW11','MIR_GENE'), #beta-TRCP2
#('23611780','miR-106b','SNAIL','MIR_GENE'), nope. means cluster + indirect: miR-106b-25 cluster may play an important role in the metastasis of human non-small cell lung cancer cells by directly suppressing the beta-TRCP2 gene expression with a consequent increase in the expression of Snail.
('23611780','miR-93','FBXW11','MIR_GENE'),
('23630358','miR-155','MSR1','MIR_GENE'), # SR-AI syn

('23643257','miR-424','FGR','MIR_GENE'),
('23643257','miR-424','MAP2K1','MIR_GENE'),
('23643257','miR-424','MAPK1','MIR_GENE'), #mitogen-activated protein kinase 1

('23667495','miR-224','DPYSL2','MIR_GENE'),
('23667495','miR-224','KRAS','MIR_GENE'),
('23667495','miR-452','DPYSL2','MIR_GENE'),
('23667495','miR-452','KRAS','MIR_GENE'),
('23667495','miR-181c','KRAS','MIR_GENE'),
('23667495','miR-340','MECP2','MIR_GENE'),
('23667495','miR-181c','MECP2','MIR_GENE'),
('23667495','miR-340','KRAS','MIR_GENE'),
('23759586','miR-34a','SIRT1','MIR_GENE'),
('23759586','miR-125b','TP53','MIR_GENE'),
('23759586','miR-125b','SIRT1','MIR_GENE'),
('23797704','miR-21','TIMP3','MIR_GENE'),
('23797704','miR-221','TIMP3','MIR_GENE'),
('23797704','miR-21','SFD','MIR_GENE'),#is TIMP3
('23797704','miR-221','SFD','MIR_GENE'),#is TIMP3
('23797704','miR-217','TIMP3','MIR_GENE'),
('23797704','miR-217','SFD','MIR_GENE'), #is TIMP3
('23797704','miR-217','SIRT1','MIR_GENE'),
('23836497','miR-20','STAT3','MIR_GENE'),
('23836497','miR-20','CCND1','MIR_GENE'),
('23836497','miR-106a','STAT3','MIR_GENE'),
('23836497','miR-106a','CCND1','MIR_GENE'),

# same genesymbols
('23846856','miR-875','PRDX3','MIR_GENE'),
('23846856','miR-875','PRX','MIR_GENE'),

('23851184','miR-200b','WNT1','MIR_GENE'),
('23851184','miR-22','WNT1','MIR_GENE'),
('23895517','mir-494','TNFSF14','MIR_GENE'),
('23895517','mir-197','TNFSF14','MIR_GENE'),

('23968734','miR-133a','PDLIM5','MIR_GENE'), # LIM
#('23968734','miR-133a','SH3 protein 1','MIR_GENE'),
('23968734','miR-133a','LASP1','MIR_GENE'), #SH3 protein 1


('24006456','miR-29b','IGF1','MIR_GENE'),
('24006456','miR-30c','IGF1','MIR_GENE'),
('24006456','miR-29b','LIF','MIR_GENE'),
('24006456','miR-30c','LIF','MIR_GENE'),
('24006456','miR-29b','PTX3','MIR_GENE'),



('24023867','miR-135a','NR3C2','MIR_GENE'),
('24023867','miR-124','NR3C2','MIR_GENE'),

('24145190','miR-203','SNAI1','GENE_MIR'),
('24145190','miR-203','CD44','MIR_GENE'), # new, not in mirtex: we found that the levels of several EMT activators and miR-203 were positively and negatively correlated with those of CD44, respectively.
('24145190','miR-203','MDU3','MIR_GENE'),
('24145190','miR-203','MIC4','MIR_GENE'),
('24145190','miR-203','MDU2','MIR_GENE'),

('24145190','miR-203','SRC','GENE_MIR'), # missing in mirtex: Finally, we discovered that c-Src kinase activity was required for the downregulation of miR-203

('24155920','miR-21','SPRY1','MIR_GENE'),
('24155920','miR-29a','MCL1','MIR_GENE'),
('24155920','miR-29b','MCL1','MIR_GENE'),

#('24219008','miR-21-5p','TGFBR3','MIR_GENE'),
('24219008','miR-21','TGFBR3','MIR_GENE'), #add
('24219008','hsa-miR-21','TGFBR3','MIR_GENE'), #add

#('24219008','miR-21-5p','PDGFD','MIR_GENE'),
('24219008','miR-21','PDGFD','MIR_GENE'),
('24219008','hsa-miR-21','PDGFD','MIR_GENE'), #add


#('24219008','miR-21-5p','PPM1L','MIR_GENE'),
('24219008','miR-21','PPM1L','MIR_GENE'),
('24219008','hsa-miR-21','PPM1L','MIR_GENE'), #add


#('24219008','miR-181a-5p','ROPN1L','MIR_GENE'),
('24219008','miR-181a','ROPN1L','MIR_GENE'),
('24219008','hsa-miR-181a','ROPN1L','MIR_GENE'),

#('24219008','miR-181a-5p','SLC37A3','MIR_GENE'),
('24219008','hsa-miR-181a','SLC37A3','MIR_GENE'),
('24219008','miR-181a','SLC37A3','MIR_GENE'),

#('24219008','miR-24-2-5p','MYC','MIR_GENE'),
('24219008','hsa-miR-24-2','MYC','MIR_GENE'),
('24219008','hsa-miR-24','MYC','MIR_GENE'),
('24219008','miR-24-2','MYC','MIR_GENE'),
('24219008','miR-24','MYC','MIR_GENE'),

#('24219008','miR-24-2-5p','KCNJ2','MIR_GENE'),
('24219008','hsa-miR-24-2','KCNJ2','MIR_GENE'),
('24219008','hsa-miR-24','KCNJ2','MIR_GENE'),
('24219008','miR-24','KCNJ2','MIR_GENE'),
('24219008','miR-24-2','KCNJ2','MIR_GENE'),


('24219349','miR-203','BMI1','MIR_GENE'),
('24220339','miR-490','FOS','MIR_GENE'),
('24223656','miR-31','RASA1','MIR_GENE'),
('24314216','miR-106','TP53','MIR_GENE'),
('24319262','miR-34a','TP53','GENE_MIR'),
('24319262','miR-145','TP53','GENE_MIR'),
('24319262','miR-155','MAF','MIR_GENE'),
('24319262','miR-34a','TWIST2','MIR_GENE'),
('24319262','miR-34a','MAF','MIR_GENE'),
('24319262','miR-145','TWIST2','MIR_GENE'),
('24319262','miR-145','MAF','MIR_GENE'),
('24330780','miR-124','FLOT1','MIR_GENE'),
('24330780','miR-124','FLOT1','MIR_GENE'),

('24376808','miR-146a','CRK','MIR_GENE'),
('24376808','miR-424','CRK','MIR_GENE'),
('24376808','miR-146a','EGFR','MIR_GENE'),
('24376808','miR-424','EGFR','MIR_GENE'),

('24376808','miR-146a','MAPK14','MIR_GENE'), #p38 / ERK
('24376808','miR-424','MAPK14','MIR_GENE'),
('24376808','miR-146a','AIMP2','MIR_GENE'),#p38 / ERK
('24376808','miR-424','AIMP2','MIR_GENE'),
('24376808','miR-146a','AHSA1','MIR_GENE'),#p38 / ERK
('24376808','miR-424','AHSA1','MIR_GENE'),
]

refDict = defaultdict(set)
for x in referenceSolution:
    refDict[x[0]].add((x[1], x[2], x[3]))


tmRemoveTMErrors = {
    ('19956414','miR-29b','MMRN1'), # ECM, extracellular matrix

('21415212','miR-486','GC'), #gastric cancer
('21415212','miR-486','HTC2'), #Array-CGH
('21415212','miR-486','EAF2'), #TRAITS
('21415212','miR-486','NF2'), #SCH cell line

    ('21703983', 'miR-632', 'PAFAH1B1'), # Notably, hsa-miR-378, hsa-miR-632, and hsa-miR-636 demonstrated particularly high discrimination between MDS and normal controls. MDS here is myelodysplastic syndromes
    ('21703983', 'hsa-miR-378', 'PAFAH1B1'),
    ('21703983', 'miR-378', 'PAFAH1B1'),
    ('21703983', 'hsa-miR-632', 'PAFAH1B1'),
    ('21703983', 'hsa-miR-636', 'PAFAH1B1'),
    ('21703983', 'miR-636', 'PAFAH1B1'),

    ('22066022', 'miR-21', 'FAM126A'), # HCC refers to hepatocellular carcinoma
    ('22066022', 'miR-21', 'ST14'), # HAI refers to histological activity index (HAI)
    ('22066022', 'miR-21', 'SPINT1'), # HAI refers to histological activity index (HAI)

    ('23643257','miR-424','FGR'), # recognizes FGR
    ('23643257','miR-424','FGFR1'),
    ('23643257','miR-424','KAL2'),
    ('23643257','miR-424','FLT2'),

    ('24330780','miR-124','TENM1'), #tumor node metastasis (TNM)
    ('24330780','miR-124','TNM'),

    ('24223656', 'miR-31', "TPT1"), # RAS p21 GTPase activating protein 1 (RASA1)  => p21
    ('24223656', 'miR-31', "CDKN1A"),
    ('24223656', 'miR-31', "H3F3AP6"),
    ('24223656', 'miR-31', "TCEAL1"),
    ('24223656', 'miR-31', "NSG1"),

    ('21609717','miR-98','MT-CO2'), # accepts COX-2
    ('21609717','miR-98','COX8A'),
    ('21609717','miR-98','CPOX'),
    ('21609717','miR-98','MT-CO2'),

    ('23527070', 'miR-21', 'SMAD5'), # SMAD2/5

    ('23190608','miR-29b','SUPT20H'), # SUPT20H has transcription-factor as syn
    ('23190608','miR-29b','SUPT20H'),

    ('18185580','miR-335', 'SUPT20H'),

    ('23113351','miR-29','RB1'), # RB1 syn: osteosarcoma
    ('23113351','miR-29a','RB1'), # RB1 syn: osteosarcoma
    ('23113351','miR-29b','RB1'), # RB1 syn: osteosarcoma

    ('22982443','miR-200c','CDH17'), # CDH17 syn for cadherin, found in E-cadherin ...
    ('20603081','miR-150','GLI2'), # THP-1 refers to cells

    ('22139444', 'miR-30c', 'NDC80'), # refers to HEC-1-B cells ...
    ('23041385', 'miR-21', 'CO'), # centenarian offspring (CO)
    ('23041385', 'miR-21', 'CALCR'), # CTR control

    ('19723773','miR-290','MEF'),#mouse embryo fibroblasts (MEF)
    ('19723773','miR-290','MEFV'),#mouse embryo fibroblasts (MEF)
    ('19723773','miR-290','ELF4'),#mouse embryo fibroblasts (MEF)

    ('19547998', 'miR-21', 'CALR'), # SSA
    ('19547998', 'miR-181b', 'CALR'),
    ('19547998', 'miR-21', 'HP'), # hyperplastic polyps
    ('19547998', 'miR-181b', 'HP'),

    ('20103675', 'miR-222', 'FAM126A'), # HCC cell lines

    ('24219349','miR-203','SP'), # side population
    ('24219349','miR-203','TFF2'), # SP

    ('21088996','miR-21','BLOC1S6'), # PDAC cells (MIA-Pa-Ca-2)
    ('21088996','miR-21','MIA'), # PDAC cells (MIA-Pa-Ca-2)
    ('21088996','miR-21','CAR2'), # PDAC cells (MIA-Pa-Ca-2)

    ('22847613','miR-130b','SLC22A3'), # epithelial-mesenchymal transition (EMT)
    ('22847613','miR-130b','ITK'), # epithelial-mesenchymal transition (EMT)

    ('22925189', 'miR-370', 'II'), #stage II <=> gene symbol
    ('22925189', 'miR-370', 'IV'), #stage IV <=> gene symbol
    ('22925189', 'miR-30a', 'II'), #stage II <=> gene symbol
    ('22925189', 'miR-30a', 'IV'), #stage IV <=> gene symbol

    ('23592910', 'miR-146a', 'IFNA1'), # TM mismatch with interferon in interfon gamma
    ('23592910', 'miR-146b', 'IFNA1'),

    ('24006456','miR-29b','INS'), # spurious hit with insulin-like growth factor
    ('24006456','miR-30c','INS'), # insulin-like

    ('20945501', 'miR-141', 'PC'), # matches PC / prostate cancer
    ('20945501', 'miR-141', 'PODXL'), # matches CRPC (castration-resitant prostate cancer)
}


# gene-mir: 36 F1: 0.88 *0.135 = 0,1188
# mir-gene: 230 F1: 0.94 *0.865 = 0,8131 
# all: F1: 0.9319

sent2rels = defaultdict(set)

allSents = sentDB.get_all_sentences()

doc2Rels = defaultdict(set)

for mirID in mmuDB.ltype2rel:

    for rel in mmuDB.ltype2rel[mirID]:

        jel = rel.toJSON()
        sent2rels[rel.assocSent].add( (rel.lid, rel.rid, rel.assocInt, rel.assocCat, rel.lPOS, rel.rPOS) )

        docID = rel.assocSent.split(".")[0]
        doc2Rels[docID].add((rel.lid, rel.rid, rel.assocInt) )


for mirID in hsaDB.ltype2rel:

    for rel in hsaDB.ltype2rel[mirID]:

        jel = rel.toJSON()
        sent2rels[rel.assocSent].add( (rel.lid, rel.rid, rel.assocInt, rel.assocCat, rel.lPOS, rel.rPOS) )
        #print(rel.assocSent, rel.lid, rel.rid, rel.assocInt, rel.assocCat, relSent)

        docID = rel.assocSent.split(".")[0]
        doc2Rels[docID].add((rel.lid, rel.rid, rel.assocInt) )

from collections import Counter

#TM, REF
elemCaseCounter = Counter()

with open("test_list.bydoc.tsv", "w") as fout:

    print("doc", "lid", "rid", "assocInt", sep="\t", file=fout)

    allDocIDs = set([x.split(".")[0] for x in allSents])

    for docID in natsorted(doc2Rels):
        for elems in doc2Rels[docID]:
            print(docID, *elems, sep="\t", file=fout)

        refOnly = refDict[docID].difference(doc2Rels[docID])
        tmOnly = doc2Rels[docID].difference(refDict[docID])

        tmOnly = [x for x in tmOnly if not (docID, x[0], x[1]) in tmRemoveTMErrors]

        correct = refDict[docID].intersection(doc2Rels[docID])

        for x in correct:
            elemCaseCounter[(True, True)] += 1

        for x in refOnly:
            elemCaseCounter[(False, True)] += 1

        for x in tmOnly:
            elemCaseCounter[(True, False)] += 1

        if len(doc2Rels[docID]) == 0 and len(refDict[docID]):
            continue
        
        if len(refOnly) == 0 and len(tmOnly) == 0:
            continue

        print(docID, len(correct), "REFONLY", refOnly)
        print(docID, len(correct), "TMONLY", tmOnly)
        print()


precision = elemCaseCounter[(True, True)] / (elemCaseCounter[(True, True)]+elemCaseCounter[(True, False)])
recall = elemCaseCounter[(True, True)] / (elemCaseCounter[(True, True)]+elemCaseCounter[(False, True)])

f1 = 2* precision * recall / (precision+recall)

#specificity = elemCaseCounter[(False, False)] / (elemCaseCounter[(True, False)] + elemCaseCounter[(False, False)])

print()
print()
print("True, True", elemCaseCounter[(True, True)])
print("TM Only", elemCaseCounter[(True, False)])
print("Ref Only", elemCaseCounter[(False, True)])
print()
print("precision", precision)
print("recall", recall)
#print("specificity", specificity)
print("f1", f1)



with open("test_list.tsv", "w") as fout:

    print("lid", "rid", "assocInt", "assocCat", "lpos", "rpos", "int_eval", "cat_eval", "sentID", "sent", sep="\t", file=fout)

    for sentID in natsorted([x for x in allSents]):
        sent = allSents[sentID]

        allElems = sent2rels.get(sentID, None)

        if allElems == None:
            allElems = [ ("", "", "", "", "", "") ]

        for lid, rid, assocInt, assocCat, lpos, rpos in allElems:

            print(lid, rid, assocInt, assocCat, lpos, rpos, "FALSE", "FALSE", sentID, sent, sep="\t", file=fout)

