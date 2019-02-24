import sys, os
from collections import Counter, OrderedDict, defaultdict

import networkx
from natsort import natsorted

from synonymes.mirnaID import miRNAPART

sys.path.insert(0, str(os.path.dirname("/mnt/d/dev/git/poreSTAT/")))


from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE

from textdb.makeNetworkView import DataBasePlotter
from utils.cytoscape_grapher import CytoscapeGrapher
from utils.tmutils import normalize_gene_names

import pandas as pd
from pandas.tools.plotting import parallel_coordinates
import matplotlib.pyplot as plt

figidx = 0



cbn2edges = {
"CV-IPN-Endothelial_cell-monocyte_interaction_1": {('CXCL16', 'CXCL16'), ('ITGAX', 'ITGAX'), ('CX3CR1', 'CX3CR1'), ('CCL2', 'PIK3CA'), ('ITGAX', 'ITGB2'), ('CXCR6', 'CXCR6'), ('SIRPA', 'ITGB2'), ('ITGB2', 'ITGB2'), ('ITGA4', 'VCAM1'), ('CCL2', 'CCL2'), ('SELE', 'SELE'), ('IL1A', 'CCR2'), ('CCL2', 'RHOA'), ('CX3CL1', 'STAT5A'), ('LTB4R', 'ITGAM'), ('IL8', 'ITGAM'), ('ITGA4', 'ITGA4'), ('CX3CR1', 'ICAM1'), ('ICAM1', 'ETS1'), ('SELP', 'SELPLG'), ('ITGA4', 'ITGB1'), ('ITGAM', 'ITGAM'), ('THY1', 'THY1'), ('CX3CL1', 'CX3CL1'), ('ITGAL', 'ITGB2'), ('ITGAM', 'ICAM1'), ('CX3CL1', 'CX3CR1'), ('ICAM1', 'ITGB2'), ('IFNG', 'CXCL16'), ('ETS1', 'CCL2'), ('ICAM1', 'ICAM1'), ('IL8', 'ITGB2'), ('CXCR6', 'CXCL16'), ('ITGB2', 'ITGAL'), ('STAT5A', 'ICAM1'), ('SELP', 'SELP'), ('ETS1', 'ETS1'), ('SELPLG', 'SELPLG'), ('ITGAL', 'ITGAL'), ('IL8', 'IL8'), ('STAT5A', 'STAT5A'), ('CCR2', 'CCR2'), ('ITGB2', 'ITGAM'), ('VCAM1', 'VCAM1'), ('CCL2', 'CCR2'), ('STAB1', 'STAB1')},
"CV-IPN-Endothelial_cell_activation_1": {('CSF3', 'CSF3'), ('SOCS2', 'STAT3'), ('CRP', 'OLR1'), ('VEGFA', 'EGR1'), ('OLR1', 'MMP14'), ('HMOX1', 'CCL2'), ('CAV1', 'NOS3'), ('IL1B', 'VCAM1'), ('MMP14', 'RAC1'), ('ATF4', 'VEGFA'), ('RHOA', 'RHOA'), ('RHOA', 'NOS3'), ('TRAF2', 'TRAF2'), ('VEGFA', 'VEGFA'), ('MMP1', 'MMP1'), ('VEGFA', 'PTK2'), ('PPARG', 'ICAM1'), ('CXCL12', 'VEGFA'), ('OLR1', 'CD40'), ('VEGFA', 'PLVAP'), ('RAC1', 'RAC1'), ('GSK3B', 'GATA6'), ('PRKCZ', 'PRKCZ'), ('PTEN', 'PTEN'), ('PPARA', 'PPARA'), ('CD40', 'CD40'), ('ETS1', 'ETS1'), ('PLCG1', 'PLCG1'), ('OLR1', 'OLR1'), ('S1PR1', 'S1PR1'), ('RELA', 'CAV1'), ('PRKACA', 'FOXO1'), ('CD40', 'MAPK8'), ('MAPK3', 'MAPK3'), ('VCAM1', 'VCAM1'), ('OLR1', 'ICAM1'), ('JUN', 'JUN'), ('CXCL12', 'EGR1'), ('TNFRSF1B', 'TNFRSF1B'), ('STAT3', 'IL8'), ('IL1B', 'SELE'), ('AKT1', 'CCL2'), ('JAK2', 'STAT3'), ('ATF4', 'IL6'), ('CSF2', 'CSF2'), ('AGTR1', 'AGTR1'), ('AGTR1', 'CCL2'), ('FOXO1', 'NOS2'), ('FOXO1', 'VCAM1'), ('CCL2', 'CCL2'), ('CD40LG', 'EGR1'), ('CDKN1A', 'EP300'), ('NR3C2', 'NOX4'), ('PRKCB', 'MMP1'), ('IL6', 'JAK2'), ('HMOX1', 'VCAM1'), ('REL', 'F3'), ('VEGFA', 'VWF'), ('TNF', 'SELE'), ('CD40', 'MAPK3'), ('CXCL2', 'CXCL2'), ('TNF', 'VWF'), ('OLR1', 'IL8'), ('STAT3', 'ICAM1'), ('CSF2', 'SELP'), ('PTK2', 'PTK2'), ('ROCK2', 'ROCK2'), ('TNF', 'MAPK3'), ('STAT1', 'ICAM1'), ('IFNG', 'CD40'), ('CXCL3', 'CXCL3'), ('IL10', 'PECAM1'), ('SCARB1', 'NOS3'), ('NOS3', 'NOS3'), ('TNFRSF1B', 'TRAF2'), ('MMP14', 'MMP14'), ('TNFRSF1A', 'TRAF2'), ('SP1', 'SP1'), ('ANGPT1', 'PRKCZ'), ('EGR1', 'EGR1'), ('RHOA', 'PTK2'), ('PPARG', 'VCAM1'), ('RAC1', 'PTK2'), ('CSF1', 'CSF1'), ('PRKCB', 'MMP3'), ('CD40', 'PLCG1'), ('REL', 'REL'), ('GATA6', 'VCAM1'), ('TRAF2', 'CHUK'), ('IL1A', 'CD40'), ('VEGFA', 'RELA'), ('GSK3B', 'GSK3B'), ('OLR1', 'CSF3'), ('RELA', 'RELA'), ('IL1B', 'VWF'), ('IL1B', 'CD40'), ('PPARG', 'PPARG'), ('ITGAV', 'CAV1'), ('EGR1', 'PECAM1'), ('ITGAV', 'ITGAV'), ('TNF', 'CCL2'), ('EIF2AK3', 'EIF2AK3'), ('TNF', 'GATA6'), ('TNF', 'ICAM1'), ('CREB1', 'CREB1'), ('FOXO1', 'FOXO1'), ('ITGAV', 'ITGB3'), ('SELE', 'SELE'), ('CD40LG', 'REL'), ('TNF', 'VCAM1'), ('S1PR1', 'MAPK3'), ('VWF', 'VWF'), ('S1PR1', 'MAPK1'), ('EP300', 'EP300'), ('ATF4', 'CCL2'), ('STAT1', 'STAT1'), ('TNF', 'CD40'), ('TNF', 'CDKN1A'), ('IFNB1', 'CD40'), ('GATA6', 'GATA6'), ('MMP14', 'RHOA'), ('HMGB1', 'AGER'), ('NOS2', 'NOS2'), ('PRKCB', 'PRKCB'), ('SPHK1', 'SPHK1'), ('F3', 'F3'), ('AGER', 'AGER'), ('ICAM1', 'ICAM1'), ('OLR1', 'CXCL3'), ('SCARB1', 'SCARB1'), ('NR3C2', 'NR3C2'), ('IL6', 'CAV1'), ('EP300', 'ETS1'), ('RELA', 'ICAM1'), ('SELP', 'SELP'), ('MMP3', 'MMP3'), ('CD40LG', 'CD40'), ('MAPK8', 'MAPK8'), ('EIF2S1', 'EIF2S1'), ('AKT1', 'AKT1'), ('IL6', 'IL6'), ('EIF2AK3', 'EIF2S1'), ('JAK2', 'STAT1'), ('GATA3', 'GATA3'), ('ATF6', 'ATF6'), ('PLVAP', 'PLVAP'), ('PECAM1', 'PECAM1'), ('TNF', 'TNF'), ('ROCK2', 'CXCL2'), ('PPARA', 'VCAM1'), ('OLR1', 'VCAM1'), ('F2', 'F2'), ('CAV1', 'CAV1'), ('ROCK2', 'IL8'), ('SRC', 'SRC'), ('CDKN1A', 'CDKN1A'), ('TRAF2', 'IKBKB'), ('TNF', 'MAPK1'), ('NOX4', 'NOX4'), ('STAT3', 'STAT3'), ('ESR1', 'ESR1'), ('OLR1', 'ROCK2'), ('JAK2', 'JAK2'), ('F2', 'JUN'), ('CREB1', 'HMOX1'), ('EGR1', 'F3'), ('CD40', 'MAPK1'), ('IL1B', 'ICAM1'), ('AGTR1', 'OLR1'), ('F2', 'PRKCZ'), ('SRC', 'JAK2'), ('CHUK', 'CHUK'), ('ERN1', 'ERN1'), ('TNFRSF1A', 'TNFRSF1A'), ('IFNG', 'JAK2'), ('ESR1', 'NOS3'), ('ATF4', 'IL8'), ('ATF4', 'ATF4'), ('TNF', 'TNFRSF1B'), ('TNF', 'TNFRSF1A'), ('ITGB3', 'ITGB3'), ('IKBKB', 'IKBKB'), ('EGR1', 'ICAM1'), ('SP1', 'NOS3'), ('ETS1', 'CCL2'), ('TNF', 'CAV1'), ('OLR1', 'CXCL2'), ('TNF', 'CHUK'), ('MAPK1', 'MAPK1'), ('CAV1', 'VCAM1'), ('STAT3', 'CCL2'), ('HMOX1', 'HMOX1'), ('IL8', 'IL8'), ('GATA3', 'VCAM1'), ('PPARG', 'SELE'), ('OLR1', 'PRKCB'), ('TNF', 'GATA3')},
"CV-IPN-Plaque_destabilization_1": {('MSR1', 'STAT1'), ('ATP2A2', 'ATP2A2'), ('IL17A', 'MMP1'), ('PPARA', 'TNF'), ('CD36', 'TLR4'), ('CD36', 'IL6'), ('CIITA', 'COL1A2'), ('TGFB1', 'COL3A1'), ('TLR2', 'IL1B'), ('TLR4', 'TLR4'), ('CD36', 'CD36'), ('IL18', 'MMP9'), ('CD36', 'TLR2'), ('MFGE8', 'IFNG'), ('TLR4', 'CXCL2'), ('IL1B', 'FGF2'), ('TIMP1', 'MMP9'), ('PPARA', 'IL6'), ('VEGFA', 'PLAU'), ('CD40', 'CEBPA'), ('TLR4', 'STAT1'), ('IL1B', 'CD40LG'), ('IL6', 'MMP9'), ('CD36', 'CXCL2'), ('PTGS2', 'MMP9'), ('MSR1', 'PLAU'), ('TLR4', 'TNF'), ('TNF', 'COL1A2'), ('IRF5', 'TNF'), ('MMP1', 'MMP1'), ('MMP9', 'MMP9'), ('HIF1A', 'VEGFA'), ('IRAK4', 'CXCL10'), ('F2', 'CXCL2'), ('CEBPA', 'MMP8'), ('IRAK4', 'NFKBIA'), ('ELK1', 'ELK1'), ('PPARA', 'F3'), ('PPARA', 'PPARA'), ('CD40', 'CD40'), ('ETS1', 'FASLG'), ('TIMP1', 'TIMP1'), ('PON2', 'PON2'), ('TLR2', 'CCL2'), ('CD36', 'CCL5'), ('CD40LG', 'SLC9A1'), ('PPARA', 'ABCA1'), ('CD40', 'VCAM1'), ('BAX', 'BAX'), ('VCAM1', 'VCAM1'), ('IRF5', 'IL6'), ('COL1A2', 'COL1A2'), ('CEBPA', 'CEBPA'), ('MMP13', 'MMP13'), ('ABCA1', 'ABCA1'), ('MSR1', 'MSR1'), ('CXCL10', 'CXCL10'), ('CD40', 'ICAM1'), ('PLAU', 'PLAU'), ('CIITA', 'COL1A1'), ('CCL2', 'CCL2'), ('IRAK4', 'TNF'), ('IL10', 'F3'), ('CCL2', 'PECAM1'), ('PON3', 'CCL2'), ('DDIT3', 'BAX'), ('TNF', 'MMP9'), ('CD40', 'SELE'), ('IFNG', 'ABCA1'), ('PON2', 'IL6'), ('TLR2', 'MYD88'), ('SLC9A1', 'IL12B'), ('MYD88', 'PTGS2'), ('TNF', 'MAP3K1'), ('TLR2', 'TNF'), ('F2', 'CCL2'), ('VCAM1', 'MMP2'), ('EGR1', 'EGR1'), ('STAT1', 'CIITA'), ('IL1B', 'MMP2'), ('CD40LG', 'MMP1'), ('TLR4', 'VCAM1'), ('IL17A', 'VCAM1'), ('PPARG', 'VCAM1'), ('CIITA', 'CIITA'), ('PECAM1', 'VCAM1'), ('IRAK4', 'IL1B'), ('TLR4', 'CCL2'), ('CD40LG', 'SELE'), ('IL1B', 'IL1B'), ('F2', 'PLAU'), ('CD36', 'MAPK9'), ('TNF', 'DDIT3'), ('CCL2', 'MMP2'), ('IFNG', 'CIITA'), ('IRAK4', 'IRAK4'), ('APOE', 'APOE'), ('PPARG', 'PPARG'), ('MAPK9', 'MAPK9'), ('IRAK4', 'CCL4'), ('CD40LG', 'PTGS2'), ('PPARA', 'CD36'), ('TLR4', 'MYD88'), ('IFNG', 'HLA-DRA'), ('TGFB1', 'COL1A2'), ('MMP8', 'MMP8'), ('F3', 'F2'), ('SELE', 'SELE'), ('FGF2', 'NOS2'), ('PPARG', 'CCL2'), ('TNF', 'VCAM1'), ('IFNG', 'COL1A2'), ('STAT1', 'STAT1'), ('CD36', 'TNF'), ('TIMP2', 'MMP2'), ('IL17A', 'MMP3'), ('TP53', 'BAX'), ('IFNG', 'COL1A1'), ('COL3A1', 'COL3A1'), ('NOS2', 'NOS2'), ('IFNG', 'IFNG'), ('F3', 'F3'), ('FGF2', 'MMP9'), ('CIITA', 'HLA-DRA'), ('TNF', 'MMP1'), ('CD36', 'CCL2'), ('ICAM1', 'ICAM1'), ('SOD1', 'IL6'), ('MAPK9', 'MSR1'), ('IKBKB', 'NFKBIA'), ('IL6', 'IL6'), ('MMP3', 'MMP3'), ('MAP3K1', 'MAP2K6'), ('CD40LG', 'CD40'), ('VEGFA', 'TIMP2'), ('DDIT3', 'DDIT3'), ('PPARA', 'ICAM1'), ('PPARA', 'CCL2'), ('MMP1', 'COL1A1'), ('PPARA', 'PTGS2'), ('FASLG', 'FASLG'), ('PPARA', 'ABCG1'), ('JAK2', 'STAT1'), ('MMP2', 'MMP2'), ('PECAM1', 'PECAM1'), ('SLC9A1', 'SLC9A1'), ('TGFB1', 'COL1A1'), ('PPARA', 'IL1B'), ('IL1B', 'MMP8'), ('TNF', 'TNF'), ('PPARA', 'VCAM1'), ('MMP1', 'COL1A2'), ('F2', 'F2'), ('TIMP3', 'TIMP3'), ('ABCG1', 'ABCG1'), ('JAK2', 'JAK2'), ('F2', 'ICAM1'), ('HLA-DRA', 'HLA-DRA'), ('IL10', 'ELK1'), ('MERTK', 'MERTK'), ('CD40', 'F3'), ('NOX4', 'NOX4'), ('CCL5', 'CCL5'), ('EGR1', 'F3'), ('SOD1', 'IL1B'), ('TGM2', 'ABCA1'), ('TLR4', 'IL1B'), ('IRAK4', 'IL6'), ('CASP1', 'CASP1'), ('CTGF', 'CTGF'), ('ERN1', 'ERN1'), ('MMP1', 'COL3A1'), ('IRAK4', 'CXCL1'), ('IRAK4', 'CCL2'), ('TIMP2', 'TIMP2'), ('IFNG', 'JAK2'), ('MYD88', 'IRF5'), ('HIF1A', 'HIF1A'), ('SOD1', 'SOD1'), ('PTGS2', 'PTGS2'), ('IL10', 'EGR1'), ('PON2', 'TNF'), ('ETS1', 'MMP1'), ('FCER1A', 'IL6'), ('TLR2', 'TLR2'), ('VEGFA', 'KDR'), ('CD40LG', 'F3'), ('SLC9A1', 'IL6'), ('IL1B', 'NFKBIA'), ('IRAK4', 'CXCL2'), ('IL1B', 'MMP9')},
"CV-IPN-Platelet_activation_1": {('F2', 'F2'), ('CXCL12', 'CXCL12'), ('P2RY1', 'P2RY1'), ('TBXA2R', 'TBXA2R'), ('TLR4', 'MYD88'), ('TLR6', 'TLR6'), ('MYD88', 'MYD88'), ('THBS1', 'CD36'), ('MAP2K4', 'MAPK9'), ('PTK2B', 'PIK3CB'), ('TLR1', 'TLR1'), ('GSN', 'GSN'), ('RAP1B', 'ITGA2B'), ('VWF', 'VWF'), ('CCL5', 'CCL5'), ('PPBP', 'PPBP'), ('P2RY12', 'P2RY12'), ('CD36', 'SYK'), ('TLR1', 'MYD88'), ('VWF', 'PRKG1'), ('ITGA2', 'PTK2B'), ('F2R', 'MAPK3'), ('CD36', 'CD36'), ('TLR4', 'TLR4'), ('PIK3CB', 'AKT1'), ('F2RL3', 'PRKCD'), ('F2RL3', 'AKT1'), ('LYN', 'LYN'), ('MAP2K4', 'MAP2K4'), ('F2R', 'F2R'), ('PRKCD', 'PRKCD'), ('PRKG1', 'PRKG1'), ('RAP1B', 'RAP1B'), ('PTGER3', 'PTGER3'), ('ITGA2B', 'ITGA2B'), ('TLR2', 'MYD88'), ('F2', 'F2RL3'), ('PRKG1', 'VASP'), ('SYK', 'SYK'), ('PRKG1', 'MAPK1'), ('PTEN', 'AKT1'), ('VEGFA', 'VEGFA'), ('F2R', 'MAPK1'), ('TLR2', 'TLR2'), ('CD40LG', 'CD40LG'), ('F2', 'F2R'), ('MAPK9', 'MAPK9'), ('PRKG1', 'RAP1B'), ('VWF', 'LYN'), ('CD36', 'VASP'), ('MAPK1', 'MAPK1'), ('ENTPD1', 'ENTPD1'), ('AKT1', 'AKT1'), ('SELP', 'SELP'), ('PF4', 'PF4'), ('F2RL3', 'F2RL3'), ('PTK2B', 'PTK2B'), ('VASP', 'VASP'), ('PIK3CB', 'PIK3CB'), ('INSR', 'AKT1'), ('IL1B', 'IL1B'), ('ITGA2B', 'ITGB3'), ('CD36', 'MAPK9'), ('FGA', 'FGA'), ('MAPK3', 'MAPK3'), ('AKT1', 'ITGA2B'), ('FGB', 'FGB')},
"CV-IPN-Smooth_muscle_cell_activation_1": {('CCND1', 'CCND1'), ('PDGFB', 'PDGFB'), ('LRP6', 'PDGFRB'), ('LRP1', 'PDGFRB'), ('PLAU', 'PTPN11'), ('CREB1', 'NR4A3'), ('TLR4', 'TLR4'), ('CREB1', 'CCND1'), ('RHOA', 'RHOA'), ('CREB1', 'PTGS2'), ('RB1', 'RB1'), ('RB1', 'E2F1'), ('CEBPD', 'CEBPD'), ('PPARA', 'CDKN2A'), ('IGF1', 'NOX4'), ('CD40LG', 'CASP1'), ('CDK4', 'E2F1'), ('PPARA', 'PPARA'), ('RAC1', 'RAC1'), ('LRP6', 'LRP6'), ('CD40', 'CD40'), ('NR4A3', 'NR4A3'), ('PDGFRA', 'PDGFRA'), ('PDGFB', 'PDGFRA'), ('LRP1', 'IRAK1'), ('MAPK3', 'MAPK3'), ('CD40', 'RELA'), ('HIF1A', 'IL8'), ('CDK4', 'CDK4'), ('HIF1A', 'NOS2'), ('AGTR1', 'AGTR1'), ('LRP1', 'RHOA'), ('PDGFRA', 'MAPK1'), ('CEBPD', 'IL6'), ('HIF1A', 'PTGS2'), ('DCBLD2', 'PDGFRB'), ('CTNNB1', 'CTNNB1'), ('IL1B', 'IRAK1'), ('HIF1A', 'IL1B'), ('LTBR', 'RELA'), ('PDGFRB', 'MAPK3'), ('MAPK3', 'CREB1'), ('E2F1', 'CCND1'), ('CCNE1', 'CCNE1'), ('MAPK1', 'CREB1'), ('IL1B', 'IL1B'), ('CTNNB1', 'CCND1'), ('AKR1B1', 'RB1'), ('GSK3B', 'GSK3B'), ('IGF1', 'RAC1'), ('PTPN11', 'PTPN11'), ('TLR4', 'MYD88'), ('CREB1', 'CREB1'), ('IRAK1', 'IRAK1'), ('MYD88', 'MYD88'), ('PDGFB', 'RAC1'), ('AKT1', 'GSK3B'), ('E2F1', 'MYC'), ('PDGFRB', 'PDGFRB'), ('PDGFA', 'PDGFRA'), ('PRKCD', 'HIF1A'), ('NOS2', 'NOS2'), ('F3', 'F3'), ('AGTR1', 'PTPN11'), ('PRKCD', 'CCND1'), ('AKR1B1', 'CDK4'), ('AKT1', 'AKT1'), ('IL6', 'IL6'), ('CD40LG', 'CD40'), ('PDGFRB', 'MAPK1'), ('TNF', 'TNF'), ('MMP2', 'MMP2'), ('LRP1', 'LRP1'), ('E2F1', 'PCNA'), ('TNF', 'AKR1B1'), ('PLAUR', 'PLAUR'), ('PCNA', 'PCNA'), ('E2F1', 'E2F1'), ('APOE', 'LRP1'), ('E2F1', 'CCNE1'), ('NOX4', 'MMP2'), ('CCR3', 'MMP2'), ('CASP1', 'CASP1'), ('CCL11', 'CCR3'), ('THBS1', 'MMP2'), ('HIF1A', 'HIF1A'), ('PTGS2', 'PTGS2'), ('GSK3B', 'CTNNB1'), ('PDGFB', 'PDGFRB'), ('PLAU', 'PLAUR'), ('IL1B', 'CEBPD'), ('RHOA', 'F3'), ('MYC', 'MYC'), ('MAPK1', 'MAPK1'), ('AKR1B1', 'AKR1B1'), ('PDGFRA', 'MAPK3'), ('IL8', 'IL8'), ('NR4A3', 'CCND1')},
"CV-IPN-Foam_cell_formation_1": {('PTK2B', 'MAPK3'), ('MSR1', 'STAT1'), ('STAT1', 'ACAT1'), ('BIRC5', 'BIRC5'), ('TNFSF15', 'CD36'), ('EDN1', 'EDNRB'), ('ADAM17', 'OLR1'), ('PRKCZ', 'VEGFA'), ('APOA1', 'APOA1'), ('SYK', 'MAPK1'), ('CD36', 'TLR4'), ('CD36', 'PTK2B'), ('TLR4', 'SYK'), ('TLR4', 'TLR4'), ('CD36', 'CD36'), ('SYK', 'MAPK3'), ('CD36', 'TLR2'), ('MAPK1', 'SREBF2'), ('FABP4', 'ABCA1'), ('PPARD', 'PPARD'), ('PTK2B', 'MAPK1'), ('VEGFA', 'VEGFA'), ('TRAF2', 'TRAF2'), ('AGTR1', 'ABCA1'), ('TLR9', 'RELA'), ('PPARG', 'ABCA1'), ('PPARA', 'CPT1B'), ('HMGB1', 'TLR4'), ('PPARA', 'PPARA'), ('RAC1', 'RAC1'), ('SREBF1', 'SREBF1'), ('CPT1A', 'CPT1A'), ('OLR1', 'OLR1'), ('PPARA', 'ABCA1'), ('SCAP', 'SREBF2'), ('MAPK3', 'MAPK3'), ('TNFSF15', 'ABCG1'), ('TLR9', 'TLR9'), ('TNFRSF1A', 'TNFRSF1A'), ('TLR9', 'MYD88'), ('CD14', 'TLR4'), ('CD14', 'CD14'), ('NCOR1', 'PPARG'), ('AGTR1', 'AGTR1'), ('FABP4', 'PPARG'), ('ABCA1', 'ABCA1'), ('MSR1', 'MSR1'), ('MAPK3', 'RAC1'), ('MAPK1', 'CDC42'), ('XDH', 'XDH'), ('IFNG', 'ABCA1'), ('EDNRB', 'EDNRB'), ('IL6', 'STAT3'), ('NCOR1', 'PPARD'), ('SPI1', 'MSR1'), ('SREBF1', 'ACAT1'), ('OLR1', 'MSR1'), ('STAT3', 'ABCA1'), ('SREBF2', 'LDLR'), ('PPARD', 'MSR1'), ('CD36', 'MAPK9'), ('MAPK9', 'MAPK9'), ('SIRT1', 'SIRT1'), ('MYD88', 'MYD88'), ('TNFSF15', 'MSR1'), ('PPARA', 'CD36'), ('SREBF2', 'ACAT1'), ('TLR4', 'IL10'), ('NCOR1', 'NCOR1'), ('IFNA1', 'MSR1'), ('TNFSF14', 'TNFSF14'), ('STAT1', 'STAT1'), ('NR1H3', 'NR1H3'), ('HMGB1', 'TLR2'), ('MYD88', 'MAPK1'), ('TNFSF15', 'SCARB1'), ('MYD88', 'MAPK3'), ('EDNRA', 'EDNRA'), ('SREBF2', 'SCARB1'), ('SREBF2', 'SREBF2'), ('NOX1', 'NOX1'), ('PPARD', 'ABCA1'), ('MAPK3', 'EGR1'), ('SCARB1', 'SCARB1'), ('NR1H3', 'ABCA1'), ('PPARA', 'CPT1A'), ('MAP3K7', 'CHUK'), ('APOA1', 'LCAT'), ('IKBKB', 'NFKBIA'), ('SIRT1', 'NR1H3'), ('CYP27A1', 'CYP27A1'), ('CDC42', 'CDC42'), ('MAP3K7', 'MAP3K7'), ('TNFRSF1A', 'TRADD'), ('MAPK1', 'EGR1'), ('PTK2B', 'PTK2B'), ('AKT1', 'AKT1'), ('CCR7', 'CCR7'), ('FABP4', 'FABP4'), ('TNFSF15', 'ABCA1'), ('PTAFR', 'PTAFR'), ('TLR4', 'IRF3'), ('PPARA', 'ABCG1'), ('LCAT', 'LCAT'), ('TRAF2', 'MAP3K7'), ('SCAP', 'SCAP'), ('ACAT1', 'ACAT1'), ('LDLR', 'LDLR'), ('MAPK3', 'CDC42'), ('TGFB1', 'NR1H3'), ('NR1H2', 'NR1H2'), ('ABCG1', 'ABCG1'), ('CPT1B', 'CPT1B'), ('SPI1', 'SPI1'), ('CHUK', 'NFKBIA'), ('STAT3', 'STAT3'), ('PPARG', 'APOA1'), ('NR1H2', 'CCR7'), ('NR3C1', 'SCAP'), ('SREBF1', 'LDLR'), ('IRF3', 'NR1H3'), ('CHUK', 'CHUK'), ('AKT1', 'MSR1'), ('TRADD', 'TRADD'), ('TNF', 'TNFRSF1A'), ('STAT3', 'XDH'), ('SOD1', 'SOD1'), ('IFNG', 'CYP27A1'), ('EDN1', 'EDNRA'), ('IKBKB', 'IKBKB'), ('TLR2', 'TLR2'), ('IFNG', 'SPI1'), ('TRADD', 'TRAF2'), ('TLR2', 'FABP4'), ('NR1H2', 'ABCA1'), ('TLR9', 'CHUK'), ('PPARG', 'XDH'), ('NR3C1', 'NR3C1'), ('MAPK1', 'MAPK1'), ('NR1H3', 'CCR7'), ('PPARG', 'CD36'), ('INSR', 'AKT1'), ('MAP3K7', 'IKBKB')}

}
genes2name = {'CXCL12': {'CXCL12'}, 'CXCL13': {'CXCL13'}, 'CXCL14': {'CXCL14'}, 'PPBP': {'PPBP'}, 'XCL1': {'XCL1'},
              'CCL2': {'CCL2'}, 'PF4': {'PF4'}, 'CXCL10': {'CXCL10'}, 'CXCL5': {'CXCL5'},
              'XCL2': {'XCL1', 'XCL2'}, 'CXCL1': {'CXCL1'}, 'CCL13': {'CCL13'}, 'CXCL11': {'CXCL11'},
              'CXCL8': {'CXCL8'}, 'CCL14': {'CCL14'}, 'CCL3': {'CCL3', 'CCL3L3'}, 'CCL5': {'CCL5'},
              'CXCL2': {'CXCL2'}, 'CCL15': {'CCL15'}, 'CXCL3': {'CXCL3'},
              'CCL21': {'CCL21C', 'CCL21A', 'CCL21', 'GM10591', 'GM13304', 'GM21541', 'CCL21B'}, 'CCL17': {'CCL17'},
              'CXCL6': {'CXCL6'}, 'CCL11': {'CCL11'}, 'CCL7': {'CCL7'}, 'CCL4': {'CCL4'}, 'CCL1': {'CCL1'},
              'CXCL16': {'CXCL16'}, 'CCL18': {'CCL18'}, 'CCL19': {'GM2564', 'CCL19'}, 'CXCL9': {'CXCL9'},
              'CCL8': {'CCL8', 'CCL12'}, 'CCL20': {'CCL20'}, 'C5': {'C5', 'HC'}, 'CCL22': {'CCL22'}, 'CCL24': {'CCL24'},
              'CX3CL1': {'CX3CL1'}, 'CCL25': {'CCL25'}, 'CCL28': {'CCL28'}, 'CCL23': {'CCL23'},
              'CCL26': {'CCL26'}, 'CCL27': {'CCL27A', 'CCL27', 'CCL27B', 'GM13306'}, 'CCL6': {'CCL6'}, 'CCL9': {'CCL9'},
              }

normGeneSymbols = normalize_gene_names(path="/mnt/d/owncloud/data/miRExplore/obodir/" + "/hgnc_no_withdrawn.syn")


cbn2restrict = OrderedDict([(

    "CV-IPN-Endothelial_cell_activation_1", {
        "sentences": "false",
        "cells": [
                {"group": "cells", "name": "endothelial cell",  "termid": "META:52"},
                {"group": "cells", "name": "HUVEC-C",  "termid": "META:02689"}
              ],
        "ncits": [
            {"group": "ncits", "name": "Vascular Cell Adhesion Protein 1", "termid": "NCIT:C48214"},
            {"group": "ncits", "name": "Intercellular Adhesion Molecule 1", "termid": "NCIT:C17304"}
        ]
        }),
    ("CV-IPN-Endothelial_cell-monocyte_interaction_1", {
        "sentences": "false",
        "cells": [{"group": "cells", "name": "foam cell", "termid": "CL:0000891"},
            {"group": "cells", "name": "monocyte", "termid": "META:148"},
            {"group": "cells","name": "endothelial cell","termid": "META:52" },
            {"group": "cells", "name": "macrophage", "termid": "META:99"},
            {"group": "cells", "name": "T cell", "termid": "META:44"}
            ],
        "go":[
            { "group": "go", "name": "inflammatory response", "termid": "GO:0006954"  },
            { "group": "go", "name": "cell proliferation", "termid": "GO:0008283"  }
        ]
    }),
    ("CV-IPN-Foam_cell_formation_1", {
        "sentences": "false",
        "cells": [
                  {"group": "cells", "name": "foam cell", "termid": "CL:0000891"},
                  {"group": "cells", "name": "smooth muscle cell", "termid": "META:83"},
                  {"group": "cells", "name": "macrophage", "termid": "META:99"}
                  ],
        "go": [{ "group": "go", "name": "apoptotic process", "termid": "GO:0006915" },
               {"group": "go", "name": "cell death", "termid": "GO:0008219"}]
    }),
    ("CV-IPN-Smooth_muscle_cell_activation_1", {
        "sentences": "false",
        "cells": [
            {"group": "cells", "name": "endothelial cell", "termid": "META:52"},
            {"group": "cells", "name": "smooth muscle cell", "termid": "META:83"},
        ]
    }),
    ("CV-IPN-Platelet_activation_1", {
            'sentences': "false",
            "cells": [
                {"group": "cells", "name": "endothelial cell", "termid": "META:52"},
                {"group": "cells", "name": "smooth muscle cell", "termid": "META:83"},
                { "group": "cells", "name": "blood cell", "termid": "META:41" },
                { "group": "cells", "name": "thromboblast", "termid": "CL:0000828" }
            ],
            "ncits": [
                { "group": "ncits", "name": "Low-Density Lipoprotein Recepter", "termid": "NCIT:C17074" },
            ]
        }),
    ("CV-IPN-Plaque_destabilization_1", {
        'sentences': "false",
        "cells": [
            {"group": "cells", "name": "endothelial cell", "termid": "META:52"},
            {"group": "cells", "name": "smooth muscle cell", "termid": "META:83"},
            {"group": "cells", "name": "blood cell", "termid": "META:41"},
            {"group": "cells", "name": "thromboblast", "termid": "CL:0000828"},
            {"group": "cells", "name": "macrophage", "termid": "META:99"},
            {"group": "cells", "name": "neutrophil", "termid": "CL:0000775"},
            {"group": "cells", "name": "leukocyte", "termid": "META:167"},
        ],
        "go": [
            {"group": "go", "name": "apoptotic process", "termid": "GO:0006915"},
            {"group": "go", "name": "autophagy", "termid": "GO:0006914"},
            {"group": "go", "name": "response to endoplasmic reticulum stress", "termid": "GO:0034976"},
        ]
    })
    ])

network2nicename = {
"CV-IPN-Plaque_destabilization_1": "Plaque destabilization",
"CV-IPN-Platelet_activation_1": "Platelet activation",
"CV-IPN-Smooth_muscle_cell_activation_1": "SMC activation",
"CV-IPN-Foam_cell_formation_1": "Foam cell formation",
"CV-IPN-Endothelial_cell-monocyte_interaction_1": "EC/MC interaction",
"CV-IPN-Endothelial_cell_activation_1": "EC activation",
}


def acceptEvidence(ev):

    return True

    if ev['data_source'] == 'miRTarBase':

        if "Weak" in ev['functional_type']:
            return False

    return True




cbn2graph = {}



totalNWGenes = set()

for cbnNW in cbn2edges:

    allNWGenes = set()
    newEdges = set()

    for edge in cbn2edges[cbnNW]:

        src = edge[0].upper()
        tgt = edge[1].upper()

        if src in normGeneSymbols:
            src = normGeneSymbols[src]

        if tgt in normGeneSymbols:
            tgt = normGeneSymbols[tgt]

        allNWGenes.add(src)
        allNWGenes.add(tgt)

        totalNWGenes.add(src)
        totalNWGenes.add(tgt)

        newEdges.add((src, tgt))

mirnaNodes = set()

stageMir2Cells = defaultdict(lambda: defaultdict(set))
stageMir2CellEvidence = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))

for cbnNW in cbn2edges:

    allNWGenes = set()
    newEdges = set()

    for edge in cbn2edges[cbnNW]:

        src = edge[0].upper()
        tgt = edge[1].upper()

        if src in normGeneSymbols:
            src = normGeneSymbols[src]

        if tgt in normGeneSymbols:
            tgt = normGeneSymbols[tgt]

    cbn2edges[cbnNW] = newEdges

    #print(cbnNW, allNWGenes)

    requestData = cbn2restrict[cbnNW]
    requestData['gene'] = list(totalNWGenes)

    requestData['disease'] = [
        {'group': 'disease', 'termid': 'DOID:1287', 'name': 'cardiovascular system disease'},
        {'group': 'disease', 'termid': 'DOID:2349', 'name': 'arteriosclerosis'}
    ]

    print(cbnNW)

    graph, nodeCounts, evcount, json = DataBasePlotter.fetchGenes(requestData, gene2name=genes2name, minPMIDEvCount=0, minTgtCount=0, acceptEv=acceptEvidence, MIRNASTRPARTS=[miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR])

    newNodes = []

    for (node, nodeAttr) in graph.nodes(data=True):

        if node in nodeCounts:
            nodeAttr['size'] = 20 + nodeCounts[node]

        if nodeAttr['color'] == 'blue':
            mirnaNodes.add(node)

        newNodes.append((node, nodeAttr))

    removeNodes = [node for node, degree in graph.degree() if degree >= 2 and node in mirnaNodes]
    #graph.remove_nodes_from(removeNodes)


    removeEdges = set()

    singleEvidenceCounter = Counter()
    edgeEvidenceCounter = Counter()

    for (edgeSrc, edgeTgt, edgeAttr) in graph.edges(data=True):

        edge = (edgeSrc, edgeTgt)
        edgeR = (edgeTgt, edgeSrc)

        evC = None

        if edge in evcount:
            evC = evcount[edge]
        elif edgeR in evcount:
            evC = evcount[edgeR]

        if evC != None:

            if len(evC) == 1:

                for x in evC:
                    singleEvidenceCounter[x] += 1

                edgeAttr['color'] = '#00ff00'

            elabel = ""
            for x in sorted([x for x in evC]):
                elabel += x[0]

            edgeAttr["label"] = elabel

            edgeEvidenceCounter[elabel] += 1

            miRName = None

            if edge[0].startswith("miR") or edge[0].startswith("let-"):
                miRName = edge[0]
            elif edge[1].startswith("miR") or edge[1].startswith("let-"):
                miRName = edge[1]


            if "celldata" in edgeAttr and miRName != None :

                for x in edgeAttr["celldata"]:
                    stageMir2Cells[cbnNW][miRName].add(x)

                for x in edgeAttr['cellEvidence']:
                    for pmid in edgeAttr['cellEvidence'][x]:
                        stageMir2CellEvidence[cbnNW][miRName][x].add(pmid)

            elif miRName == None:
                print("INVALID MIRNA:", edge, file=sys.stderr)

        else:
            print(edge, "not in evcount")

    print(singleEvidenceCounter)
    print(edgeEvidenceCounter)

    for edge in removeEdges:
        graph.remove_edge(edge[0], edge[1])

    graph.remove_nodes_from(networkx.isolates(graph))

    for edge in cbn2edges[cbnNW]:
        if edge[0] == edge[1]:
            continue

        graph.add_edge(edge[0], edge[1], {'color': '#aa0000'})

    cbn2graph[cbnNW] = graph


### which cell types are most common per stage?
print()
print("cells per stage")
print()
print()
print()

importantCellTypesPerMirna = defaultdict(lambda: Counter())

for stage in stageMir2Cells:

    cellCounter = Counter()

    stageMirnaCellPairs = Counter()
    mirnaCellPairs = defaultdict(lambda: Counter())
    stageCellCount = Counter()

    for mirna in stageMir2Cells[stage]:


        for cell in stageMir2Cells[stage][mirna]:
            cellCounter[cell] += 1

        mirnaStageCells = [x for x in stageMir2Cells[stage][mirna]]

        for i in range(0, len(mirnaStageCells)):
            for j in range(i+1, len(mirnaStageCells)):

                if mirnaStageCells[i][0].startswith("CVCL") or mirnaStageCells[j][0].startswith("CVCL"):
                    continue

                cellpair = tuple(sorted((mirnaStageCells[i],mirnaStageCells[j])))

                mirnaCellPairs[mirna][cellpair] += 1
                stageMirnaCellPairs[cellpair] += 1

                stageCellCount[cellpair[0]] += 1
                stageCellCount[cellpair[1]] += 1

                importantCellTypesPerMirna[mirna][cellpair] += 1


    mostCommonCellPairs = []
    for (mpair, count) in stageMirnaCellPairs.most_common(): #20
        mostCommonCellPairs.append(mpair)


    edge2support = defaultdict(set)

    for mirna in mirnaCellPairs:
        for cellpair in mirnaCellPairs[mirna]:
            if cellpair in mostCommonCellPairs:

                edge2support[cellpair].add(mirna)

                print(stage, mirna, cellpair[0], cellpair[1], mirnaCellPairs[mirna][cellpair], stageMirnaCellPairs[cellpair], stageMir2CellEvidence[stage][mirna].get(cellpair[0]),stageMir2CellEvidence[stage][mirna].get(cellpair[1]) )

    cellgraph = networkx.Graph()

    allnodes = set()
    for edge in edge2support:
        allnodes.add(edge[0])
        allnodes.add(edge[1])

    for node in allnodes:
        cellgraph.add_node(node[1] + " ("+node[0]+")", size=20 + stageCellCount[node])


    cellCommunicatorDF = DataFrame()
    cellCommunicatorDF.addColumns(["miRNA", "cells"])

    mirna2cells = defaultdict(set)

    for edge in edge2support:
        cellgraph.add_edge(
            edge[0][1] + " (" + edge[0][0] + ")",
            edge[1][1] + " (" + edge[1][0] + ")",
            label=", ".join(edge2support.get(edge, [])))

        mirnas = edge2support.get(edge, [])

        for mirna in mirnas:
            mirna2cells[mirna].add(edge[0][1] + " (" + edge[0][0] + ")")
            mirna2cells[mirna].add(edge[1][1] + " (" + edge[1][0] + ")")


    cells2mirnas = defaultdict(set)
    for mirna in mirna2cells:
        cells = tuple(sorted(mirna2cells[mirna]))

        cells2mirnas[cells].add(mirna)


    for cells in cells2mirnas:

        ncells = []
        for cell in cells:
            ncells.append( ""+cell+"}" )

        rowdict = {
            'miRNA': "\\makecell[l]{" + "\\\\".join(sorted(cells2mirnas[cells])) + "}",
            'cells': "\\makecell[l]{" + "\\\\".join(cells) + "}"
        }
        cellCommunicatorDF.addRow(DataRow.fromDict(rowdict))


    cellCommunicatorDF.export("/mnt/d/yanc_network/stage_cells_cliques"+stage+".latex", ExportTYPE.LATEX)
    cellCommunicatorDF.export("/mnt/d/yanc_network/stage_cells_cliques"+stage+".tsv", ExportTYPE.TSV)


    CytoscapeGrapher.showGraph(cellgraph, '/mnt/d/yanc_network/', name="stage_cells_" + stage)
    figidx = CytoscapeGrapher.plotNXGraph(cellgraph, stage, ["/mnt/d/yanc_network/stage_cells_" + stage.replace(" ", "_") + ".png", "/mnt/d/yanc_network/stage_cells_" + stage.replace(" ", "_") + ".pdf"], figidx)

    print()
    print()
    print()
    print()
    print()

    plotKeys = []
    plotValues = []
    for (cell, count) in cellCounter.most_common(20):
        print(stage, cell[0], cell[1], count)
        plotKeys.append(cell[1] + " ("+cell[0]+")")
        plotValues.append(count)


    plt.figure(figidx)
    figidx+= 1
    plt.bar(plotKeys, plotValues)
    plt.xticks(rotation=90)
    plt.tight_layout()

    plt.savefig("/mnt/d/yanc_network/cells_per_stage_" + stage + ".png",bbox_inches='tight')
    plt.savefig("/mnt/d/yanc_network/cells_per_stage_" + stage + ".pdf",bbox_inches='tight')

    plt.show()
    print()
    print()
    print()
    print()
    print()

for mirna in importantCellTypesPerMirna:

    for (cp, count) in importantCellTypesPerMirna[mirna].most_common(10):
        print(mirna, cp[0], cp[1], count)

print()
print()
print()
print()
print()

### which miRNAs play a role in all networks?

### are there miRNAs which orchestrate the transition from one state to the other? which ones are unique to each stage? only for miRNAs connected to genes present in both stages

makeStory = [
    ['CV-IPN-Endothelial_cell_activation_1','CV-IPN-Endothelial_cell-monocyte_interaction_1','CV-IPN-Foam_cell_formation_1','CV-IPN-Smooth_muscle_cell_activation_1','CV-IPN-Platelet_activation_1',
     'CV-IPN-Plaque_destabilization_1'
     ]
]

storyTransitions = []


for storyline in makeStory:
    for i in range(0, len(storyline)-1):
        storyTransitions.append( (storyline[i], storyline[i+1]) )


def getMirsForGene(gene, graph):

    alledges = graph.edges(gene)

    targetMirs = set()

    for u,v in alledges:

        if u == gene:
            if v.startswith("miR"):
               targetMirs.add(v)

        elif v == gene:
            if u.startswith("miR"):
                targetMirs.add(u)

    return targetMirs

outfile = open("/mnt/d/yanc_network/" + "mirs_per_stage.tsv", 'w')

print("Condition1", "Condition2", "Gene", "miRNAs Only Before", "miRNAs Only After", "Common Mirs", sep="\t")
print("Condition1", "Condition2", "Gene", "miRNAs Only Before", "miRNAs Only After", "Common Mirs", sep="\t", file=outfile)
for storyTransition in storyTransitions:

    graphFrom = cbn2graph[storyTransition[0]]
    graphTo = cbn2graph[storyTransition[1]]

    fromGenes = set()
    toGenes = set()

    for (node, nodeAttr) in graphFrom.nodes(data=True):
        if not node.startswith("miR"):
            fromGenes.add(node)

    for (node, nodeAttr) in graphTo.nodes(data=True):
        if not node.startswith("miR"):
            toGenes.add(node)

    commonGenes = fromGenes.intersection(toGenes)

    genesWithDiff = {}

    for cgene in commonGenes:

        fromMirs = getMirsForGene(cgene, graphFrom)
        toMirs = getMirsForGene(cgene, graphTo)

        diffMirs = set()
        toDiffMirs = set()
        fromDiffMirs = set()
        commonMirs = set()
        for x in fromMirs.union(toMirs):
            if not x in fromMirs and x in toMirs:
                toDiffMirs.add(x)
            elif x in fromMirs and not x in toMirs:
                fromDiffMirs.add(x)
            elif x in fromMirs and x in toMirs:
                commonMirs.add(x)


        if len(toDiffMirs) > 0 or len(fromDiffMirs) > 0:
            genesWithDiff[cgene] = (fromDiffMirs, toDiffMirs, commonMirs)


    for cgene in genesWithDiff:
        mirset = genesWithDiff[cgene]
        print(storyTransition[0], storyTransition[1], cgene, ",".join(sorted(mirset[0])), ",".join(sorted(mirset[1])),",".join(sorted(mirset[2])), sep="\t")
        print(storyTransition[0], storyTransition[1], cgene, ",".join(sorted(mirset[0])), ",".join(sorted(mirset[1])),",".join(sorted(mirset[2])), sep="\t", file=outfile)

outfile.close()


latexStageMirsFile = open("/mnt/d/yanc_network/cbn_stages_most_reg.txt", 'w')

mirSByCBN = defaultdict(set)

mir2networks = defaultdict(set)
interactionsPerNetworkAndMirna = defaultdict(lambda: dict())

overallMostRegulatingMIRs = Counter()

for cbnNW in cbn2restrict:

    graph = cbn2graph[cbnNW]

    for (node, nodeAttr) in graph.nodes(data=True):

        if 'color' in nodeAttr and nodeAttr['color'] == 'blue':

            otherGraphHasNode = False

            for oNW in cbn2graph:
                if oNW == cbnNW:
                    continue

                ograph = cbn2graph[oNW]

                if ograph.has_node(node):
                    otherGraphHasNode = True
                    break

            if not otherGraphHasNode:
                nodeAttr['color'] = '#d942f4'


    intNodes = [(node,degree) for node, degree in graph.degree() if node in mirnaNodes]
    intNodes = sorted(intNodes, key=lambda x: x[1], reverse=True)


    for node in intNodes:
        if node[0].startswith("miR"):
            mirSByCBN[cbnNW].add(node[0])
            mir2networks[node[0]].add(cbnNW)

            overallMostRegulatingMIRs[node[0]] += node[1]

            interactionsPerNetworkAndMirna[cbnNW][node[0]] = node[1]


    print(cbnNW)
    cbnNWMirs = []
    for i in range(0, min(10, len(intNodes))):
        print(i, intNodes[i])
        cbnNWMirs.append(  intNodes[i][0] + " ("+str(intNodes[i][1])+")"  )

    latexStageMirsFile.write(cbnNW + "&" + ", ".join(cbnNWMirs) + "\n")

    print()



    print()

    mygraph = CytoscapeGrapher.showGraph(graph, location='/mnt/d/yanc_network/', name=cbnNW)

latexStageMirsFile.close()

print("Overall most regulating")
print(overallMostRegulatingMIRs.most_common(10))
print(", ".join([x[0] + " (" + str(x[1]) + ")" for x in overallMostRegulatingMIRs.most_common(10)]))
print()
print()


mirna2nwcount = Counter()
for mirna in mir2networks:
    mirna2nwcount[mirna] = len(mir2networks[mirna])

print("most common mirnas")
mcMirnas = set()
for x in mirna2nwcount.most_common(10):
    print(x)
    mcMirnas.add(x[0])

print(mcMirnas)

mirna2stagecount = defaultdict(list)
stages = [x for x in cbn2restrict]
for mirna in natsorted(mcMirnas):

    mirnarow = [mirna]
    for stage in stages:
        icount = interactionsPerNetworkAndMirna[stage].get(mirna, 0)
        mirnarow.append(icount)

        print(stage, mirna, icount, sep="\t")


    mirna2stagecount[mirna] = mirnarow




stageDF = pd.DataFrame.from_dict(mirna2stagecount, orient='index', columns=["miRNA"] + [network2nicename[x] for x in stages])

ax = parallel_coordinates(stageDF, 'miRNA', colormap=plt.get_cmap("Set2"))
plt.xticks(rotation=90)

lgd = plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=5)

plt.tight_layout()

plt.savefig("/mnt/d/yanc_network/mirna_stage_parallel.png",bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.savefig("/mnt/d/yanc_network/mirna_stage_parallel.pdf",bbox_extra_artists=(lgd,), bbox_inches='tight')

plt.show()


print()
print()


allSet = None

for nw in mirSByCBN:

    if allSet == None:
        allSet = mirSByCBN[nw]
    else:
        allSet = allSet.intersection(mirSByCBN[nw])

print("miRNAs in all nws")
print(allSet)


import matplotlib.pyplot as plt

for stages in makeStory:

    mergedGraph = cbn2graph[stages[0]]

    for i in range(1, len(stages)):
        mergedGraph = networkx.compose(mergedGraph, cbn2graph[stages[i]])


    pos = networkx.spring_layout(mergedGraph)

    for stage in stages:
        plt.figure(figidx, figsize=(50,30))
        figidx += 1

        networkGraph = cbn2graph[stage]

        edges = networkGraph.edges()

        colors = []
        for u,v in edges:
            elem = networkGraph[u][v]

            if not 'color' in elem:
                elem['color'] = "#0000FF"

            colors.append(elem['color'])


        nodes = networkGraph.nodes()
        nodeColors = []

        mirNodes = 0
        geneNodes = 0

        for x in nodes:
            if any([x.lower().startswith(y) for y in ['mir', 'let']]):
                nodeColors.append('blue')
                mirNodes += 1
            else:
                nodeColors.append('green')
                geneNodes += 1


        lineWidth = 0.5

        if len(nodes) < 200 and len(edges) < 200:
            print("change linewidth", stage)
            lineWidth = 3

        else:
            print("linewidth unchanged", stage)
            lineWidth=0.75

        networkx.draw(networkGraph, pos, font_size=25, with_labels=False, node_color=nodeColors, edges=edges, edge_color=colors, node_size=1200, width=lineWidth, font_weight='bold', dpi=1000)
        for p in pos:  # raise text positions
            clist = list(pos[p])
            clist[1] = clist[1] + 0.02
            pos[p] = tuple(clist)

        networkx.draw_networkx_labels(networkGraph, pos, font_weight='bold', font_size=25)

        plt.suptitle("{title}: {mirc} miRs, {genec} genes, {intc} interactions".format(title=stage, mirc=mirNodes, genec=geneNodes, intc=len(colors)))

        plt.savefig("/mnt/d/yanc_network/" + stage.replace(" ", "_") + ".png")
        plt.savefig("/mnt/d/yanc_network/" + stage.replace(" ", "_") + ".pdf")



