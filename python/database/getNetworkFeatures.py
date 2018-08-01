import json
import sys
import platform
import os
import datetime

class ExpressionNetwork:

    MIRNA = 'mirna'
    PROTEIN_CODING = 'protein_coding'
    MM10 = 'Network.MM10'
    HG38 = 'Network.HG38'
    LOG2FC = 'log2fc'
    RIGHT = 'right'
    LEFT = 'left'
    RIGHT_gene = 'right_gene'
    LEFT_gene = 'left_gene'
    INF = 1000000

    ### files and folders

    BIOCLIENT = 'Linux-4.4.140-94.42-default-x86_64-with-SuSE-12-x86_64'
    MARKUS = 'Linux-4.4.0-17713-Microsoft-x86_64-with-Ubuntu-18.04-bionic'
    dict_species = {
        MM10: ['Mouse', 'Mus musculus', MM10],
        HG38: ['Human', 'Homo sapiens', HG38]
    }

    def __init__(self):
        if sys.platform == 'win32':
            self.ABS_ROOT = 'E:\\masterpraktikum\\DiffExp\\integration\\'
        elif platform.platform() == ExpressionNetwork.BIOCLIENT:
            self.ABS_ROOT = '/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/data/integration/'
        elif platform.platform() == ExpressionNetwork.MARKUS:
            self.ABS_ROOT = '/home/mjoppich/ownCloud/data/miRExplore/obodir/expression/'
        else:
            pass

        self.dict_expression = {
            ExpressionNetwork.MM10: {
                'mirTrap mir103': self.ABS_ROOT + 'expression/M103/mirtrap_mouse_103.txt_confident.json',
                'mirTrap let7': self.ABS_ROOT + 'expression/Mlet7/mirtrap_mouse_let7.txt_confident.json'
            },
            ExpressionNetwork.HG38: {
                'mirTrap mir103': self.ABS_ROOT + 'expression/H103/mirtrap_human_103.txt_confident.json',
                'mirTrap let7': self.ABS_ROOT + 'expression/Hlet7/mirtrap_human_let7.txt_confident.json',
                'SRR20546 control/IL1a': self.ABS_ROOT + 'expression/SRR20546_control_IL1a/ena_SRR20546_control_IL1a.txt_confident.json',
                'SRR20546 control/both': self.ABS_ROOT + 'expression/SRR20546_control_both/ena_SRR20546_control_both.txt_confident.json',
                'SRR20546 control/PDGF1': self.ABS_ROOT + 'expression/SRR20546_control_PDGF1/ena_SRR20546_control_PDGF1.txt_confident.json'
                }
        }

        self.dict_coexpression = {
            ExpressionNetwork.MM10: {},
            ExpressionNetwork.HG38: {
                'coexpression SRR20546': self.ABS_ROOT + 'coexpression/ena.txt.json_protein_coding_coexpression_scoring.txt'}
            }

        self.dict_lncDetails = {
            ExpressionNetwork.MM10: self.ABS_ROOT + 'neighbors/neighbor_protein_coding_mm10_index.json',
            ExpressionNetwork.HG38: self.ABS_ROOT + 'neighbors/neighbor_protein_coding_hg38_index.json'
        }

        self.dict_anno = {
            ExpressionNetwork.MM10: {
                ExpressionNetwork.PROTEIN_CODING:self.ABS_ROOT + 'anno/mm10_protein_coding.json_new',
                ExpressionNetwork.MIRNA:self.ABS_ROOT + 'anno/mm10_miRNA.json_new'
            },
            ExpressionNetwork.HG38: {
                ExpressionNetwork.PROTEIN_CODING:self.ABS_ROOT + 'anno/hg38_protein_coding.json_new',
                ExpressionNetwork.MIRNA:self.ABS_ROOT + 'anno/hg38_miRNA.json_new'
            }
        }
    ########################################## Aux methods #######################################################

    # private method
    def __loadJson(self,file):
        try:
            dict = {}
            with open(file, 'r') as f:
                dict = json.load(f)
            return dict
        except:
            print("file not found! " + file)

    # private method
    def __getSpecies(self,name):
        species_id = ''
        if name in ExpressionNetwork.dict_species:
            species_id = name
        else:
            for key in ExpressionNetwork.dict_species:
                if name in ExpressionNetwork.dict_species[key]:
                    species_id = key
                    break
        return species_id

    ########################################## LNC details #######################################################
    # Get Expression
    # Get Neighbors/Othologs/Hubs info

    # return json with the expression per gene (ENS/LNC)
    # @gene (can be either ensembl ID or LNC  id)
    # @species - Network.MM10/Network.HG38 or aliases
    def getExpressionOrthoHubs(self,gene, species):

        # get species
        species_id = self.__getSpecies(species)
        if len(species_id) == 0:  # species not found
            return None

        # iterate over expression dict
        dict = {}
        dict[gene] = {}
        dict[gene]['expression'] = {}
        for exp_alias in self.dict_expression[species_id]:
            exp_dict = self.__loadJson(self.dict_expression[species_id][exp_alias])
            if exp_dict:

                for expGene in exp_dict:

                    if expGene.startswith(gene):

                        for exp_name in exp_dict[expGene]:  # gene diff exp in the experiment
                            if len(exp_dict[expGene]) > 1:  # check if one or more experiments in the dictionary
                                dict[gene]['expression'][exp_name] = exp_dict[expGene][exp_name]
                            else:
                                dict[gene]['expression'][exp_alias] = exp_dict[expGene][exp_name]

        # get hubs information
        dict_index = self.__loadJson(self.dict_lncDetails[species_id])

        for elem in dict_index:

            if elem.startswith(gene):

                dict_hubs_ortho = self.__loadJson(self.ABS_ROOT + "neighbors/" + dict_index[elem])
                dict[gene].update(dict_hubs_ortho[elem])

        return dict

    ########################################## Network connections #######################################################

    def __splitGeneByType(self, geneList, species):

        dict_ens = self.__loadJson(self.dict_anno[species][ExpressionNetwork.PROTEIN_CODING])

        allGenes = [x[:x.index(".")] for x in dict_ens]

        proteinCodingList = [item for item in geneList if item in allGenes]
        lncList = [item for item in geneList if item.startswith('LNC')]
        mirnaList = list(set(geneList)-set(proteinCodingList)-set(lncList))

        return (lncList,proteinCodingList,mirnaList)


    
    def __getCooexp(self,lnc,protein_coding,species):

        edges_aux = []
        for file in self.dict_coexpression[species]:
            skipHeader = 0
            with open(self.dict_coexpression[species][file],'r') as f:
                for line in f:
                    if skipHeader == 0:
                        skipHeader = 1
                        continue
                    else:
                        arr = line.replace("\n","").split("\t")
                        score = float(arr[0])
                        lnc_id = arr[1]
                        ens_id = arr[5]

                        if "." in ens_id:
                            ens_id = ens_id[:ens_id.index(".")]

                        if lnc_id in lnc and ens_id in protein_coding:
                            # add tuple (LNCID, ENSID, 'coexpression SRR20546', 19)
                            edges_aux.append((lnc_id,ens_id,file,score))

        return edges_aux



    def __getMirTrap(self, lncList, mirnaList, species):

        mirna_dict = {
            ExpressionNetwork.MM10:{
                'miR-103':{
                    'alias':['ENSMUSG00000065563','MI0000587','ENSMUSG00000065553','MI0000587','mmu-miR-103', "miR-103", "MI0007262"],
                    'source':[self.dict_expression[ExpressionNetwork.MM10]['mirTrap mir103']]
                },
                'let-7':{
                    'alias':['ENSMUSG00000105621','Mirlet7f-1','Mirlet7i','ENSMUSG00000065406','Mirlet7d','ENSMUSG00000065453'],
                    'source':[self.dict_expression[ExpressionNetwork.MM10]['mirTrap let7']]
                },
            },
            ExpressionNetwork.HG38:{
                'miR-103':{
                    'alias':['hsa-miR-103'],
                    'source':[self.dict_expression[ExpressionNetwork.HG38]['mirTrap mir103']]
                },
                'let-7':{
                    'alias':['hsa-let-7'],
                    'source':[self.dict_expression[ExpressionNetwork.HG38]['mirTrap let7']]
                }
            }
        }
        edges_mirtrap = []
        for mirna in mirnaList:
            id = ''
            if mirna in mirna_dict[species]:
                id = mirna
            else:
                for m in mirna_dict[species]:
                    if mirna in mirna_dict[species][m]['alias']:
                        id = m
                        break
            if len(id) != 0: # no miRNA found
                for file in mirna_dict[species][id]['source']:
                    dict = self.__loadJson(file)
                    for lnc in lncList:
                        if lnc in dict:
                            for exp in dict[lnc]:
                                edges_mirtrap.append((lnc,mirna,exp,dict[lnc][exp][ExpressionNetwork.LOG2FC]))

        return edges_mirtrap



    def __getNeighbors(self,lncList,species):

        aux_list = []
        edges_neigh = []
        protein_coding_list = []

        dict_index = self.__loadJson(self.dict_lncDetails[species])

        aux_lnc = [] # remove lnc from the list which are not in dictionary
        for lnc in lncList:
            if lnc in dict_index:
                aux_lnc.append(lnc)

        if len(aux_lnc)>0:

            lncList = aux_lnc
            print(lncList)
            dict = self.__loadJson(self.ABS_ROOT + "neighbors/" + dict_index[lncList[0]])
            sub = lncList
            while True:
                for lnc in sub:
                    if lnc not in aux_list and lnc in dict:
                        aux_list.append(lnc)

                        if dict[lnc][ExpressionNetwork.RIGHT] < ExpressionNetwork.INF: # check if dist to right gene ist smaller than 1000000
                            edges_neigh.append((lnc,dict[lnc][ExpressionNetwork.RIGHT_gene],ExpressionNetwork.RIGHT_gene,dict[lnc][ExpressionNetwork.RIGHT]))
                            protein_coding_list.append(dict[lnc][ExpressionNetwork.RIGHT_gene])
                        if dict[lnc][ExpressionNetwork.LEFT] > (-1)*ExpressionNetwork.INF: # check if dist to left gene ist greater than -1000000
                            edges_neigh.append((lnc,dict[lnc][ExpressionNetwork.LEFT_gene],ExpressionNetwork.LEFT_gene,dict[lnc][ExpressionNetwork.LEFT]))
                            protein_coding_list.append(dict[lnc][ExpressionNetwork.LEFT_gene])

                if len(aux_list) == len(lncList):
                    break
                else:
                    sub = list(set(lncList)-set(aux_lnc))
                    if len(sub) > 0:
                        dict = self.__loadJson(self.ABS_ROOT + "neighbors/" + dict_index[sub[0]])
                    else:
                        break


        return (edges_neigh, protein_coding_list)



    # get neighobrs (LNC, ENS, right_gene, dist)
    # get coexp (LNC, ENS, coexp, score)
    # get miRNA - lnc interaction from mirtrap (LNC, mirRNA, mirTrap exp, log2FC)
    def getEdgesFeature(self,geneList,species):

        edges = []

        # get species
        species_id = self.__getSpecies(species)
        if len(species_id) == 0:  # species not found
            return None

        # split by gene type (LNC, protein coding, microRNA)
        (lnc,protein_coding,miRNA) = self.__splitGeneByType(geneList, species_id)

        # get neighbors lnc
        (edges_neighbors,proteinSet) = self.__getNeighbors(lnc,species_id)
        edges.extend(edges_neighbors)
        # add the neighbor proteins to the protein_coding list - compute for those also the coexp
        protein_coding.extend(proteinSet)
        protein_coding = list(set(protein_coding))

        # get cooexpr
        edges.extend(self.__getCooexp(lnc,protein_coding,species_id))

        # check if mirtrap data
        edges.extend(self.__getMirTrap(lnc,miRNA,species_id))

        retedges = []
        for elem in edges:

            newelem = {
                "source": elem[0],
                "target": elem[1]
            }

            for i in range(2, len(elem), 2):
                newelem[elem[i]] = elem[i+1]

            if "." in newelem["source"]:
                newelem["source"] = newelem["source"][:newelem["source"].index(".")]

            if "." in newelem["target"]:
                newelem["target"] = newelem["target"][:newelem["target"].index(".")]

            retedges.append(newelem)


        return retedges


if __name__ == "__main__":


    netw = ExpressionNetwork()

    # call expression
    print(datetime.datetime.utcnow())
    #print(netw.getExpressionOrthoHubs("ENSMUSG00000012848","Mouse"))
    print(datetime.datetime.utcnow())
    #print(netw.getExpressionOrthoHubs("ENSG00000178015.4","Human"))
    print(datetime.datetime.utcnow())
    #print(netw.getExpressionOrthoHubs("LNC_GE_hg38_00044140","Human"))
    print(datetime.datetime.utcnow())
    print(netw.getExpressionOrthoHubs("ENSMUSG00000012848","Mouse"))



    # call edges features in network

    print(datetime.datetime.utcnow())
    print(netw.getEdgesFeature(['LNC_GE_hg38_00125866','ENSG00000163735','LNC_GE_hg38_00056231',
                               "LNC_GE_hg38_00111342",
                               "LNC_GE_hg38_00025954",
                                "LNC_GE_hg38_00106589",
                                "LNC_GE_hg38_00097816",
                                "LNC_GE_hg38_00151536",
                                "LNC_GE_hg38_00058485",
                                "ENSG00000147206"
                                    ],'Human'))

    print(netw.getEdgesFeature(['LNC_GE_hg38_00016618', 'LNC_GE_hg38_00018475', 'LNC_GE_hg38_00075136', 'LNC_GE_hg38_00074648', 'ENSG00000181374', 'LNC_GE_hg38_00018474', 'MI0007262', 'LNC_GE_hg38_00016546'], "human"))
    #
    print(datetime.datetime.utcnow())
    print(netw.getEdgesFeature(['MI0000587','LNC_GE_mm10_00514627'],'Mouse'))    #
    print(datetime.datetime.utcnow())

    print("bla")
    print(netw.getEdgesFeature(['LNC_GE_hg38_00016618', 'LNC_GE_hg38_00018475', 'LNC_GE_hg38_00075136', 'LNC_GE_hg38_00074648', 'ENSG00000181374', 'LNC_GE_hg38_00018474', "MI0007262", 'miR-103', 'LNC_GE_hg38_00016546'],'Human'))    #

