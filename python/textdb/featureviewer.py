import sys
import os
import json

#from textdb.rfamDB import RFamDB
from mongodb.pymongodb import MongoDB
from textdb.rfamDB import RFamDB


class FeatureViewer:

    #def __init__(self, my_gene, my_mirna, basedir, org, rfamDB=None):
    def __init__(self, org, basedir=None, rfamDB = None):

        self.rfamDB = rfamDB
        self.org = org

        self.mongo = None

        if basedir == None:
            print("FeatureViewer:", "db based")


            self.mongo = MongoDB('minglerna')



        else:
            print("FeatureViewer:", "file based")

            interaction_file = os.path.join(basedir, "mm10_interactions_allDBs_and_pc.json")
            with open(interaction_file) as interaction_input:
                self.df_interactions = json.load(interaction_input)

            gene_file = os.path.join(basedir, "mm10_primary_assembly_and_lncRNA.json")
            with open(gene_file) as gene_input:
                self.df_features = json.load(gene_input)

    def queryDBForGeneID(self, my_collection, geneID):

        collection = self.mongo.getCollection(my_collection)
        result = self.mongo.getJsonObjectForGeneID(collection, geneID)
        json_result = json.loads(result)

        return json_result

    def getFeatures(self, gene, mirna):

        valid_gene = 1
        gene_type = None
        if 'ENS' in gene:
            gene_type = 'gene'

        elif 'LNC' in gene:
            gene_type = 'lncrna'
        else:
            valid_gene = 0

        container = {}
        json_res = []

        if (valid_gene == 1):

            #print("gene: " + self.gene)
            #print("mirna: " + str(self.mirna))

            container['gene_id'] = gene
            container['gene_type'] = gene_type

            (geneInfo, container['transcript_list']) = self.getGeneAnnotations(gene, gene_type, mirna)

            for elemID in geneInfo:
                container[elemID] = geneInfo[elemID]

            if 'chr' in container and self.rfamDB != None:
                rfams = self.rfamDB.get_entries(self.org, container['chr'], container['gene_start'], container['gene_stop'], container['strand'])
                container['rfams'] = rfams

            json_res.append(container)
        else:
            print("ERROR: Gene ID not recognized " + gene)
        return json_res

    def getGeneAnnotations(self, gene, my_type, mirna):

        if self.mongo:

            self.df_features = self.queryDBForGeneID('assembly', gene)
            self.df_interactions = self.queryDBForGeneID('interactions', gene)


        json_transcripts = []
        gene_length = 0

        geneInfo = {}

        for geneElem in self.df_features:

            gene_id = geneElem["gene_id"]

            compGeneID = gene_id

            if "." in compGeneID:
                compGeneID = compGeneID[0:compGeneID.index(".")]

            if compGeneID == gene:

                gene_start = int(geneElem['start'])
                gene_stop = int(geneElem['stop'])
                gene_length = gene_stop - gene_start

                geneInfo['gene_length'] = gene_length
                geneInfo['gene_start'] = gene_start
                geneInfo['gene_stop'] = gene_stop
                geneInfo['strand'] = geneElem['strand']
                geneInfo['chr'] = geneElem['chr']

                for transcript in geneElem['transcripts_list']:
                    transcript_id = transcript['transcript_id']
                    transcript_start = int(transcript['start'])
                    transcript_stop = int(transcript['stop'])
                    transcript_length = transcript_stop - transcript_start

                    container_transcript = {}
                    container_transcript['transcript_id'] = transcript_id
                    container_transcript['start'] = transcript_start - gene_start
                    container_transcript['stop'] = transcript_stop - transcript_start
                    container_transcript['coding_list'] = []
                    container_transcript['mirnas'] = {}
                    container_transcript['other_mirnas'] = {}

                    if (my_type == 'gene'):
                        if transcript['UTR']:
                            for utr in transcript['UTR']:
                                utr_assembly_start = int(utr['start'])
                                utr_assembly_stop = int(utr['stop'])
                                utr_length = utr_assembly_stop - utr_assembly_start
                                utr_start = utr_assembly_start - transcript_start
                                utr_stop = utr_start + utr_length
                                container_transcript['coding_list'].append({'x': utr_start, 'y': utr_stop, 'type': 'UTR'})
                        if transcript['CDS']:
                            for cds in transcript['CDS']:
                                cds_assembly_start = int(cds['start'])
                                cds_assembly_stop = int(cds['stop'])
                                cds_length = cds_assembly_stop - cds_assembly_start
                                cds_start = cds_assembly_start - transcript_start
                                cds_stop = cds_start + cds_length
                                container_transcript['coding_list'].append({'x': cds_start, 'y': cds_stop, 'type': 'CDS'})
                    elif transcript['exon_list']:
                        for exons in transcript['exon_list']:
                            exon_assembly_start = int(exons['start'])
                            exon_assembly_stop = int(exons['stop'])
                            exox_length = exon_assembly_stop - exon_assembly_start
                            exon_start = exon_assembly_start - transcript_start
                            exon_stop = exon_start + exox_length
                            container_transcript['coding_list'].append({'x': exon_start, 'y': exon_stop, 'type': 'exon'})
                    else:
                        print("ERROR: No UTR, CDS nor Exons found")

                    (container_transcript['mirnas'], container_transcript['other_mirnas']) = self.getTranscriptInteractions(mirna)

                    json_transcripts.append(container_transcript)
        return (geneInfo, json_transcripts)


    def getTranscriptInteractions(self, mirna):
        interaction_list = {}
        other_interaction_list = {}

        for interaction in self.df_interactions:




            for transcript in interaction['transcript_list']:
                for inter in transcript['interaction_list']:
                    interaction_mirna = inter['mirna'] 
                    if interaction_mirna == mirna:
                        interaction_list[interaction_mirna] = []
                        for alignment in inter['alignment_list']:
                            x = alignment['lnc_start']
                            y = alignment['lnc_end']
                            interaction_list[interaction_mirna].append({'x': x, 'y': y})
                    else:
                        other_interaction_list[interaction_mirna] = []
                        for alignment in inter['alignment_list']:
                            x = alignment['lnc_start']
                            y = alignment['lnc_end']
                            other_interaction_list[interaction_mirna].append({'x': x, 'y': y})  
        return (interaction_list, other_interaction_list)

if __name__ == '__main__':

    if len(sys.argv) >= 3:
        input_gene = sys.argv[1]
        input_mirna = sys.argv[2]
    else:
        input_gene = "ENSMUSG00000045382"
        input_mirna = "" 

    rfDB = RFamDB.loadFromFile('/mnt/c/ownCloud/data/miRExplore/textmine/aggregated_pmid/rfam.regions.mirexplore')

    fv = FeatureViewer( "mmu", "/home/mjoppich/ownCloud/data/miRExplore/obodir/", rfamDB=rfDB)
    #fv = FeatureViewer(input_gene, input_mirna, "/home/mjoppich/ownCloud/data/miRExplore/obodir/", "mmu")

    json_result = fv.getFeatures(input_gene, input_mirna,)

    print(json.dumps(json_result, indent=4))

# python3 featureviewer.py ENSMUSG00000056476.13 mmu-let-7g-5p
