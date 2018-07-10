import sys
import os
import json
import pprint
from pprint import pprint
from collections import defaultdict

from textdb.rfamDB import RFamDB


class FeatureViewer:

    """
    Given gene_id (and mirna_id), extracts gene features and mirna interactions
    input: unique gene ID (LNC_GE_... / ENSG...), optional miRNA ID (mmu-let-7g-5p)
    output: json containing gene features and mirna / other mirnas interactions
        [
            {
                "gene_id": "LNC_GE_mm10_00500789",
                "gene_type": "lncrna",
                "gene_length": 938,
                "transcript_list": [
                    {
                        "transcript_id": "LNC_TR_mm10_00607615",
                        "start": 0,
                        "stop": 938,
                        "coding_list": [
                            {
                                "x": 0,
                                "y": 737,
                                "type": "exon"
                            },
                            {
                                "x": 934,
                                "y": 938,
                                "type": "exon"
                            }
                        ],
                        "mirnas": {
                            "mmu-let-7g-5p": [
                                {
                                    "x": "307",
                                    "y": "328"
                                },
                                {
                                    "x": "307",
                                    "y": "328"
                                }
                            ]
                        },
                        "other_mirnas": {
                            "mmu-let-7g-3p": [
                                {
                                    "x": "307",
                                    "y": "328"
                                }
                            ],
                            "mmu-let-7g-3ppp": [
                                {
                                    "x": "307",
                                    "y": "328"
                                }
                            ]
                        }
                    }
                ]
            }
        ]
    """

    def __init__(self, basedir, org, rfamDB=None):

        #my_dir = "/home/mjoppich/ownCloud/data/miRExplore/obodir/"
        #my_dir = "./mm10"

        self.rfamDB = rfamDB
        self.org = org

        #interaction_file = os.path.join(my_dir, "mm10_interactions_allDBs.json")
        interaction_file = os.path.join(basedir, "mm10_interactions_allDBs_and_pc.json")
        with open(interaction_file) as interaction_input:
            self.df_interactions = json.load(interaction_input)

        gene_file = os.path.join(basedir, "mm10_primary_assembly_and_lncRNA.json")
        with open(gene_file) as gene_input:
            self.df_genes = json.load(gene_input)

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

            print("gene: "+gene)
            print("mirna: " + str(mirna))

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

    def getGeneAnnotations(self, my_gene, type, mirna_id):
        json_transcripts = []
        gene_length = 0

        geneInfo = {}

        for gene in self.df_genes:
            gene_id = gene["gene_id"]

            compGeneID = gene_id

            if "." in compGeneID:
                compGeneID = compGeneID[0:compGeneID.index(".")]

            if compGeneID == my_gene:
                gene_start = int(gene['start']) 
                gene_stop = int(gene['stop'])  
                gene_length = gene_stop - gene_start

                geneInfo['gene_length'] = gene_length
                geneInfo['gene_start'] = gene_start
                geneInfo['gene_stop'] = gene_stop
                geneInfo['strand'] = gene['strand']
                geneInfo['chr'] = gene['chr']


                for transcript in gene['transcripts_list']:
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
                    if (type == 'gene'):
                        if transcript['UTR']:
                            # print("found UTR for " + transcript_id)
                            for utr in transcript['UTR']:
                                utr_assembly_start = int(utr['start'])
                                utr_assembly_stop = int(utr['stop'])
                                utr_length = utr_assembly_stop - utr_assembly_start
                                utr_start = utr_assembly_start - transcript_start
                                utr_stop = utr_start + utr_length
                                container_transcript['coding_list'].append({'x': utr_start, 'y': utr_stop, 'type': 'UTR'})
                        if transcript['CDS']:
                            # print("found CDS for " + transcript_id)
                            for cds in transcript['CDS']:
                                cds_assembly_start = int(cds['start'])
                                cds_assembly_stop = int(cds['stop'])
                                cds_length = cds_assembly_stop - cds_assembly_start
                                cds_start = cds_assembly_start - transcript_start
                                cds_stop = cds_start + cds_length
                                container_transcript['coding_list'].append({'x': cds_start, 'y': cds_stop, 'type': 'CDS'})
                    elif transcript['exon_list']:
                        #print("found exons for " + transcript_id)
                        for exons in transcript['exon_list']:
                            exon_assembly_start = int(exons['start'])
                            exon_assembly_stop = int(exons['stop'])
                            exox_length = exon_assembly_stop - exon_assembly_start
                            exon_start = exon_assembly_start - transcript_start
                            exon_stop = exon_start + exox_length
                            container_transcript['coding_list'].append({'x': exon_start, 'y': exon_stop, 'type': 'exon'})
                    else:
                        print("ERROR: No UTR, CDS nor Exons found")

                    if (gene_id and transcript_id):
                        (container_transcript['mirnas'], container_transcript['other_mirnas']) = self.getTranscriptInteractions(transcript_id, gene_id, mirna_id)
                    else:
                        print("ERROR: Invalid gene_id or transcript_id (None Type)")
                    json_transcripts.append(container_transcript)
        return (geneInfo, json_transcripts)


    def getTranscriptInteractions(self, my_transcript, my_gene, my_mirna):
        interaction_list = {}
        other_interaction_list = {}
        for interaction in self.df_interactions:
            if interaction['gene_id'] == my_gene:
                transcripts_list = interaction['transcript_list']
                for transcript in transcripts_list:
                    transcript_id = transcript['transcript_id']
                    if transcript_id == my_transcript:
                        inter_list = transcript['interaction_list']
                        for inter in inter_list:
                            interaction_mirna = inter['mirna'] 
                            if interaction_mirna == my_mirna:
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

            #                     print("No interaction found for " + my_mirna + " and transcript " + my_transcript)
            #         else:
            #             print("No interactions found for transcript " + my_transcript)
            # else:
            #     print("No interactions found for gene " + my_gene)

        return (interaction_list, other_interaction_list)

if __name__ == '__main__':

    if len(sys.argv) >= 3:
        input_gene = sys.argv[1]
        input_mirna = sys.argv[2]
    else:
        input_gene = "ENSMUSG00000045382"
        input_mirna = "" 

    rfDB = RFamDB.loadFromFile('/mnt/c/ownCloud/data/miRExplore/textmine/aggregated_pmid/rfam.regions.mirexplore')

    fv = FeatureViewer("/home/mjoppich/ownCloud/data/miRExplore/obodir/", "mmu", rfamDB=rfDB)

    json_result = fv.getFeatures(input_gene, input_mirna)
    '''
        do stuff with json json_result
    '''
    print(json.dumps(json_result, indent=4))

# python3 featureviewer.py ENSMUSG00000056476.13 mmu-let-7g-5p
