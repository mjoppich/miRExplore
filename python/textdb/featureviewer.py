import sys
import os
import json
import pprint
from pprint import pprint
from collections import defaultdict


class FeatureViewer:

    def __init__(self):

        my_dir = "/home/mjoppich/ownCloud/data/miRExplore/obodir/map_data"

        interaction_file = os.path.join(my_dir, "miranda_interactions.json")
        with open(interaction_file) as interaction_input:
            self.df_interactions = json.load(interaction_input)

        mapping_file = os.path.join(my_dir, "mapping_keys.json")
        with open(mapping_file) as mapping_input:
            self.df_mapping = json.load(mapping_input)

        gene_file = os.path.join(my_dir, "mm10_primary_assembly_and_lncRNA.json")
        with open(gene_file) as gene_input:
            self.df_genes = json.load(gene_input)



    def get_features(self, my_gene, my_mirna):

        valid_gene = 1
        gene_id = my_gene
        if 'ENS' in gene_id:
            gene_type = 'gene'
            gene_id = my_gene
        else:
            gene_type = 'lncrna'
            gene_id = self.getMappedToUniqueID(gene_id)
            if (gene_id):
                print(gene_id)
            else:
                valid_gene = 0
                print("ERROR: No mapping found for " + gene_id)


        container = {}
        json_string = []

        if (valid_gene == 1):
            container['gene_id'] = gene_id
            container['gene_type'] = gene_type
            container['gene_length'] = '2000'
            container['transcript_list'] = self.getGeneAnnotations(gene_id, gene_type, my_mirna)
            json_string.append(container)
            print(json.dumps(json_string, indent=4))
        else:
            print("ERROR: Gene ID not recognized " + gene_id)

    def getTranscriptInteractions(self, my_transcript, my_gene, my_mirna):
        interaction_list = {}
        for interaction in self.df_interactions:
            if interaction["gene"] == my_gene:
                transcripts_list = interaction["transcript_list"]
                for transcript in transcripts_list:
                    transcript_id = transcript["transcript_id"]
                    if transcript_id == my_transcript:
                        inter_list = transcript["interaction_list"]
                        for inter in inter_list:
                            mirna_interaction = {}
                            if inter["mirna"] == my_mirna:  # remove if you want all mirnas
                                interaction_list[inter['mirna']] = {'mirna': inter['mirna'], 'x': inter['lnc_start'], 'y': inter['lnc_end']}
                            else:
                                print("ERROR: No interaction found for " + my_mirna + " and transcript " + my_transcript)
                    else:
                        print("ERROR: No interactions found for transcript " + my_transcript)
            else:
                print("ERROR: No interactions found for gene " + my_gene)

        return interaction_list


    def getMappedToUniqueID(self, my_id):
        '''
            input: noncode lncRNA transcript ID: NONMMU...
            output: unique lncRNA transcript ID: LNC_TR_mm10_...
        '''
        for key, value in self.df_mapping.items():
            if value["alias"] == my_id:
                return key
        return None


    def getMappedFromUniqueID(self, my_id):
        '''
            input: unique lncRNA transcript ID: LNC_TR_mm10_...
            output: noncode lncRNA transcript ID: NONMMU...
        '''
        for key, value in self.df_mapping.items():
            if key == my_id:
                return value["alias"]
        return None


    def getGeneAnnotations(self, my_gene, type, mirna_id):
        json_transcripts = []

        for gene in self.df_genes:
            gene_id = gene["gene_id"]
            if (gene_id == my_gene):
                print(gene_id)
                gene_start = int(gene['start'])  # 92103000
                gene_stop = int(gene['stop'])  # 92105221
                for transcript in gene['transcripts_list']:
                    transcript_id = transcript['transcript_id']
                    transcript_start = int(transcript['start'])  # 92103000
                    transcript_stop = int(transcript['stop'])  # 92103102
                    transcript_length = transcript_stop - transcript_start  # 102

                    container_transcript = {}
                    container_transcript['transcript_id'] = transcript_id
                    container_transcript['start'] = transcript_start - gene_start
                    container_transcript['stop'] = transcript_stop - transcript_start
                    container_transcript['coding_list'] = []
                    container_transcript['mirnas'] = []
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
                        # getTranscriptInteractionsForGene(transcript_id, gene_id)
                        if transcript['CDS']:
                            # print("found CDS for " + transcript_id)
                            for cds in transcript['CDS']:
                                print(cds.get('start') + "  stop " + cds['stop'])
                                cds_assembly_start = int(cds['start'])
                                cds_assembly_stop = int(cds['stop'])
                                cds_length = cds_assembly_stop - cds_assembly_start
                                cds_start = cds_assembly_start - transcript_start
                                cds_stop = cds_start + cds_length
                                container_transcript['coding_list'].append({'x': cds_start, 'y': cds_stop, 'type': 'CDS'})
                    # container_transcript['mirnas'] = getTranscriptInteractions(transcript_id, my_gene)
                    elif transcript['exon_list']:
                        # print("found exons for " + transcript_id)
                        for exons in transcript['exon_list']:
                            exon_assembly_start = int(exons['start'])
                            exon_assembly_stop = int(exons['stop'])
                            exox_length = exon_assembly_stop - exon_assembly_start
                            exon_start = exon_assembly_start - transcript_start
                            exon_stop = exon_start + exox_length
                            container_transcript['coding_list'].append({'x': exon_start, 'y': exon_stop, 'type': 'exon'})
                        # map back to old IDs (remove once Miranda is parsed)
                        transcript_id = self.getMappedFromUniqueID(transcript_id)
                        gene_id = self.getMappedFromUniqueID(my_gene)
                    else:
                        print("ERROR: No UTR, CDS nor Exons found")

                    if (gene_id and transcript_id):
                        container_transcript['mirnas'] = self.getTranscriptInteractions(transcript_id, gene_id, mirna_id)
                    else:
                        print("ERROR: Invalid gene_id or transcript_id (None Type)")
                    json_transcripts.append(container_transcript)
        return json_transcripts


if __name__ == '__main__':


    fv = FeatureViewer()

    #fv.get_features('NONMMUG078714.1', 'mmu-let-7g-5p')
    fv.get_features('ENSMUSG00000079800.2', None)




# python3 featureviewer.py "NONMMUG078714.1" "mmu-let-7g-5p"

# LNC_GE_mm10_00566701 --> NONMMUG078714.1
# LNC_TR_mm10_00704988 --> NONMMUT125591.1
# LNC_TR_mm10_00704989 --> NONMMUT125592.1



