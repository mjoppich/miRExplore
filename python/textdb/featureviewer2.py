import json
import pprint
from pprint import pprint
from collections import defaultdict

mapping_file = "mapping_keys.json"
with open(mapping_file) as mapping_input:    
	df_mapping = json.load(mapping_input)

noncode_file = "mm10_primary_assembly_and_lncRNA.json"
with open(noncode_file) as noncode_input:    
    df_noncode = json.load(noncode_input)

miranda_file = "miranda_test.tsv"
miranda_input = open(miranda_file, "r")
next(miranda_input) # skip header

def getUniqueID(map_id):
	'''
		input: noncode lncRNA transcript ID: NONMMU... 
		output: unique lncRNA transcript ID: LNC_TR_mm10_...
	'''
	for key, value in df_mapping.items():
		if value["alias"] == map_id:
			return key
	return None

def getLncRNAAnnotations(tr_id):
	for entry in df_noncode:
		gene_id = entry['gene_id']
		for transcript in entry['transcripts_list']:
			transcript_id = transcript['transcript_id']
			if transcript_id == tr_id: 
				transcript_start = int(transcript['start'])
				transcript_stop = int(transcript['stop'])
				transcript_length = transcript_stop - transcript_start
				exon_list = []
				for exons in transcript['exon_list']:
					exon_assembly_start = int(exons['start'])
					exon_assembly_stop = int(exons['stop'])
					exox_length = exon_assembly_stop - exon_assembly_start
					exon_start = exon_assembly_start - transcript_start
					exon_stop = exon_start + exox_length
					exon_list.append({'x':exon_start, 'y':exon_stop})
				return(exon_list, gene_id, transcript_length)
	return None


interactions = {}
for line in miranda_input:
	line_array = line.split("\t")
	mirna = line_array[0]
	lncrna = line_array[1] # noncode lncRNA transcript ID: NONMMU... 
	binding_start = int(line_array[6])
	binding_stop = int(line_array[7])
	interaction_id = mirna + "@@" + lncrna
	if interaction_id in interactions.keys():
		interactions[interaction_id].append({'x':binding_start,'y':binding_stop})
	else:
		interactions[interaction_id] = [{'x':binding_start,'y':binding_stop}]		

json_string = []
for interaction_id, interaction_sites in interactions.items():
	interaction_array = interaction_id.split("@@")
	
	mirna = interaction_array[0]
	lncrna = interaction_array[1]

	unique_transcript_id = getUniqueID(lncrna)
	if unique_transcript_id:
		(exon_list, gene_id, transcript_length) = getLncRNAAnnotations(unique_transcript_id)
		if(exon_list and gene_id and transcript_length):

			container = {}
			container[mirna] = interaction_sites
			container[lncrna] = exon_list
			json_string.append(container)		

		else:
			print("no exons found for " + lncrna + "\t" + unique_transcript_id)
	else:
		print("no mapping if found for " + lncrna)
	
			
print(json.dumps(json_string, indent=4))
#my_json = json.dumps(json_string)
#pprint(my_json)
