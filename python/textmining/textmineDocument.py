import argparse
import string

import ahocorasick

from synonymes.SynonymFile import Synfile
from utils.idutils import ParseObject, dataDir

parser = argparse.ArgumentParser(description='Textmine documents')
parser.add_argument('-i', '--input', type=str, help='inputfile', required=True)
parser.add_argument('-s', '--synonyms', nargs='+', type=str, help='syn files', required=True)
parser.add_argument('-c', '--characters', type=str, required=False, default='-()[]{}')

#args = parser.parse_args()

args = ParseObject()
args.synonyms = [ dataDir + "/miRExplore/textmine/synonyms/mirna_mimat.syn"]
args.characters = ' -()[]{}=!"§$%&/=?+*\'#-,.;:'

synfileMap = {}
mapSynfile = {}
for file in args.synonyms:
    idx = len(synfileMap)

    synFile = Synfile(file)

    synfileMap[idx] = synFile
    mapSynfile[synFile] = idx



A = ahocorasick.Automaton()


for idx in synfileMap:

    synFile = synfileMap[idx]

    for synonym in synFile:

        for synWord in synonym.syns:
            A.add_word( synWord , (idx, synWord, synonym) )
            A.add_word( synWord.upper() , (idx, synWord.upper(), synonym) )



A.make_automaton()

allowedWordChars = set(string.ascii_letters).union({str(x) for x in range(0,10)}).union({x for x in args.characters})

testSentence = "hsa-miR-155 is possibly contained in this sentence in contrast to HSA-MIR-155."

for end_index, (insert_order, synWord, original_value) in A.iter(testSentence):
    start_index = end_index - len(synWord) + 1

    testIdx = end_index+1
    if testIdx < len(testSentence):
        testChar = testSentence[testIdx]
        accept =  testChar in args.characters
    else:
        accept = True

    print((start_index, end_index, accept, (insert_order, synWord, original_value)))

