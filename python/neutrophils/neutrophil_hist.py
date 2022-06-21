import argparse
import json
import sys, os
from lxml import etree
from xml.dom import minidom

from neutrophils.allowedSentences import AllowedSentences

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")


import requests

from textmining.SentenceID import SentenceID

parser = argparse.ArgumentParser(description='db query', add_help=False)
parser.add_argument('-o', '--output', type=argparse.FileType("w"), help='outfile', default=None, required=False)

args = parser.parse_args()

fetchAddData = False

def makeSubstr(sent, s, e):

    if len(sent) <= e:
        return ""

    return sent[s:e]



for maxSentDist in [0,1,2,3,4,5]:

    query = {"elements": [{'group': 'NEUTROPHIL', 'name': 'PMN', 'termid': 'PMN'}, {'group': 'NEUTROPHIL', 'name': 'neutrophils', 'termid': 'neutrophils'}], "sentences": str(fetchAddData), "obolevel": 1, "messenger_obolevel": 1, "sentence_distance": maxSentDist}
    r = requests.post("http://localhost:65522/query", data=json.dumps(query))

    res = json.loads(r.content.decode())

    sep = "\t"

    allPMIDs = set()
    allPMC = set()

    print("received", len(res['rels']), "relations", file=sys.stderr)

    for rel in res['rels']:

        for ev in rel['evidences']:

            docid = ev['docid']

            if docid == None:
                continue

            if docid.startswith("PMC"):
                allPMC.add(docid.replace('PMC', ""))
            else:
                allPMIDs.add(docid)

    elemcount = 0

    relEvCount = 0
    evCount = 0
    evPMID = set()
    evSent = set()
    evPMC = set()

    allowedSentsGen = AllowedSentences(maxSentDist)

    for rel in res['rels']:

        for ev in rel['evidences']:

            docid = ev['docid']

            sentid = ev['rel_sentence']
            sentence = ev.get("sentence", "")

            aSent = sentid.split(".")
            aSentNum = int(aSent[-1])

            aSent = aSent[0:2]

            # remove citations/references
            if aSent[0].startswith("PMC") and aSent[1] == '4':
                continue

            if len(sentence) > 0:
                if "///" in sentence:
                    continue

            allowedSentIDs = allowedSentsGen.getDistanceDictBySentID(sentid)

            verbdir = ev['rel_direction_verb']

            lid = ev['lid']
            loid = ev['lontid']

            rid = ev['rid']
            roid = ev['rontid']

            trusts = (ev['trust']['verb'], ev['trust']['stack'], ev['trust']['relex'], ev['trust']['conj'])

            # if len([x for x in trusts[0:3] if x > 0]) == 0:
            #    continue

            textLeft = ""
            textRight = ""

            if len(sentence) > 0:
                textLeft = sentence[ev['lpos'][0]:ev['lpos'][1]]
                textRight = sentence[ev['rpos'][0]:ev['rpos'][1]]

            effect = res['pmidinfo']['categories'].get(docid, [None])
            message = res['pmidinfo']['messengers'].get(docid, [None])

            effectWords = {}
            messageWords = {}

            effectDistances = {}
            messageDistances = {}

            foundEffects = set()
            if effect != None:

                for termEffect in effect:
                    for x in termEffect['evidences']:
                        if x[0] == sentid:
                            effectElem = (termEffect['termid'], termEffect['termname'], x[0], x[1], x[2])
                            foundEffects.add(effectElem)

                            if effectElem in messageDistances:
                                if allowedSentIDs[x[0]] < effectDistances[effectElem]:
                                    effectDistances[effectElem] = allowedSentIDs[x[0]]

                                    if x[0] == sentid:
                                        effectWords[effectElem] = makeSubstr(sentence, x[1], x[2])
                            else:
                                effectDistances[effectElem] = allowedSentIDs[x[0]]

                                if x[0] == sentid:
                                    effectWords[effectElem] = makeSubstr(sentence, x[1], x[2])

                if len(foundEffects) == 0:
                    for termEffect in effect:
                        for x in termEffect['evidences']:
                            if x[0] in allowedSentIDs:

                                effectElem = (termEffect['termid'], termEffect['termname'], x[0], x[1], x[2])
                                foundEffects.add(effectElem)

                                if effectElem in messageDistances:

                                    if allowedSentIDs[x[0]] < effectDistances[effectElem]:
                                        effectDistances[effectElem] = allowedSentIDs[x[0]]

                                        if x[0] == sentid:
                                            effectWords[effectElem] = makeSubstr(sentence, x[1], x[2])


                                else:
                                    effectDistances[effectElem] = allowedSentIDs[x[0]]

                                    if x[0] == sentid:
                                        effectWords[effectElem] = makeSubstr(sentence, x[1], x[2])

            foundMessages = set()
            if message != None:

                for termMessage in message:

                    for x in termMessage['evidences']:
                        if x[0] == sentid:
                            messageElem = (termMessage['termid'], termMessage['termname'], x[0], x[1], x[2])
                            foundMessages.add(messageElem)

                            if messageElem in messageDistances:
                                if messageElem in messageDistances:
                                    if allowedSentIDs[x[0]] < messageDistances[messageElem]:
                                        messageDistances[messageElem] = allowedSentIDs[x[0]]

                                        if x[0] == sentid:
                                            messageWords[messageElem] = makeSubstr(sentence, x[1], x[2])

                            else:
                                messageDistances[messageElem] = allowedSentIDs[x[0]]

                                if x[0] == sentid:
                                    messageWords[messageElem] = makeSubstr(sentence, x[1], x[2])

                if len(foundMessages) == 0:
                    for termMessage in message:
                        for x in termMessage['evidences']:
                            if x[0] in allowedSentIDs:

                                messageElem = (termMessage['termid'], termMessage['termname'], x[0], x[1], x[2])
                                foundMessages.add(messageElem)

                                if messageElem in messageDistances:
                                    if allowedSentIDs[x[0]] < messageDistances[messageElem]:
                                        messageDistances[messageElem] = allowedSentIDs[x[0]]

                                        if x[0] == sentid:
                                            messageWords[messageElem] = makeSubstr(sentence, x[1], x[2])

                                else:
                                    messageDistances[messageElem] = allowedSentIDs[x[0]]
                                    if x[0] == sentid:
                                        messageWords[messageElem] = makeSubstr(sentence, x[1], x[2])

            trustStr = sep.join([str(x) for x in trusts])


            # if len(foundEffects) == 0 or len(foundMessages) == 0:
            #    print("INC EFF MESS", len(foundEffects), len(foundMessages))

            relEvCount += 1

            for messageElem in foundMessages:
                for effectElem in foundEffects:
                    mDist = messageDistances[messageElem]
                    eDist = effectDistances[effectElem]

                    messageWord = messageWords.get(messageElem, "")
                    effectWord = effectWords.get(effectElem, "")

                    if messageElem[2] == effectElem[2]:
                        messageInterval = (messageElem[3], messageElem[4])
                        effectInterval = (effectElem[3], effectElem[4])

                        overlap = max([messageInterval[0], effectInterval[0]]) <= min(
                            [messageInterval[1], effectInterval[1]])

                        if overlap:
                            continue

                    assert(mDist <= maxSentDist)
                    assert(eDist <= maxSentDist)


                    evCount += 1


                    if docid.startswith("PMC"):
                        evPMC.add(docid)
                    else:
                        evPMID.add(docid)

                    evSent.add(sentid)


    print("Distance", maxSentDist)
    print("relEvCount", relEvCount)
    print("evCount", evCount)
    print("evPMID", len(evPMID))
    print("evPMC", len(evPMC))
    print("evSents", len(evSent))
    print("PMIDs", len(allPMIDs))
    print("PMCs", len(allPMC))
    print("\n\n\n")
# Save the file


#print(elemcount, file=sys.stderr)