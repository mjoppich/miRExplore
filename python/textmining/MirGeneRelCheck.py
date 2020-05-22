from collections import defaultdict
import regex as re

#nlp = spacy.load("/mnt/d/spacy/models/en_ner_bionlp13cg_md-0.2.4/en_ner_bionlp13cg_md/en_ner_bionlp13cg_md-0.2.4")


class SentenceRelationChecker:

    def __init__(self, nlp):
        self.relCheck = MirGeneRelCheck()
        self.nlp = nlp

    def fixPos(self, pos, sentence):
        posElems = [(idx, i, ord(i)) for idx, i in enumerate(sentence) if ord(i) > 127 and idx < pos]

        if len(posElems) > 0:
            pos -= len(posElems)

        return pos


    def __offset_positions(self, entDict, value):
        entDict["entity_location"] = (entDict["entity_location"][0] + value, entDict["entity_location"][1] + value)
        return entDict

    def processLocation(self, sent, e1p, e2p, verbose=False, relClassifier=None):

        osent = sent
        sent = sent.lstrip()

        if len(osent) != len(sent):
            diff = len(osent) - len(sent)
            e1p = self.__offset_positions(e1p, -diff)
            e2p = self.__offset_positions(e2p, -diff)

        sent = sent.rstrip()

        if not sent.endswith("."):
            sent += "."

        doc = self.nlp(sent)

        e1Token = None
        e2Token = None

        tOvl = 0
        for t in doc:
            tInt = (t.idx, t.idx + len(t))
            eInt = e1p["entity_location"]

            if eInt[0] < tInt[0]:
                tmp = tInt
                tInt = eInt
                eInt = tmp

            if eInt[0] < tInt[1]:
                ovl = tInt[1] - eInt[0]

                if ovl > tOvl:
                    e1Token = t

        tOvl = 0
        for t in doc:
            tInt = (t.idx, t.idx + len(t))
            eInt = e2p["entity_location"]

            if eInt[0] < tInt[0]:
                tmp = tInt
                tInt = eInt
                eInt = tmp

            if eInt[0] < tInt[1]:
                ovl = tInt[1] - eInt[0]

                if ovl > tOvl:
                    e2Token = t

        if e1Token == None or e2Token == None:
            return None, None

        if e1p["entity_type"] == "gene":
            geneWord = e1Token
            mirWord = e2Token

        else:
            geneWord = e2Token
            mirWord = e1Token

        if verbose:
            print(e1Token, e1Token.i, e2Token, e2Token.i)
            print("(u\"{}\", {}, {})".format(doc.text_with_ws, e1Token.i, e2Token.i))
        
            print("geneword", geneWord, "mirword", mirWord)

        checkResults = self.relCheck.checkRelation(doc, mirWord, geneWord, verbose=verbose, relClassifier=relClassifier)

        if verbose:
            print(checkResults)

        acceptInteraction = checkResults[0]

        ent1 = doc.text_with_ws[e1p["entity_location"][0]:e1p["entity_location"][1]]
        ent2 = doc.text_with_ws[e2p["entity_location"][0]:e2p["entity_location"][1]]

        aSent = sent  # doc.text_with_ws

        esentence = "".join([aSent[:e1p["entity_location"][0]],
                             "<{}>".format(e1p["entity_type_token"]), ent1, "</{}>".format(e1p["entity_type_token"]),
                             aSent[e1p["entity_location"][1]:e2p["entity_location"][0]],
                             "<{}>".format(e2p["entity_type_token"]), ent2, "</{}>".format(e2p["entity_type_token"]),
                             aSent[e2p["entity_location"][1]:]])

        if verbose:
            print(" hit sentence", doc.text_with_ws)
            print("full sentence", esentence)

        retDict = {"accept_relation": acceptInteraction, "check_results": checkResults[1],"tagged_sent": esentence}

        return retDict


    def check_sentence(self, sentence, e1, e2, fix_special_chars=True, checkSigPahtway=True, verbose=False, relClassifier=None):
        """

        :param sentence: a sentence to process
        :param e1: entity 1 dict: ('entType': [mirna|gene], 'entLocation': (startpos, endpos), 'entTypeToken: "e1|e2")
        :param e2: entity 1 dict: ('entType': [mirna|gene], 'entLocation': (startpos, endpos), 'entTypeToken: "e1|e2")
        :return:
        """

        assert (all([x in e1 for x in ["entity_type", "entity_location", "entity_type_token"]]))
        assert (all([x in e2 for x in ["entity_type", "entity_location", "entity_type_token"]]))

        if fix_special_chars:
            e1Pos = (self.fixPos(e1["entity_location"][0], sentence), self.fixPos(e1["entity_location"][1], sentence))
            e2Pos = (self.fixPos(e2["entity_location"][0], sentence), self.fixPos(e2["entity_location"][1], sentence))

        else:
            e1Pos = e1["entity_location"]
            e2Pos = e2["entity_location"]

        e1["entity_location"] = tuple(e1Pos)
        e2["entity_location"] = tuple(e2Pos)

        if not e1Pos[0] < e2Pos[0]:
            tmpE = e1
            e1 = e2
            e2 = tmpE

        ent1 = sentence[e1["entity_location"][0]:e1["entity_location"][1]]
        ent2 = sentence[e2["entity_location"][0]:e2["entity_location"][1]]

        e1['entity_text'] = ent1
        e2['entity_text'] = ent2

        if e1["entity_type"] == "mirna":
            if ent1.startswith("miR") and len(ent1) > 4 and ent1[3] == " ":
                ent1 = list(ent1)
                sentence = list(sentence)

                ent1[3] = '-'
                sentence[e1["entity_location"][0]+3] = '-'

                ent1 = "".join(ent1)
                sentence = "".join(sentence)

                if verbose:
                    print("fixed sentence", sentence)

        elif e2["entity_type"] == "mirna":
            if ent2.startswith("miR") and len(ent2) > 4 and ent2[3] == " ":
                ent2 = list(ent2)
                sentence = list(sentence)

                ent2[3] = '-'
                sentence[e2["entity_location"][0]+3] = '-'

                ent2 = "".join(ent2)
                sentence = "".join(sentence)

                if verbose:
                    print("fixed sentence", sentence)


        assert (e1["entity_type_token"] != e2["entity_type_token"])
        assert (e1["entity_location"][0] < e2["entity_location"][0])

        """
        if e1[4] == "e1":
            ent1 = "REGULATOR"
        else:
            ent1 = "REGULATEE"

        if e2[4] == "e1":
            ent2 = "REGULATOR"
        else:
            ent2 = "REGULATEE"
        """

        for entity in [ent1, ent2]:
            if verbose:
                print("add special case", entity)
            self.nlp.tokenizer.add_special_case(entity, [{'ORTH': entity, 'TAG': 'NNP'}])


        """
        Here we split the sentence - maybe also useful for and, concatenated sentences!
        """
        splitSigns = [": ", "; "]
        accSemics = set()

        foundSplitSigns = [x for x in splitSigns if x in sentence]

        if len(foundSplitSigns) > 0:

            semics = [x for x in re.finditer("({})[A-Za-z0-9]".format("|".join(splitSigns)), sentence)]
            for x in semics:
                xreg = sentence[max(0, x.regs[0][0] - 10):min(x.regs[0][1] + 10, len(sentence))]
                if "(" in xreg and ")" in xreg:
                    continue

                accSemics.add(x.regs[0][0])


        sentFrags = []

        if len(accSemics) >= 0:
            lastSplit = 0

            for x in sorted(accSemics):
                splitAt = x#x.regs[0][0]
                sentFrags.append((sentence[lastSplit:splitAt], lastSplit, splitAt))
                lastSplit = splitAt + 1

            if lastSplit < len(sentence):
                sentFrags.append((sentence[lastSplit:], lastSplit, len(sentence)))

        else:
            if verbose:
                print("no acc semics")

        commonFragFound = False

        for sidx, x in enumerate(sentFrags):
            e1Inside = x[1] <= e1['entity_location'][0] <= e1['entity_location'][1] <= x[2]
            e2Inside = x[1] <= e2['entity_location'][0] <= e2['entity_location'][1] <= x[2]

            if e1Inside and e2Inside:

                e1Token = None
                e2Token = None

                e1p = e1.copy()
                e2p = e2.copy()

                e1p = self.__offset_positions(e1p, -x[1])
                e2p = self.__offset_positions(e2p, -x[1])


                pDict = self.processLocation(x[0], e1p, e2p, verbose=verbose, relClassifier=relClassifier)

                if pDict['tagged_sent'] != None and pDict['accept_relation'] != None:
                    commonFragFound = True
                    fullsentence = "; ".join([x[0] for x in sentFrags[:sidx]] + [pDict['tagged_sent']] + [x[0] for x in sentFrags[sidx + 1:]])
                    pDict['full_sentence'] = fullsentence

                    pDict["entity1"] = e1p
                    pDict["entity2"] = e2p

                    break

        if not commonFragFound:

            pDict = self.processLocation(sentence, e1, e2, verbose=verbose, relClassifier=relClassifier)
            acceptInteraction = pDict['accept_relation']
            fullsentence = pDict['tagged_sent']
            pDict['full_sentence'] = fullsentence

            pDict["entity1"] = e1
            pDict["entity2"] = e2

            if not ": " in foundSplitSigns and len(sentFrags) > 2:
                pDict['accept_relation'] = False
                pDict["sentfrags_check"] = False


        return pDict


class MirGeneRelCheck:

    def __init__(self):

        self.surGeneContexts = [
            (re.compile('(%s){e<=3}' % "signaling pathway"), 30),
            (re.compile('(%s){e<=3}' % "signaling circuit"), 22),
            (re.compile('(%s){e<=3}' % "driven pathogenic signaling"), 30),
            (re.compile('(%s){e<=3}' % "receptor signaling"), 30),
            (re.compile('(%s){e<=2}' % "pathway"), 15),
            (re.compile('(%s){e<=1}' % "generation"), 12),
            #(re.compile('(%s){e<=1}' % "-mediated"), 15),
            (re.compile('(%s){e<=1}' % "-sensitive"), 15),
            (re.compile('(%s){e<=1}' % "-driven"), 13),
            (re.compile('(%s){e<=1}' % "adenomas"), 15),
            (re.compile('(%s){e<=1}' % "transgenic"), 15),
            (re.compile('(%s){e<=1}' % "signaling through"), 22),
            (re.compile('(%s){e<=1}' % "signaling"), 10),
            (re.compile('(%s){e<=1}' % "knockout"), 15),
            (re.compile('(%s){e<=1}' % "treatment"), 10),
            (re.compile('(%s){e<=1}' % "strains"), 20), # was 15
            (re.compile('(%s){e<=1}' % "treated"), 15),
            (re.compile('(%s){e<=1}' % "peptides"), 10),
            (re.compile('(%s){e<=1}' % "-positive"), 15),
            (re.compile('(%s){e<=1}' % "-negative"), 15),
            (re.compile('(%s){e<=1}' % "-producing"), 12),
            (re.compile('(%s){e<=1}' % "stressed.{1-7}cells"), 30),
            (re.compile('(%s){e<=2}' % "expressing vector"), 25),
            (re.compile('(%s){e<=1}' % "\(-/-\)"), 10),
            (re.compile('(%s){e<=1}' % "\(-/+\)"), 10),
            (re.compile('(%s){e<=1}' % "\(+/-\)"), 10),
            (re.compile('(%s){e<=1}' % "\(+/+\)"), 10),
            (re.compile('(%s){e<=1}' % "copy number"), 20),
            (re.compile('(%s){e<=0}' % " ratio"), 6),
            (re.compile('(%s){e<=0}' % "rs[0-9]+"), 10),
            (re.compile('(%s){e<=0}' % "family member"), 16),
        ]

        self.surMiRContexts = [
            (re.compile('(%s){e<=2}' % "knockout"), 15),
            (re.compile('(%s){e<=2}' % "treatment"), 15),
            (re.compile('(%s){e<=1}' % "seed site"), 15),
            (re.compile('(%s){e<=2}' % "luciferase reporter"), 30),
            (re.compile('(%s){e<=2}' % "transfected"), 30),
            (re.compile('(%s){e<=3}' % "family microRNAs"), 20),
            (re.compile('(%s){e<=0}' % "family"), 8),
            (re.compile('(%s){e<=0}' % "miRNA family"), 14),
            (re.compile('(%s){e<=0}' % "-targeted"), 9),

        ]

        self.surPreGeneContexts = [
            (re.compile('(%s){e<=2}' % "analysis of active"), 30),
            (re.compile('(%s){e<=1}' % "downstream .* target"), 35),
            (re.compile('(%s){e<=1}' % "downstream .* molecule"), 35),
            
        ]

        self.surPreMiRContexts = [
            (re.compile('(%s){e<=0}' % "near"), 10),
            (re.compile('(%s){e<=0}' % "but not"), 8),
            (re.compile('(%s){e<=1}' % "gene locus of"), 25),
            #(re.compile('(%s){e<=2}' % "binding sites"), 17),
        ]

        pass

    def __findByDep(self, token, dep):
        tks = []
        for x in token.children:
            if x.dep_ in dep:
                tks.append(token)
        
        return tks

    def _deepCheck(self, doc, mirword, geneword, verbose):

        sentenceCompartments = []
        thisSubtree = [t for t in doc]

        allProcessedTokens = set()

        for token in doc:


            createCompartment = False

            CONJDEP=["cconj", "xconj", "conj", "ccomp", "parataxis", "advcl", "xcomp"]

            if token.pos_ in ["VERB", "AUX"] and token.dep_ in CONJDEP: # remove "xcomp": appears to act as ... #"acl:relcl"

                #createCompartment = createCompartment or True

                # do not split: by, including

                newSubtree = [x for x in token.subtree]
                newSubtree = sorted(newSubtree, key=lambda x: x.idx)

                subsentCheck = True#any(["subj" in x.dep_ for x in newSubtree])

                
                for i in range(0, min(len(newSubtree), 2)):

                    if newSubtree[i].pos_ in ["NOUN", "PROPN"]:
                        if verbose:
                            print("deepcheck connection rule noun found", token, newSubtree[i], newSubtree)
                        break

                    if str(newSubtree[i]).lower() in ["by", "including", "because", "through", "of", "via", "to"]:
                        
                        if verbose:
                            print("deepcheck connection rule", newSubtree)

                        subsentCheck = False
                        break

                if token.dep_ in ["conj"]:
                    #should only be applied to "and" conjunction
                    allToks = [(t.i, t) for t in doc]


                    dobjs = self.__findByDep(token.head, ["dobj"])

                    if len(dobjs) > 0:

                        allelems = [x for x in dobjs]
                        for x in dobjs:
                            chds = self.__followChild(x, ["det", "amod"])

                            for c in chds:
                                if not c in allelems:
                                    allelems.append(c)

                        #print("deep check deps", allelems)
                        if mirword in allelems or geneword in allelems:
                            subsentCheck = False # 

                    else:
                        #subsentCheck = False # do not split!
                        pass
                    

                    for tidx, t in allToks:
                        if tidx < token.i and token.i -3 < tidx:
                            if t.dep_ in ["cc"] and str(t) in ["and"]:
                                
                                leftTree = [x for x in token.lefts]

                                #if verbose:
                                #    print("deepcheck left rule", token, leftTree)

                                nounFound = False
                                for x in leftTree:
                                    if x.pos_ in ["NOUN", "PROPN"]:
                                        nounFound = True
                                        break

                                #if not nounFound:
                                #    subsentCheck = False

                                break

                    

                #subsentCheck = subsentCheck and nounFound
                createCompartment = createCompartment or subsentCheck



            elif token.dep_ in ["amod"]:
                newSubtree = [x for x in token.subtree]
                newSubtree = sorted(newSubtree, key=lambda x: x.idx)

                amodCheck = False
                for i in range(0,min(len(newSubtree),3)):
                    if str(newSubtree[i]).lower() in ["while", "whereas", "thereby", "resulting", "suggestive", "whereby"]:
                        amodCheck = True
                        break

                createCompartment = createCompartment or amodCheck

            elif token.pos_ in ["VERB", "AUX"] and token.dep_ in ["acl:relcl"]:
                newSubtree = [x for x in token.subtree]
                newSubtree = sorted(newSubtree, key=lambda x: x.idx)

                aclclauseCheck = False
                for i in range(0,min(len(newSubtree),3)):
                    if str(newSubtree[i]).lower() in ["whereby"]:
                        aclclauseCheck = True
                        break

                #not working!
                createCompartment = createCompartment or aclclauseCheck

            elif token.dep_ in ["conj"]:
                newSubtree = [x for x in token.subtree]
                newSubtree = sorted(newSubtree, key=lambda x: x.idx)

                conjCheck = False
                for i in range(0,min(len(newSubtree),3)):
                    if str(newSubtree[i]).lower() in ["actually"]:
                        conjCheck = True
                        break

                createCompartment = createCompartment or conjCheck


            if createCompartment:
                newSubtree = [x for x in token.subtree]
                hasSubj = len([x for x in newSubtree if "subj" in x.dep_]) > 0

                #if len(newSubtree) < 7:
                #    continue

                """
                prechecks completed
                """

                thisSubtree = newSubtree

                #thisSubtree must have a subj?

                thisAncestor = [a for a in token.ancestors]

                ancestorSubtree = []
                for anc in thisAncestor:
                    ast = [x for x in anc.subtree if not (x in allProcessedTokens or x in thisSubtree)]

                    for x in ast:
                        allProcessedTokens.add(x)

                    ancestorSubtree += ast

                ancestorSubtree = [x for x in sorted(ancestorSubtree, key=lambda x: x.idx)]

                #print("atree", ancestorSubtree)
                wasAdded = False
                for x in ancestorSubtree:
                    if str(x) in [";"]:
                        lVerb = "VERB" in [y.pos_ for y in ancestorSubtree if y.i < x.i] or "AUX" in [y.pos_ for y in ancestorSubtree if y.i < x.i]
                        rVerb = "VERB" in [y.pos_ for y in ancestorSubtree if y.i > x.i] or "AUX" in [y.pos_ for y in ancestorSubtree if y.i > x.i]

                        if lVerb and rVerb:
                            wasAdded = True

                            sentenceCompartments.append([y for y in ancestorSubtree if y.i < x.i])
                            sentenceCompartments.append([y for y in ancestorSubtree if y.i > x.i])

                            if verbose:
                                print("compartment split by verb")

                if not wasAdded:
                    sentenceCompartments.append(ancestorSubtree)



        if len(thisSubtree) > 0:
            sentenceCompartments.append(thisSubtree)

        remainingTokens = []
        for x in doc:
            wasFound = False
            
            for comp in sentenceCompartments:
                if x in comp:
                    wasFound = True
                    break
            
            if not wasFound:
                remainingTokens.append(x)

        if len(remainingTokens) > 0:
            sentenceCompartments.append(remainingTokens)

        splitPositions = []
        finalCompartments = []

        #remove pos_ PUNCT,  CCONJ from tail
        #remove pos_ SCONJ from head
        for scIdx, compart in enumerate(sentenceCompartments):
            origcompart = compart

            othercomparts = []
            for comps in sentenceCompartments[scIdx+1:]:
                othercomparts += comps

            compart = [x for x in compart if not x in othercomparts]


            kpart = []
            for x in compart:
                if len(kpart) == 0 and x.pos_ in ["PUNCT", "SCONJ", "CCONJ"]:
                    continue

                kpart.append(x)

            kpart = reversed(kpart)
            compart = []
            for x in kpart:
                if len(compart) == 0 and x.pos_ in ["PUNCT", "SCONJ", "CCONJ"]:
                    continue

                compart.append(x)
            compart = list(reversed(compart))

            if verbose:
                addinfo = ()
                if len(compart) > 0:
                    addinfo = (compart[0].idx, doc.text_with_ws[compart[0].idx:])
                print(scIdx, compart, len(compart), *addinfo)

            if len(compart) > 0:
                splitPositions.append(compart[0].idx)
                finalCompartments.append(compart)
            else:
                if verbose:
                    print("empty compartment", origcompart)

        return finalCompartments, splitPositions

    def checkCompartments(self, doc, mirword, geneword, verbose=False):

        compartments, splitPositions = self._deepCheck(doc, mirword, geneword, verbose)

        compCheck = False
        for comp in compartments:
            if mirword in comp and geneword in comp:
                compCheck = True

        return compCheck

    def getCompartment(self, doc, mirword, geneword, verbose=False):

        compartments, splitPositions = self._deepCheck(doc, mirword, geneword, verbose)

        compCheck = False
        for comp in compartments:
            if mirword in comp and geneword in comp:
                return comp

        return None


    def checkSurContext(self, doc, mirword, geneword, verbose=False):

        """
        BLA-mediated blubb
        BLA-signaling pathway
        """

        def checkPostContext(doc, word, contexts):

            for surContext, preview in contexts:

                geneEnd = word.idx+1
                maxEnd = min(len(doc.text_with_ws), geneEnd + len(word) + preview)
                contextStr = doc.text_with_ws[geneEnd:maxEnd]

                mRes = surContext.search(contextStr)

                #print(contextStr, surContext)

                if mRes != None:
                    
                    if verbose:
                        print("Context Fail", mRes, surContext)
                    return False

            return True

        if not checkPostContext(doc, geneword, self.surGeneContexts):
            return False

        if not checkPostContext(doc, mirword, self.surMiRContexts):
            return False


        for surContext, preview in self.surPreGeneContexts:

            geneEnd = geneword.idx
            minStart = max(0, geneEnd-1-preview)
            contextStr = doc.text_with_ws[minStart:geneEnd]

            mRes = surContext.search(contextStr)


            if mRes != None:
                if verbose:
                    print("Context Fail pre", mRes, surContext)
                return False

        for surContext, preview in self.surPreMiRContexts:

            mirEnd = mirword.idx
            minStart = max(0, mirEnd-1-preview)
            contextStr = doc.text_with_ws[minStart:mirEnd]

            mRes = surContext.search(contextStr)

            if mRes != None:
                if verbose:
                    print("Context Fail pre", mRes, surContext)
                return False

        """
        BLA, BLI and BLUBB cells
        """
        commonEdges = self.__getConjuncts(doc, verbose)

        for cidx, cname in enumerate([x for x in commonEdges]):

            edges = sorted(commonEdges[cname], key=lambda x: x.i)
            edgeStr = [str(x).lower() for x in edges]

            if geneword in edges:

                for eidx, edge in enumerate(edges):
                    if str(edge) in ["cells",]:

                        if eidx > 0 and str(edges[eidx-1]) in ["tested", "arrested"]:
                            if verbose:
                                print("Context Fail for cells ignored",edges[eidx-1])
                            continue
                        else:
                            if verbose:
                                print("Context Fail for cells")

                            return False

        """
        GENE-miR-126-GENE signaling circuit
        """
        if mirword.idx == geneword.idx:
            if verbose:
                print("Context Fail for same entity")
            return False

        """
        miR-659 binding-site of GRN
        """

        for t in doc:

            if t.pos_ in ["NOUN", "PROPN"]:

                allLefts = [x for x in t.lefts]
                allRights = [x for x in t.rights]

                if (mirword in allLefts and geneword in allRights) or (mirword in allRights and geneword in allLefts):
                    if verbose:
                        print("Context Fail for Noun Rule")
                    #return False
    
        """
        all gene tokes
        """

        allGeneTokens = []
        if geneword.dep_ in ["nummod"]:

            allGeneTokens.append(geneword)
            allGeneTokens.append(geneword.head)

            for child in geneword.head.children:
                if child.dep_ in ["nummod", "compound"] and not child in allGeneTokens:
                    allGeneTokens.append(child)

            allGeneTokens = sorted(allGeneTokens, key=lambda x: x.i)

            #print("all gene tokens", allGeneTokens)

            for gt in allGeneTokens:

                for gtc in gt.children:

                    if gtc.dep_ in ["conj"]:
                        
                        gtcStr = str(gtc)

                        for surContext, preview in self.surGeneContexts:
                            mRes = surContext.search(gtcStr)

                            if mRes != None:
                                if verbose:
                                    print("Context Fail GTC", mRes, surContext)
                                return False





        return True


    def checkRelation(self, doc, mirword, geneword, verbose=False, relClassifier=None):

        singleResults = {}

        """
        Create single results        
        """

        conjResult, conjs = self.checkCommonConj(doc, mirword, geneword, verbose)
        singleResults["conj"] = conjResult

        sdpPass, passive, negated = self.checkSDP(doc, mirword, geneword, verbose)
        singleResults["sdp"] = sdpPass

        compPass = self.checkCompartments(doc, mirword, geneword, verbose)
        singleResults["compartment"] = compPass

        sigPathway = self.checkSurContext(doc, mirword, geneword, verbose)
        singleResults["context"] = sigPathway

        """
        Create overall result        
        """
        trues = 0
        falses = 0
        for x in singleResults:
            if singleResults[x]:
                trues += 1
            else:
                falses += 1

        #if trues > falses:
        if falses == 0:
            overallResult = True
        else:
            overallResult = False

        """
        ADD INFORMATION
        """
        singleResults["passive"] = passive
        singleResults["negation"] = singleResults.get("negation", False) or negated

        """
        CLASSIFY RELATION
        """
        
        if relClassifier != None:
            singleResults["classification"] = relClassifier(doc, mirword, geneword, self, verbose=verbose)


        return overallResult, singleResults

    def checkPassive(self, doc, mirword, geneword, verbose):

        lowerI = mirword.i if mirword.i < geneword.i else geneword.i
        upperI = mirword.i if mirword.i > geneword.i else geneword.i

        compartments, splitPositions = self._deepCheck(doc, mirword, geneword, verbose)

        passiveSeen = False

        for comp in compartments:

            for tkn in comp:
                if "pass" in tkn.dep_:
                    if verbose:
                        print(tkn, tkn.dep_, tkn.pos_)
                    passiveSeen = True

        return passiveSeen and ("subj" in geneword.dep_ or "subj" in mirword.dep_)

        

    def checkSDP(self, doc, ent1, ent2, verbose):

        if ent1.i < ent2.i:
            sdpRes = self._get_sdp_path(doc, ent1.i, ent2.i)
        else:
            sdpRes = self._get_sdp_path(doc, ent2.i, ent1.i)

        sdpPass = True
        passive = False
        negated = False
        subjCount = 0
        verbCount = 0

        if verbose:
            print("SDP sent:", doc.text_with_ws)
            print("SDP     :", sdpRes)

        for sdp in sdpRes:
            """
            REASONING: VERB and AUX indicate a VERB, AUX="is"
            
            """

            if sdp[1] in ["VERB", "AUX"] and not sdp[2] in ["acl"]:
                verbCount += 1

            if sdp[2] in ['nsubjpass']:
                passive=True

            if sdp[2] in ['nsubj', 'nsubjpass']:
                subjCount += 1

            if sdp[2] in ['neg']:
                negated = True

        for idx in range(0, len(sdpRes)-1):

            if sdpRes[idx][1] in ["VERB"] and sdpRes[idx+1][1] in ["NOUN"] and str(sdpRes[idx+1][0]) in ["pathway", "pathways"]:
                if verbose:
                    print("SDP VERB NOUN PATHWAY")
                sdpPass = False

        if len(sdpRes) >= 2:

            #[(Shh, 'PROPN', 'compound', False), (N-myc, 'PROPN', 'nsubj', False), (induced, 'VERB', 'ROOT', False), (expression, 'NOUN', 'dobj', False), (miR-17/92, 'PROPN', 'compound', False)]
            #In CGNPs, the <e2>Shh</e2> effector N-myc, but not Gli1, induced <e1>miR-17/92</e1> expression.
            #should be false!

            targetWord = sdpRes[0]
            nWord = sdpRes[1]

            if targetWord[1] == "PROPN" and not targetWord[2] in ["conj"] and nWord[1] == "PROPN" and "nsubj" in nWord[2]:
                if verbose:
                    print("SDP PROPN SUBJ rule")
                sdpPass = False


        #if subjCount != 1:
        #    sdpPass = False

        #if verbCount == 0 and len(sdp) > 3:
        #    sdpPass = False

        if len(sdpRes) == 1:
            sdpPass = False

        return sdpPass, passive, negated

    def __followChildSel(self, tk, deps, fwDeps, bwDeps, verbose=False):

        tks = set()
        tks.add(tk)
        for c in tk.children:
            if c.dep_ in deps:
                if verbose:
                    print("fc add children for", tk)
                tks  = tks.union(self.__followChildSel(c, deps, fwDeps, bwDeps))

            elif c.dep_ in fwDeps and c.idx > tk.idx:
                if verbose:
                    print("fc add children fw for", tk)
                tks  = tks.union(self.__followChildSel(c, deps, fwDeps, bwDeps))
            elif c.dep_ in bwDeps and c.idx < tk.idx:
                if verbose:
                    print("fc add children bw for", tk)
                tks  = tks.union(self.__followChildSel(c, deps, fwDeps, bwDeps))

        return tks

    def __followHeadSel(self, tk, deps, fwDeps, bwDeps, verbose=False):

        tks = set()
        tks.add(tk)
        for c in [tk.head]:
            if tk.dep_ in deps:
                if verbose:
                    print("fc add children for", tk)
                tks  = tks.union(self.__followHeadSel(c, deps, fwDeps, bwDeps))

            elif tk.dep_ in fwDeps and c.idx > tk.idx:
                if verbose:
                    print("fc add children fw for", tk)
                tks  = tks.union(self.__followHeadSel(c, deps, fwDeps, bwDeps))
            elif tk.dep_ in bwDeps and c.idx < tk.idx:
                if verbose:
                    print("fc add children bw for", tk)
                tks  = tks.union(self.__followHeadSel(c, deps, fwDeps, bwDeps))

        return tks

    def __followChild(self, tk, deps, verbose=False):

        tks = set()
        tks.add(tk)
        for c in tk.children:
            if c.dep_ in deps:
                if verbose:
                    print("fc add children for", tk)
                tks  = tks.union(self.__followChild(c, deps))

        return tks

    def __getConjuncts(self, doc, verbose):
        commonEdges = defaultdict(set)

        def addElemsToCommon(telems):
            added = False
            for x in commonEdges:
                if any([y in commonEdges[x] for y in telems]):
                    for y in telems:
                        commonEdges[x].add(y)

                        """
                        Association studies using all common variants detected in the 3' UTR of BACE1 and the miR-29 gene cluster did not identify an association with AD risk.

                        UTR conj cluster

                        also adds : 3' UTR of BACE1     the miR-29 gene
                        """
                        added = True

                    if added:
                        break

                if added:
                    break

            if not added:
                for y in telems:
                    commonEdges[t].add(y)

        for t in doc:
            #print(t.idx, t.text, t.dep_, t.pos_, t.head.text, t.conjuncts)

            if len(t.conjuncts) > 0:
                telems = list(t.conjuncts) + [t]

                if verbose:
                    print("Conjuncts", telems)

                idx2t = {e.i: e for e in doc}
                for e in t.conjuncts:
                    n = idx2t.get(e.i+1, None)

                    if verbose:
                        print("neighbour", e, n)

                    if n.pos_ in ["PUNCT"]:
                        n = idx2t.get(e.i+2, None)

                    if verbose:
                        print("neighbour", e, n)

                    if n != None:
                        if n.head == e.head and n.dep_ in ["dep"]:
                            telems.append(n)


                if verbose:
                    print("Conjunction before fc", t)

                #childTokens = self.__followChild(t, ["case", "compound", "amod", "nmod", "dep", "appos", "acl", "dobj", "nummod"], verbose=verbose)
                childTokens = self.__followChildSel(t, ["case", "amod", "nmod", "dep", "appos", "acl", "dobj", "nummod"], ["compound"], [], verbose=verbose)
                
                for c in childTokens:
                    telems.append(c)

                    #if childDep in ["compound", "amod", "nmod", "appos"]:
                    #    telems.append(child)

                addElemsToCommon(telems)
                continue

            if t.dep_ in ["appos"]:
                telems = [t.left_edge, t.right_edge]
                addElemsToCommon(telems)

        return commonEdges

    def checkCommonConj(self, doc, mirword, geneword, verbose):

        commonEdges = self.__getConjuncts(doc, verbose)

        # this is meant to add description of the conj as well.
        #commonEdges[cname] = commonEdges[cname].union(subtreeEdges)

        if verbose:
            print("Conjunctions")
            for cname in commonEdges:
                celes = sorted(commonEdges[cname])
                print(cname, celes)

        for cname in commonEdges:

            celes = sorted(commonEdges[cname], key=lambda x: x.i)
            if mirword in celes and geneword in celes:
                
                rems = list(set(celes).difference(set([mirword, geneword])))
                rems = sorted(rems, key=lambda x: x.idx)

                if len(celes) > 2 and str(rems[0]) in ["between"]:
                    if verbose:
                        print("Conjunction case rule", rems)
                    continue

                mirPos = celes.index(mirword)
                if mirPos != None:
                    oElems = celes[max(mirPos-2, 0): mirPos]

                    for oElem in oElems:
                        if str(oElem) in ["target"]:
                            if verbose:
                                print("Conjunction target rule", celes)
                            continue

                aChildren = self.__followHeadSel(mirword, ["appos", "amod"], [], [])
                
                if verbose:
                    print("achildren", aChildren)
                if geneword in aChildren:
                    continue

                return False, celes

            elif geneword in celes:

                for t in celes[-2:]:
                    if str(t).lower() in ["pathways", "pathway"]:
                        return False, celes


            

        return True, None

    def _get_sdp_path(self, doc, subj, obj):
        lca_matrix = doc.get_lca_matrix()
        lca = lca_matrix[subj, obj]

        #    for x in [(t.idx, t.text, t.dep_, t.pos_, t.head.text) for t in doc]:

        current_node = doc[subj]
        subj_path = [(current_node, current_node.pos_, current_node.dep_, False)]
        if lca != -1:
            if lca != subj:
                while current_node.head.i != lca:

                    if current_node.i == current_node.head.i:
                        break

                    current_node = current_node.head
                    subj_path.append((current_node, current_node.pos_, current_node.dep_, subj_path[-1][0] in current_node.conjuncts))
                subj_path.append((current_node.head, current_node.head.pos_, current_node.head.dep_, subj_path[-1][0] in current_node.head.conjuncts))
                current_node = doc[obj]
        obj_path = [(current_node, current_node.pos_, current_node.dep_, False)]
        if lca != -1:
            if lca != obj:
                while current_node.head.i != lca:
                    if current_node.i == current_node.head.i:
                        break

                    current_node = current_node.head
                    obj_path.append((current_node, current_node.pos_, current_node.dep_, obj_path[-1][0] in current_node.conjuncts))
                obj_path.append((current_node.head, current_node.head.pos_, current_node.head.dep_, obj_path[-1][0] in current_node.head.conjuncts))

        intElem = (subj_path[-1][0], subj_path[-1][1], subj_path[-1][2], subj_path[-1][3] or obj_path[-1][3])

        retpath = subj_path[:-1] + [intElem] + obj_path[::-1][1:]
        return retpath