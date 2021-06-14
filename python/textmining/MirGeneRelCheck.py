from collections import defaultdict, Counter
import regex as re
import os,sys
import copy
from intervaltree import Interval, IntervalTree

#nlp = spacy.load("/mnt/d/spacy/models/en_ner_bionlp13cg_md-0.2.4/en_ner_bionlp13cg_md/en_ner_bionlp13cg_md-0.2.4")

class SentenceRelationClassifier:

    def __init__(self, relationCSV, active_checks = None):

        print("Setting relfile", relationCSV, file=sys.stderr)
        self.relfile = relationCSV

        allStems, allRels, stem2class = self.loadStems()
        self.allStems = allStems
        self.allRels = allRels
        self.stem2class = stem2class


        self.dir2opp = {"MIR_GENE": "GENE_MIR", "GENE_MIR": "MIR_GENE"}
        self.idir2opp = {"POS": "NEG", "NEG": "POS"}
        self.rdir2opp = {"DOWN": "UP", "UP": "DOWN"}

        self.major_checks = ["compartment", "between", "counts","return"]

        self.checks = [
            'compartment gene mir pos corr',
            'compartment gene mir neg corr',
            'compartment gene mir',
            'compartment mir gene pos corr',
            'compartment mir gene neg corr',
            'compartment mir gene',

            'between gene mir pos corr',
            'between gene mir inv corr', 
            'between gene mir neg corr',
            'between gene mir reg corr',
            'between gene mir sites',

            "between gene mir",
            "counts between gene mir",

            'between mir gene pos corr',
            'between mir gene neg corr',
            'between mir gene reg corr',
            'between mir gene',
            'between mir gene rb',
            'between mir gene nrb',

            'counts opp',
            'counts between equal',
            'counts between alternating',
            'counts between',

            "return mediated",
            "return m g by",
            "return"
            
        ]

        self.active_checks = None

        if active_checks is None:
            self.active_checks = [x for x in self.checks]+[x for x in self.major_checks]
        else:
            self.active_checks = [x for x in active_checks if x in self.checks or x in self.major_checks]

    def loadStems(self):
    
        #print("Loading stems")
        allRels = []
        stem2class = {}

        renameClass = {
            'INT': 'NEU',
            "BOOST": "POS",
            "COMBINE": "NEU",
            "CHANGE": "NEU",
        }

        if not os.path.exists(self.relfile):
            print("Error Loading allrels file!", self.relfile, sys.stderr)
            exit(-1)

        with open(self.relfile) as fin:

            for line in fin:
                aline = line.strip().split(",", 1)
                aregSyns = aline[1].split(":")
                regSyns = []
                if len(aregSyns) > 1:
                    regSyns = aregSyns[1].split("|")
                
                aline[1] = aline[1].split(":")[0]           

                regClass = aline[0]
                regStem = aline[1]

                if regClass in renameClass:
                    regClass = renameClass[regClass]
                    
                stem2class[regStem] = regClass
                    
                for regSyn in regSyns:
                    if not regSyn in stem2class:
                        stem2class[regSyn] = regClass


            allRels.append((regClass, regStem))

        allStems = [x for x in stem2class]
    
        return allStems, allRels, stem2class


    def getTokenDist(self, tk1, tk2):
        return abs(tk1.i-tk2.i)

    def countStems(self, tokens, geneword, mirword, keepNEU=False):
        stemTypes2Count = Counter()
        mirgene2stem = defaultdict(lambda: Counter())
        stemEvs = []

        for tidx,x in enumerate(tokens):
            tknStem, matchStem, wasFound = self.getstem(x)

            if wasFound and (tknStem != "NEU" or keepNEU):
                
                if tidx >= len(tokens)-2 and str(tokens[tidx-1]) == "to":
                    continue
                
                geneDist = self.getTokenDist(x, geneword)
                mirDist = self.getTokenDist(x, mirword)
                            
                if geneDist < mirDist and geneDist < 3:
                    mirgene2stem["gene"][tknStem] += 1
                elif mirDist < geneDist and mirDist < 3:
                    mirgene2stem["mir"][tknStem] += 1

            if wasFound and (tknStem != "NEU" or keepNEU):
                stemTypes2Count[tknStem] += 1
                stemEvs.append((x, tknStem))

        if len(stemTypes2Count) == 0:
            if "regulated" in str(geneword):
                stemTypes2Count["NEU"] += 1

        return stemTypes2Count,mirgene2stem,stemEvs

    def getRelsForTokens(self, tokens, desc,verbose):
        rels = []
        verbSet = set(["VERB"])
        for x in tokens:
            
            tokenStemClass, s, _ = self.getstem(x)
            rels.append((s, None, tokenStemClass, x.pos_, x))
            if verbose:
                print(desc, x, s, tokenStemClass, x.pos_, x)
                        
        hasVerb = len(list(set([x[3] for x in rels]) & verbSet)) > 0
        
        if hasVerb:
            rels = [x for x in rels if x[3] in verbSet]
        
        return rels

    def evalCounts(self, dirReg, mirPos, mirNeg, genePos, geneNeg):

        invcor = False
        regdir = None

        #if dirReg == 0:
        if (mirPos > 0 and geneNeg > 0) or (genePos > 0 and mirNeg > 0):
            return "DOWN", True

        #else:
        if ((mirPos > 0 and genePos > 0) and(mirNeg == 0 and geneNeg == 0)) or ((mirNeg > 0 and geneNeg >0) and(mirPos == 0 and genePos == 0)):
            return "UP", False

        return None, None

    def scoreWord(self, tkn):
        
        jtkstem, s, _ = self.getstem(tkn)
                            
        if jtkstem == None:

            if jtk.startswith(("loss", "decreas", "low")):
                return 0,1

            elif jtk.startswith(("increas", "high", "posit")):
                return 1,0

        else:
            if jtkstem == "POS":
                return 1,0
            elif jtkstem == "NEG":
                return 0,1
            
        return 0,0


    def getTextBetween(self, tks, lword, rword, dist=3):
        tokenBetween = [x for x in tks if lword.i-dist < x.i < rword.i+dist]
        
        lit = lword
        while lit.pos_ in ["nmod", "conj"]:
            if not lit in tokenBetween:
                tokenBetween.append(lit)
                
            lit = lit.head
            
        rit = rword
        while rit.pos_ in ["nmod", "conj"]:
            if not rit in tokenBetween:
                tokenBetween.append(rit)
                
            rit = rit.head


        tokenBetween = list(sorted(tokenBetween, key=lambda x: x.i))
        tokenBetweenStrict = list(sorted([x for x in tks if lword.i < x.i < rword.i], key=lambda x: x.i))
        textBetween = ' '.join(" ".join([str(x) for x in tokenBetween]).split())
        textAfter = ' '.join(" ".join([str(x) for x in [y for y in rword.doc if rword.i <= y.i < rword.i+5]]).split())
        textBefore = ' '.join(" ".join([str(x) for x in [y for y in rword.doc if lword.i-5 <= y.i < lword.i]]).split())

        return tokenBetween, tokenBetweenStrict, textBetween, textBefore, textAfter

    def getstem(self, tkn):
        tknstr = str(tkn).lower()
        for sidx, s in enumerate(self.stem2class):
            if tknstr.startswith(s):
                
                if tknstr in ["enhancer","inhibitor"]:
                    continue

                tokenStemClass = self.stem2class[s]
                
                if tokenStemClass == "NEU":
                    for i in range(1,2):

                        if tkn.i-i<0:
                            continue

                        #tknPre = str(tkn.doc[tkn.i-i])
                        #if tknPre.startswith(("high", "elevat", "increas")):
                        #    tokenStemClass = "POS"

                        #if tknPre.startswith(("low", "decreas")):
                        #    tokenStemClass = "NEG"
                                                    
                return tokenStemClass, s, True
            
            elif tknstr.endswith(("-stressed")):
                return "POS", s, True
            
        return "NEU", str(tkn), False

    def allowed_check(self, checkname, majorname):

        assert( checkname in self.checks )

        if majorname != None:
            #if not majorname in self.major_checks:
            #    print(majorname)

            assert( majorname in self.major_checks)

            return majorname in self.active_checks

        return checkname in self.active_checks


    def classify(self, doc, mirword, geneword, relchecker, verbose):
        #isPassive = relchecker.checkPassive(doc, mirword, geneword, verbose=False)
        _, isPassive, _ = relchecker.checkSDP(doc, mirword, geneword, verbose=False)

        compartment = relchecker.getCompartment(doc, mirword, geneword, verbose=False)
        
        if verbose:
            print("passive?", isPassive)
            
            
        sdpRes = relchecker._get_sdp_path(doc, mirword.i, geneword.i)
            
        stemTypes2Count = Counter()
        
        for x in sdpRes:
            tknStem, matchStem, _ = self.getstem(x[0])
            
            if tknStem != "NEU":
                stemTypes2Count[tknStem] += 1

        if verbose:        
            print(stemTypes2Count)
            
        if len(stemTypes2Count) > 0:
            regStem = stemTypes2Count.most_common(1)[0][0]
            
            if regStem == "NEG":
                regStem = "DOWN"
            elif regStem == "POS":
                regStem = "UP"
                
            if stemTypes2Count["POS"] == stemTypes2Count["NEG"]:
                regStem = "DOWN"
                
            idir = "MIR_GENE"
            if isPassive and mirword.i < geneword.i:
                idir = "GENE_MIR"
                
            #return {"regulation_dir": regStem, "interaction_dir": idir, "passive": isPassive, "reg_detect": "sdp counts"}
            
        targetWords = ["associated with","regulator of", " regulation of", "direct regulator of", "direct regulator for", "binding to", "binds to", "bind to", "bound to", "interacting with","targeting of"]
        targetWordsBy = [ "regulation by", "regulated by", "targets of", "target of", "effect of", "effects of", "effect by", "effects by", "recognized by"]
            
        if compartment != None:
            tokenBetween = compartment
            textBetween = ' '.join(" ".join([str(x) for x in tokenBetween]).split())

            if verbose:
                print("Compartment Text", textBetween)
            
            if geneword.i < mirword.i:

                textLeft = ' '.join(" ".join([str(x) for x in [x for x in compartment if x.i < geneword.i]]).split())
                textBetweenStrict = ' '.join(" ".join([str(x) for x in [x for x in compartment if geneword.i < x.i < mirword.i]]).split())
                textRight = ' '.join(" ".join([str(x) for x in [x for x in compartment if mirword.i < x.i]]).split())
                
                if self.allowed_check("compartment gene mir pos corr", "compartment") and any([x in textBetween for x in ["positive correlation"]]):
                    return {"regulation_dir": "UP", "interaction_dir": "GENE_MIR", "passive": isPassive, "reg_detect": "compartment gene mir pos corr", "reg_detect_major": "compartment"}
                
                if self.allowed_check("compartment gene mir neg corr", "compartment") and any([x in textBetween for x in ["inverse correlation", "inversely correlated", "inversely regulated"]]):
                    return {"regulation_dir": "DOWN", "interaction_dir": "GENE_MIR", "passive": isPassive, "reg_detect": "compartment gene mir neg corr", "reg_detect_major": "compartment"}
                
                if self.allowed_check("compartment gene mir", "compartment"):
                    regDetect = "compartment gene mir"

                    for x in targetWords:
                        if x in textBetween and not x in textRight:
                            if x in textLeft and not x in textBetweenStrict and not textRight.startswith("promoter"):
                                return {"regulation_dir": "NEU", "interaction_dir": "MIR_GENE", "passive": isPassive, "reg_detect": regDetect, "reg_detect_major": "compartment"}
                            else:
                                return {"regulation_dir": "NEU", "interaction_dir": "GENE_MIR", "passive": isPassive, "reg_detect": regDetect, "reg_detect_major": "compartment"}

                            
                if self.allowed_check("compartment gene mir", "compartment"):
                    regDetect = "compartment gene mir by"

                    for x in targetWordsBy:
                        if x in textBetween and not x in textRight:
                            if x in textLeft and not x in textBetweenStrict:                               
                                return {"regulation_dir": "NEU", "interaction_dir": "GENE_MIR", "passive": isPassive, "reg_detect": regDetect, "reg_detect_major": "compartment"}   
                            else:
                                return {"regulation_dir": "NEU", "interaction_dir": "MIR_GENE", "passive": isPassive, "reg_detect": regDetect, "reg_detect_major": "compartment"}   
                
            elif mirword.i < geneword.i:

                textLeft = ' '.join(" ".join([str(x) for x in [x for x in compartment if x.i < mirword.i]]).split())
                textBetweenStrict = ' '.join(" ".join([str(x) for x in [x for x in compartment if mirword.i < x.i < geneword.i]]).split())
                textRight = ' '.join(" ".join([str(x) for x in [x for x in compartment if geneword.i < x.i]]).split())
                
                if self.allowed_check("compartment mir gene pos corr", "compartment") and any([x in textBetween for x in ["positive correlation"]]):
                    return {"regulation_dir": "UP", "interaction_dir": "MIR_GENE", "passive": isPassive, "reg_detect": "compartment mir gene pos corr", "reg_detect_major": "compartment"}
                
                if self.allowed_check("compartment mir gene neg corr", "compartment") and any([x in textBetween for x in ["inverse correlation", "inversely correlated", "inversely regulated"]]):
                    return {"regulation_dir": "DOWN", "interaction_dir": "MIR_GENE", "passive": isPassive, "reg_detect": "compartment mir gene neg corr", "reg_detect_major": "compartment"}
                
                if self.allowed_check("compartment mir gene", "compartment"):
                    regDetect = "compartment mir gene"
                    #"regulator of"
                    for x in targetWords:
                        if x in textBetween and not x in textRight:
                            if x in textLeft and not x in textBetweenStrict:
                                return {"regulation_dir": "NEU", "interaction_dir": "GENE_MIR", "passive": isPassive, "reg_detect": regDetect, "reg_detect_major": "compartment"}
                            else:
                                return {"regulation_dir": "NEU", "interaction_dir": "MIR_GENE", "passive": isPassive, "reg_detect": regDetect, "reg_detect_major": "compartment"}

                            
                if self.allowed_check("compartment mir gene", "compartment"):
                    regDetect = "compartment mir gene by"
                    #"regulation by"
                    for x in targetWordsBy:
                        if x in textBetween and not x in textRight:
                            if x in textLeft and not x in textBetweenStrict:                               
                                return {"regulation_dir": "NEU", "interaction_dir": "MIR_GENE", "passive": isPassive, "reg_detect": regDetect, "reg_detect_major": "compartment"}   
                            else:
                                return {"regulation_dir": "NEU", "interaction_dir": "GENE_MIR", "passive": isPassive, "reg_detect": regDetect, "reg_detect_major": "compartment"}   

           
        if compartment != None:
            tokSrc = compartment
            tokenDist = len(doc)
        else:
            tokSrc = doc
            tokenDist = 3
        
        if geneword.i < mirword.i:

            tokenBetween, tokenBetweenStrict, textBetween, textBefore, textAfter = self.getTextBetween(tokSrc, geneword, mirword, dist=tokenDist)

            textLeft = ' '.join(" ".join([str(x) for x in [x for x in tokSrc if x.i < geneword.i]]).split())
            textBetweenStrict = ' '.join(" ".join([str(x) for x in [x for x in tokSrc if geneword.i < x.i < mirword.i]]).split())
            textRight = ' '.join(" ".join([str(x) for x in [x for x in tokSrc if mirword.i < x.i]]).split())

            if verbose: 
                print(tokSrc)  
                print("TextBetween:", textBetween)
            

            if self.allowed_check("between gene mir pos corr", "between") and any([x in textBetween for x in ["positive correlat", "positively correlat"]]):
                return {"regulation_dir": "UP", "interaction_dir": "GENE_MIR", "passive": isPassive, "reg_detect": "between gene mir pos corr", "reg_detect_major": "between"}
            
            if self.allowed_check("between gene mir inv corr", "between") and any([x in textBetween for x in ["inversely correlated to", "inversely regulated to", "inversely related to", "inverse relation to", "inverse relationship to", "negative correlat", "negatively correlat"]]):
                return {"regulation_dir": "DOWN", "interaction_dir": "MIR_GENE", "passive": isPassive, "reg_detect": "between gene mir neg corr", "reg_detect_major": "between"}
            
            if self.allowed_check("between gene mir neg corr", "between") and any([x in textBetween for x in ["inverse correlation", "inversely correlated", "inversely regulated", "inverse relation", "inversely related"]]):
                return {"regulation_dir": "DOWN", "interaction_dir": "GENE_MIR", "passive": isPassive, "reg_detect": "between gene mir neg corr", "reg_detect_major": "between"}
            

            if self.allowed_check("between gene mir reg corr", "between"):
                regDetect = "between gene mir reg corr"
                #"regulation by"
                for x in ["negative regulation of"]:
                    if x in textBetween and not x in textRight:
                        if x in textLeft and not x in textBetweenStrict:                               
                            return {"regulation_dir": "DOWN", "interaction_dir": "MIR_GENE", "passive": isPassive, "reg_detect": regDetect, "reg_detect_major": "between"}
                        else:
                            return {"regulation_dir": "DOWN", "interaction_dir": "GENE_MIR", "passive": isPassive, "reg_detect": regDetect, "reg_detect_major": "between"}


            if self.allowed_check("between gene mir sites", "between") and any([x in textBefore for x in ["binding site", "binding sites", "binding to", "binding directly to"]]):
                return {"regulation_dir": "NEU", "interaction_dir": "MIR_GENE", "passive": isPassive, "reg_detect": "between gene mir sites", "reg_detect_major": "between"}
            
            
    
            if self.allowed_check("between gene mir", "between") and any([x in textBetweenStrict for x in targetWords]) or any([x in textBetweenStrict for x in ["target", "targeted by", "binding efficiency", "binding site for", "recognizing"]]):
                return {"regulation_dir": "NEU", "interaction_dir": "MIR_GENE", "passive": isPassive, "reg_detect": "between gene mir", "reg_detect_major": "between"}
            

            #stemTypes2Count = Counter()
            #for x in tokenBetween:
            #    tknStem, matchStem = getstem(x)
            #    if tknStem != "NEU":
            #        stemTypes2Count[tknStem] += 1
                    
            if False and self.allowed_check("counts between gene mir", "counts") and len(stemTypes2Count) > 0:
                regStem = stemTypes2Count.most_common(1)[0][0]
                if regStem == "NEG":
                    regStem = "DOWN"
                elif regStem == "POS":
                    regStem = "UP"
                    
                if stemTypes2Count["POS"] == stemTypes2Count["NEG"]:
                    regStem = "DOWN"

                idir = "MIR_GENE"
                if isPassive and geneword.i < mirword.i:
                    idir = "GENE_MIR"

                return {"regulation_dir": regStem, "interaction_dir": idir, "passive": isPassive, "reg_detect": "counts between gene mir", "reg_detect_major": "counts"}
            
        else:
            #mirword.i < geneword.i
            tokenBetween, tokenBetweenStrict, textBetween, textBefore, textAfter = self.getTextBetween(tokSrc, mirword, geneword, dist=tokenDist)

            textLeft = ' '.join(" ".join([str(x) for x in [x for x in tokSrc if x.i < mirword.i]]).split())
            textBetweenStrict = ' '.join(" ".join([str(x) for x in [x for x in tokSrc if mirword.i < x.i < geneword.i]]).split())
            textRight = ' '.join(" ".join([str(x) for x in [x for x in tokSrc if geneword.i < x.i]]).split())

            if verbose: 
                print(tokSrc)  
                print("TextBetween:", textBetween)
            
            if self.allowed_check("between mir gene pos corr", "between") and any([x in textBetween for x in ["positive correlation"]]) and stemTypes2Count.get("NEG", 0) == 0:
                return {"regulation_dir": "UP", "interaction_dir": "MIR_GENE", "passive": isPassive, "reg_detect": "between mir gene pos corr", "reg_detect_major": "between"}
            
            if self.allowed_check("between mir gene reg corr", "between") and  any([x in textBetween for x in ["negative regulation of", "inhibition of"]]):

                if isPassive:
                    return {"regulation_dir": "DOWN", "interaction_dir": "GENE_MIR", "passive": isPassive, "reg_detect": "between mir gene reg corr by", "reg_detect_major": "between"}
                else:
                    return {"regulation_dir": "DOWN", "interaction_dir": "MIR_GENE", "passive": isPassive, "reg_detect": "between mir gene reg corr", "reg_detect_major": "between"}
            
            if self.allowed_check("between mir gene neg corr", "between") and  any([x in textBetween for x in ["inverse correlation", "inversely correlated", "inversely regulated", "inverse relation", "inversely related"]]):
                return {"regulation_dir": "DOWN", "interaction_dir": "MIR_GENE", "passive": isPassive, "reg_detect": "between mir gene neg corr", "reg_detect_major": "between"}
            
            if self.allowed_check("between mir gene", "between") and  any([x in textBetween for x in ["directly inhibits", "directly represses", "translational inhibition" ]]):
                return {"regulation_dir": "DOWN", "interaction_dir": "MIR_GENE", "passive": isPassive, "reg_detect": "between mir gene", "reg_detect_major": "between"}
            
            if self.allowed_check("between mir gene", "between") and  any([x in textBetweenStrict for x in targetWords + ["by targeting", "by regulating", "directly regulated", "its target gene", "binding site"]]): #, "regulation of"
                return {"regulation_dir": "NEU", "interaction_dir": "MIR_GENE", "passive": isPassive, "reg_detect": "between mir gene", "reg_detect_major": "between"}
            
            if self.allowed_check("between mir gene rb", "between") and  any([x in textBetweenStrict for x in ["regulated by", "targeted by"]]):
                return {"regulation_dir": "NEU", "interaction_dir": "GENE_MIR", "passive": isPassive, "reg_detect": "between mir gene rb", "reg_detect_major": "between"}

            if self.allowed_check("between mir gene nrb", "between") and  any([x in textBetween for x in ["negatively regulated by"]]):
                return {"regulation_dir": "DOWN", "interaction_dir": "GENE_MIR", "passive": isPassive, "reg_detect": "between mir gene nrb", "reg_detect_major": "between"}
            
        tBtwnRes = None
            
        if len(tokenBetween) > 0:
            if verbose:
                print("infer direction from token set", tokenBetween)
            tBtwnRes = self.inferDirectionFromTokenSet(doc, tokenBetween, geneword, mirword, isPassive, verbose)

        if tBtwnRes != None:
            if verbose:
                print("Return After Infer")
                print(tBtwnRes)
            return tBtwnRes
        
            


        """
        FIRST GUESS
        """

        madeGuess = False
        regDir = "MIR_GENE"
        regStem = "NEU"

        if mirword.i < geneword.i:
            
            if str(geneword).endswith( ("-regulating") ):
                regDir = "MIR_GENE"
                regStem = "NEU"
                madeGuess = True

            elif str(geneword).endswith( ("dependent", "independent") ):
                regDir = "GENE_MIR"
                regStem = "NEU"
                madeGuess = True

            elif (len(doc) > geneword.i+1) and any([str(doc[geneword.i+j]) in ["target"] for j in [1]]):
                regDir = "MIR_GENE"
                regStem = "NEU"
                madeGuess = True


        elif geneword.i < mirword.i:

            if str(mirword).endswith( ("-regulating") ):
                regDir = "GENE_MIR"
                regStem = "NEU"
                madeGuess = True

            elif str(mirword).endswith( ("dependent", "independent") ):
                regDir = "MIR_GENE"
                regStem = "NEU"
                madeGuess = True

            elif str(geneword).endswith(("-regulating", "-targeting")):
                regStem = "NEU"
                regDir = "MIR_GENE"
                madeGuess = True

            elif str(geneword).endswith(("dependent", "induced")):
                regDir = "GENE_MIR"
                regStem = "NEU"
                madeGuess = True

            elif (len(doc) > mirword.i+1) and any([str(doc[mirword.i+j]) in ["target", "targets"] for j in [1]]):
                regDir = "MIR_GENE"
                regStem = "NEU"
                madeGuess = True


        if madeGuess:

            if self.allowed_check("return", "return"):
                return {"regulation_dir": regStem, "interaction_dir": regDir, "passive": isPassive, "reg_detect": "return", "reg_detect_major": "return"}

        else:
        
            regeluateWords = ["mediated",]
            
            if mirword.i < geneword.i and "mediated" in str(mirword):
                regStem = "DOWN"
                
            if self.allowed_check("return mediated", "return") and  geneword.i < mirword.i and ("mediated" in str(mirword) or any([str(doc[mirword.i-j]) in ["mediated",] for j in [1]])):
                regStem = "NEU"
                regDir = "MIR_GENE"
                
                return {"regulation_dir": regStem, "interaction_dir": regDir, "passive": isPassive, "reg_detect": "return mediated", "reg_detect_major": "return"}
                
            if self.allowed_check("return m g by", "return") and  mirword.i<geneword.i and any([str(doc[geneword.i-j]) in ["by", "with"] for j in [1,2,3]]):
                regStem = "NEU"
                regDir = "MIR_GENE"
                
                return {"regulation_dir": regStem, "interaction_dir": regDir, "passive": isPassive, "reg_detect": "return m g by", "reg_detect_major": "return"}
                
            if geneword.i < mirword.i:
                regStem = "NEU"
                regDir = "GENE_MIR"
            
            if geneword.i < mirword.i and any([str(doc[mirword.i-j]) in ["by", "with"] for j in [1,2,3]]):
                regStem = "NEU"
                regDir = "MIR_GENE"
                
            if mirword.i < geneword.i and str(doc[mirword.i+1]) == "target":
                regStem = "NEU"
                regDir = "MIR_GENE"
        
        #if mirword.i < geneword.i:
        #    regDir = "MIR_GENE"
        #else:
        #    regDir = "GENE_MIR"
        
        #if isPassive:
        #    if mirword.i < geneword.i:
        #        regDir = "GENE_MIR"
        #    else:
        #        regDir = "MIR_GENE"
            
        if self.allowed_check("return", "return"):
            return {"regulation_dir": regStem, "interaction_dir": regDir, "passive": isPassive, "reg_detect": "return", "reg_detect_major": "return"}

        else:
            return {"regulation_dir": "DOWN", "interaction_dir": "MIR_GENE", "passive": isPassive, "reg_detect": "static", "reg_detect_major": "static"}


    def inferDirectionFromTokenSet(self, doc, inTokens, geneword, mirword, isPassive, verbose):

        if geneword.i < mirword.i:
            tokenBetween, tokenBetweenStrict, textBetween, textBefore, textAfter = self.getTextBetween(inTokens, geneword, mirword)
        else: #mirword.i < geneword.i
            tokenBetween, tokenBetweenStrict, textBetween, textBefore, textAfter = self.getTextBetween(inTokens, mirword, geneword)

        stemTypes2Count,mirgene2stem, tknEvs = self.countStems(inTokens, geneword, mirword)

        if verbose:
            print("stemTypes2Count:", stemTypes2Count)
            print("mirgene2stem:", mirgene2stem)
                        
        if len(stemTypes2Count) > 0:
            regStem = stemTypes2Count.most_common(1)[0][0]
            
            if self.allowed_check("counts opp", "counts") and  len(stemTypes2Count) >= 2 and len(mirgene2stem["gene"]) >0 and len(mirgene2stem["mir"])>0 and mirgene2stem["gene"].most_common(1)[0][0] == self.idir2opp[mirgene2stem["mir"].most_common(1)[0][0]]:
                if verbose:
                    print(mirgene2stem)
                if geneword.i < mirword.i:
                    if mirgene2stem["gene"].most_common(1)[0][0] == "POS":
                        regDir = "MIR_GENE"
                        regStem = "DOWN"
                    else:
                        regDir = "MIR_GENE"
                        regStem = "DOWN"
                else:
                                            
                    if mirgene2stem["mir"].most_common(1)[0][0] == "NEG":
                        regDir = "MIR_GENE"
                        regStem = "DOWN"
                                                    
                    else:
                        regDir = "MIR_GENE"
                        regStem = "DOWN"
                        
                    if "after" in textBetween:
                        regDir = self.dir2opp[regDir]
                        
                if isPassive:
                    regDir = self.dir2opp[regDir]
                    regStem = self.rdir2opp[regStem]
                    
                return {"regulation_dir": regStem, "interaction_dir": regDir, "passive": isPassive, "reg_detect": "counts opp", "reg_detect_major": "counts"}
                

        stemTypes2Count, mirgene2stem, tknEvs = self.countStems(tokenBetween, geneword, mirword, keepNEU=True)
        if verbose:
            print("counts between 1 stemTypes2Count", stemTypes2Count)
            print("counts between 1 stemTypes2Count", mirgene2stem)
            print("counts between 1 stemTypes2Count", tknEvs)

        if len(stemTypes2Count) > 0:

            regStem = "NEU" if len(stemTypes2Count) == 0 else stemTypes2Count.most_common(1)[0][0]
            
            if self.allowed_check("counts between equal", "counts") and stemTypes2Count["POS"] == stemTypes2Count["NEG"] == 1 and stemTypes2Count["NEU"] <= 3:
                regDir = "MIR_GENE"
                regStem = "DOWN"
                
                #if any(x in textBetween for x in ["after"]):
                #    regDir = self.dir2opp[regDir]
                #for (tkn, tkndir) in tknEvs:
                #    print(tkn,str(doc[tkn.i - 1]), str(doc[tkn.i + 1]))

                if any(str(doc[tkn.i - 1]) in ["by", "after"] for (tkn, tkndir) in tknEvs):

                    if mirword.i < geneword.i:
                        regDir = "MIR_GENE"
                    else:
                        regDir = "GENE_MIR"

                if any(str(doc[tkn.i + 1]) in ["following", "after"] for (tkn, tkndir) in tknEvs):

                    if mirword.i < geneword.i:
                        regDir = "GENE_MIR"
                    else:
                        regDir = "MIR_GENE"
                
                return {"regulation_dir": regStem, "interaction_dir": regDir, "passive": isPassive, "reg_detect": "counts between equal", "reg_detect_major": "counts"}
            
            if self.allowed_check("counts between alternating", "counts") and  stemTypes2Count["POS"] == 2 and stemTypes2Count["NEG"] == 1 and stemTypes2Count["NEU"] <= 2:
                regDir = "MIR_GENE"
                regStem = "DOWN"
                
                return {"regulation_dir": regStem, "interaction_dir": regDir, "passive": isPassive, "reg_detect": "counts between alternating", "reg_detect_major": "counts"}
            
            if regStem == "NEU":
                if stemTypes2Count["POS"] > 0 and stemTypes2Count["POS"] > stemTypes2Count["NEG"]:
                    regStem = "POS"
                    
                elif stemTypes2Count["NEG"] > 0 and stemTypes2Count["NEG"] > stemTypes2Count["POS"]:
                    regStem = "NEG"
                
            
            #regStem = "NEU"

            if self.allowed_check("counts between", "counts"):

                if mirword.i < geneword.i:
                    regDir = "MIR_GENE"
                else:
                    regDir = "GENE_MIR"
                    
                if mirword.i < geneword.i and "-stressed" in textBetween:
                    regDir = "GENE_MIR"
                    
                if mirword.i < geneword.i and str(doc[geneword.i-1]) in ["by", "with"]:
                    regDir = "GENE_MIR"

                    stemTypes2Count,mirgene2stem, tknEvs = self.countStems(tokenBetweenStrict, geneword, mirword)

                    if len(stemTypes2Count) > 0:
                        regStem = stemTypes2Count.most_common(1)[0][0]
                    #else:
                    #    regStem = "NEU"
                    if verbose:
                        print("counts between 2", stemTypes2Count, regStem)
                        
                elif geneword.i < mirword.i and any([str(doc[mirword.i-j]) in ["by", "with", "via", "through"] for j in [1,2]]):
                    regDir = "MIR_GENE"
                    #print("counts between 3 (before)", stemTypes2Count, regStem)

                    stemTypes2Count,mirgene2stem, tknEvs = self.countStems(tokenBetweenStrict, geneword, mirword)
                    if len(stemTypes2Count) > 0:
                        regStem = stemTypes2Count.most_common(1)[0][0]
                    #else:
                    #    regStem = "NEU"
                    if verbose:
                        print("counts between 3", stemTypes2Count, regStem)
                        
                elif isPassive and geneword.i < mirword.i:
                    regDir = self.dir2opp[regDir]

                elif len(stemTypes2Count) == 1 and "NEU" in stemTypes2Count:
                    
                    for (tkn, tknDir) in tknEvs:

                        if str(doc[tkn.i + 1]) in ["with"]:
                            #regStem = "NEU"
                            regDir = "MIR_GENE"

                
                lBound = mirword.i if mirword.i < geneword.i else geneword.i
                rBound = mirword.i if mirword.i > geneword.i else geneword.i

                betweenTokens = [x for x in tknEvs if lBound < x[0].i < rBound]
                betweenTokenTypes = Counter([x[1] for x in betweenTokens])

                if verbose:
                    print("Between Tokens", betweenTokens)

                if len(betweenTokenTypes) == 1 and 'NEG' in betweenTokenTypes:
                    if len(betweenTokens) % 2 == 0:
                        regDir = "MIR_GENE"
                        regStem = "NEG"

                        if verbose:
                            print("Between Tokens Switch")

                
                if len(betweenTokens) > 0:
                    
                    rTkn = sorted(betweenTokens, key=lambda x: rBound-x[0].i)[0]

                    #for tkn, tknStem in betweenTokens:
                    if str(doc[rTkn[0].i+1]) in ["by", "on"]:

                        if verbose:
                            print("visited switch dir")

                        if mirword.i < geneword.i:
                            regDir = "GENE_MIR"
                        
                        #break

                if stemTypes2Count["POS"] > 0 and stemTypes2Count["POS"] == stemTypes2Count["NEG"]:
                    regStem = "NEU"                       
                    
                #if geneword.i < mirword.i:
                #    if str(geneword).endswith("-regulating"):
                #        regStem = "NEU"

                if mirword.i < geneword.i:
                    
                    if str(geneword).endswith( ("-regulating") ):
                        regDir = "MIR_GENE"
                        regStem = "NEU"
                    elif str(geneword).endswith( ("dependent", "independent") ) or ((len(doc) > geneword.i+1) and str(doc[geneword.i+1]) in ["dependent", "independent"]):
                        regDir = "GENE_MIR"
                        regStem = "NEU"

                    elif (len(doc) > geneword.i+1) and any([str(doc[geneword.i+j]) in ["target"] for j in [1]]):
                        regDir = "GENE_MIR"
                        regStem = "NEU"

                    elif (len(doc) > geneword.i+1) and any([str(doc[mirword.i+j]) in ["target", "targets"] for j in [1]]):
                        regDir = "MIR_GENE"
                        regStem = "NEU"


                elif geneword.i < mirword.i:

                    if str(mirword).endswith( ("-regulating") ):
                        regDir = "GENE_MIR"
                        regStem = "NEU"
                    elif str(mirword).endswith( ("dependent", "independent") ):
                        regDir = "MIR_GENE"
                        regStem = "NEU"

                    elif str(geneword).endswith(("-regulating", "-targeting")):
                        regStem = "NEU"
                        regDir = "MIR_GENE"
                    elif str(geneword).endswith(("dependent", "induced")):
                        regDir = "GENE_MIR"
                        regStem = "NEU"

                    elif (len(doc) > mirword.i+1) and any([str(doc[mirword.i+j]) in ["target", "targets"] for j in [1]]):
                        regDir = "MIR_GENE"
                        regStem = "NEU"



                        
                if mirword.i < geneword.i:
                    if "after the addition" in textBetween:
                        regStem = "DOWN"
                        regDir = "GENE_MIR"
                    
                    elif any([str(doc[geneword.i-j]) in ["after", "following"] for j in [1,2,3]]):
                        regStem = "DOWN"
                        regDir = "GENE_MIR"
                        
                    
                if regStem == "NEG":
                    regStem = "DOWN"
                elif regStem == "POS":
                    regStem = "UP"

                """
                if isPassive:

                    if verbose:
                        print("Passive", geneword.dep_, geneword.head, geneword.head.dep_)
                        print("Passive", mirword.dep_, mirword.head, mirword.head.dep_)

                    if geneword.i < mirword.i:
                        if geneword.dep_ in ["compound"] and geneword.head.dep_ in ["nsubjpass"]:
                            regDir = self.dir2opp[regDir]

                    else:
                        if mirword.dep_ in ["compound"] and mirword.head.dep_ in ["nsubjpass"]:
                            regDir = self.dir2opp[regDir]
                """


                if verbose:
                    print("counts between final", isPassive)
                    print([str(x) for x in tokenBetween])
                    print(stemTypes2Count)
                    print(tknEvs)

                return {"regulation_dir": regStem, "interaction_dir": regDir, "passive": isPassive, "reg_detect": "counts between", "reg_detect_major": "counts"}

        return None


class SentenceRelationChecker:

    def __init__(self, nlp, nlp_ent=None, active_checks=None):
        self.relCheck = MirGeneRelCheck(nlp_ent=nlp_ent)
        self.nlp = nlp
        self.nlp_ent = nlp_ent

        if active_checks is None:
            active_checks = self.relCheck.checks_available

        self.active_checks = [x for x in active_checks if x in self.relCheck.checks_available]



    def fixPos(self, pos, sentence):
        posElems = [(idx, chr(i), i) for idx, i in enumerate(sentence) if i > 127 and idx < pos]

        if len(posElems) > 0:
            pos -= len(posElems)

        return pos


    def __offset_positions(self, entDict, value):
        entDict["entity_location"] = (entDict["entity_location"][0] + value, entDict["entity_location"][1] + value)
        return entDict

    def __decode_bytes(self, byteSent):
        posMapping = {}
        decStr = ""
        xidx = 0
        while xidx < len(byteSent):
            x = byteSent[xidx]
                
            if x > 127:
                
                #print(xidx, "{0:b}".format(x))

                utfLen = 1

                # refer to https://en.wikipedia.org/wiki/UTF-8 for the byte shape to length
                if (x & 0b11100000) == 0b11000000:
                    # this and next pos
                    utfLen = 2
                elif (0b11110000 & x) == 0b11100000:
                    # this and next pos
                    utfLen = 3
                elif (0b11111000 & x) == 0b11110000:
                    # this and next pos
                    utfLen = 4
                elif (0b11111100 & x) == 0b11111000:
                    # this and next pos
                    utfLen = 5
                elif (0b11111110 & x) == 0b11111100:
                    # this and next pos
                    utfLen = 6
                    
                for y in range(xidx, xidx+utfLen):
                    posMapping[y] = len(decStr)

                utfBytes = byteSent[xidx:xidx+utfLen]
                #print(utfLen, utfBytes)
                utfStr = utfBytes.decode()
                decStr += utfStr
                
                xidx += utfLen

            else:
                posMapping[xidx] = len(decStr)
                decStr += chr(x)
                
                xidx += 1

        return decStr, posMapping


    def decodeSentence(self, byteSent, e1p, e2p):

        decStr, posM = self.__decode_bytes(byteSent)

        newE1pEntL = (
            posM.get(e1p["entity_location"][0], e1p["entity_location"][0]),
            posM.get(e1p["entity_location"][1], e1p["entity_location"][1])
        )

        newE2pEntL = (
            posM.get(e2p["entity_location"][0], e2p["entity_location"][0]),
            posM.get(e2p["entity_location"][1], e2p["entity_location"][1])
        )

        e1pC = copy.deepcopy(e1p)
        e2pC = copy.deepcopy(e2p)

        e1pC["entity_location"] = newE1pEntL
        e2pC["entity_location"] = newE2pEntL

        e1pC["entity_text"] = decStr[newE1pEntL[0]:newE1pEntL[1]]
        e2pC["entity_text"] = decStr[newE2pEntL[0]:newE2pEntL[1]]

        return decStr, e1pC, e2pC


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
            print(sent, file=sys.stderr)
            print(e1p, file=sys.stderr)
            print(e2p, file=sys.stderr)
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

        allActiveResults = [checkResults[1][x] for x in self.active_checks]
        acceptInteraction = all(allActiveResults)

        if len(self.active_checks) == 0:
            acceptInteraction = True

        #acceptInteraction = checkResults[0]

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
        assert(fix_special_chars == False)

        assert (all([x in e1 for x in ["entity_type", "entity_location", "entity_type_token"]]))
        assert (all([x in e2 for x in ["entity_type", "entity_location", "entity_type_token"]]))

        if type(sentence) == bytes:
            origSent = sentence
            origE1 = e1
            origE2 = e2
            
            sentence, e1, e2 = self.decodeSentence(origSent, e1, e2)
        
        else:
            if fix_special_chars:
                e1Pos = (self.fixPos(e1["entity_location"][0], sentence), self.fixPos(e1["entity_location"][1], sentence))
                e2Pos = (self.fixPos(e2["entity_location"][0], sentence), self.fixPos(e2["entity_location"][1], sentence))

            else:
                e1Pos = e1["entity_location"]
                e2Pos = e2["entity_location"]

            e1["entity_location"] = tuple(e1Pos)
            e2["entity_location"] = tuple(e2Pos)

        if not e1["entity_location"][0] < e2["entity_location"][0]:
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

        if e1["entity_type"] == "mirna" and not e1["entity_text"].upper().startswith(("MIR", "MICRO", "HSA", "MMU", "LET")):
            print(e1, file=sys.stderr)
            print(e2, file=sys.stderr)

        if e2["entity_type"] == "mirna" and not e2["entity_text"].upper().startswith(("MIR", "MICRO", "HSA", "MMU", "LET")):
            print(e1, file=sys.stderr)
            print(e2, file=sys.stderr)

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
            self.nlp.tokenizer.add_special_case(entity, [{'ORTH': entity}]) #, 'POS': 'PROPN'

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

    def __init__(self, nlp_ent=None):

        self.nlp_ent = nlp_ent

        #(re.compile('(%s){e<=1}' % "-mediated"), 15),
        self.surGeneContexts = [
            (re.compile('(%s){e<=3}' % "signaling pathway"), 30),
            (re.compile('(%s){e<=3}' % "signaling circuit"), 22),
            (re.compile('(%s){e<=3}' % "driven pathogenic signaling"), 30),
            (re.compile('(%s){e<=3}' % "receptor signaling"), 30),
            #p53 and ERK1/2 pathways
            (re.compile('(%s){e<=2}' % "pathway"), 20),
            (re.compile('(%s){e<=1}' % "generation"), 12),
            (re.compile('(%s){e<=1}' % "-sensitive"), 15),
            (re.compile('(%s){e<=1}' % "-driven"), 13),
            (re.compile('(%s){e<=1}' % "adenomas"), 15),
            (re.compile('(%s){e<=1}' % "transgenic"), 15),
            (re.compile('(%s){e<=1}' % "genotype"), 15),
            (re.compile('(%s){e<=1}' % "differentiation"), 20), #Th1/Th2 differentiation
            #(re.compile('(%s){e<=1}' % "cells"), 11), #Th1/Th2 cells
            (re.compile('(%s){e<=1}' % "signaling through"), 22),
            (re.compile('(%s){e<=1}' % "signaling"), 10),
            (re.compile('(%s){e<=1}' % "knockout"), 15),
            (re.compile('(%s){e<=1}' % "treatment"), 10),
            (re.compile('(%s){e<=1}' % "strains"), 20), # was 15
            (re.compile('(%s){e<=1}' % "treated"), 15),
            (re.compile('(%s){e<=1}' % "peptides"), 10),
            (re.compile('(%s){e<=1}' % "-positive"), 15),
            (re.compile('(%s){e<=1}' % "-negative"), 10),
            (re.compile('(%s){e<=1}' % "-producing"), 12),
            (re.compile('(%s){e<=1}' % "stressed.{1-7}cells"), 30),
            (re.compile('(%s){e<=2}' % "expressing vector"), 25),
            (re.compile('(%s){e<=1}' % "\(-/-\)"), 10),
            (re.compile('(%s){e<=1}' % "\(-/\+\)"), 10),
            (re.compile('(%s){e<=1}' % "\(\+/-\)"), 10),
            (re.compile('(%s){e<=1}' % "\(\+/\+\)"), 10),
            (re.compile('(%s){e<=1}' % "copy number"), 20),
            (re.compile('(%s){e<=0}' % " ratio"), 6),
            (re.compile('(%s){e<=0}' % "rs[0-9]+"), 10),
            (re.compile('(%s){e<=0}' % "family member"), 16),
            # new mirtex
            (re.compile('(%s){e<=0}' % "mutation"), 16),
            (re.compile('(%s){e<=0}' % "status"), 10),
            (re.compile('(%s){e<=0}' % "\+"), 1),
            (re.compile('(%s){e<=0}' % "\-\W"), 2),
            (re.compile('(%s){e<=0}' % "cells"), 10),
            (re.compile('(%s){e<=2}' % "phosphorylation"), 30),
            (re.compile('(%s){e<=0}' % "methylation"), 30),

            (re.compile('(%s){e<=0}' % "administration"), 20),
            (re.compile('(%s){e<=0}' % "axis"), 10),
            (re.compile('(%s){e<=0}' % "siRNA"), 10),
            (re.compile('(%s){e<=1}' % "patient"), 12),
            (re.compile('(%s){e<=0}' % "fraction"), 10),
            (re.compile('(%s){e<=0}' % "-triggered"), 16),
            (re.compile('(%s){e<=0}' % "-enriched"), 15),
            (re.compile('(%s){e<=0}' % "-inducing"), 15),
            (re.compile('(%s){e<=0}' % "-purified"), 15),

            (re.compile('(%s){e<=0}' % "components"), 16),
            (re.compile('(%s){e<=0}' % "deposition"), 16),
            (re.compile('(%s){e<=0}' % "stage"), 10),
            (re.compile('(%s){e<=0}' % "mutants"), 10),

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
            # new mirtex
            (re.compile('(%s){e<=1}' % "-antagon"), 12),
            (re.compile('(%s){e<=1}' % "cluster"), 15),
            (re.compile('(%s){e<=0}' % "cells"), 20),
            (re.compile('(%s){e<=0}' % "-modulated"), 15),
            (re.compile('(%s){e<=0}' % "axis"), 15),
            (re.compile('(%s){e<=0}' % "pathway"), 15),

        ]

        self.surPreGeneContexts = [
            (re.compile('(%s){e<=2}' % "analysis of active"), 30),
            (re.compile('(%s){e<=1}' % "downstream .* target"), 35),
            (re.compile('(%s){e<=1}' % "downstream .* molecule"), 35),
            # new mirtex
            (re.compile('(%s){e<=1}' % "intron of"), 15),
            (re.compile('(%s){e<=1}' % "experienced"), 18), #disease/condition
            (re.compile('(%s){e<=0}' % "healthy"), 18), #disease/condition

            (re.compile('(%s){e<=1}' % "acetylated"), 18),
            (re.compile('(%s){e<=1}' % "phosphorylated"), 22),
            (re.compile('(%s){e<=1}' % "methylated"), 18),

            (re.compile('(%s){e<=2}' % "phosphorylation of"), 30),
            (re.compile('(%s){e<=0}' % "methylation of"), 30),
            (re.compile('(%s){e<=0}' % "acetylation of"), 30),

            (re.compile('(%s){e<=1}' % "soluble"), 18),
        ]

        self.surPreMiRContexts = [
            (re.compile('(%s){e<=0}' % "near"), 10),
            (re.compile('(%s){e<=0}' % "but not"), 8),
            (re.compile('(%s){e<=1}' % "gene locus of"), 25),
            (re.compile('(%s){e<=1}' % "antagomir to"), 18),
            (re.compile('(%s){e<=1}' % "antagomir of"), 18),
            # new mirtex
            (re.compile('(%s){e<=0}' % "anti"), 6),
            (re.compile('(%s){e<=1}' % "acetylated"), 18),
            (re.compile('(%s){e<=1}' % "phosphorylated"), 22),
            (re.compile('(%s){e<=1}' % "methylated"), 18),
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
                elif token.dep_ in ["advcl"]:

                    keywordFound = False
                    for i in range(0, min(len(newSubtree), 3)):
                        if newSubtree[i].text in ["importantly", "significantly"]:
                            keywordFound = True
                            break

                    if keywordFound:
                        if verbose:
                            print("DeepCheck advcl keyword found")
                        subsentCheck = False

                #subsentCheck = subsentCheck and nounFound
                createCompartment = createCompartment or subsentCheck



            elif token.dep_ in ["amod"]:
                newSubtree = [x for x in token.subtree]
                newSubtree = sorted(newSubtree, key=lambda x: x.idx)
                connectionWords = ["while", "whereas", "thereby", "resulting", "suggestive", "whereby"]
                amodCheck = False
                for i in range(0,min(len(newSubtree),3)):
                    if str(newSubtree[i]).lower() in connectionWords:
                        amodCheck = True

                        if verbose:
                            print("amodCheck is true")
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

                if verbose:
                    print("DeepCheck Conj Subtree")
                    print(newSubtree)

                conjCheck = False
                for i in range(0,min(len(newSubtree),3)):

                    if i == 0 and newSubtree[i].pos_ in ["VERB"]:
                        if verbose:
                            print("deepcheck conj verb")
                        conjCheck = False
                        createCompartment = False

                    if str(newSubtree[i]).lower() in ["actually"]:
                        conjCheck = True
                        break

                createCompartment = createCompartment or conjCheck

            if token.dep_ in ["conj"]:
                newSubtree = [x for x in token.subtree]
                newSubtree = sorted(newSubtree, key=lambda x: x.idx)

                conjCheck = False
                for i in range(0,min(len(newSubtree),1)):

                    if i == 0 and newSubtree[i].pos_ in ["VERB"]:
                        if verbose:
                            print("deepcheck conj verb")
                        conjCheck = False
                        createCompartment = False
                        break

                                 


            if createCompartment:
                newSubtree = [x for x in token.subtree]
                hasSubj = len([x for x in newSubtree if "subj" in x.dep_]) > 0
                #print("creating subtree at", str(token), token.dep_)
                #print("creating subtree", newSubtree)

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
        # apparantly this was not needed!

        if False:
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

            # bubble around actual gene token
            allGeneTokens = sorted(allGeneTokens, key=lambda x: x.i)

            #print("all gene tokens", allGeneTokens)

            for gt in allGeneTokens:

                for gtc in gt.children:

                    # if listing of genes
                    if gtc.dep_ in ["conj"]:
                        
                        gtcStr = str(gtc)

                        for surContext, preview in self.surGeneContexts:
                            mRes = surContext.search(gtcStr)

                            if mRes != None:
                                if verbose:
                                    print("Context Fail GTC", mRes, surContext)
                                return False





        return True

    def check_entity(self, doc, mirword, geneword, verbose):

        if self.nlp_ent != None:
            entDoc = self.nlp_ent(doc.text_with_ws)
    
            enttree = IntervalTree()
            for ent in entDoc.ents:
                enttree.addi(ent.start_char, ent.end_char, (ent.label_, ent.text))

            possTokenEnts = enttree[geneword.idx:(geneword.idx+len(geneword.text))]

            if verbose:
                print("Recognized overlapping entities:", possTokenEnts)

            ents = set([x.data[0] for x in possTokenEnts])

            if 'CELL' in ents:
                if verbose:
                    print("Gene reject:", ents)                
                return False

            for x in possTokenEnts:
                if x.data[0] in ['ORGANISM']:

                    if x.begin == geneword.idx:

                        if verbose:                    
                            print("Organism reject:", ents, geneword, doc.text_with_ws)

                        return False
            
            
            ## SDP check
            return True


            enttree = IntervalTree()

            for ent in entDoc:
                enttree.addi(ent.idx, ent.idx+len(ent.text), (ent.ent_type_, ent.text,ent))

            possGeneTokens = sorted(enttree[geneword.idx:(geneword.idx+len(geneword.text))], key=lambda x: x.begin)
            possmiRTokens = sorted(enttree[mirword.idx:(mirword.idx+len(geneword.text))], key=lambda x: x.begin)


            if len(possGeneTokens) > 0 and len(possmiRTokens) > 0:

                pMir = possmiRTokens[0].data[2]
                pGene = possGeneTokens[0].data[2]

                if verbose:
                    print(possGeneTokens, possmiRTokens)
                    print(pMir, pGene)

                if pMir.i < pGene.i:
                    sdpRes = self._get_sdp_path(entDoc, pMir.i, pGene.i)
                else:
                    sdpRes = self._get_sdp_path(entDoc, pGene.i, pMir.i)

                if len(sdpRes) > 3:
                    if verbose:
                        print("try context sdp",len(sdpRes)-3, len(sdpRes)-1)
                        print(sdpRes)

                    for i in range(len(sdpRes)-3, len(sdpRes)-1):
                        if sdpRes[i][0].ent_type_ in ["GENE_OR_GENE_PRODUCT"] and not sdpRes[i+1][0].dep_ in ["cc", "conj"]:
                            if verbose:
                                print("Context SDP Check")
                                print(sdpRes[i])
                            return False

            return True
                
                
        return True

    @property
    def checks_available(self):
        return ["conj", "sdp", "compartment", "context", "entity"]


    def checkRelation(self, doc, mirword, geneword, verbose=False, relClassifier=None):

        singleResults = {}

        """
        Create single results        
        """

        conjResult, conjs = self.checkCommonConj(doc, mirword, geneword, verbose)
        singleResults["conj"] = conjResult
        conjResult, conjs = self.checkCommonConj(doc, mirword, geneword, verbose, "nmod")
        singleResults["conj_nmod"] = conjResult

        sdpPass, passive, negated = self.checkSDP(doc, mirword, geneword, verbose)
        singleResults["sdp"] = sdpPass

        compPass = self.checkCompartments(doc, mirword, geneword, verbose)
        singleResults["compartment"] = compPass

        sigPathway = self.checkSurContext(doc, mirword, geneword, verbose)
        singleResults["context"] = sigPathway

        entityCheck = self.check_entity(doc, mirword, geneword, verbose)
        singleResults["entity"] = entityCheck

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

            if (sdpRes[idx][1] in ["VERB"] or sdpRes[idx][2] in ["acl:relcl"]) and str(sdpRes[idx][0]) in ["embedded", "located", "contained", "present"]:
                if verbose:
                    print("SDP CONTAIN/EMBED/LOCATE RULE")
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

    def __getConjuncts(self, doc, verbose, conjLabel="conj"):
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

            if conjLabel in ["conj"]:
                tokenConjs = t.conjuncts
            else:
                tokenConjs = []
                if t.dep_ in [conjLabel]:
                    tokenConjs.append(t.head)

            if len(tokenConjs) > 0:
                telems = list(tokenConjs) + [t]

                if verbose:
                    print("Conjuncts", telems)

                idx2t = {e.i: e for e in doc}
                for e in tokenConjs:
                    n = idx2t.get(e.i+1, None)

                    if verbose:
                        print("neighbour", e, n)

                    if n != None:
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
                #childTokens = self.__followChildSel(t, ["case", "amod", "nmod", "dep", "appos", "acl", "dobj", "nummod"], ["compound"], [], verbose=verbose)
                # there is really no reason why compound should only look forward?!
                childTokens = self.__followChildSel(t, ["case", "amod", "dep", "appos", "acl", "dobj", "nummod","compound","nmod"], [], [], verbose=verbose) # "nmod",
                
                for c in childTokens:
                    telems.append(c)

                    #if childDep in ["compound", "amod", "nmod", "appos"]:
                    #    telems.append(child)

                addElemsToCommon(telems)
                continue

            if t.dep_ in ["appos"]:
                telems = [t.left_edge, t.right_edge]
                addElemsToCommon(telems)

        retEdges = defaultdict(list)
        for x in commonEdges:
            retEdges[x] = sorted(commonEdges[x], key=lambda x: x.i)

        return retEdges

    def checkCommonConj(self, doc, mirword, geneword, verbose, conjLabel="conj"):

        commonEdges = self.__getConjuncts(doc, verbose, conjLabel)

        # this is meant to add description of the conj as well.
        #commonEdges[cname] = commonEdges[cname].union(subtreeEdges)

                # now check whether the conjuncts are connected?
        # bla as well as blubb, bla including blubb, ...

        ceNames = sorted([x for x in commonEdges], key=lambda x: commonEdges[x][0].i)
        newConjuncts = {}

        for i in range(0, len(ceNames)-1):

            if not ceNames[i] in newConjuncts:
                newConjuncts[ceNames[i]] = commonEdges[ceNames[i]]

            lEdges = newConjuncts[ceNames[i]]
            rEdges = commonEdges[ceNames[i+1]]
            
            textBetween = str(doc[ lEdges[-1].i : rEdges[0].i ])

            if verbose:
                print("Conunct merge check:", ceNames[i], ceNames[i+1])
                print(textBetween)

            if len(textBetween) > 30:
                continue

            if any([y in textBetween for y in ["as well as", "including", "and"]]):
                if verbose:
                    print("Double Conjunction Found")

                newConjuncts[ceNames[i+1]] = lEdges + rEdges

            elif len(textBetween.replace(" ", "")) == 0:

                if verbose:
                    print("Double Conunction: Empty text bewteen")


                connectionFound = False
                for x in lEdges:
                    if x.head in rEdges and x.dep_ in ["compound"]:
                        connectionFound = True
                        if verbose:
                            print(x, "in other set r")

                for x in rEdges:
                    if x.head in lEdges and x.dep_ in ["compound"]:
                        connectionFound = True
                        if verbose:
                            print(x, "in other set l")

                if connectionFound:
                    if verbose:
                        print("Double Conjunction Found")

                    newConjuncts[ceNames[i+1]] = lEdges + rEdges


        for x in newConjuncts:
            if x in commonEdges and commonEdges[x] == newConjuncts[x]:
                continue

            else:
                commonEdges[(x, "conj")] = newConjuncts[x]
            


        if verbose:
            print("Conjunctions (", conjLabel, ")")
            for cname in commonEdges:
                celes = sorted(commonEdges[cname])
                print(cname, celes)



        for cname in commonEdges:

            celes = sorted(commonEdges[cname], key=lambda x: x.i)
            if mirword in celes and geneword in celes:

                if mirword.i < geneword.i:
                    tokenBetween = doc[ mirword.i : geneword.i ]
                else:
                    tokenBetween = doc[ geneword.i : mirword.i ]

                conjTargetRule = False
                for tkn in tokenBetween:

                    if str(tkn) in ["target", "targets"]:
                        if str(doc[tkn.i+1]) == "of":
                            conjTargetRule = True
                            break
                            
                        if str(doc[tkn.i-1]) == "its":
                            conjTargetRule = True
                            break
                    
                if conjTargetRule:
                    if verbose:
                        print("CONJ TARGET OF RULE", tokenBetween)
                    continue

                
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

                aChildren = self.__followHeadSel(mirword, ["appos", "amod"], [], []) # dep
                
                if verbose:
                    print("achildren", aChildren)
                if geneword in aChildren:

                    if verbose:
                        print("achildren saved", aChildren)
                        exit()
                    continue

                if verbose:
                    print("achildren False")

                return False, celes

            elif geneword in celes:

                for t in celes[-2:]:
                    if str(t).lower() in ["pathways", "pathway"]:

                        if verbose:
                            print("Pathway False")

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