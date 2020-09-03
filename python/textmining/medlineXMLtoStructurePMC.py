import glob
import os, sys
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")


from collections import defaultdict
from lxml import etree

from utils.idutils import eprint
import logging
import nltk.data

from utils.parallel import MapReduce

logger = logging.getLogger('convertJatsToText')


class PubmedJournal:
    def __init__(self, journal, isoabbrev):
        self.journal = journal
        self.isoabbrev = isoabbrev


class PubmedEntry:
    
    def __init__(self, pubmedId):

        self.pmid = pubmedId
        self.pmc = None
        self.created = None
        self.journal = None

        self.title = None
        self.abstract = None
        self.sections = None

        self.pub_types = None
        self.cites = None

        self.authors = []
        self.doi = None

    def getID(self):
        return self.pmid
            

    def _makeSentences(self, content, tokenizer):

        returns = []

        if type(content) == str:
            returns = tokenizer.tokenize(content)
        else:
            for x in content:
                sents = tokenizer.tokenize(x)
                returns += sents

        #returns = [x + "." for x in returns]

        return returns

    def to_sentences(self, tokenizer):

        finalSents = []

        abstracts = [self.abstract[x] for x in self.abstract]
        abstractSents = self._makeSentences(abstracts, tokenizer)

        sections = [self.sections[x] for x in self.sections]
        textSents = self._makeSentences(sections, tokenizer)

        titleSents = self._prepareSentences(str(self.pmid), 1, [self.title])
        for x in titleSents:
            finalSents.append(x)

        if len(abstractSents) > 0:
            abstractSents = self._prepareSentences(str(self.pmid), 2, abstractSents)

            for x in abstractSents:
                finalSents.append(x)

        if len(textSents) > 0:
            textSents = self._prepareSentences(str(self.pmid), 3, textSents)

            for x in textSents:
                finalSents.append(x)

        return finalSents

    def _prepareSentences(self, articleName, module, sents):

        iSent = 1
        outContent = []

        for x in sents:

            if x == None:
                continue

            x = x.strip()
            #x = x.strip(',.;')

            if len(x) > 0:

                content = str(articleName) + "." + str(module) + "." + str(iSent) + "\t" + str(x)
                outContent.append(content)
                iSent += 1

        return outContent


    def printIgnoration(self):

        for (section, cnt) in self.ignoredSections.most_common():
            print(str(section) + " " + str(cnt))

    @classmethod
    def get_node_with_attrib(cls, node, path, attrib, attribvalue, default=None):
        try:
            values = [x for x in node.findall(path)]

            for idNode in values:
                if not attrib in idNode.attrib:
                    continue

                idType = idNode.attrib[attrib]

                if idType == attribvalue:
                    return idNode

            return default
        except:
            return default

    @classmethod
    def get_nodes_with_attrib(cls, node, path, attrib, attribvalue, default=None):
        try:
            values = [x for x in node.findall(path)]
            keepNodes = []

            for idNode in values:
                if not attrib in idNode.attrib:
                    continue

                idType = idNode.attrib[attrib]

                if idType == attribvalue:
                    keepNodes.append(idNode)

            return keepNodes
        except:
            return default

    @classmethod
    def get_node(cls, node, path, default=None):
        try:
            value = node.find(path)
            return value
        except:
            return default

    @classmethod
    def get_value_from_node(cls, node, path, default=None):
        try:
            if path != None:
                valueElem = node.find(path)
                value = valueElem.text
                return value
            else:
                return node.text
        except:
            return default

    @classmethod
    def get_inner_text_from_node(cls, node, default=[]):
        if node == None:
            return default
        texts = [x.strip() for x in node.itertext() if len(x.strip()) > 0]

        if len(texts) == 0:
            return default
        elif len(texts) == 1:
            return texts[0]
        else:
            return texts

    @classmethod
    def get_inner_text_from_path(cls, node, path, default=None):
        fnode = cls.get_node(node, path, None)
        textReturn =  cls.get_inner_text_from_node(fnode, default=default)

        if type(textReturn) == list:
            textReturn = " ".join(textReturn)
        
        return textReturn

    @classmethod
    def get_nodes(cls, node, path):
        try:
            return [x for x in node.find(path)]
        except:
            return []

    @classmethod
    def _find_doi(cls, node):
        if node == None:
            return None

        idnodes = cls.get_nodes(node, 'PubmedData/ArticleIdList')

        for idnode in idnodes:
            if not 'IdType' in idnode.attrib:
                continue

            idType = idnode.attrib['IdType']

            if idType == 'doi':
                return idnode.text

        return None

    @classmethod
    def _find_authors(cls, node):
        if node == None:
            return None

        authNodes = cls.get_nodes_with_attrib(node, 'front/article-meta/contrib-group//contrib',"contrib-type", "author")

        allAuthors = []

        for authorNode in authNodes:

            lastName = cls.get_value_from_node(authorNode, 'name/given-names', '')
            foreName = cls.get_value_from_node(authorNode, 'name/surname', '')
            initials = ''
            affiliation = ''

            if lastName == '' and foreName == '':
                continue 

            allAuthors.append( (lastName, foreName, initials, affiliation) )

        return allAuthors

    @classmethod
    def get_node_values(cls, node):
        if node == None:
            return None

        childNodes = [x for x in node]

        allValues = []

        for child in childNodes:

            value = cls.get_value_from_node(child, None, None)

            if value != None:
                allValues.append(value)

        return allValues

    def cleanAbstract(self, abstract):

        for x in ['ABSTRACT TRUNCATED AT 250 WORDS', 'ABSTRACT TRUNCATED AT 400 WORDS', 'ABSTRACT TRUNCATED']:

            if abstract.endswith(x):
                return abstract.replace(x, '')

        return abstract

    @classmethod
    def get_abstract(cls, node):

        # TODO care for structued abstracts https://www.nlm.nih.gov/bsd/licensee/elements_descriptions.html#publicationtypelist
        structuredAbstract = {}

        abstractTexts = cls.get_nodes(node, '//abstract')

        for atext in abstractTexts:

            label = atext.attrib.get('id', 'p1').upper()
            text = "".join([x for x in atext.itertext()])#atext.text

            if text != None:
                text = cls.removeLinebreaks(text)
                structuredAbstract[label] = text


        return structuredAbstract

    @classmethod
    def get_text(cls, node):

        # TODO care for structued abstracts https://www.nlm.nih.gov/bsd/licensee/elements_descriptions.html#publicationtypelist
        structuredText = {}

        abstractTexts = cls.get_nodes(node, '/body')

        for secText in abstractTexts:

            secTitle = cls.get_inner_text_from_path(secText, "title")

            

            if secTitle == None:
                secTitle = secText.attrib.get('id', 'general').upper()

            text = "".join([x for x in secText.itertext()])#atext.text

            if text != None:
                text = cls.removeLinebreaks(text)
                structuredText[secTitle] = text


        return structuredText

    @classmethod
    def get_literature(cls, node):

        if node == None:
            return None

        foundNodes = node.findall('back/ref-list//element-citation/pub-id')

        allReturnValues = []

        for subNode in foundNodes:
            PMID = subNode.text

            if PMID != None and PMID != '':
                allReturnValues.append( PMID )

        return allReturnValues

    @classmethod
    def removeLinebreaks(cls, text):

        if text == None or type(text) != str:
            return text

        text = text.replace('\n', ' ').replace('\r', '')
        return text

    @classmethod
    def fromXMLNode(cls, node):

        pmidIDNode = cls.get_node_with_attrib(node, 'front/article-meta/article-id', 'pub-id-type', 'pmid')
        pmcIDNode = cls.get_node_with_attrib(node, 'front/article-meta/article-id', 'pub-id-type', 'pmc')

        pmid = cls.get_value_from_node(pmidIDNode, None)
        pmc = cls.get_value_from_node(pmcIDNode, None)

        print(pmid)
        print(pmc)

        if pmc == None:
            raise ValueError("could not find PMC ID")

        pmc = "PMC{}".format(pmc)

        journal_title = cls.get_inner_text_from_path(node, 'front/journal-meta/journal-title-group/journal-title')
        journal_abbrev_title_node = cls.get_node_with_attrib(node, 'front/journal-meta/journal-id', 'journal-id-type', 'iso-abbrev')
        journal_abbrev_title = cls.get_value_from_node(journal_abbrev_title_node, None)


        journal_doi_node = cls.get_node_with_attrib(node, 'front/article-meta/article-id', 'pub-id-type', 'doi')
        journal_doi = cls.get_value_from_node(journal_doi_node, None)


        print(journal_title)
        print(journal_abbrev_title)
        print(journal_doi)

        title = cls.get_inner_text_from_path(node, 'front/article-meta/title-group/article-title')
        date_published_node = cls.get_node_with_attrib(node, 'front/article-meta/pub-date', "pub-type", "nihms-submitted")

        if date_published_node == None:
            date_published_node = cls.get_node_with_attrib(node, 'front/article-meta/pub-date', "pub-type", "pmc-release")

        date_day = cls.get_value_from_node(date_published_node, "day", "1")
        date_month = cls.get_value_from_node(date_published_node, "month", "1")
        date_year = cls.get_value_from_node(date_published_node, "year", "1")

        print(title)

        abstract = cls.get_abstract(node)
        text = cls.get_text(node)

        authors = cls._find_authors( node )

        publicationTypes = node.find('.').attrib.get("article-type", "unknown")
        citedLiterature = cls.get_literature(node)

        pubmed = PubmedEntry(pmc)
        pubmed.pmc = pmid
        pubmed.created = (date_year, date_month, date_day)
        pubmed.journal = (journal_title, journal_abbrev_title)

        pubmed.title = cls.removeLinebreaks(title)
        pubmed.abstract = abstract
        pubmed.sections = text

        pubmed.authors = authors
        pubmed.doi = journal_doi

        pubmed.pub_types = [publicationTypes]
        pubmed.cites = citedLiterature

        return pubmed

class PubmedXMLParser:

    def __init__(self):
        self.tree = None

    def remove_namespace(self, tree):
        """
        Strip namespace from parsed XML, assuming there's only one namespace per node
        """
        if tree == None:
            return

        for node in tree.iter():
            try:
                has_namespace = node.tag.startswith('{')
            except AttributeError:
                continue  # node.tag is not a string (node is a comment or similar)
            if has_namespace:
                node.tag = node.tag.split('}', 1)[1]

    def parseXML(self, path):

        self.tree = None

        try:
            self.tree = etree.parse(path)
        except:
            try:
                self.tree = etree.fromstring(path)
            except Exception as e:
                eprint("Unable to load graph:", str(e))
                raise
        if '.nxml' in path:
            self.remove_namespace(self.tree)  # strip namespace for

        return self.tree

class PubmedArticleIterator:

    def __init__(self, parser):
        self.parser = parser

    def __iter__(self):

        if self.parser == None:
            return self

        return self.parser.tree.findall('article').__iter__()

    def __next__(self):
        raise StopIteration()


"""

python3 ~/python/miRExplore/python/textmining/downloadPubmedAbstracts.py
python3 medlineXMLtoStructure.py
python3 removeDuplicateSentences.py

"""



if __name__ == '__main__':

    #nltk.data.path.append("/mnt/d/dev/nltk_data/")

    tokenizer_loc = 'tokenizers/punkt/english.pickle'
    tokenizer = nltk.data.load(tokenizer_loc)

    storagePath = '/mnt/f/dev/data/pmcxml/raw/'

    allXMLFiles = [y for x in os.walk(storagePath) for y in glob.glob(os.path.join(x[0], '*.xml'))]
    #allXMLFiles = ['/mnt/f/dev/data/pmcxml/raw/PMC0012XXXXX/PMC1249490.xml', '/mnt/f/dev/data/pmcxml/raw/PMC0012XXXXX/PMC1249491.xml', '/mnt/f/dev/data/pmcxml/raw/PMC0012XXXXX/PMC1249508.xml', '/mnt/f/dev/data/pmcxml/raw/PMC0012XXXXX/PMC1266050.xml', '/mnt/f/dev/data/pmcxml/raw/PMC0012XXXXX/PMC1266051.xml']
    #allXMLFiles = ["/mnt/f/dev/data/pmcxml/raw/PMC0014XXXXX/PMC1455165.xml"]

    print("Found", len(allXMLFiles), "files")

    print("Going through", len(allXMLFiles), " files.")

    def senteniceFile(filenames, env):


        for filename in filenames:
            print(filename)
            storagePath = os.path.dirname(filename) + "/"

            basefile = os.path.basename(filename)
            sentfile = basefile.replace(".xml", ".sent")
            titlefile = basefile.replace(".xml", ".title")
            authorfile = basefile.replace(".xml", ".author")
            citationfile = basefile.replace(".xml", ".citation")
            datefile = basefile.replace(".xml", ".date")
            typefile = basefile.replace(".xml", ".pubtype")
            pmidfile = basefile.replace(".xml", ".pmid")

            pmid2title = {}
            pmid2authors = defaultdict(set)
            pmid2citations = defaultdict(set)

            with open(storagePath + sentfile, 'w') as outfile, open(storagePath + datefile, 'w') as outdate, open(storagePath + typefile, "w") as outtype, open(storagePath + pmidfile, "w") as outpmid:

                pubmedParser = PubmedXMLParser()
                pubmedParser.parseXML(filename)

                for elem in [pubmedParser.tree]:

                    try:
                        entry = PubmedEntry.fromXMLNode(elem)

                        if entry == None:
                            continue

                        sents = entry.to_sentences(tokenizer)

                        for x in sents:
                            outfile.write(x + "\n")

                        pmidID = entry.getID()

                        if entry.created != None:
                            print(pmidID, "\t".join([str(x) for x in entry.created]), sep="\t", file=outdate)

                        if entry.pub_types != None:
                            for ept in entry.pub_types:
                                print(pmidID, ept, sep="\t", file=outtype)

                        if entry.pmc != None:
                            print(pmidID, entry.pmc, sep="\t", file=outpmid)

                        if entry.title != None:
                            pmid2title[pmidID] = entry.title

                        if entry.authors != None and len(entry.authors) > 0:
                            for author in entry.authors: #first, initials, last
                                pmid2authors[pmidID].add( (author[1], author[2], author[0]) )

                        if entry.cites != None and len(entry.cites) > 0:
                            for cite in entry.cites:

                                try:
                                    val = int(cite)
                                    pmid2citations[pmidID].add( val )
                                except:
                                    continue


                    except:

                        eprint("Exception", sentfile)
                        exit(-1)
                        try:

                            pmid = PubmedEntry.fromXMLNode(elem)
                            eprint(pmid)

                        except:
                            pass

                        continue

            with open(storagePath + titlefile, 'w') as outfile:

                print(titlefile)

                for pmid in pmid2title:
                    title = pmid2title[pmid]
                    if title == None or len(title) == 0:
                        continue

                    outfile.write(str(pmid) + "\t" + str(title) + "\n")

            with open(storagePath + authorfile, 'w') as outfile:

                print(authorfile)

                for pmid in pmid2authors:
                    authors = pmid2authors[pmid]

                    if authors == None or len(authors) == 0:
                        continue

                    for author in authors:

                        first = author[0] if author[0] != None else ''
                        initials = author[1] if author[1] != None else ''
                        last = author[2] if author[2] != None else ''

                        outfile.write(str(pmid) + "\t" + "\t".join([first, initials, last]) + "\n")

            with open(storagePath + citationfile, 'w') as outfile:

                print(citationfile)

                for pmid in pmid2citations:
                    citations = pmid2citations[pmid]

                    if citations == None or len(citations) == 0:
                        continue

                    for quote in citations:

                        outfile.write(str(pmid) + "\t" + str(quote) + "\n")


    ll = MapReduce(8)
    result = ll.exec( allXMLFiles, senteniceFile, None, 1, None)

    print("Done")
