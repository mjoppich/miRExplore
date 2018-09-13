import glob
import os, sys
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")


from collections import defaultdict
from lxml import etree

from utils.idutils import eprint
import logging

from utils.parallel import MapReduce

logger = logging.getLogger('convertJatsToText')

import nltk.data

class PubmedJournal:
    def __init__(self, journal, isoabbrev):
        self.journal = journal
        self.isoabbrev = isoabbrev


class PubmedEntry:



    def __init__(self, pubmedId):

        self.pmid = pubmedId
        self.created = None
        self.journal = None

        self.title = None
        self.abstract = None
        self.abstract_text = None

        self.pub_types = None
        self.cites = None

        self.authors = []
        self.doi = None

    def getID(self):

        try:
            val = int(self.pmid)
            return val
        except:
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

        titleSents = self._prepareSentences(str(self.pmid), 1, [self.title])
        for x in titleSents:
            finalSents.append(x)

        if len(abstractSents) > 0:
            abstractSents = self._prepareSentences(str(self.pmid), 2, abstractSents)

            for x in abstractSents:
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
        return cls.get_inner_text_from_node(fnode, default=default)

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

        authNodes = cls.get_nodes(node, 'AuthorList')

        allAuthors = []

        for authorNode in authNodes:
            if 'ValidYN' in authorNode.attrib:
                validAuthor = authorNode.attrib['ValidYN'] == 'Y'

                if not validAuthor:
                    continue

            lastName = cls.get_value_from_node(authorNode, 'LastName', '')
            foreName = cls.get_value_from_node(authorNode, 'ForeName', '')
            initials = cls.get_value_from_node(authorNode, 'Initials', '')
            affiliation = cls.get_inner_text_from_path(authorNode, 'AffiliationInfo', '')

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

        abstractTexts = cls.get_nodes(node, 'Abstract')
        if abstractTexts == None or len(abstractTexts) == 0:
            abstractTexts = cls.get_nodes(node, 'OtherAbstract')

        for atext in abstractTexts:

            if atext.tag != 'AbstractText':
                continue

            label = atext.attrib.get('Label', 'GENERAL').upper()
            text = "".join([x for x in atext.itertext()])#atext.text

            if text != None:
                text = cls.removeLinebreaks(text)
                structuredAbstract[label] = text


        return structuredAbstract

    @classmethod
    def get_literature(cls, node):

        if node == None:
            return None

        foundNodes = cls.get_nodes(node, 'MedlineCitation/CommentsCorrectionsList')

        allReturnValues = []

        for subNode in foundNodes:
            if not 'RefType' in subNode.attrib:
                continue

            validAuthor = subNode.attrib['RefType'] == 'Cites'

            if not validAuthor:
                continue

            PMID = cls.get_value_from_node(subNode, 'PMID')
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


        pmid = cls.get_value_from_node(node, 'MedlineCitation/PMID')

        date_created = cls.get_inner_text_from_path(node, 'MedlineCitation/DateCreated')

        articleNode = cls.get_node(node, 'MedlineCitation/Article')

        if articleNode == None:
            return None

        journal_title = cls.get_inner_text_from_path(articleNode, 'Journal/Title')
        journal_abbrev_title = cls.get_inner_text_from_path(articleNode, 'Journal/ISOAbbreviation')

        title = cls.get_inner_text_from_path(articleNode, 'ArticleTitle')
        if title != None and type(title) != str:
            title = "".join(title)

        if title != None and len(title) > 1 and title[0] == '[':
            if title[len(title)-1] == ']':
                title = title[1:len(title)-1]
            elif title[len(title)-2] == ']':
                title = title[1:len(title)-2] + title[len(title)-1]

        abstract = cls.get_abstract(articleNode)

        doi = cls._find_doi( node )
        authors = cls._find_authors( articleNode )

        publicationTypes = cls.get_node_values(cls.get_node(articleNode, 'PublicationTypeList'))
        citedLiterature = cls.get_literature(node)

        pubmed = PubmedEntry(pmid)
        pubmed.created = date_created
        pubmed.journal = (journal_title, journal_abbrev_title)

        pubmed.title = cls.removeLinebreaks(title)
        pubmed.abstract = abstract
        pubmed.abstract_text = abstract

        pubmed.authors = authors
        pubmed.doi = doi

        pubmed.pub_types = publicationTypes
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

        return self.parser.tree.findall('PubmedArticle').__iter__()

    def __next__(self):
        raise StopIteration()


if __name__ == '__main__':
    tokenizer_loc = 'tokenizers/punkt/english.pickle'
    tokenizer = nltk.data.load(tokenizer_loc)

    storagePath = '/mnt/raidtmpbio/joppich/pmid/'
    baseFileName = 'pubmed18n'

    allXMLFiles = glob.glob(storagePath+baseFileName+'*.xml.gz')


    startFrom = 0
    endOn = 2000


    allfiles = []
    for filename in allXMLFiles:
        basefile = os.path.basename(filename)
        basefile = basefile.split('.')[0]
        basefile = basefile.replace(baseFileName, '')
        number = int(basefile)

        if startFrom <= number and number <= endOn:

            allfiles.append(filename)

    print(len(allfiles), allfiles)

    def senteniceFile(filenames, env):


        for filename in filenames:
            print(filename)

            basefile = os.path.basename(filename)
            sentfile = basefile.replace(".xml.gz", ".sent")

            titlefile = basefile.replace(".xml.gz", ".title")
            authorfile = basefile.replace(".xml.gz", ".author")
            citationfile = basefile.replace(".xml.gz", ".citation")

            pmid2title = {}
            pmid2authors = defaultdict(set)
            pmid2citations = defaultdict(set)

            with open(storagePath + sentfile, 'w') as outfile:

                pubmedParser = PubmedXMLParser()
                pubmedParser.parseXML(filename)

                for elem in PubmedArticleIterator(pubmedParser):

                    try:
                        entry = PubmedEntry.fromXMLNode(elem)
                        sents = entry.to_sentences(tokenizer)

                        for x in sents:
                            outfile.write(x + "\n")

                        pmidID = entry.getID()

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
                        try:

                            pmid = elem.find('MedlineCitation/PMID').text
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


    ll = MapReduce(6)
    result = ll.exec( allfiles, senteniceFile, None, 1, None)

    print("Done")
