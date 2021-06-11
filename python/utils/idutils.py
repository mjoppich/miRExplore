import codecs
import re
import os
from collections import defaultdict
import sys

dataDir = "/mnt/t/owncloud/data/"
#dataDir = "/home/users/joppich/ownCloud/data/"

miRExploreDir = dataDir + "miRExplore/"

class ParseObject(object):
    pass


def isNumber(sText):
    """

    :param sText: a text string
    :return: returns TRUE if sText is a number
    """
    if len(sText) == 0:
        return False

    if str(sText[0]).isdigit():

        try:
            iValue = int(sText)
            return True

        except:

            try:
                fValue = float(sText)
                return True
            except:

                return False

    return False

ltype2label = {
    "gene with protein product": 'GENE_PROTEIN_CODING',
    "pseudogene": 'GENE_PSEUDO',
    "RNA, transfer": 'GENE_TRNA',
    "RNA, micro": 'GENE_MIRNA'
}

aminoAcids = ['Phe', 'Leu', 'Ser', 'Tyr', 'Cys', 'Trp', 'Leu', 'Pro', 'His', 'Gln','Arg', 'Ile', 'Met', 'Thr', 'Asn', 'Lsy', 'Ser', 'Arg', 'Val', 'Ala', 'Asp', 'Glu', 'Gly']

def loadExludeWords(common=True, generic=True, syngrep=True, disease=True, taxnames=True, cell_co=True, manual=True):

    exclWords = defaultdict(set)


    if common:
        """
        loads 10.000 most common english words
        """
        with open(dataDir + "/miRExplore/textmine/excludes/exclude_words.common.syn") as infile:
            for line in infile:
                line = line.strip()

                exclWords['common'].add(line)
                exclWords['common'].add(line.upper())


    if manual:
        """
        manually curated words to remove to given ontologies!
        """
        with open(dataDir + "/miRExplore/textmine/excludes/manual_curated.syn") as infile:
            for line in infile:
                line = line.strip()
                exclWords['manual'].add(line)

    if generic:
        """
        loads words that have been removed manually/curated
        """
        with open(dataDir + "/miRExplore/textmine/excludes/exclude_words.generic.syn") as infile:
            for line in infile:
                line = line.strip()

                exclWords['generic'].add(line)
                exclWords['generic'].add(line.upper())
    
    if generic and syngrep:
        """
        loads words that have been removed manually/curated
        """
        with open(dataDir + "/miRExplore/textmine/excludes/exclude_words.syngrep.syn") as infile:
            for line in infile:
                line = line.strip()

                exclWords['syngrep'].add(line)
                exclWords['syngrep'].add(line.upper())


    if disease:
        """
        loads diseases from pharmgkb?
        """
        with open(dataDir + "/miRExplore/textmine/excludes/exclude_words.disease.syn") as infile:
            for line in infile:
                line = line.strip()

                exclWords['disease'].add(line)
                exclWords['disease'].add(line.upper())


    if taxnames:
        """
        loads species names
        """
        with open(dataDir + "/miRExplore/textmine/excludes/exclude_words.names.syn") as infile:
            for line in infile:
                line = line.strip()

                exclWords['taxnames'].add(line)
                #exclWords['taxnames'].add(line.upper())


    if cell_co:
        """
        loads entities from go cellular component to remove any conflicts with cell types
        """
        with open(dataDir + "/miRExplore/textmine/excludes/exclude_words.cell_co.syn") as infile:
            for line in infile:
                line = line.strip()

                exclWords['cell_co'].add(line)
                exclWords['cell_co'].add(line.upper())

    return exclWords


def makeDBGeneID( symbol ):

    sym = symbol.upper()
    sym = sym.replace(':', '_')

    return sym

speciesName2TaxID = {
    'Homo sapiens': 9606,
    'Mus musculus': 10116,
}


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def mirtarbase_function_label(function):

    label = function.replace('(', '').replace(')', '')
    label = label.replace(' ', '_').replace('-', '_')
    label = label.replace('\"', '').replace('\'', '')
    label = label.replace('/','').replace('\\', '')
    label = re.sub('[^A-Za-z0-9_]', '', label)

    if label[-1] == '_':
        label = label[:-1]

    label = label.upper()

    return label

def mirtarbase_exp_type(expType):

    label = mirtarbase_function_label(expType)
    return label

def printToFile(content, filename, codec='latin1'):

    with codecs.open(filename, 'wb', codec) as outfile:
        for elem in content:
            outfile.write((str(elem)+os.linesep).encode(codec, 'ignore').decode(codec, 'ignore'))


def containsDigitAndChars( someStr ):

    hasDigit = False
    hasChar = False

    for char in someStr:
        if char.isdigit():
            hasDigit = True

        if char.isalpha():
            hasChar = True

    return hasDigit and hasChar