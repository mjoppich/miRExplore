import codecs
import re
import os
from collections import defaultdict

dataDir = "/mnt/c/ownCloud/data/"
dataDir = "/home/users/joppich/ownCloud/data/"

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

def loadExludeWords():

    exclWords = defaultdict(set)

    with open(dataDir + "/miRExplore/textmine/excludes/exclude_words.generic.syn") as infile:
        for line in infile:
            line = line.strip()

            exclWords['generic'].add(line)
            exclWords['generic'].add(line.upper())

    with open(dataDir + "/miRExplore/textmine/excludes/exclude_words.disease.syn") as infile:
        for line in infile:
            line = line.strip()

            exclWords['disease'].add(line)
            exclWords['disease'].add(line.upper())

    with open(dataDir + "/miRExplore/textmine/excludes/exclude_words.names.syn") as infile:
        for line in infile:
            line = line.strip()

            exclWords['names'].add(line)
            exclWords['names'].add(line.upper())

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