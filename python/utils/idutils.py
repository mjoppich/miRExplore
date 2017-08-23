import re
import os


dataDir = "/mnt/c/ownCloud/data/"

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

def printToFile(content, filename):

    with open(filename, 'w') as outfile:

        for elem in content:

            outfile.write(str(elem) + os.linesep)


def containsDigitAndChars( someStr ):

    hasDigit = False
    hasChar = False

    for char in someStr:
        if char.isdigit():
            hasDigit = True

        if char.isalpha():
            hasChar = True

    return hasDigit and hasChar