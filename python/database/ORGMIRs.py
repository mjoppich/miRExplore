from collections import defaultdict

from porestat.utils.DataFrame import DataFrame


class ORGMIRDB:

    def __init__(self, filename):


        self.orgmir2mimat = defaultdict(set)

        self.mimat2mi = {}
        self.mimat2orgmi = {}
        self.mimat2orgmir = {}

        self.mi2orgmi = defaultdict(set)
        self.mi2orgmir = defaultdict(set)
        self.orgmi2mi = defaultdict(set)
        self.orgmir2mi = defaultdict(set)

        allRelations = DataFrame.parseFromFile(filename, bConvertTextToNumber=False)

        for row in allRelations:

            mimat = row['MIMAT']
            orgmir = row['ORGMIR']
            mi = row['MI']
            orgmi = row['ORGMI']

            self.mimat2mi[mimat] = mi
            self.mimat2orgmir[mimat] = orgmir
            self.mimat2orgmi[mimat] = orgmi
            self.orgmir2mimat[orgmir].add(mimat)

            self.mi2orgmi[mi].add(orgmi)
            self.mi2orgmir[mi].add(orgmir)
            self.orgmir2mi[orgmir].add(mi)
            self.orgmi2mi[orgmi].add(mi)