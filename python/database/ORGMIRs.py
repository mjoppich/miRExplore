from collections import defaultdict

from porestat.utils.DataFrame import DataFrame


class ORGMIRDB:

    def __init__(self, filename):


        self.orgmir2mimat = defaultdict(set)
        self.mimat2orgmir = {}
        self.mimat2mi = {}
        self.mi2orgmir = defaultdict(set)
        self.orgmir2mi = defaultdict(set)

        allRelations = DataFrame.parseFromFile(filename, bConvertTextToNumber=False)

        for row in allRelations:

            mimat = row['MIMAT']
            orgmir = row['ORGMIR']
            mi = row['MI']

            self.mimat2orgmir[mimat] = orgmir
            self.orgmir2mimat[orgmir].add(mimat)

            self.mimat2mi[mimat] = mi
            self.mi2orgmir[mi].add(orgmir)
            self.orgmir2mi[orgmir].add(mi)
