from porestat.utils.DataFrame import DataFrame

class MIRFamily:

    def __init__(self, familyID, familyName):
        self.familyID = familyID
        self.familyName = familyName

        self.childMIMATs = set()

    def addMIMAT(self, mimat):
        self.childMIMATs.add(mimat)


    def __str__(self):
        return self.familyID + "\t" + self.familyName + "\t" + str(len(self.childMIMATs))

class MIRFamilyDB:

    def __init__(self, filename):

        self.families = []
        self.iterIdx = 0

        if not filename is None:


            with open(filename, 'r') as infile:

                currentMIRFamily = None
                allFamilies = {}

                for line in infile:

                    if line.startswith('//'):
                        continue

                    aline = [x.strip() for x in line.split(' ') if len(x) > 0]


                    if aline[0] == 'AC':
                        currentMIRFamily = MIRFamily(aline[1], None)
                        self.families.append( currentMIRFamily )
                    elif aline[0] == 'ID':
                        currentMIRFamily.familyName = aline[1]
                    elif aline[0] == 'MI':
                        currentMIRFamily.addMIMAT( (aline[1], aline[2]) )

    def __iter__(self):
        self.iterIdx = 0
        return self

    def __next__(self):

        if self.iterIdx == len(self.families):
            raise StopIteration

        elem = self.families[self.iterIdx]
        self.iterIdx += 1

        return elem


    def __getitem__(self, item):

        if type(item) == int:
            return self.families[item]

        for family in self:
            if family.familyName == item:
                return family

        raise KeyError("does not contain item: " + str(item))

    def findMI(self, MIID):

        for family in self:

            for (mi, miname) in family.childMIMATs:
                if mi == MIID:
                    return family

        return None