import os
from collections import defaultdict
import intervaltree

class RFamEntry:

    # RF00001 CM001002.2      102127887       102128005       94.50   9.8e-16 1       119     0       full

    def __init__(self):

        self.rfam_id=None

        self.orgid=None
        self.chr = None

        self.rfam_start=None
        self.rfam_end = None
        self.strand = None

        self.bit_score = None
        self.evalue_score = None

        self.cm_start = None
        self.cm_end = None
        self.truncated = None

        self.type = None

    def __str__(self):

        return str(self.__dict__)


class RFamOrgDB:

    def __init__(self):

        self.accid2info = {}


    @classmethod
    def loadFromFile(cls, filepath):

        ret = RFamOrgDB()

        with open(filepath, 'r') as fin:

            for line in fin:

                line = line.strip().split()

                #CM000663.2      hsa chr1 1

                print(line)

                accid = line[0]
                orgid = line[1]
                chrid = line[2]

                try:
                    chrnum = int(line[3])
                except:
                    chrnum = line[3]


                ret.accid2info[accid] = {
                    'acc': accid,
                    'org': orgid,
                    'chr': chrid,
                    'chrnum': chrnum
                }

        return ret


class RFamDB:

    def __init__(self):

        self.all_entries = defaultdict(lambda: intervaltree.IntervalTree())


    def add_rfam_entry(self, entry):

        assert(isinstance(entry, RFamEntry))

        entryID = (entry.orgid, entry.chr)


        self.all_entries[entryID].addi(entry.rfam_start, entry.rfam_end, entry)


    def get_entries(self, org, chr, start, stop):

        if not (org, chr) in self.all_entries:
            return None

        it = self.all_entries[(org, chr)]

        return it[start:stop]




    @classmethod
    def loadFromFile(cls, filepath):

        file_base = os.path.basename(filepath)
        fileDir = os.path.dirname(filepath)

        accDB = RFamOrgDB.loadFromFile(fileDir + "/seqaccid_org_chr")

        ret = RFamDB()


        with open(filepath, 'r') as fin:

            for line in fin:

                line = line.strip().split()

                #RF00001 CM001002.2      102127887       102128005       94.50   9.8e-16 1       119     0       full

                rfamid = line[0]
                accid = line[1]

                newaccid = accDB.accid2info.get(accid, None)


                startpos = int(line[2])
                endpos = int(line[3])

                bitscore = float(line[4])
                evaluescore = float(line[5])

                cmstart = int(line[6])
                cmend = int(line[7])

                cmtruncated = int(line[8])

                rfam_type = line[9]

                rfentry = RFamEntry()
                rfentry.rfam_id = rfamid
                rfentry.orgid = accid

                if newaccid != None:

                    rfentry.orgid = newaccid['org']
                    rfentry.chr = newaccid['chr']



                if startpos < endpos:
                    rfentry.strand = '+'

                    rfentry.rfam_start = startpos
                    rfentry.rfam_end = endpos

                else:
                    rfentry.strand = '-'

                    rfentry.rfam_start = endpos
                    rfentry.rfam_end = startpos

                rfentry.bit_score = bitscore
                rfentry.evalue_score = evaluescore

                rfentry.cm_start = cmstart
                rfentry.cm_end = cmend

                rfentry.truncated = cmtruncated
                rfentry.type = rfam_type

                #print(rfentry)


                ret.add_rfam_entry(rfentry)



        return ret



if __name__ == '__main__':

    rfDB = RFamDB.loadFromFile('/mnt/c/ownCloud/data/miRExplore/textmine/aggregated_pmid/rfam.regions.mirexplore')

    allEntries = rfDB.get_entries('hsa', 'chr19', 9050000, 9150000)

    for x in allEntries:
        print(x, x.data)


