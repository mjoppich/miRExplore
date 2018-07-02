from collections import defaultdict

import os

from utils.DataFrame import DataFrame

from textdb.AbstractDBClasses import DataBaseRel, DataBaseDescriptor
from typing import Any


class MirandaRel(DataBaseRel):

    def __init__(self, lent, rent, datasource, transcript, align_score, energy, mirna_start, mirna_end, 
                    lnc_start, lnc_end, align_length, mirna_identity, lncrna_identity, mirna_alignment, 
                    alignment, lncrna_alignment):

        self.lent = lent
        self.rent = rent
        self.data_source = datasource

        self.transcript = transcript
        self.align_score = align_score
        self.energy = energy
        self.mirna_start = mirna_start
        self.mirna_end = mirna_end
        self.lnc_start = lnc_start
        self.lnc_end = lnc_end
        self.align_length = align_length
        self.mirna_identity = mirna_identity
        self.lncrna_identity = lncrna_identity
        self.mirna_alignment = mirna_alignment
        self.alignment = alignment
        self.lncrna_alignment = lncrna_alignment

    @property
    def lid(self):
        return self.lent[0]

    @property
    def rid(self):
        return self.rent[0]

    @property
    def ltype(self):
        return self.lent[1]

    @property
    def rtype(self):
        return self.rent[1]

    @property
    def docid(self):
        return -1

    def toJSON(self):

        retJSON= {
            'ltype': self.ltype,
            'rtype': self.rtype,
            'rid': self.rid,
            'lid': self.lid,
            'data_source': self.data_source,

            'transcript': self.transcript,
            'align_score': self.align_score,
            'energy': self.energy,
            'mirna_start': self.mirna_start,
            'mirna_end': self.mirna_end,
            'lnc_start': self.lnc_start,
            'lnc_end': self.lnc_end,
            'align_length': self.align_length,
            'mirna_identity': self.mirna_identity,
            'lncrna_identity': self.lncrna_identity,
            'mirna_alignment': self.mirna_alignment,
            'alignment': self.alignment,
            'lncrna_alignment': self.lncrna_alignment
        }

        return retJSON

class MirandaRelDB(DataBaseDescriptor):


    def __init__(self, ltype, rtype):

        super().__init__()

        self.ltyped = ltype
        self.rtyped = rtype

    @property
    def ltype(self):
        return self.ltyped

    @property
    def rtype(self):
        return self.rtyped

    @classmethod
    def loadFromFile(cls, filepath, ltype='mirna', rtype='lncrna'):

        ret = MirandaRelDB(ltype, rtype)
        file_base = os.path.basename(filepath)

        mirandaEvidences = DataFrame.parseFromFile(filepath,bConvertTextToNumber=False)

        for mirtEntry in mirandaEvidences:

            lid = mirtEntry['Name_miRNA']
            rid = mirtEntry['Name_gene']
            org, lid = cls.harmonizeMIRNA(lid)

            print(org, lid, rid)

            if not org in ['hsa', 'mmu']:
                continue

            transcript = mirtEntry['Name_transcript']
            align_score = mirtEntry['align_score']
            energy = mirtEntry['energy']
            mirna_start = mirtEntry['mirna_start']
            mirna_end = mirtEntry['mirna_end']
            lnc_start = mirtEntry['lnc_start']
            lnc_end = mirtEntry['lnc_end']
            align_length = mirtEntry['align_len']
            mirna_identity = mirtEntry['mirna_iden']
            lncrna_identity = mirtEntry['lncrna_iden']
            mirna_alignment = mirtEntry['mirna_alignment']
            alignment = mirtEntry['alignment']
            lncrna_alignment = mirtEntry['lncrna_alignment']
            dataSource = 'miranda'

            relations = set([
                MirandaRel((lid, ltype), (rid, rtype), dataSource, transcript, align_score, energy, mirna_start, mirna_end, 
                    lnc_start, lnc_end, align_length, mirna_identity, lncrna_identity, mirna_alignment, alignment, lncrna_alignment)
                ])
         

            for rel in relations:

                ret.ltype2rel[lid].add(rel)
                ret.rtype2rel[rid].add(rel)

            ret.all_ltypes.add(lid)
            ret.all_rtypes.add(rid)

        return ret
