from collections import defaultdict

import os

from synonymes.SynonymFile import Synfile
from textdb.AbstractDBClasses import DataBaseRel, DataBaseDescriptor
from utils.DataFrame import DataFrame


class MiRNADiseaseEntry(DataBaseRel):

    #MiRNADiseaseEntry(("", ltype), (rid, rtype), diseaseDOID, diseaseDescr, diseasePMID, diseaseEffect, diseaseMeasure, orgs)

    def __init__(self, lent, rent, diseaseDOID, diseaseDescr, diseasePMID, diseaseEffect, diseaseMeasure, orgs, datasource, dataid):

        self.lent = lent
        self.rent = rent

        self.doid = diseaseDOID
        self.disease = diseaseDescr

        self.pmid = diseasePMID
        self.effect = diseaseEffect
        self.measurement = diseaseMeasure

        self.organism = tuple(orgs)
        self.data_id = dataid
        self.data_source = datasource

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
        return self.pmid

    def toJSON(self):

        retJSON= {
            'ltype': self.ltype,
            'rtype': self.rtype,
            'rid': self.rid,
            'lid': self.lid,

            'doid': self.doid,
            'descr': self.disease,


            'exp_support': self.measurement,
            'effect': self.effect,

            'data_source': self.data_source,
            'data_id': self.data_id
        }

        if self.docid != None:
            retJSON['docid']= self.pmid

        if self.organism != None:
            retJSON['orgs'] = tuple(self.organism)

        return retJSON

class MiRNADiseaseDB(DataBaseDescriptor):


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
    def loadFromFile(cls, filepath, ltype='gene', rtype='mirna',normGeneSymbols=None):


        ret = MiRNADiseaseDB(ltype, rtype)
        file_base = os.path.basename(filepath)

        mir2disease = DataFrame.parseFromFile(filepath, bConvertTextToNumber=False)

        datasource = "mir2disease"


        # mirna   disease effect  measurement     year    title   pmid    doid
        # hsa-let-7f-2    kidney cancer   up-regulated    microarray      2007    Micro-RNA profiling in kidney and bladder cancers.      17826655        DOID:263

        for idx,disentry in enumerate(mir2disease):

            rid = disentry['mirna']

            diseaseDescr = disentry['disease']
            diseaseDOID = disentry['doid']

            diseasePMID = disentry['pmid']
            diseaseEffect = disentry['effect']
            diseaseMeasure = disentry['measurement']

            org, rid = cls.harmonizeMIRNA(rid)

            if not org in ['hsa', 'mmu']:
                continue


            orgs = set()
            orgs.add(org)

            dataid = file_base + "_" + str(idx)

            relations = set([
                MiRNADiseaseEntry(("", ltype), (rid, rtype), diseaseDOID, diseaseDescr, diseasePMID, diseaseEffect, diseaseMeasure, orgs, datasource, dataid)
                ])

            for rel in relations:

                ret.ltype2rel[""].add(rel)
                ret.rtype2rel[rid].add(rel)

            ret.all_ltypes.add("")
            ret.all_rtypes.add(rid)


        return ret









