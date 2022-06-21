from synonymes.SynonymFile import Synfile

synFile = Synfile("/mnt/d/mgi.syn")

def makeSynStats(synfile):
    synCounts = []
    for syn in synfile:
        synCounts.append(len(syn))

    print("Synonym location", synfile.location)
    print("Synonym count", len(synCounts))
    print("Avg. Synonyms with 0", sum(synCounts) / len(synCounts))
    print("Avg. Synonyms without 0", sum([x for x in synCounts if x > 0]) / len(synCounts))
    print("Total synonyms", sum(synCounts))

synFile = Synfile("/mnt/d/owncloud/data/miRExplore/textmine/synonyms/mgi.syn")
makeSynStats(synFile)

synFile = Synfile("/mnt/d/owncloud/data/miRExplore/textmine/synonyms/hgnc.syn")
makeSynStats(synFile)