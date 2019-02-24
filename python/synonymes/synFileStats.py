from synonymes.SynonymFile import Synfile

synFile = Synfile("/mnt/d/mgi.syn")
synFile = Synfile("/mnt/d/owncloud/data/miRExplore/textmine/synonyms/hgnc.syn")



synCounts = []
for syn in synFile:

    synCounts.append(len(syn))


print("Synonym location", synFile.location)
print("Synonym count", len(synCounts))
print("Avg. Synonyms with 0", sum(synCounts)/len(synCounts))
print("Avg. Synonyms without 0", sum([x for x in synCounts if x > 0])/len(synCounts))