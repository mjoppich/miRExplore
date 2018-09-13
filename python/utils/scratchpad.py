import spacy
from spacy import displacy

nlp = spacy.load('en')
doc = nlp("".join([u'Microvesicles released by apoptotic human neutrophils suppress proliferation and IL-2/IL-2 receptor expression of resting T helper cells.', u'Macrophages are influenced by chemokines from neutrophils']))

alldeps = [(t.idx, t.text, t.dep_, t.pos_, t.head.text) for t in doc]

for t in doc:
    print(t.idx, t.text, t.dep_, t.pos_, t.head.text, [x for x in t.conjuncts], [x for x in t.children])

displacy.serve(doc, style='dep', port=5005)