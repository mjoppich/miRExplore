import plac
import random
from pathlib import Path
import spacy


# training data
from spacy import displacy

TRAIN_DATA = [
    ("They trade mortgage-backed securities.", {
        'heads': [1, 1, 4, 4, 5, 1, 1],
        'deps': ['nsubj', 'ROOT', 'compound', 'punct', 'nmod', 'dobj', 'punct']
    }),
    ("I like London and Berlin.", {
        'heads': [1, 1, 1, 2, 2, 1],
        'deps': ['nsubj', 'ROOT', 'dobj', 'cc', 'conj', 'punct']
    })
]


@plac.annotations(
    model=("Model name. Defaults to blank 'en' model.", "option", "m", str),
    output_dir=("Optional output directory", "option", "o", Path),
    n_iter=("Number of training iterations", "option", "n", int))
def main(model=None, output_dir=None, n_iter=10):


    nlp = spacy.load('en')  # create blank Language class
    print("Created blank 'en' model")

    # test the trained model
    test_text = "I like securities."
    #test_text= u"At the same time, inhibition of miR-210 significantly reduced the expression of pro-inflammatory cytokines (TNF-α, IL-1β, and IL-6) and chemokines (CCL2 and CCL3), but had no significant effect on anti-inflammatory factors (TGF-β and IL-10)."
    test_text= u"At the same time, inhibition of miR-210 significantly reduced the expression of proinflammatory cytokines and chemokines, but had no significant effect on anti-inflammatory factors."
    test_text2 = u"We conclude that miR-223 directly targets the chemoattractants CXCL2, CCL3, and IL-6 in myeloid cells."

    allsents = [test_text, test_text2]
    alldocs = []

    for sent in allsents:
        doc = nlp(sent)
        alldocs.append(doc)

        alldeps = [(t.text, t.dep_, t.pos_, t.head.text) for t in doc]
        print('Dependencies', alldeps)



        for x in alldeps:
            print(x)

    displacy.serve(alldocs, style='dep', port=5005)


if __name__ == '__main__':
    plac.call(main)

    # expected result:
    # [
    #   ('I', 'nsubj', 'like'),
    #   ('like', 'ROOT', 'like'),
    #   ('securities', 'dobj', 'like'),
    #   ('.', 'punct', 'like')
    # ]