from cube.api import Cube

cube=Cube(verbose=True)

cube.load("en", local_models_repository="/mnt/d/nlpcube/")

text="One potential microRNA that regulates Bcan is miR-9 and overexpression of miR-9 can partly rescue the effects of Dicer1 deletion on the MG phenotype."


sentences=cube(text)

for sentence in sentences:
  for entry in sentence:
    print(str(entry.index)+"\t"+entry.word+"\t"+entry.lemma+"\t"+entry.upos+"\t"+entry.xpos+"\t"+entry.attrs+"\t"+str(entry.head)+"\t"+str(entry.label)+"\t"+entry.space_after)
  print("")