import argparse

import sys

import glob
import os

import spacy

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")


from collections import Counter, defaultdict

import nltk

from porestat.utils.DataFrame import DataFrame
import re

from database.ORGMIRs import ORGMIRDB
from synonymes.SynfileMap import SynfileMap
from synonymes.SynonymFile import Synfile, AssocSynfile
from synonymes.mirnaID import miRNA, miRNAPART
from textmining.SentenceDB import SentenceDB, RegPos
from textmining.SyngrepHitFile import SyngrepHitFile


from utils.parallel import MapReduce
from enum import Enum




if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='aggregate tm results', add_help=False)
    parser.add_argument('-s', '--sentdir', type=str, help='where are the sentences?', required=True)
    parser.add_argument('-r', '--resultdir', type=str, help='where are all the index-files?', required=True)



    args = parser.parse_args()


    allfiles = glob.glob(args.sentdir + "/*.sent")
    allfileIDs = [os.path.basename(x).replace(".sent", "") for x in allfiles]
    allfileIDs = sorted(allfileIDs, reverse=True)

    # allfileIDs = [894]
    threads = 4

    if __debug__:
        threads = 4
        sys.stderr.write("Running on threads:" + str(threads) + "\n")

    sys.stderr.write("Debug Mode? " + str(__debug__) + " and threads " + str(threads) + "\n")


    def analyseFile(splitFileIDs, env):

        fileCoocs = []

        for splitFileID in splitFileIDs:

            origFile = args.sentdir + "/" + splitFileID + ".sent"
            resFile = args.resultdir + "/" + splitFileID + ".sent"

            print(resFile)

            with open(origFile, 'r') as fin, open(resFile, 'w') as fout:

                for line in fin:

                    line = line.strip()

                    if line[-1] == '.':
                        line = line[0:-1]

                    fout.write(line + "\n")


    ll = MapReduce(threads)
    result = ll.exec(allfileIDs, analyseFile, None, 1, None)