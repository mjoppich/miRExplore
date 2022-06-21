import os, sys

from collections import Counter


import argparse


if __name__ == '__main__':

    
    parser = argparse.ArgumentParser(description='Convert Medline XML to miRExplore base files')
    parser.add_argument('-i', '--index', type=argparse.FileType("r"), nargs="+", required=True, help="input ontology file")
    parser.add_argument('-m', '--min-count', type=int, default=1000, help="if above count, print")
    args = parser.parse_args()


    occCounter = Counter()
    for inFile in args.index:

        print(inFile.name, file=sys.stderr)

        for line in inFile:

            line = line.split("\t")
            occCounter[line[2]] += 1


    for x in sorted([x for x in occCounter], key=lambda x: occCounter[x]):
        if occCounter[x] > args.min_count:
            print(x, occCounter[x], file=sys.stderr)
            print(x)
            