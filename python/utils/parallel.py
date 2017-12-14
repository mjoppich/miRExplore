import queue
import threading
import time
import sys
import os
import multiprocessing
import pickle

from pathos import multiprocessing as mp

import itertools

class MapReduce:

    def __init__(self, procs = 4):
        self.pool = None
        self.nprocs = procs

    def exec(self, oIterable, oFunc, sEnvironment, chunkSize = 1, pReduceFunc = None):

        self.pool = mp.ProcessPool(nodes=self.nprocs)

        allResults = []

        resultObj = None

        for x in self.chunkIterable(oIterable, chunkSize):
            allResults.append( self.pool.apipe( oFunc, x, sEnvironment ) )

        while len(allResults) > 0:

            i=0
            while i < len(allResults):

                if allResults[i].ready():

                    result = allResults[i].get()

                    if pReduceFunc != None:

                        resultObj = pReduceFunc(resultObj, result, sEnvironment)

                    else:

                        if resultObj == None:
                            resultObj = []

                        resultObj.append(result)

                    del allResults[i]

                else:
                    i += 1

            time.sleep(0.5)

        self.pool.join()
        self.pool.close()

        return resultObj



    @classmethod
    def chunkIterable(cls, iterable, size):

        it = iter(iterable)
        item = list(itertools.islice(it, size))

        while item:
            yield item
            item = list(itertools.islice(it, size))