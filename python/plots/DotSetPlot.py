import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches


class DotSetPlot:

    def __init__(self):

        pass


    def plot(self, labels, data, order=None, sortData=True, max=-1, label2orderspan=None, numbers=None):

        if order == None:
            order = [x for x in labels]

        if label2orderspan == None:
            label2orderspan = {}
            for x in order:
                label2orderspan[x] = 3

        assert(len(order) == len(labels))

        allElements = [x for x in data]

        if sortData:
            allElements = sorted(allElements, key=lambda x: sum([len(data[x][y]) for y in data[x]]), reverse=True)

        if max != -1 and len(allElements) > max:
            allElements = allElements[0:max]

        allElements = allElements[::-1]

        maxCols = sum([len(labels[x]) for x in labels])

        fig = plt.figure(figsize=(10,16))

        prevAx = None
        allAx = []

        orderSpan = 3
        numberSpan = 1

        maxSpan = 0

        if numbers != None:
            maxSpan += len(numbers) * numberSpan

        maxSpan += sum([label2orderspan[x] for x in order])

        if numbers != None:
            for nidx, colName in enumerate(numbers):

                #ax1 = plt.subplot(1, len(labels)+len(numbers), len(allAx) + 1, sharey=allAx[0] if len(allAx) > 0 else None)
                ax1 = plt.subplot2grid((1, maxSpan), (0,nidx*numberSpan), colspan=numberSpan, sharey=allAx[0] if len(allAx) > 0 else None )

                for i in range(0, len(allElements)):

                    if i % 2 == 0:
                        rect = patches.Rectangle((-0.50, i + 0.5), len([colName]), 1, linewidth=2, edgecolor='none',
                                                 zorder=-1, facecolor='gainsboro')
                        ax1.add_patch(rect)


                for idx, num in enumerate(reversed(numbers[colName])):
                    plt.text(0, -0.3+idx, str(num), horizontalalignment='center')

                plt.setp(ax1, yticks=range(0, len(allElements)), yticklabels=allElements)
                ax1.set_xticks(range(0, len([colName])))
                ax1.set_xticklabels([colName], rotation=90, ha="right")

                ax1.get_yaxis().set_visible(len(allAx) == 0)

                ax1.set_xlim([-0.1,  0.1])
                ax1.set_ylim([-0.1, len(allElements) - 0.5])
                ax1.tick_params(axis=u'both', which=u'both', length=0)
                ax1.spines["top"].set_visible(False)
                ax1.spines["right"].set_visible(False)
                ax1.spines["left"].set_visible(False)
                ax1.spines["bottom"].set_visible(False)

                allAx.append(ax1)


        curColIdx = 0

        if numbers != None:
            curColIdx = len(numbers)*numberSpan


        for oidx, catLabel in enumerate(order):

            if isinstance(order, dict):
                catLabels = order[catLabel]
            else:
                catLabels = list(labels[catLabel])

            print(catLabels)

            #ax1 = plt.subplot(1, len(labels)+len(numbers), len(allAx)+1, sharey=allAx[0] if len(allAx) > 0 else None)
            ax1 = plt.subplot2grid((1, maxSpan), (0, curColIdx), colspan=label2orderspan[catLabel],
                                   sharey=allAx[0] if len(allAx) > 0 else None)
            curColIdx += label2orderspan[catLabel]

            for i in range(0, len(allElements)):

                if i % 2 == 0:
                    rect = patches.Rectangle((-0.50, i+0.5), len(catLabels), 1, linewidth=2, edgecolor='none', zorder=-1, facecolor='gainsboro')
                    ax1.add_patch(rect)


            arraySize = (len(allElements), len(catLabels))
            npArray = np.zeros(arraySize, dtype=int)

            print(arraySize)

            layoutX = []
            layoutY = []

            for i in range(arraySize[0]):
                elemData = data[allElements[i]][catLabel]
                for j in range(arraySize[1]):

                    if catLabels[j] in elemData:
                        npArray[i,j] = 1
                    else:
                        npArray[i, j] = 0

                    layoutX.append(j)
                    layoutY.append(i)


            ax1.scatter(layoutX, layoutY, c="darkgrey", s=50)

            dotIdx = np.nonzero(npArray>=1)

            allX = []
            allY = []

            for ze in zip(dotIdx[0], dotIdx[1]):
                allX.append(ze[1])
                allY.append(ze[0])

            ax1.scatter(allX, allY, c="black", s=50)

            print(allX)
            print(allY)
            #ax1.xticks(np.arange(len(allElements)), allElements)

            plt.setp(ax1, yticks=range(0, len(allElements)), yticklabels=allElements)
            plt.setp(ax1, xticks=range(0, len(catLabels)), xticklabels=catLabels)

            ax1.set_xticks(range(0, len(catLabels)))
            ax1.set_xticklabels(catLabels, rotation=90,)


            ax1.get_yaxis().set_visible(len(allAx)==0)

            ax1.set_xlim([-0.5, len(catLabels)-0.5])
            ax1.set_ylim([-0.5, len(allElements)-0.5])
            ax1.tick_params(axis=u'both', which=u'both', length=0)
            ax1.spines["top"].set_visible(False)
            ax1.spines["right"].set_visible(False)
            ax1.spines["left"].set_visible(False)
            ax1.spines["bottom"].set_visible(False)

            allAx.append(ax1)

        plt.tight_layout()




if __name__ == '__main__':


    labels = [
        ["A", "B", "C"],
        ["D", "E", "F"]
    ]

    data = {
        'miR-9': [["A", "C"], ["E"]],
        'miR-91': [["B", "C"], ["D", "E"]],
        'miR-11': [[], ["D"]]
    }

    DotSetPlot().plot(labels, data, max=2)