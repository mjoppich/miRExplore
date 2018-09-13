import codecs
from collections import defaultdict
from enum import Enum

import os
from jinja2 import Template
from openpyxl import Workbook

__author__ = 'joppich'

import argparse

import operator
from copy import deepcopy

def toFloat(sText, nDefault=None):
    """

    :param sText: a text string
    :param nDefault: a default value to return if sText is not a number
    :return: int value of sText or None if not a number
    """
    if type(sText) is int:
        return sText

    if type(sText) is float:
        return int(sText)

    if (len(sText) == 0):
        return nDefault

    if str(sText[0]).isdigit():

        try:
            iValue = float(sText)
            return iValue

        except:

            return nDefault

    return nDefault


def toInt(sText, nDefault=None):
    """

    :param sText: a text string
    :param nDefault: a default value to return if sText is not a number
    :return: int value of sText or None if not a number
    """
    if type(sText) is int:
        return sText

    if type(sText) is float:
        return int(sText)

    if (len(sText) == 0):
        return nDefault

    if str(sText[0]).isdigit():

        try:
            iValue = int(sText)
            return iValue

        except:

            return nDefault

    return nDefault

def isNumber(sText):
    """

    :param sText: a text string
    :return: returns TRUE if sText is a number
    """
    if len(sText) == 0:
        return False

    if str(sText[0]).isdigit():

        try:
            iValue = int(sText)
            return True

        except:

            try:
                fValue = float(sText)
                return True
            except:

                return False

    return False


def toNumber(sText, nDefault=None):
    """

    :param sText: a text string
    :return: returns TRUE if sText is a number
    """
    oValue = toInt(sText, None)

    if oValue is None:
        return toFloat(sText, nDefault)

    return oValue

class DataFrameException(Exception):

    def __init__(self, msg):
        super(DataFrameException, self).__init__()

        self.msg = msg

    def __str__(self):
        return self.msg

class DataRowException(Exception):

    def __init__(self, msg):
        super(DataRowException, self).__init__()

        self.msg = msg

class ExportTYPEAction(argparse.Action):

    def __init__(self, option_strings, dest, nargs=None, const=None, default=None, type=None, choices=None, required=False, help=None, metavar=None):
        super(ExportTYPEAction, self).__init__(option_strings, dest, nargs, const, default, type, choices, required, help, metavar)

        self.help = 'Sets type of file export and must be one of {m}'.format(m=', '.join([str(x.value) for x in ExportTYPE]))

    def __call__(self, parser, args, values, option_string=None):

        try:
            eVal = ExportTYPE[values.upper()]
            args.__dict__[ self.dest ] = eVal

        except:

            raise argparse.ArgumentError(None, 'ExportTYPE can not be {n}, '
                                               'it must be one of {m}'.format(n=values,
                                                                              m=', '.join([str(x.value) for x in ExportTYPE])))


class ExportTYPE(Enum):
    CSV=0
    TSV=1
    XLSX=2
    HTML=3
    HTML_STRING=4


class DataSeries:
    """

    A simple list for data

    """

    def __init__(self, data):
        self.data = data
        self.dataIdx = 0

    def __iter__(self):
        self.dataIdx = 0
        return self

    def __next__(self):
        try:
            item = self.__getitem__(self.dataIdx)
        except (IndexError, DataRowException):
            raise StopIteration()

        self.dataIdx += 1

        return item

    def __getitem__(self, item):
        return self.data[item]

    def __setitem__(self, key, value):

        if type(self.data) == tuple:
            oldData = list(self.data)
            oldData[key] = value
            self.data = tuple(oldData)
        else:
            self.data[key] = value

    def append_data(self, value):
        self.data = [x for x in self.data] + [value]

    def remove_data(self, idx):
        self.data = [ self.data[i] for i in range(0, len(self.data)) if i != idx ]

    def to_tuple(self):
        return tuple(self.data)

    def to_list(self):
        return list(self.data)

    def to_set(self):
        return set(self.data)

    def to_scalar(self):
        if len(self.data) > 1:
            raise DataFrameException("Data is not scalar: " + str(self.data))

        return self.data[0]

class DataColumnAccess:
    """
    A DataSeries with column names
    """

    def __init__(self, column2idx):


        self.column2idx = column2idx
        self.idx2default = []

    def getColumnIndex(self, oColumn):

        if oColumn in self.column2idx:
            return self.column2idx[oColumn]

        else:

            try:
                idx = int(oColumn)

                if idx < 0 or idx > len(self.column2idx):
                    raise DataRowException("Invalid column number: " + str(oColumn))

                return idx

            except:
                raise DataFrameException("Invalid column: " + str(oColumn) + "\n\n Available columns: " + str([x for x in self.column2idx]))




    def copyHeader(self):
        return deepcopy(self.column2idx)

    def getHeader(self):

        vReturn = [''] * len(self.column2idx)

        for x in self.column2idx:
            iIdx = self.column2idx[x]
            vReturn[iIdx] = x

        return vReturn

    def columnExists(self, name):
        return name in self.column2idx

    def addColumn(self, name):
        """

        :param name: name of new column (must not be used already)
        :param default: default value to be added to all existing rows
        :return: -1 if not succeeded (name already used), column index otherwise
        """

        if name in self.column2idx:
            raise DataFrameException("Column already exists: " + name)

        newColIdx = len(self.column2idx)
        self.column2idx[name] = newColIdx

        return newColIdx

    def removeColumn(self, name):

        colIdx = self.getColumnIndex(name)
        DataSeries.remove_data(self, colIdx)



class DefaultColumns():

    def __init__(self, **kwargs):

        self.hasDefault = False

        if 'global_default' in kwargs:
            self.hasDefault = True
            self.idx2default = defaultdict(lambda x: kwargs['global_default'])
        else:
            self.idx2default = {}

    def add_default(self, colIdx, default):
        self.idx2default[colIdx] = default

    def get_default(self, colIdx):

        if not colIdx in self.idx2default and not self.hasDefault:
            raise DataFrameException("No default value for column: " + str(colIdx))

        return self.idx2default[colIdx]

class DefaultDataColumnAccess(DefaultColumns, DataColumnAccess):

    def __init__(self, dHeader, **kwargs):

        DataColumnAccess.__init__(self, dHeader)
        DefaultColumns.__init__(self, **kwargs)

    def addColumn(self, name, default=None):

        newColIdx = DataColumnAccess.addColumn(self, name)
        DefaultColumns.add_default(self, newColIdx, default)

        return newColIdx


    def addColumns(self, names, default=None):

        for x in names:
            if self.columnExists(x):
                raise DataFrameException("Column already exists: " + x)

        for x in names:
            self.addColumn(x, default)


class DataRow(DefaultDataColumnAccess, DataSeries):

    def __init__(self, dHeader, data = [], **kwargs):

        DefaultDataColumnAccess.__init__(self, dHeader, **kwargs)
        DataSeries.__init__(self, data)

    def __getitem__(self, item):
        """

        :param item: the key to look for
        :return: if @item is column name, return column name. if @item is int, the item-th element is return.
        """
        try:
            idx = self.getColumnIndex(item)
            return self.data[idx]
        except:
            raise DataRowException("Column not found: " + str(item))

    def addColumn(self, name, default=None):
        DefaultDataColumnAccess.addColumn(self, name, default)
        DataSeries.append_data(self, default)

    def to_tuple(self):
        return DataSeries.to_tuple(self)

    def to_pairs(self):
        return [ (x, self.__getitem__(x)) for x in self.column2idx ]

    @classmethod
    def fromDict(cls, dictionary):

        allitems = list(dictionary.items())

        dHeader = {}
        velements = []

        for i in range(0, len(allitems)):
            dHeader[ allitems[i][0] ] = i
            velements.append(allitems[i][1])

        return DataRow(dHeader, data=velements)


class DataFrame(DataSeries, DefaultDataColumnAccess):
    """

        behave like a dataseries, but have columns

    """

    def __init__(self, default=None):

        DefaultDataColumnAccess.__init__(self, {}, global_default=default)# for a col in row
        DataSeries.__init__(self, []) # for each row


    def __getitem__(self, item):

        if type(item) == int:
            return DataRow(self.column2idx, DataSeries.__getitem__(self, item))
        elif type(item) == str or (item in self.column2idx):

            seriesData = [None] * len(self.data)

            idx = self.getColumnIndex(item)
            for i in range(0, len(self.data)):
                val = self.data[i][idx]
                seriesData[i] = val

            return DataSeries(seriesData)

        elif type(item) == tuple or type(item) == list:

            if len(item) == 1:
                return DataRow(self.column2idx, DataSeries.__getitem__(self, item[0]))

            elif len(item) == 2:
                rowData = DataSeries.__getitem__(self, item[0])
                row = DataRow(self.column2idx, rowData)

                return row[item[1]]

        elif type(item) == slice:

            elemIdx = range(item.start, item.stop, item.step)
            return [DataRow(self.column2idx, DataSeries.__getitem__(self, x)) for x in elemIdx]

    def __setitem__(self, key, value):

        if type(key) == int:
            DataSeries.__setitem__(self, key, value)
        elif type(key) == tuple or type(key) == list:

            if len(key) == 1:
                DataSeries.__setitem__(self, key, value)

            elif len(key) == 2:
                rowData = DataSeries.__getitem__(self, key[0])
                row = DataRow(self.column2idx, rowData)
                row[key[1]] = value
                row_tuple = row.to_tuple()

                DataSeries.__setitem__(self, key[0], row_tuple)

        elif type(key) == slice:

            elemIdx = range(key.start, key.stop, key.step)

            if len(elemIdx) == len(value):
                for x in elemIdx:
                    DataSeries.__setitem__(x, value[x])

    def merge(self, other):

        if not type(other)==DataFrame:
            raise DataFrameException("Can only merge two dataframes")

        owncols = set(self.column2idx)
        othercols = set(other.column2idx)

        for x in owncols:
            if x not in othercols:
                raise DataFrameException("Missing col in other: " + str(x))

        for row in other:
            self.addRow( row )

    def addColumn(self, name, default=None):
        newColIdx = DefaultDataColumnAccess.addColumn(self, name, default)

        for i in range(0, len(self)):
            oldData = list(self.data[i])
            oldData.append(default)

            self.data[i] = tuple(oldData)

        return newColIdx


    def __len__(self):
        return len(self.data)

    def toRowDict(self, oTuple):

        if len(oTuple) != len(self.column2idx):
            return None

        dReturn = {}

        for x in self.column2idx:
            dReturn[x] = oTuple[self.column2idx[x]]

        return dReturn


    def toDataRow(self, idxcol = None, datacol = None):

        if datacol == None:
            datacol = len(self.column2idx)-1

        if idxcol == None:
            idxcol = [x for x in range(0, len(self.data))]
        else:
            idxcol = [x[idxcol] for x in self.data]

        dataDict = {}

        for i in range(0, len(self.data)):

            key = idxcol[i]
            value = self.data[i][datacol]

            dataDict[key] = value

        return DataRow.fromDict(dataDict)




    def applyByRow(self, oColumn, oFunc):

        iColumnIndex = self.getColumnIndex(oColumn)

        for i in range(0, len(self.data)):
            vElemLine = list(self.data[i])

            self[i, iColumnIndex] = oFunc(vElemLine)

    def applyToRow(self, oFunc):
        """

        :param oFunc: must return tuple/list of length of header
        :return:
        """

        for i in range(0, len(self.data)):
            vElemLine = list(self.data[i])
            vReturnLine = oFunc(vElemLine)

            self[i] = vReturnLine

    def getColumn(self, oColumn):

        vReturn = []

        iColumnIndex = self.getColumnIndex(oColumn)

        for i in range(0, len(self.data)):
            vReturn.append(self[i, iColumnIndex])

        return vReturn

    def setColumn(self, oColumn, vNewValues):

        if len(vNewValues) != len(self.data):
            raise ValueError("elements have different lengths")

        iColumnIndex = self.getColumnIndex(oColumn)

        for i in range(0, len(self.data)):
            self.data[i, iColumnIndex] = vNewValues[i]

    def addRow(self, row):

        if not (type(row) == DataRow):
            raise DataFrameException( 'Trying to insert invalid row type: ' + str(type(row)))

        rowheaders = set([x for x in row.column2idx])
        selfheaders = set([x for x in self.column2idx])
        setDifference = rowheaders.difference(selfheaders)

        if len(setDifference) > 0:
            raise DataFrameException( 'Row does not contain all needed headers: ' + str(selfheaders))

        newrow = [None] * len(self.column2idx)
        for x in self.column2idx:
            newrow[ self.column2idx[x] ] = row[x]

        self.data.append(tuple(newrow))

    def getRow(self, oColumn, oValue, oDefaultValue = None):

        iColumnIndex = self.getColumnIndex(oColumn)

        for i in range(0, len(self.data)):

            if self.data[i][iColumnIndex] == oValue:
                return DataRow(self.data[i], self.column2idx)

        return oDefaultValue

    def getRows(self, oColumn, vValues):

        iColumnIndex = self.getColumnIndex(oColumn)

        vReturn = []
        for i in range(0, len(self.data)):

            if self.data[i][iColumnIndex] in vValues:
                vReturn.append( DataRow(self.data[i], self.column2idx) )

        return vReturn


    def findRow(self, oColumn, oValue, oDefaultValue=None):

        iColumnIndex = self.getColumnIndex(oColumn)

        for i in range(0, len(self.data)):

            if self.data[i][iColumnIndex] == oValue:
                return self[i]

        return oDefaultValue

    def findRows(self, oColumn, vValues):

        iColumnIndex = self.getColumnIndex(oColumn)

        vReturn = []

        for i in range(0, len(self.data)):

            if self.data[i][iColumnIndex] in vValues:
                vReturn.append(self.data[i])

        return vReturn

    def __str__(self):

        sortedHeader = sorted(self.column2idx.items(), key=operator.itemgetter(1))

        vHeader = [str(x[0]) for x in sortedHeader]

        sStr = "\t".join(vHeader)

        for oLine in self.data:
            sStr += "\n"
            sStr += "\t".join([str(x) for x in oLine])

        return sStr

    def _makeStr(self, sep='\t'):

        sortedHeader = sorted(self.column2idx.items(), key=operator.itemgetter(1))

        vHeader = [str(x[0]) for x in sortedHeader]

        sStr = sep.join(vHeader)

        for oLine in self.data:
            sStr += "\n"
            sStr += sep.join([str(x) for x in oLine])

        return sStr

    def _writeToFile(self, content, filename):

        with open(filename, 'w') as file:
            file.write(content)

    def _makeXLSX(self, outFile):
        wb = Workbook()
        # grab the active worksheet
        ws = wb.active

        # Data can be assigned directly to cells
        sortedHeader = sorted(self.column2idx.items(), key=operator.itemgetter(1))
        vHeader = [str(x[0]) for x in sortedHeader]

        ws.append(vHeader)

        for oLine in self.data:
            ws.append([str(x) for x in oLine])

        # Save the file
        wb.save( outFile )

    def _makeHTMLString(self):

        headpart = """
                <link rel="stylesheet" href="https://cdn.datatables.net/1.10.15/css/jquery.dataTables.min.css">
                <script src="https://code.jquery.com/jquery-1.12.4.js"></script>
                <script src="https://cdn.datatables.net/1.10.15/js/jquery.dataTables.min.js"></script>
        """

        bodypart = """
        <table id="example" class="display" cellspacing="0" width="100%">
                <thead>
                <tr>
                {% for column in columns %}
                    <th>{{column}}</th>
                    {% endfor %}
                </tr>
                </thead>

                <tbody>
                {%- for row in rows %}
                <tr>
                    {% for idx in indices %}
                    <td>{{ row[idx] }}</td>
                    {%- endfor -%}
                </tr>
                {% endfor -%}
                </tbody>

                <tfoot>
                <tr>
                {% for column in columns %}
                    <th>{{column}}</th>
                    {% endfor %}
                </tr>
                </tfoot>

                </table>

                <script>
                $(document).ready(function() {
                    $('#example tfoot th').each( function () {
                        var title = $(this).text();
                        $(this).html( '<input type="text" placeholder="Search '+title+'" />' );
                    } );

                    // DataTable
                    var table = $('#example').DataTable({
                                                        "columnDefs": [
                                                            { "type": "numeric-comma" }
                                                        ]
                                                    } );

                    // Apply the search
                    table.columns().every( function () {
                        var that = this;

                        $( 'input', this.footer() ).on( 'keyup change', function () {
                            if ( that.search() !== this.value ) {
                                that
                                    .search( this.value )
                                    .draw();
                            }
                        } );
                    } );

                } );
                </script>
        """

        sortedHeader = sorted(self.column2idx.items(), key=operator.itemgetter(1))
        vHeader = [str(x[0]) for x in sortedHeader]
        vIndices = [x[1] for x in sortedHeader]

        jinjaTemplate = Template(bodypart)
        output = jinjaTemplate.render(rows=self.data, indices=vIndices, columns=vHeader)

        return (headpart, output)


    def _makeHTML(self, outFile):

        (headpart, bodypart) = self._makeHTMLString()

        htmlfile="""

        <html>
            <head>
        """ + headpart + """
            </head>
            <body>
        """ + bodypart + """
            </body>
        </html>
        """

        with open(outFile, 'w') as outHtml:
            outHtml.write(htmlfile)







    def export(self, outFile, exType=ExportTYPE.TSV):

        outputText = None

        if exType == ExportTYPE.XLSX and outFile != None:
            self._makeXLSX(outFile)
            return

        if exType == ExportTYPE.HTML and outFile != None:
            self._makeHTML(outFile)
            return

        if exType == ExportTYPE.HTML_STRING:
            return self._makeHTMLString()

        if exType == ExportTYPE.TSV:
            outputText = self._makeStr('\t')
        elif exType == ExportTYPE.CSV:
            outputText = self._makeStr(';')
        elif exType == None:
            outputText = self._makeStr('\t')

        if outFile == None:
            print(outputText)
            return
        else:
            self._writeToFile( outputText, outFile )
            return

        raise DataFrameException('Invalid export type {n} with output file {m}'.format(n=str(exType), m=str(outFile)))


    @classmethod
    def createHeader(cls, oHeaderFrom, cDelim):

        aHeaderLine = None

        if type(oHeaderFrom) == str:
            aHeaderLine = oHeaderFrom.strip().split(cDelim)

        if type(oHeaderFrom) == tuple or type(oHeaderFrom) == list:
            aHeaderLine = oHeaderFrom

        if aHeaderLine == None:
            raise Exception("invalid headerfrom type: must be str, tuple or list")

        dDesc2Idx = {}
        dIdx2Desc = {}

        iIndex = 0
        # oHeader is header
        for sElem in aHeaderLine:
            dDesc2Idx[sElem] = iIndex
            dIdx2Desc[iIndex] = sElem

            iIndex += 1

        return (dDesc2Idx, dIdx2Desc)

    @classmethod
    def createLineTuple(cls, sLine, oHeader, cDelim, vNumberHeader=None, oEmptyValue=None):

        sLine = sLine.strip()

        aLine = sLine.split(cDelim)

        vLine = list(aLine)

        if vNumberHeader is None:
            vNumberHeader = [False] * len(oHeader)

        for e in range(0, len(vNumberHeader)):
            if vNumberHeader[e]:
                oRes = toNumber(vLine[e], vLine[e])
                vLine[e] = oRes

                # if oRes == None:
                #    print("some error " + str(aLine))

        while len(vLine) < len(oHeader):
            vLine.append(oEmptyValue)

        return tuple(vLine)

    @classmethod
    def parseFromFile(cls, sFileName, oHeader=None, cDelim='\t', bConvertTextToNumber=True, encoding="utf-8"):

        if not os.path.isfile(sFileName):
            print("error loading file: " + sFileName)
            print("file does not exist")

            return None

        oNewDataFrame = DataFrame()

        vLines = []

        with codecs.open(sFileName, mode='r', encoding=encoding) as fin:
            vLines = fin.readlines()

        iStartLine = 0

        if oHeader is None:

            # retrieve header from first line
            iStartLine += 1

            oHeaderRes = cls.createHeader(vLines[0], cDelim)
            oHeader = oHeaderRes[0]

        else:

            oHeaderRes = cls.createHeader(oHeader, cDelim)
            oHeader = oHeaderRes[0]

        oNewDataFrame.column2idx = oHeader

        vNumberHeader = [False] * len(oHeader)

        for i in range(iStartLine, len(vLines)):

            sLine = vLines[i]

            if i == iStartLine and bConvertTextToNumber:

                sLine = vLines[i].strip()
                aLine = sLine.split(cDelim)
                vLine = list(aLine)

                for e in range(0, len(vLine)):
                    if isNumber(vLine[e]):
                        vNumberHeader[e] = True

            lineTuple = cls.createLineTuple(sLine, oHeader, cDelim, vNumberHeader, None)

            oNewDataFrame.data.append(lineTuple)

        return oNewDataFrame