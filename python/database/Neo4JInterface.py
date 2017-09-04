import neo4j.exceptions as neo4je
from neo4j.v1 import GraphDatabase, basic_auth

from utils.idutils import eprint


class neo4jInterface:

    def __init__(self, simulate=True, printQueries=False):
        self.dbQueries = 0
        self.simulateDB = simulate
        self.printQueries = printQueries

        if simulate == False:
            self.driver = GraphDatabase.driver("bolt://localhost:7687", auth=basic_auth("neo4j", "neo4j2"))
            self.session = self.driver.session()
        else:
            self.session = None
            self.driver = None

    def close(self):

        if self.session:
            self.session.close()


    def createIndex(self, label, property):
        createIdx = 'CREATE INDEX ON :{label}({property})'.format(label=label, property=property)
        self.runInDatabase(createIdx)

    def createUniquenessConstraint(self, label, property):
        createConstraint = 'CREATE CONSTRAINT ON (a:{label}) ASSERT a.{property} IS UNIQUE'.format(label=label, property=property)
        self.runInDatabase(createConstraint)

    def dropIndex(self, label, property):
        dropIdx = 'DROP INDEX ON :{label}({property})'.format(label=label, property=property)
        self.runInDatabase(dropIdx)

    def dropUniquenessConstraint(self, label, property):
        dropContraint = 'DROP CONSTRAINT ON (a:{label}) ASSERT a.{property} IS UNIQUE'.format(label=label, property=property)
        self.runInDatabase(dropContraint)

    def deleteNode(self, labels, props, nodename='n'):
        slabels = self.makeLabels(labels)
        propclause = self.makeWherePropsMatch(nodename, props)

        delString = "match ({nodename}{labels}) {propclause} DELETE {nodename}".format(nodename=nodename, labels=slabels, propclause=propclause)
        self.runInDatabase(delString)

    def deleteRelationship(self, selectorLeft, selLeftLabels, propsLeft, selectorRight, selRightLabels, propsRight, labels, props, relname='r'):

        matchRel = self.__createMatchRelationShip(selectorLeft, selLeftLabels, propsLeft, selectorRight, selRightLabels, propsRight, labels, props, relname)
        matchRel += ' DELETE {relname}'.format(relname=relname)

        self.runInDatabase(matchRel)

    def createNode(self, labels, props, nodename='n'):

        slabels = self.makeLabels(labels, sep=':')
        propclause = self.makePropsCreate(props)

        insertString = "create ({nodename}{labels} {propclause}) return {nodename}".format(nodename=nodename,
                                                                                           labels=slabels,
                                                                                           propclause=propclause)
        self.runInDatabase(insertString)

    def matchNodes(self, labels, props, nodename='n'):
        slabels = self.makeLabels(labels)
        propclause = self.makeWherePropsMatch(nodename, props)

        delString = "match ({nodename}{labels}) {propclause} return {nodename}".format(nodename=nodename, labels=slabels, propclause=propclause)
        return self.runInDatabase(delString)

    def makePropsMerge(self, props):

        outStr = ""

        if props == None or len(props) == 0:
            return outStr

        parts = []
        for key in props:

            value = props[key]
            if type(value) == int:
                valueStr = str(value)
            else:
                valueStr = "\'" + str(value) + "\'"

            partStr = "{key}: {value}".format(key=str(key), value=valueStr)

        return "{" + ", ".join(parts) + "}"

    def mergeNode(self, labels, props, addprops):

        nodename = 'n'

        slabels = self.makeLabels(labels)
        propclause = self.makePropsCreate(props)
        newpropsclause = self.makeWherePropsMatch(nodename, addprops, returnPrefix='', junctor=', ')

        delString = "MERGE ({nodename}{labels} {propclause}) ON MATCH SET {newpropsclause} return {nodename}".format(nodename=nodename, labels=slabels, propclause=propclause, newpropsclause=newpropsclause)
        return self.runInDatabase(delString)

    def nodeExists(self, labels, props, nodename='n'):

        slabels = self.makeLabels(labels, sep=':')
        propclause = self.makeWherePropsMatch(nodename, props)

        testString = "match ({nodename}{labels}) {propclause} return {nodename}".format(nodename=nodename, labels=slabels, propclause=propclause)

        res = self.runInDatabase(testString)

        if res == None:
            return False

        allVals = [x for x in res]

        if len(allVals) == 0:
            return False

        return True

    def createNodeIfNotExists(self, labels, props, nodename='n', labelsCheck = None, propsCheck = None):

        # TODO this can be done much nicer -> merge n:labelscheck {propsCheck} set props return n

        if labelsCheck!= None or propsCheck != None:
            if not self.nodeExists(labelsCheck, propsCheck, nodename):
                self.createNode(labels, props, nodename)
        else:
            slabels = self.makeLabels(labels, sep=':')
            propclause = self.makePropsCreate(props)

            insertString = "MERGE ({nodename}{labels} {propclause}) return {nodename}".format(nodename=nodename,
                                                                                               labels=slabels,
                                                                                               propclause=propclause)
            self.runInDatabase(insertString)


    def runInDatabase(self, query):
        self.dbQueries += 1

        if self.dbQueries % 10000 == 0:
            print(self.dbQueries)

        if self.printQueries:
            print(query)

        if self.simulateDB:
            pass
        else:

            try:

                returnVal = self.session.run( query )
                return returnVal

            except neo4je.ClientError as e:
                eprint(e)
                exit(-1)

    def makeLabels(self, labels, sep='|'):

        labelStr = ""

        if labels == None:
            return labelStr

        if len(labels) > 0:
            labelStr = ":" + sep.join(labels)

        return labelStr

    def makeWherePropsMatch(self, selector, props, returnPrefix='WHERE ', junctor=' AND '):

        whereProps = []

        if props != None:
            for key in props:

                value = props[key]
                if value == None:
                    continue

                if type(value) == int:
                    valueStr = str(value)
                else:
                    valueStr = str(value)

                    if "\"" in valueStr:
                        valueStr = valueStr.replace("\"", "\'")

                    valueStr = "\"" + valueStr + "\""

                whereProps.append(  selector + "." + key + "="+valueStr )

        if len(whereProps) == 0:
            return ""

        return returnPrefix + junctor.join(whereProps)

    def makePropsCreate(self, props, nodename=''):

        if len(nodename) > 0:
            if nodename[len(nodename)-1] != '.':
                nodename += "."

        propElems = []

        if props != None:
            for key in props:

                value = props[key]
                if value == None:
                    continue

                if type(value) == int:
                    valueStr = str(value)
                elif type(value) == list:
                    valueStr = str(value)
                elif type(value) == set:
                    valueStr = str(list(value))
                else:
                    valueStr = "\"" + str(value) + "\""

                propElems.append( nodename + str(key) + ": " + valueStr )

        if len(propElems) == 0:
            return ""

        propStr = "{" + ", ".join(propElems) + "}"
        return propStr

    def __createMatchRelationShip(self, selectorLeft, selLeftLabels, propsLeft, selectorRight, selRightLabels, propsRight, labels, props, relname):

        leftLabelStr = self.makeLabels(selLeftLabels)
        rightLabelStr = self.makeLabels(selRightLabels)

        matchPropsLeft = self.makeWherePropsMatch(selectorLeft, propsLeft, returnPrefix='')
        matchPropsRight = self.makeWherePropsMatch(selectorRight, propsRight, returnPrefix='')

        relLabelStr = self.makeLabels(labels)
        mathPropsRel = self.makeWherePropsMatch(relname, props, returnPrefix='')

        # MATCH (n:MIRNA)-[r:MATURE_OF]-(m:MIRNA_PRE)
        relMatch = "MATCH ({selLeft}{selLeftLabels})-[{relname}{rellabels}]-({selRight}{selRightLabels}) ".format(
            selLeft=selectorLeft, selRight=selectorRight, selLeftLabels=leftLabelStr, selRightLabels=rightLabelStr,
            relname=relname, rellabels=relLabelStr)

        if len(matchPropsLeft) > 0 or len(mathPropsRel) > 0 or len(matchPropsRight) > 0:

            middleNode = ""
            if len(matchPropsLeft) > 0 and len(matchPropsRight) > 0:
                middleNode = " AND "

            middleRel = ""
            if len(matchPropsRight) > 0 and len(mathPropsRel) > 0:
                middleRel = " AND "

            whereClause = " WHERE {leftCond}{middleNode}{rightCond}{middleRel}{relCond} ".format(
                leftCond=matchPropsLeft, middleNode=middleNode, middleRel=middleRel, rightCond=matchPropsRight,
                relCond=mathPropsRel)
            relMatch += whereClause

        return relMatch

    def relationshipExists(self, selectorLeft, selLeftLabels, propsLeft, selectorRight, selRightLabels, propsRight, labels, props, relname='r'):

        relCreate = self.__createMatchRelationShip(selectorLeft, selLeftLabels, propsLeft, selectorRight, selRightLabels, propsRight, labels, props, relname)

        relCreate += " return {selLeft},{relname},{selRight}".format(selLeft = selectorLeft, selRight = selectorRight, relname=relname)

        res = self.runInDatabase(relCreate)

        if res == None or res._keys == None or len(res._keys) == 0:
            return False

        return True


    def createRelationshipIfNotExists(self, selectorLeft, selLeftLabels, propsLeft, selectorRight, selRightLabels, propsRight, labels, props, relname='r'):

        if not self.relationshipExists(selectorLeft, selLeftLabels, propsLeft, selectorRight, selRightLabels, propsRight, labels, props, relname):
            self.createRelationship(selectorLeft, selLeftLabels, propsLeft, selectorRight, selRightLabels, propsRight, labels, props, relname)

    def createRelationship(self, selectorLeft, selLeftLabels, propsLeft, selectorRight, selRightLabels, propsRight, labels, props, relname='r'):

        leftLabelStr = self.makeLabels(selLeftLabels)
        rightLabelStr = self.makeLabels(selRightLabels)


        relCreate = "MATCH ({selLeft}{selLeftLabels}), ({selRight}{selRightLabels}) ".format(selLeft=selectorLeft, selRight=selectorRight, selLeftLabels=leftLabelStr, selRightLabels=rightLabelStr)

        matchPropsLeft = self.makeWherePropsMatch(selectorLeft, propsLeft, returnPrefix='')
        matchPropsRight = self.makeWherePropsMatch(selectorRight, propsRight, returnPrefix='')

        if len(matchPropsLeft) > 0 or len(matchPropsRight) > 0:

            middle = ""
            if len(matchPropsLeft) > 0 and len(matchPropsRight) > 0:
                middle = " AND "

            whereClause = " WHERE {leftCond}{middle}{rightCond} ".format(leftCond=matchPropsLeft, middle=middle, rightCond=matchPropsRight)
            relCreate += whereClause

        relLabelStr = self.makeLabels(labels)
        relPropStr = self.makePropsCreate(props)

        relCreate += " CREATE ({selLeft})-[{relname}{rellabels} {relProps}]->({selRight}) return {relname}".format(selLeft = selectorLeft, selRight = selectorRight, relname=relname, rellabels=relLabelStr, relProps=relPropStr)

        retVal = self.runInDatabase(relCreate)

        return retVal

    def createNodeKeyConstraint(self, label, properties, nodeName='n'):

        eprint("WILL ONLY WORK WITH ENTERPRISE EDITION!")
        return None

        labelStr = self.makeLabels(label)
        propStr = ", ".join([nodeName + "." + x for x in properties])

        createConstraint = "CREATE CONSTRAINT ON (n{label}) ASSERT ({propstr}) IS NODE KEY".format(label=labelStr, propstr = propStr)

        return self.runInDatabase(createConstraint)
