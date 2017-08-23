from neo4j.v1 import GraphDatabase, basic_auth

class neo4jInterface:

    def __init__(self, simulate=True, printQueries=True):
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

    def createNode(self, labels, props, nodename='n'):

        slabels = self.makeLabels(labels)
        propclause = self.makePropsCreate(props)

        insertString = "create ({nodename}{labels} {propclause}) return {nodename}".format(nodename=nodename,
                                                                                           labels=slabels,
                                                                                           propclause=propclause)
        self.runInDatabase(insertString)

    def nodeExists(self, labels, props, nodename='n'):

        slabels = self.makeLabels(labels)
        propclause = self.makeWherePropsMatch(nodename, props)

        testString = "match ({nodename}{labels}) {propclause} return n".format(nodename=nodename, labels=slabels, propclause=propclause)

        res = self.runInDatabase(testString)

        if res == None or res._keys == None or len(res._keys) == 0:
            return False

        return True

    def createNodeIfNotExists(self, labels, props, nodename='n', labelsCheck = None, propsCheck = None):

        if labelsCheck == None:
            labelsCheck = labels

        if propsCheck == None:
            propsCheck = props

        if not self.nodeExists(labelsCheck, propsCheck, nodename):
            self.createNode(labels, props, nodename)


    def runInDatabase(self, query):
        self.dbQueries += 1

        if self.dbQueries % 10000 == 0:
            print(self.dbQueries)

        if self.simulateDB:
            if self.printQueries:
                print(query)
        else:
            returnVal = self.session.run( query )
            return returnVal

    def makeLabels(self, labels):

        labelStr = ""
        if len(labels) > 0:
            labelStr = ":" + ":".join(labels)

        return labelStr

    def makeWherePropsMatch(self, selector, props):

        whereProps = []

        if props != None:
            for key in props:

                value = props[key]
                if value == None:
                    continue

                if type(value) == int:
                    valueStr = str(value)
                else:
                    valueStr = "\"" + str(value) + "\""

                whereProps.append(  selector + "." + key + "="+valueStr )

        if len(whereProps) == 0:
            return ""

        return " AND ".join(whereProps)

    def makePropsCreate(self, props):

        propElems = []

        if props != None:
            for key in props:

                value = props[key]
                if value == None:
                    continue

                if type(value) == int:
                    valueStr = str(value)
                else:
                    valueStr = "\"" + str(value) + "\""

                propElems.append( str(key) + ": " + valueStr )

        if len(propElems) == 0:
            return ""

        propStr = "{" + ", ".join(propElems) + "}"
        return propStr

    def relationshipExists(self, selectorLeft, selLeftLabels, propsLeft, selectorRight, selRightLabels, propsRight, labels, props, relname='r'):

        leftLabelStr = self.makeLabels(selLeftLabels)
        rightLabelStr = self.makeLabels(selRightLabels)

        matchPropsLeft = self.makeWherePropsMatch(selectorLeft, propsLeft)
        matchPropsRight = self.makeWherePropsMatch(selectorRight, propsRight)

        relLabelStr = self.makeLabels(labels)
        mathPropsRel = self.makeWherePropsMatch(props)

        # MATCH (n:MIRNA)-[r:MATURE_OF]-(m:MIRNA_PRE) return n,r,m LIMIT 25

        relCreate = "MATCH ({selLeft}{selLeftLabels})-[{relname}{rellabels}]-({selRight}{selRightLabels}) ".format(selLeft=selectorLeft, selRight=selectorRight, selLeftLabels=leftLabelStr, selRightLabels=rightLabelStr, relname=relname, rellabels = relLabelStr)

        if len(matchPropsLeft) > 0 or len(mathPropsRel) > 0 or len(matchPropsRight) > 0:

            middleNode = ""
            if len(matchPropsLeft) > 0 and len(matchPropsRight) > 0:
                middleNode = " AND "

            middleRel = ""
            if len(matchPropsRight) > 0 and len(mathPropsRel) > 0:
                middleRel = " AND "

            whereClause = " WHERE {leftCond}{middleNode}{rightCond}{middleRel}{relCond} ".format(leftCond=matchPropsLeft, middleNode=middleNode, middleRel=middleRel, rightCond=matchPropsRight, relCond=mathPropsRel)
            relCreate += whereClause

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

        matchPropsLeft = self.makeWherePropsMatch(selectorLeft, propsLeft)
        matchPropsRight = self.makeWherePropsMatch(selectorRight, propsRight)

        if len(matchPropsLeft) > 0 or len(matchPropsRight) > 0:

            middle = ""
            if len(matchPropsLeft) > 0 and len(matchPropsRight) > 0:
                middle = " AND "

            whereClause = " WHERE {leftCond}{middle}{rightCond} ".format(leftCond=matchPropsLeft, middle=middle, rightCond=matchPropsRight)
            relCreate += whereClause

        relLabelStr = self.makeLabels(labels)
        relPropStr = self.makePropsCreate(props)

        relCreate += " CREATE ({selLeft})-[{relname}{rellabels} {relProps}]->({selRight}) return r".format(selLeft = selectorLeft, selRight = selectorRight, relname=relname, rellabels=relLabelStr, relProps=relPropStr)

        self.runInDatabase(relCreate)