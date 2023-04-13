import pandas as pd
from itertools import combinations
import multiprocessing as mp
class ppiDataset:

    def __init__(self, filename,proteinLabels:list=[] ,**readerKwargs):

        self.data: pd.DataFrame = pd.read_csv(filename, compression='gzip', **readerKwargs)
        self.proteinLabels = proteinLabels
        self.ppis = set()
    
    def getPPIs(self, isCorum: bool = False) -> set:
        
        data = self.data.copy()

        if isCorum:
            
            def combinationsOfProteins(complx):
                
                if ';' in list(complx['subunits(Gene name)']) :
                    complx['proteinTuple'] = list(combinations(
                        complx['subunits(Gene name)'].split(';'), 2))
                    
                    return complx

            data  = data.apply(lambda complx: combinationsOfProteins(complx), axis=1)
            ppiList = list(data.dropna()['proteinTuple'])
            ppiSet = {item for sublist in ppiList for item in sublist}
            

        else:
            data['proteinTuple'] = list(zip(data[self.proteinLabels[0]], data[self.proteinLabels[1]])) 
            ppiSet = set(data['proteinTuple'])

        self.ppis = ppiSet

        return ppiSet


#Deprecated
class TreeNode:

    def __init__(self, value, children:set=set()):

        self.value = value
        self.children = children
    
    def addChild(self, childNode):
        
        self.children.add(childNode)
    
    def removeChild(self, childNode):

        self.children = [child for child in self.children if child is not childNode]

    def transverse(self):

        nodesToVisit=[self]

        while len(nodesToVisit) > 0:

            currentNode=nodesToVisit.pop()
            print(currentNode.value)
            nodesToVisit += list(currentNode.children) # convert the set into a list so that the trasnverse works 
        
    def getNode(self, nodeValue):
        
        if self.value == nodeValue:
            return self
        
       
        elif(len(self.children)):
            

            for node in self.children:

                hasNode = node.getNode(nodeValue)
                if hasNode:
                    return hasNode
                
        return None
    def getChildrenValue(self):
        return {child.value for child in self.children}
    
    def getNodeFirstLayer(self, nodeValue):

        for child in self.children:

            if child.value==nodeValue:

                return child
