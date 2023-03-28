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
