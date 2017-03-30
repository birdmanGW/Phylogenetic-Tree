from numpy import *
data = loadtxt('dmatrix.txt')
data = data[:400,:400]
ddict = {} 
dist = 0

class Node:
## Class to create and form nodes, storing their NAME, DISTANCE, PARENT / also appends children if called any
    def __init__(self, name, distance, parent):
        self.name = name
        self.distance = distance
        self.children = []  
        self.parent = parent
      
    def addChild(self, child):
        self.children.append(child)
        
    def isLeaf(self):
        return len(self.children) == 0
    
    def newickOutput(self):
        if self.isLeaf():
            return "%s:%f" % (self.name, self.distance)
        else:
            childStrings = [x.newickOutput() for x in self.children]
            inside = ",".join(childStrings)
            return "(%s):%f" % (inside, self.distance)    

    def displayNode(self, depth=0):
        if self.name == None:
            print '\t' * depth,':  {0}'.format(self.distance)
        else:
            print '\t' * depth, '{0}:{1}'.format(self.name,self.distance)       
            
    def displayTree(self, depth=0):
        if self.isLeaf():
            self.displayNode(depth)
        else:
            self.displayNode(depth)
            for child in self.children:
                child.displayTree(depth+1) 
                              
def splitTopComma(input_string):
## Parses newickString to break down node names and distances
    ans = []
    inParent = cumsum([x == ")" for x in input_string]) - cumsum([x == "(" for x in input_string])
    isComma = array([x == "," for x in input_string])
    topCommas = where((inParent == 0) * isComma)[0]
    start = 0
    for x in topCommas:
        ans.append(input_string[start:x])
        start = x+1
    ans.append(input_string[start:])
    return ans
def newickNode(newickString, parent=None):
## Reads newickString, breaks down accordingly calling in splittopCommas, and creates node in Node class
    ## if not a leaf...
    if newickString[0] == '(':
        ## find the last colon
        colonLoc = newickString.rfind(':')
        ## get distance to parent by taking to the end of the string
        ## and convert to number
        parentD = float(newickString[colonLoc+1:])
        newickString1 = newickString[1:colonLoc-1]
        ## make a new empty node with that distance and parent as given by the function
        newNode = Node(None, parentD, parent)
        ## for each substring below it,
        for childString in splitTopComma(newickString1):
            ## make a child node with the same function, but with newNode as parent
            childNode = newickNode(childString, newNode)
            ## then add the child node to our new Node
            newNode.addChild(childNode)
    ## else it's a leaf
    else:
        ## and parse the leaf
        splitString = newickString.split(':')
        nodeName = splitString[0]
        nodeDistance = float(splitString[1])
        ## make a node
        newNode = Node(nodeName, nodeDistance, parent) 
    ## return node
    return newNode

def distance2leaves(node, ans, dist):
## Searches childDescendents for ancestral relation, and records corresponding distances
    ##clear previous ans
    ans = []
    #iterate through nodes children
    for child in node.children:
        ##If child is leaf
        if child.isLeaf():
            ##create new distance with added distance to child's parent, 'node'
            dist2 = dist + child.distance + node.distance
            ##store child in Dict
            childDescendants = [[child.name, dist2]]
        ##Child is 'None'
        else:
            ##append dist with node.distance
            dist = dist + node.distance
            ##Add dist with dist2 to form entire length distance for child to ancestral parent 'node'
            childDescendants = distance2leaves(child, ans, dist)
        ##append children to ans
        ans += childDescendants
    return ans
def storeDescendants(node, descendDict, dist):
## Stores the descendents and adds them to descendDict
        ## if node is a leaf
        if node.isLeaf():
            ## make a list with the singleton of the node.name
            ans = [[node.name, 0.]]
            ## add as entry to dictionary
            descendDict[node] = ans
            ##recreate ans with new child.distance dist
            ans = [[node.name, dist]]
        else:
            ## make empty list
            ans = []
            ## for each child in the node
            for child in node.children:
                ##Assign child.distance to dist
                dist = child.distance 
                ## store the descendants for that child
                childDescendants = storeDescendants(child, descendDict, dist)
                ## and join the list to answer
                ans += childDescendants
            ## put answer in dictionary
            descendDict[node] = ans
            ##create blank dist variable
            dist = 0
            ##create new ans with distance to parents in 'distance2leaves
            ans = distance2leaves(node, ans, dist)
        ## return answer
        return ans    
def calcDist(x, yo, y, dist):
# Distance creation for original dMatrix
    ##Search through parent of column leaf...
    for z in ddict[y.parent]:
        ##If parent descendants contain the row leaf...
        if z[0] == x.name:
            ##Search through parent of column leaf again...
            for q in ddict[y.parent]:
                ##If parent descendants contain the original column leaf...
                if q[0] == yo.name:
                    dist = float(z[1] + q[1])
                    return dist
        else:
            continue ##next z in y.parent descendants
    ##Change y to y's parent
    y = y.parent
    ##Send parent back into calcDist
    dist = calcDist(x, yo, y, dist)
    ##Return final calculated dist for output
    return dist 
def dMatrix(dmatrix):  ## Uses leafNames, if creating matrix move leafNames list above this function**
    # enumeration through row and column does not work here
    for r, x in enumerate(leafNames):
            for c, y in enumerate(leafNames):
                    ##If leaf is finding distance to self...
                    if x == y:
                        dist = 0.
                    else:
                        ##Search through leaves parents for relationships...
                            dist = calcDist(x, y, y, dist)
                    dmatrix[r,c] = dist 
                    
leafNames = ['Z%d' % x for x in range(data.shape[0])]
## ADD NEWICK FORMAT CODE HERE 
#inputTree = ''
#intree = newickNode(inputTree)  
#storeDescendants(intree, ddict, dist)
cherryNodes = []

def initializeNodes(leafNames):
    ## function makes empty trees with the names in the list.
    ans =[]
    for name in leafNames:
        ans.append(Node(name, 0., None))
    return ans                   
def getQMatrix(dMatrix, leafNodes):
## function creates a Qmatrix such that
## Q[i,j] = (n-2)Dij - sumi - sumj
## input dMatrix is a distance matrix
## output: qMatrix as above
    n = len(dMatrix)
    rowsum = sum(dMatrix, axis=0)
    ans = (n-2)*dMatrix - rowsum[:,newaxis] - rowsum[newaxis,:]
    ## Iteration prints here...
    print n
    for i in range(len(leafNodes)):
        for j in range(len(leafNodes)):
            if i == j:
                ans[i,j] = 0
    return ans
def findMinimum(qMatrix, leafNodes):
## find the indicies for the min value in the upper triangle
## input: a matrix,
## output: (i,j) where row i, column j contains minimum element.
    minValue = 10000
    for i in range(len(leafNodes)):
        for j in range(len(leafNodes)):
                if qMatrix[i,j] < minValue:
                    minValue = qMatrix[i,j]
                    f = i
                    g = j
    return f, g   
def computeLeafDistances(dMatrix, f, g):
## job of this function is to find the distance from cherry leaves
## to new parent leaf
    n = len(dMatrix)
    rowsum = sum(dMatrix, axis=0)
    if n-2 == 0:
        dist_f = .5*dMatrix[f,g]
        dist_g = .5*dMatrix[f,g]
    else:
        dist_f = .5*dMatrix[f,g] + (rowsum[f] - rowsum[g]) / (2*(n-2))
        dist_g = dMatrix[f,g] - dist_f
    return dist_f, dist_g
def joinTwoNodes(leafNodes, f, g, dist_f, dist_g):
## job of function is to take two elements from the list of nodes as given by the indicies, and join them as a cherry to a new node.
        f = leafNodes.pop(f)
        g = leafNodes.pop(g-1)
        cherryNode = Node(None, 0, None) ## make empty parent
        cherryNodes.append(cherryNode)
        if f in cherryNodes:
            #f is cherryNode, add child and distance to self
            f.distance = dist_f
            f.parent = cherryNode
            cherryNode.addChild(f)
            if g in cherryNodes:
                #g is cherryNode, add child to self
                g.distance = dist_g
                g.parent = cherryNode
                cherryNode.addChild(g) 
            else:
                #g not cherryNode, create g node
                leafNode_g = Node(g.name, dist_g, cherryNode)
                cherryNode.addChild(leafNode_g) 
        elif g in cherryNodes:
            #g is cherryNode, add child to self
            g.distance = dist_g
            g.parent = cherryNode
            cherryNode.addChild(g)
            #f not cherryNode, create f node
            leafNode_f = Node(f.name, dist_f, cherryNode)
            cherryNode.addChild(leafNode_f)
        else:
            # create f_node
            node_f = Node(f.name, dist_f, cherryNode)
            # create g_node
            node_g = Node(g.name, dist_g, cherryNode)
            # add cherryNode child f
            cherryNode.addChild(node_f)
            # add cherryNode child g
            cherryNode.addChild(node_g)
        leafNodes.append(cherryNode)
def updateDistance(dmatrix, f, g):
## update distanceMatrix
## given a distance matrix and two nodes to join,
## makes the new distance matrix with the new node at the end.
    N = len(dmatrix)
    ## make a new list of indices
    indexList = range(N)
    if N == 2:
        indexList.pop(f)
        indexList.pop(g-1)
    elif g > f:
        indexList.pop(g)
        indexList.pop(f)
    else:
        indexList.pop(f)
        indexList.pop(g)  
    ## find distances for popped rows
    di = dmatrix[f]
    dj = dmatrix[g]
    ## and their pairwise distance
    dij = dmatrix[f,g]
    ## make a new matrix, one smaller
    dmatrix2 = zeros((N-1, N-1))
    ## assign value from old matrix over indexList
    dmatrix2[:N-2, :N-2] =  dmatrix[indexList, :][:, indexList]
    ## compute distance to new node, over indexList
    dux = 0.5*(di + dj - dij)[indexList]
    ## assign to last row and column
    dmatrix2[N-2, :N-2] = dux
    dmatrix2[:N-2, N-2] = dux
    return dmatrix2
    
def neighborJoining(dMatrix, leafNames):
    ## build a list of nodes
    leafNodes = initializeNodes(leafNames)
    ## while our distance matrix is big enough...
    while len(dMatrix) > 1:
        ## make a qmatrix
        qMatrix = getQMatrix(dMatrix, leafNodes)
        ## find upperT minimum
        f, g = findMinimum(qMatrix, leafNodes)
        ## find dist_f, dist_g
        dist_f, dist_g = computeLeafDistances(dMatrix, f, g)
        ## make new parent
        joinTwoNodes(leafNodes, f, g, dist_f, dist_g)
        ## update distances
        dMatrix = updateDistance(dMatrix, f, g)
    return leafNodes[0]  

## Contains the root
Tree = neighborJoining(data, leafNames)  

## Output in newickFormat
outputNewick = Tree.newickOutput()
print "outputTree =", outputNewick
## Output in visual Tree
outputTree = newickNode(outputNewick)
outputTree.displayTree()
