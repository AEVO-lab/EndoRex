import sys
from collections import defaultdict
import argparse


#example command line
#python endorex.py -s "(A, B)r;" -g "((((A__0__gA1, A__1__gA2)w, A__0__gA3)x, A__0__gA4)y, B__1__gB1)z;"


parser = argparse.ArgumentParser(description='EndoRex : Endosymbiotic Reconciliation Software.')
parser.add_argument('--species_tree', '-s', required = True, dest = 'snewick', help = '''
                        Species tree in newick format.Each leaf must have a distinct name.''')
parser.add_argument('--gene_tree', '-g', required = True, dest = 'gnewick', help = '''
                        Gene tree in newick format.  Each leaf must be labeled with 3 parameters 
                        separated by a double underscore __, as follows : 
                        [name]__[species]__[Location]
                        where 
                        [name] is an arbitrary identifier of the gene
                        [species] is the name of the species containing the gene and must be 
                                  a leaf of the species tree, 
                        [location] is either 0 or 1''')
                        
parser.add_argument('--dupcost', '-d', default = 1, dest = 'dupcost', type=float, help = 'Duplication cost')
parser.add_argument('--losscost', '-l', default = 1, dest = 'losscost', type=float, help = 'Loss cost')
parser.add_argument('--transfercost01', '-f0', default = 1, dest = 'transfercost01', type=float, help = 'EGT Transfer cost from 0 to 1')
parser.add_argument('--transfercost10', '-f1', default = 1, dest = 'transfercost10', type=float, help = 'EGT Transfer cost from 1 to 0')
parser.add_argument('--transpocost01', '-p0', default = 1, dest = 'transpocost01', type=float, help = 'EGT Transposition cost from 0 to 1')
parser.add_argument('--transpocost10', '-p1', default = 1, dest = 'transpocost10', type=float, help = 'EGT Transposition cost from 1 to 0')

                        
args = parser.parse_args()


dupcost = args.dupcost 
losscost = args.losscost
transfercosts = [args.transfercost01, args.transfercost10]
transpocosts = [args.transpocost01, args.transpocost10]



class Node:
    
    def __init__(self, content=None, parent=None, children=None, ID=None):
        self.content = content
        self.parent = parent
        self.children = children or []

    
    def __str__(self):
        return f'Cell({str(self.content)})'
    
    def __repr__(self):
        return self.content
    
    def __iter__(self):  # this works, somehow!
        for child in self.children:
            if child:
                for c in child:
                    yield c
        yield self
            
    def add_child(self, child=None):
        if child is None:
            child = Node()
        child.parent=self
        self.children.append(child)
        
    def add_children(self, *children):
        for child in children:
            self.add_child(child)
    
    def is_leaf(self):
        return not bool(self.children)
    
    def is_root(self):
        return self.parent is None
    
    def to_newick(self):
        if self.is_leaf():
            return self.content
        else:
            ls = [node.to_newick() for node in self.children]
            s = ','.join(ls)
            return '(' + s + ')' + self.content + (';' if self.is_root() else '')
            
            


    @staticmethod
    def from_newick(newick):
        stack = []
        i = 0
        while i < len(newick) - 1:
            if newick[i] == '(':  # beginning node
                stack.append(Node())
                i += 1
            elif newick[i] == ')':  # end node
                content = ''
                i += 1
                while newick[i] not in ',();':
                    content += newick[i]
                    i += 1

                stack[-1].content = content.strip()
                if len(stack) > 1:  
                    stack[-2].add_child(stack[-1])
                    stack.pop()
                else:
                    break
            elif newick[i] not in ',();':
                content = ''
                while newick[i] not in ',();':
                    content += newick[i]
                    i += 1
                child = Node(content = content.strip())
                stack[-1].add_child(child)
            else:
                i += 1
        root=stack[0]
        
        return root

    
    def get_clone(self):
        node = Node()
        node.content = self.content 
        
        for c in self.children:
            cclone = c.clone()
            node.add_child(cclone)
        
        return node
        
    #return lca of self and x.  Take time O(h), with h the height of the tree
    def get_lca_with(self, x):
        ancestors = {self}
        y = self
        while not y.is_root():
            y = y.parent
            ancestors.add(y)
        z = x 
        while z not in ancestors:
            z = z.parent
        return z
        
    #counts number of nodes between low_node and up_node, including both
    @staticmethod
    def get_path_length(low_node, up_node):
        cpt = 1
        while low_node != up_node:
            low_node = low_node.parent
            cpt += 1
        return cpt


class DELCostInfo:
    
    #self.D is a 3d dictionary where, for x a node of genetree, bx in {0,1} and evtype in Spe, Dup, Egtf
    #self.D[x, bx, ev]['cost'] is the min cost possible 
    #self.D[x, bx, ev]['origins'] is an array of pairs (bl, br), each pair meaning that it is possible to get an opt solution 
    #                             for x, bx and ev such that children are assigned bl and br, respectively
    
    def __init__(self, gene_tree, species_tree, lca_map, leaf_bx):
        self.D = defaultdict(dict)
        self.gene_tree = gene_tree 
        self.species_tree = species_tree
        self.lca_map = lca_map
        self.leaf_bx = leaf_bx
        
    def set_cost_info(self, x, bx, evtype, cost, origins):
        self.D[x, bx, evtype] = {'cost' : cost, 'origins' : origins}
        
    def get_cost(self, x, bx, evtype):
        return self.D[x, bx, evtype]['cost']
        
    def get_min_cost(self, x, bx):
        if x.is_leaf():
            return self.get_cost(x, bx, 'Spe')
            
        return min(self.get_cost(x, bx, 'Spe'), self.get_cost(x, bx, 'Dup'), self.get_cost(x, bx, 'Egtf'))
        
    def get_best_event(self, x, bx):
        if x.is_leaf():
            return 'Spe'
           
        bev = 'Spe'
        if self.get_cost(x, bx, 'Dup') < self.get_cost(x, bx, bev):
            bev = 'Dup'
        if self.get_cost(x, bx, 'Egtf') < self.get_cost(x, bx, bev):
            bev = 'Egtf'
        return bev
        
    def get_reconciliation(self):
        #performs backtracking on the gene tree
        #TODO : add unary nodes for Egtr
        x = self.gene_tree
        
        bev = {}
        bev[0] = self.get_best_event(x, 0)
        bev[1] = self.get_best_event(x, 1)
            
        
        bx = 0 if self.get_min_cost(x, 0) < self.get_min_cost(x, 1) else 1
        ev = bev[bx]
        
        return self.get_reconciliation_rec(x, bx, ev) 
        
        
    #For backtracking.  We make a copy of self.genetree along the way and return it.
    #Parent knows which bx and event to use on x, here we must figure what to call children on
    def get_reconciliation_rec(self, x, bx, ev):
        node = Node()   #will end up being a copy of x
        
        if x.is_leaf():
            node.content = x.content
            return node
        else:
            node.content = ev + "__" + str(bx)
            
            origin = self.D[x, bx, ev]['origins'][0]    #only return first solution 
            
            b = {}
            bnodes = {}
            
            for i in range(len(x.children)):
                b[i] = origin[i]
                bev = self.get_best_event(x.children[i], b[i])
                bnodes[i] = self.get_reconciliation_rec(x.children[i], b[i], bev)
            
            
            #a bunch of checks for transpositions
            if ev == 'Spe' or ev == 'Dup':
                for i in range(len(x.children)):
                    if b[i] != bx:
                        child_tr = Node(content = 'Egtp')
                        node.add_child(child_tr)
                        child_tr.add_child(bnodes[i])
                    else:
                        node.add_child(bnodes[i])
            else:   #Egtf
                if (b[0] != bx and b[1] != bx) or (b[0] == bx and b[1] == bx):
                    child_tr = Node(content = 'Egtp')
                    node.add_child(child_tr)
                    child_tr.add_child(bnodes[0])
                    node.add_child(bnodes[1])
                else:   #exactly one is different --> no transpo needed
                    node.add_child(bnodes[0])
                    node.add_child(bnodes[1])
                    
            
            return node
        
        
        
    
def parse_genetree_info(genetree, speciestree):
    spmap = {}
    for s in speciestree:
        if s.is_leaf():
            spmap[s.content] = s
    
    lca_map = {}
    leaf_bx = {}
    for x in genetree:
        if x.is_leaf():
            pz = x.content.split('__')
            lca_map[x] = spmap[pz[0]]
            
            leaf_bx[x] = int(pz[1])
        else:
            lca_map[x] = lca_map[x.children[0]].get_lca_with(lca_map[x.children[1]])
            

    return (lca_map, leaf_bx)



def compute_delrecon(gene_tree, species_tree, dupcost, losscost, trfcosts, trpcosts):
    
    (lca_map, leaf_bx) = parse_genetree_info(gene_tree, species_tree)
    
    delinfo_before = DELCostInfo(gene_tree, species_tree, lca_map, leaf_bx)
        
    delinfo = compute_delrecon_rec(gene_tree, species_tree, lca_map, leaf_bx, delinfo_before, dupcost, losscost, trfcosts, trpcosts)
    
    return delinfo.get_reconciliation()
    
                
def compute_delrecon_rec(x, species_tree, lca_map, leaf_bx, delinfo, dupcost, losscost, trfcosts, trpcosts):
    #ASSUMPTION : trpcosts has two indices [0] and [1], one for each possible switch source
    #ASSUMPTION : trpcosts[i] <= trfcosts[i] + losscost (ie we never replace a transposition by a transfer + loss)
    if x.is_leaf():
        #leaves must have the given b mapping, are speciation and have no origins
        delinfo.set_cost_info(x, leaf_bx[x], 'Spe', 0, [])
        delinfo.set_cost_info(x, 1 - leaf_bx[x], 'Spe', float('inf'), [])
    else:
        l_x = 0
        for c in x.children:
            compute_delrecon_rec(c, species_tree, lca_map, leaf_bx, delinfo, dupcost, losscost, trfcosts, trpcosts)
            l_x += Node.get_path_length(lca_map[c], lca_map[x])
        
        
        
        canBeSpec = True
        for c in x.children:
            if lca_map[x] == lca_map[c]:
                canBeSpec = False
        
        
        for bx in [0, 1]:
            cost_spe = losscost * (l_x - 4)
            cost_dup = dupcost + losscost * (l_x - 2)
            cost_trf = trfcosts[bx] + losscost * (l_x - 2)
            
            origins_spe_dup = {}
            
            #this just implements recurrence as seen in paper
            #spec and dup are essentially the same.  Origin remembers bx of children for backtracking.
            for i in range(len(x.children)):
                c = x.children[i]
                mincost_c_bx = delinfo.get_min_cost(c, bx)
                mincost_c_antibx = delinfo.get_min_cost(c, 1 - bx)
                
                if mincost_c_bx <= trpcosts[bx] + mincost_c_antibx:
                    cost_spe += mincost_c_bx
                    cost_dup += mincost_c_bx
                    origins_spe_dup[i] = [bx]
                    if mincost_c_bx == trpcosts[bx] + mincost_c_antibx: #if equality, two origins are possible
                        origins_spe_dup[i].append(1 - bx)
                else:
                    cost_spe += trpcosts[bx] + mincost_c_antibx
                    cost_dup += trpcosts[bx] + mincost_c_antibx
                    origins_spe_dup[i] = [1 - bx]

            #all combinations of origins are possible at the children
            origin_combos_spe_dup = [ [x, y] for x in origins_spe_dup[0] for y in origins_spe_dup[1]]
            
                
            if not canBeSpec:
                cost_spe = float('inf')
                
            delinfo.set_cost_info(x, bx, 'Spe', cost_spe, origin_combos_spe_dup)
            delinfo.set_cost_info(x, bx, 'Dup', cost_dup, origin_combos_spe_dup)
                
            #now EGTF nodes...
            cost_egtf = trfcosts[bx] + losscost * (l_x - 2)
            
            c1 = x.children[0]
            c2 = x.children[1]
            mincost_c1_bx = delinfo.get_min_cost(c1, bx)
            mincost_c2_bx = delinfo.get_min_cost(c2, bx)
            mincost_c1_antibx = delinfo.get_min_cost(c1, 1 - bx)
            mincost_c2_antibx = delinfo.get_min_cost(c2, 1 - bx)
            
            min_auxcost = min(mincost_c1_bx + mincost_c2_antibx, \
                              mincost_c1_antibx + mincost_c2_bx, \
                              trpcosts[1 - bx] + mincost_c1_bx + mincost_c2_bx, \
                              trpcosts[bx] + mincost_c1_antibx + mincost_c2_antibx)
                              
            
            cost_egtf += min_auxcost
            
            origin_combos_etrf = []
            
            if mincost_c1_bx + mincost_c2_antibx == min_auxcost:
                origin_combos_etrf.append( [bx, 1 - bx] )
            if mincost_c1_antibx + mincost_c2_bx == min_auxcost:
                origin_combos_etrf.append( [1 - bx, bx] )
            if trpcosts[1 - bx] + mincost_c1_bx + mincost_c2_bx == min_auxcost:
                origin_combos_etrf.append( [bx, bx] )
            if trpcosts[bx] + mincost_c1_antibx + mincost_c2_antibx == min_auxcost:
                origin_combos_etrf.append( [1 - bx, 1 - bx] )
            
            delinfo.set_cost_info(x, bx, 'Egtf', cost_egtf, origin_combos_etrf)
                
    return delinfo        
    
    
#snewick = '(A, B)r;'
#gnewick = '((((A__0__gA1, A__1__gA2)w, A__0__gA3)x, A__0__gA4)y, B__1__bB1)z;'



species_tree = Node.from_newick(args.snewick)
gene_tree = Node.from_newick(args.gnewick)
   

reconciled_tree = compute_delrecon(gene_tree, species_tree, dupcost, losscost, transfercosts, transpocosts)

print(reconciled_tree.to_newick())
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
