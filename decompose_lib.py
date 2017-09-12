# uym2 addeddecompose_by_diameter(
# June 2017
# utils for tree decomposition


from dendropy import Tree, Node
try:
    from queue import Queue # python 3
except:
    from Queue import Queue # python 2
#from tree import PhylogeneticTree
from sepp import get_logger

_LOG = get_logger(__name__)

def decompose_by_diameter(a_tree,strategy,max_size=None,min_size=None,max_diam=None):
    def __ini_record__():
        for node in a_tree.postorder_node_iter():
               node.marked = False
        a_tree.seed_node.marked = True
        for node in a_tree.postorder_node_iter():
               __updateNode__(node)
    
    #def __find_midpoint_edge__(t):
    def __find_midpoint_edge__(seed_node):
        u = seed_node.bestLCA.anchor
        d = 0
        while (d+u.edge_length < seed_node.diameter/2):
            d += u.edge_length
            u = u.parent_node
        return u.edge
    
    #def __find_centroid_edge__(t):
    def __find_centroid_edge__(seed_node):
        u = seed_node
        product = -1
        acc_nleaf = 0

        while u and not u.is_leaf():
            max_child = None
            max_child_nleaf = 0
            for ch in u.child_node_iter():
                if ch.marked:
                    continue
                if ch.nleaf >= max_child_nleaf:
                    max_child_nleaf = ch.nleaf
                    max_child = ch
            acc_nleaf += (u.nleaf-max_child.nleaf)
            new_product = max_child.nleaf * acc_nleaf
            if new_product < product:
                break
            product = new_product
            u = max_child

        return u.edge

    def __bisect__(t,e):
#        e = __find_centroid_edge__(t)
        
        u = e.tail_node
        v = e.head_node

        u.remove_child(v)
        t1 = Tree(seed_node = v)

        if u.num_child_nodes() == 1:
            p = u.parent_node
            v = u.child_nodes()[0]
            l_v = v.edge_length
            u.remove_child(v)
            if p is None: # u is the seed_node; this means the tree runs out of all but one side
                t.seed_node = v
                return t,t1
            l_u = u.edge_length
            p.remove_child(u)
            p.add_child(v)
            v.edge_length = l_u+l_v
            u = p

        while u is not None:
            __updateNode__(u)
            u = u.parent_node

        return t,t1

    def __clean_up__(t):
        for node in t.postorder_node_iter():
            delattr(node,"nleaf")
            delattr(node,"anchor")
#            delattr(node,"maxheight")
            delattr(node,"maxdepth")
            delattr(node,"diameter")
#            delattr(node,"topo_diam")
            delattr(node,"bestLCA")

    def __updateNode__(node):
        if node.is_leaf():
            node.anchor = node
#            node.maxheight = 0
            node.maxdepth = 0
            node.diameter = 0
#            node.topo_diam = 0
            node.bestLCA = node
            node.nleaf = 1
            return

#        n1 = -1
#        n2 = -1
        d1 = -1
        d2 = -1
        anchor1 = None
        anchor2 = None
        node.diameter = 0
#        node.topo_diam = 0
        node.bestLCA = None
        node.nleaf = 0

        for ch in node.child_node_iter():
               if ch.marked:
                   continue 
               node.nleaf += ch.nleaf
#               n = ch.maxheight + 1
               d = ch.maxdepth + ch.edge_length
#               if n > n1:
#                   n2 = n1
#                   n1 = n
#                   anchor2 = anchor1
#                   anchor1 = ch.anchor
#               elif n > n2:
#                   n2 = n
#                   anchor2 = ch.anchor
               if d > d1:
                   d2 = d1
                   d1 = d
                   anchor2 = anchor1
                   anchor1 = ch.anchor
               elif d > d2:
                   d2 = d
                   anchor2 = ch.anchor
               if ch.diameter > node.diameter:
                   node.diameter = ch.diameter
                   node.bestLCA = ch.bestLCA
#               node.diameter = max(ch.diameter,node.diameter)

#        node.diameter = max(d1+d2, node.diameter)
        node.maxdepth = max(d1,0)
#        node.maxheight = n1
        node.anchor = anchor1
        if d1+d2 > node.diameter:
            node.diameter = d1+d2
            node.bestLCA = node

    #def __get_breaking_edge__(t,edge_type):
    def __get_breaking_edge__(seed_node,edge_type):
        if seed_node.nleaf <= max_size and seed_node.diameter <= max_diam:
            return None
        if edge_type == 'midpoint':
            e = __find_midpoint_edge__(seed_node)
        elif edge_type == 'centroid':
            e = __find_centroid_edge__(seed_node)
        else:
            _LOG.warning("Invalid decomposition type! Please use either 'midpoint' or 'centroid'")
            return None

        n = e.head_node.nleaf
        if (n < min_size) or (seed_node.nleaf - n) < min_size:
            return None
        return e

    def __check_stop__(t):
        return ( (t.seed_node.nleaf <= max_size and t.seed_node.diameter <= max_diam) or
                 (t.seed_node.nleaf//2 < min_size) )     

    #def __break_by_MP_centroid__(t):
    def __break_by_MP_centroid__(seed_node):
        e = __get_breaking_edge__(seed_node,'midpoint')
        if e is None:
#            print("Midpoint failed. Trying centroid decomposition...")
            e = __get_breaking_edge__(seed_node,'centroid')
#        else:
#            print("Successfully splitted by midpoint")
        return e

    #def __break(t):
    def __break(seed_node):
        if strategy == "centroid":
            return __get_breaking_edge__(seed_node,'centroid')
        elif strategy == "midpoint":
            return __break_by_MP_centroid__(seed_node)
        else:
            raise Exception("strategy not valid: %s" %strategy)

    _LOG.debug("Starting brlen decomposition ...")
    #tqueue = Queue()
    tstk = []
    __ini_record__()
    min_size = min_size if min_size else 0
    max_size = max_size if max_size else a_tree.seed_node.nleaf
    max_diam = max_diam if max_diam else a_tree.seed_node.diameter


    r = 0
    

    #e = __break(a_tree)
    e = __break(a_tree.seed_node)

    if e is None:
        #__clean_up__(a_tree)
        #return [(a_tree,"r0d1"]
        return [(a_tree.seed_node,"r0d1")]
        
    treeMap = {} 
    e.head_node.marked = True
    u = e.head_node.parent_node
    #print(a_tree.seed_node.nleaf)
    
    while u is not None:
        __updateNode__(u)
        if u.marked:
            break
        u = u.parent_node
    
    #print(a_tree.seed_node.nleaf)
    #print(e.head_node.nleaf)
     
    #t1,t2 = __bisect__(a_tree,e)
    e1 = __break(a_tree.seed_node)
    e2 = __break(e.head_node)
    
    #print(e1.head_node is a_tree.seed_node)
     

    if e1 is None:
         #__clean_up__(t1)
         #treeMap.append((t1,'r'+str(r)+'d1'))
         #treeMap.append((a_tree.seed_node,'r'+str(r)+'d1'))
         name = 'r' + str(r) + 'd1'
         a_tree.seed_node.groupName = name
         treeMap[name] = a_tree.seed_node
    else:
        tstk.append((a_tree.seed_node,e1))
    
    if e2 is None:
         #__clean_up__(t2)            
         #treeMap.append((t2,'r'+str(r)+'d2'))
         #treeMap.append((e.head_node,'r'+str(r)+'d2'))
         name = 'r' + str(r) + 'd2'
         e.head_node.groupName = name
         treeMap[name] = e.head_node
    else:
        tstk.append((e.head_node,e2))
    r = r+1
    
    #tstk.append((a_tree,e))
    
    while len(tstk):
        #t,e = tstk.pop()
        node,e = tstk.pop()
        #print(node is a_tree.seed_node)
        e.head_node.marked = True
        u = e.tail_node
        while u is not None:
            __updateNode__(u)
            if u.marked:
                break
            u = u.parent_node
         
        #t1,t2 = __bisect__(t,e)
        e1 = __break(node)
        e2 = __break(e.head_node)
        if e2 is None:
             #__clean_up__(t2)            
             #treeMap.append((e.head_node,'r'+str(r)+'d1'))
             name = 'r' + str(r) + 'd1'
             e.head_node.name = name
             treeMap[name] = e.head_node
        #else:
        #    tstk.append((t2,e2))
        if e1 is None:
             #__clean_up__(t1)
             #treeMap.append((node,'r'+str(r)+'d2'))
             name = 'r' + str(r) + 'd2'
             node.name = name
             treeMap[name] = node
        #else:
        #    tstk.append((t1,e1))
        
        if e1 is not None:
            tstk.append((node,e1))
        if e2 is not None:
            tstk.append((e.head_node,e2))
            
        r = r+1
    
    return treeMap

def place_group_onto_tree(a_tree,grouping):
    treeMap = {}
    for node in a_tree.postorder_node_iter():
        if node.is_leaf():
            node.name = grouping[node.taxon.label]
            node.nleaf = 1
        else:
            node.name = None
            node.nleaf = 0
            children_names = {}
            for ch in node.child_node_iter():
                if not ch.marked and ch.name is not None::
                      children_names[ch.name] = [ch] if ch.name not in children_names else children_names[ch.name].append(ch)
            if len(children_names.keys()) == 1:
                node.name = children_names.keys()[0]
            else:
                for name in children_names:
                    if len(children_names[name]) > 1:
                        node.name = name
                        node.marked = True
                        treeMap[name] = node
                    else:
                        ch = children_names[name][0]
                        ch.marked = True 
                        treeMap[name] = ch
        # compute nleaf
           for ch in node.child_node_iter():
               if not ch.marked:
                   node.nleaf += ch.nleaf
    treeMap[a_tree.seed_node.name] = a_tree.seed_node
    return treeMap                   


def compute_group_distance_matrix(a_tree,treeMap):
# This function should be call AFTER the tree was decomposed and annotated (see decompose_by_diameter(...))
    def __compute_sumIn():
    # Assume that a_tree is already decomposed and annotated with sub-groups
    # Compute the sum of the distances to the leaves inside each group to its marked (root) node 
        for node in a_tree.postorder_node_iter():
            node.sumIn = 0
            for ch in node.child_node_iter():
                if ch.marked:
                    continue
                node.sumIn += (ch.nleaf*ch.edge_length + ch.sumIn)

    def __compute_ingroup_distance(node):
    # Compute the average distance from a node to all the leaves in its group
    # Assume that sumIn was computed for each of the nodes
        v = node
        SUM = node.sumIn
        acc_edge = 0
        
        while not v.marked:
            u = v.parent_node
            acc_edge += v.edge_length
            s = 0
            for ch in u.child_node_iter():
                if not ch.marked and ch is not v:
                    s += (ch.sumIn + (ch.edge_length + acc_edge)*ch.nleaf)
            SUM += s        
            v = u

        return SUM/v.nleaf

    def __preprocess():
        __compute_sumIn()
        for name in treeMap:
            v = treeMap[name]
            u = v.parent_node
                
            if u is not None:
                u.distance = __compute_ingroup_distance(u)
            v.distance = __compute_ingroup_distance(v)

    def __compute_group_distance(A,B):
        if A == B:
            return 0

        nodeA = treeMap[A]
        nodeB = treeMap[B]

        u = nodeA.parent_node
        lastGroup = A
        lastDistance = 0
        d  = nodeA.edge_length
        mapping = {u:d}

        while u is not None:
            if u is nodeB:
                v = treeMap[lastGroup]
                return nodeA.distance + lastDistance + v.edge_length +  v.parent_node.distance
            mapping[u] = d
            if u.marked:
                lastGroup = u.name
                lastDistance = d
            d += u.edge_length if u.edge_length else 0
            u = u.parent_node

        u = nodeB.parent_node
        lastGroup = B
        lastDistance = 0
        d = nodeB.edge_length

        while u is not None:
            if u is nodeA:
                v = treeMap[lastGroup]
                return nodeB.distance + lastDistance + v.edge_length + v.parent_node.distance
            if u in mapping:
                return nodeA.distance + nodeB.distance + mapping[u] + d
            if u.marked:
                lastGroup = u.name
                lastDistance = d    
            d += u.edge_length
            u = u.parent_node
            
        return -1 # should not happen unless something went wrong
    
    __preprocess()   
    D = {}
    # Main part: compute distances for all pairs of groups
    groups = treeMap.keys()
    for i,A in enumerate(groups):
       for B in groups[i+1:]:
           D[(A,B)] = __compute_group_distance(A,B)
           
    return D  
