import numpy as np

def get_index(length,x):
        index_y=x%length
        index_x=(x-index_y)/length
        return int(index_x),index_y

def built_tree(edge_set,start,tree_set):
        rep_edge_set=[[] for c in range(len(edge_set))]
        for c in range(len(edge_set)):
            for d in range(len(edge_set[c])):
                rep_edge_set[c].append(edge_set[c][d])
        for c in rep_edge_set[start]:
                tree_set.append([start,c])
                idx=rep_edge_set[c].index(start)
                print(rep_edge_set[c][idx],c,idx)
                print(edge_set)                
                del rep_edge_set[c][idx]
                print(edge_set)
                if len(rep_edge_set[c])>0:
                        tree_set=built_tree(rep_edge_set,c,tree_set)
        return tree_set



def deter_circle(edge_mapping,index_x,index_y):
        if edge_mapping[index_x]==[]:
                return True
        else:
                tree_set=built_tree(edge_mapping,index_x,[])
                node_set=[]
                for [x,y] in tree_set:
                        node_set.append(x)
                        node_set.append(y)
        if index_y not in node_set:
                return True
        else:
                return False





def min_spanning_tree(nodes_name,d_array):
        #d_array is a list of distance matrix
        #set distance from itself to itself as a big number, say 999
        map_to_name=[]
        for c in nodes_name:
            map_to_name.append(c)
        for c in range(len(d_array)):
            if d_array[c]==0:
                d_array[c]=999
        array_sorted=np.argsort(d_array)
        print(array_sorted)
        idx=0
        edge_mapping=[[] for i in range(len(nodes_name))]
        num_of_edge=1
        edge_list=[]
        edge_list_index=[]
        while num_of_edge<len(nodes_name):
                index_x,index_y=get_index(len(nodes_name),array_sorted[idx])
                if index_x < index_y:
                    if deter_circle(edge_mapping,index_x,index_y)==True:
                        edge_mapping[index_x].append(index_y)
                        edge_mapping[index_y].append(index_x)
                        edge_list.append([map_to_name[index_x],map_to_name[index_y]])
                        edge_list_index.append([index_x,index_y])
                        num_of_edge=num_of_edge+1
                        print(index_x,index_y)
                        print(edge_mapping)                        
                idx=idx+1
        sum_of_distance=0
        for [x,y] in edge_list_index:
            index=x*len(nodes_name)+y
            sum_of_distance=sum_of_distance+d_array[index]
        return edge_list,sum_of_distance

