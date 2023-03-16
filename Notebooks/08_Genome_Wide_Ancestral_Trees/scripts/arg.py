
from math import exp, log
#import networkx.algorithms.non_randomness

import json

import random
# print("Setting random seed")
# random.seed(3)
# numpy.random.seed(3)

from random import random, shuffle, choices
from functools import partial
from copy import deepcopy, copy
# import networkx as nx
from functools import reduce
from numpy.random import exponential, sample
import numpy as np
from time import sleep
# import matplotlib.pyplot as plt
from collections import defaultdict


class Lineage(object):
    """
    Lineages are the edges of our graph
    """
    def __init__(self, lineageid=None, down=None, up=None, intervals=None):
        self.lineageid = lineageid
        self.down = down # node at bottom of edte
        self.up = up # node at top of edge
        self.intervals = intervals

    def __hash__(self):
        return self.lineageid

    def __eq__(self, other):
        return hasattr(other, 'lineageid') and self.lineageid == other.lineageid  

    def __repr__(self):
        return f'{self.lineageid}:{self.intervals}'

    def get_dict(self):
        d = self.__dict__.copy()
        if d['up']:
            d['up'] = d['up'].nodeid
        d['down'] = d['down'].nodeid
        return d

    def toJSON(self):
        return json.dumps(self, default=self.get_dict, sort_keys=True, indent=4)

    def __deepcopy__(self, memo):
        if self in memo:
            return memo.get(self)
        dup = type(self)(
            lineageid=deepcopy(self.lineageid, memo),
            # down=deepcopy(self.down, memo),
            # up=deepcopy(self.up, memo),
            intervals=deepcopy(self.intervals, memo)
            )
        memo[self] = dup
        dup.down=deepcopy(self.down, memo)
        dup.up=deepcopy(self.up, memo)
        return dup

class Node():
    """
    Leaf, Coalscent and Recombination are the nodes of the graph
    """
    def __hash__(self):
        return self.nodeid

    def __eq__(self, other):
        return hasattr(other, 'nodeid') and self.nodeid == other.nodeid  

    def __repr__(self):
        return f'{self.nodeid}'

class Leaf(Node):

    def __init__(self, nodeid=None, height=None, parent=None, intervals=[(0, 1)], xpos=None):
        self.nodeid = nodeid
        self.height = height
        self.intervals = intervals
        self.parent = parent
        self.xpos = xpos

    def get_dict(self):
        d = self.__dict__.copy()
        d['parent'] = d['parent'].lineageid
        return d

    def toJSON(self):
        return json.dumps(self, default=self.get_dict, sort_keys=True, indent=4)

    def __deepcopy__(self, memo):
        if self in memo:
            return memo.get(self)
        dup = type(self)(
            nodeid=deepcopy(self.nodeid, memo),
            height=deepcopy(self.height, memo),
            # parent=deepcopy(self.parent, memo),
            intervals=deepcopy(self.intervals, memo),
            xpos=deepcopy(self.xpos, memo)
            )
        memo[self] = dup
        dup.parent = deepcopy(self.parent, memo)
        return dup            

class Coalescent(Node):

    def __init__(self, nodeid=None, height=None, children=None, parent=None, xpos=None):
        self.nodeid = nodeid
        self.height = height
        self.children = children
        self.parent = parent
        self.xpos = xpos

    def get_dict(self):
        d = self.__dict__.copy()
        d['children'] = [c.lineageid for c in d['children']]
        d['parent'] = d['parent'].lineageid
        return d

    def toJSON(self):
        return json.dumps(self, default=self.get_dict, sort_keys=True, indent=4)

    def __deepcopy__(self, memo):
        if self in memo:
            return memo.get(self)
        dup = type(self)(
            nodeid=deepcopy(self.nodeid, memo),
            height=deepcopy(self.height, memo),
            # children=deepcopy(self.children, memo),
            # parent=deepcopy(self.parent, memo),
            xpos=deepcopy(self.xpos, memo)
            )
        memo[self] = dup
        dup.children=[deepcopy(c, memo) for c in self.children]
        dup.parent=deepcopy(self.parent, memo)
        return dup            

class Recombination(Node):

    def __init__(self, nodeid=None, height=None, child=None, recomb_point=None, left_parent=None, right_parent=None, xpos=None):
        self.nodeid = nodeid
        self.height = height
        self.child = child
        self.left_parent = left_parent
        self.right_parent = right_parent
        self.recomb_point = recomb_point
        self.xpos = xpos

    def get_dict(self):
        d = self.__dict__.copy()
        d['left_parent'] = d['left_parent'].lineageid
        d['right_parent'] = d['right_parent'].lineageid
        d['child'] = d['child'].lineageid
        return d

    def toJSON(self):
        return json.dumps(self, default=self.get_dict, sort_keys=True, indent=4)

    def __deepcopy__(self, memo):
        if self in memo:
            return memo.get(self)
        dup = type(self)(
            nodeid=deepcopy(self.nodeid, memo),
            height=deepcopy(self.height, memo),
            # child=deepcopy(self.child, memo),
            recomb_point=deepcopy(self.recomb_point, memo),
            # left_parent=deepcopy(self.left_parent, memo),
            # right_parent=deepcopy(self.right_parent, memo),
            xpos=deepcopy(self.xpos, memo)
            )  
        memo[self] = dup
        dup.child=deepcopy(self.child, memo)
        dup.left_parent=deepcopy(self.left_parent, memo)
        dup.right_parent=deepcopy(self.right_parent, memo)
        return dup                  

def flatten(list_of_tps):
 
    return reduce(lambda ls, ival: ls + list(ival), list_of_tps, [])


def unflatten(list_of_endpoints):

    return [ [list_of_endpoints[i], list_of_endpoints[i + 1]]
          for i in range(0, len(list_of_endpoints) - 1, 2)]

def merge(query, annot, op):

    a_endpoints = flatten(query)
    b_endpoints = flatten(annot)

    assert a_endpoints == sorted(a_endpoints), "not sorted or non-overlaping"
    assert b_endpoints == sorted(b_endpoints), "not sorted or non-overlaping"


    sentinel = max(a_endpoints[-1], b_endpoints[-1]) + 1
    a_endpoints += [sentinel]
    b_endpoints += [sentinel]

    a_index = 0
    b_index = 0

    res = []

    scan = min(a_endpoints[0], b_endpoints[0])
    while scan < sentinel:
        in_a = not ((scan < a_endpoints[a_index]) ^ (a_index % 2))
        in_b = not ((scan < b_endpoints[b_index]) ^ (b_index % 2))
        in_res = op(in_a, in_b)

        if in_res ^ (len(res) % 2):
            res += [scan]
        if scan == a_endpoints[a_index]: 
            a_index += 1
        if scan == b_endpoints[b_index]: 
            b_index += 1
        scan = min(a_endpoints[a_index], b_endpoints[b_index])

    return unflatten(res)

def interval_diff(a, b):
    if not (a and b):
        return a and a or b
    return merge(a, b, lambda in_a, in_b: in_a and not in_b)

def interval_union(a, b):
    if not (a and b):
        return []
    return merge(a, b, lambda in_a, in_b: in_a or in_b)

def interval_intersect(a, b):
    if not (a and b):
        return []
    return merge(a, b, lambda in_a, in_b: in_a and in_b)
    
def interval_sum(intervals):
    return sum(e - s for (s, e) in intervals)

def interval_span(intervals):
    starts, ends = zip(*intervals)
    return max(ends) - min(starts)

def interval_split(intervals, pos):
    left, right = list(), list()
    for s, e in intervals:
        if pos >= e:
            left.append((s, e))
        elif pos >= s and pos < e:
            left.append((s, pos))
            right.append((pos, e))
        else:
            right.append((s, e))
    return left, right

def interval_any_shared_borders(a, b):
    flat = flatten(a) + flatten(b)
    return len(flat) > len(set(flat))

def x_positions_traverse(node, offset):
    """
    Recursive function for adding x positions
    """
    offset -= 1
    if type(node) is Leaf:
        return 1 # offset
        # do nothing
    if type(node) is Recombination:
        # add xpos to children
        node.child.down.xpos = node.xpos
        offset = x_positions_traverse(node.child.down, offset)
        # add xpos to node if both parents are positioned
        if node.left_parent.up.xpos is not None and node.right_parent.up.xpos is not None:
            left, right = sorted([node.left_parent.up.xpos, node.right_parent.up.xpos])
            x = left + (right - left) / 2 
            node.xpos = x
    if type(node) is Coalescent:
        # add xpos to children
        node.children[0].down.xpos = node.xpos - offset
        node.children[1].down.xpos = node.xpos + offset
        for child in node.children:
            x_positions_traverse(child.down, offset)

####################################################

def segments_crossing(line1, line2):
    xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
    ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)
    if div == 0:
        return False
    #    raise Exception('lines do not intersect')

    d = (det(*line1), det(*line2))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div

    # return x, y

    return \
        min(line1[0][0], line1[1][0]) < x < max(line1[0][0], line1[1][0]) and \
        min(line2[0][0], line2[1][0]) < x < max(line2[0][0], line2[1][0]) and \
        min(line1[0][1], line1[1][1]) < y < max(line1[0][1], line1[1][1]) and \
        min(line2[0][1], line2[1][1]) < y < max(line2[0][1], line2[1][1])


def crossing(a, b):
    return segments_crossing(
        ((a.down.xpos, a.down.height), (a.up.xpos, a.up.height)), 
        ((b.down.xpos, b.down.height), (b.up.xpos, b.up.height))
        )

def get_all_crossing_pairs(lineages):
    crossing_pairs = []
    for i in range(len(lineages)):
        for j in range(i+1, len(lineages)):
            if crossing(lineages[i], lineages[j]):
                crossing_pairs.append([i, j])
    return crossing_pairs

def reduce_crossovers(nodes):

    lineages = get_parent_lineages(nodes)
    lineages = [x for x in lineages if x.up is not None]
    crossing_pairs = get_all_crossing_pairs(lineages)

    for i, j in crossing_pairs:
        # swap downs
        lineages[i].down.xpos, lineages[j].down.xpos = lineages[j].down.xpos, lineages[i].down.xpos
        new_crossing_pairs = get_all_crossing_pairs(lineages)
        if not len(new_crossing_pairs) < len(crossing_pairs):
            # swap downs back
            lineages[i].down.xpos, lineages[j].down.xpos = lineages[j].down.xpos, lineages[i].down.xpos
        else:
            crossing_pairs = new_crossing_pairs

        # swap ups
        lineages[i].up.xpos, lineages[j].up.xpos = lineages[j].up.xpos, lineages[i].up.xpos
        new_crossing_pairs = get_all_crossing_pairs(lineages)
        if not len(new_crossing_pairs) < len(crossing_pairs):
            # swap ups back
            lineages[i].up.xpos, lineages[j].up.xpos = lineages[j].up.xpos, lineages[i].up.xpos
        else:
            crossing_pairs = new_crossing_pairs


        leaf_idx = [i for i, n in enumerate(nodes) if type(n) is Leaf]

        if type(lineages[i].down) is Leaf and type(lineages[j].down) is not Leaf:
            # try to find new position for leaf i

            for k in leaf_idx:
                # swap
                lineages[i].down.xpos, lineages[k].down.xpos = lineages[k].down.xpos, lineages[i].down.xpos

                new_crossing_pairs = get_all_crossing_pairs(lineages)
                if not len(new_crossing_pairs) < len(crossing_pairs):
                    # swap back
                    lineages[i].down.xpos, lineages[k].down.xpos = lineages[k].down.xpos, lineages[i].down.xpos
                else:
                    break


        if type(lineages[i].down) is Leaf and type(lineages[j].down) is not Leaf:
            # try to find new position for leaf j
            for k in leaf_idx:
                # swap
                lineages[j].down.xpos, lineages[k].down.xpos = lineages[k].down.xpos, lineages[j].down.xpos

                new_crossing_pairs = get_all_crossing_pairs(lineages)
                if not len(new_crossing_pairs) < len(crossing_pairs):
                    # swap back
                    lineages[j].down.xpos, lineages[k].down.xpos = lineages[k].down.xpos, lineages[j].down.xpos
                else:
                    break


    # TODO: if cross over involves an internal and an external branch, see if
    # you can move the leaf to a new position that removes the crossover



    # leaf_idx = [i for i, n in enumerate(nodes) if type(n) is Leaf]
    # suffled_leaf_idx = leaf_idx[:]
    # shuffle(suffled_leaf_idx)
    # for i, j in sorted(zip(leaf_idx, suffled_leaf_idx)):
    #     print(i, j)
    #     nodes[i].xpos = nodes[j].xpos

    # new_crossing_pairs = get_all_crossing_pairs(lineages)
    # if not len(new_crossing_pairs) < len(crossing_pairs):
    #     for i, j in sorted(zip(leaf_idx, suffled_leaf_idx), reverse=True):
    #         nodes[i].xpos = nodes[j].xpos




def branch_length(lineage):
    return ((lineage.down.xpos - lineage.up.xpos)**2 + (lineage.down.height - lineage.up.height)**2)**0.5

def reduce_total_branch_length(nodes, shift):
    lineages = get_parent_lineages(nodes)
    lineages = [x for x in lineages if x.up is not None]

    cur_tot = sum(branch_length(l) for l in lineages)
    for i in range(len(lineages)):
        orig_xpos = lineages[i].up.xpos
        lineages[i].up.xpos += shift
        tot = sum(branch_length(l) for l in lineages)
        if tot < cur_tot:
            cur_tot = tot
        else:
            lineages[i].up.xpos = orig_xpos
        lineages[i].up.xpos -= shift
        tot = sum(branch_length(l) for l in lineages)
        if tot < cur_tot:
            cur_tot = tot
        else:
            lineages[i].up.xpos = orig_xpos

def redistribute_leaves(nodes):

    leaves = [n for n in nodes if type(n) is Leaf]
    # min_x = min(n.xpos for n in leaves)
    # max_x = max(n.xpos for n in leaves)
    leaves = sorted(leaves, key=lambda x: x.xpos)
    min_x = 0
    max_x = 1
    for xpos, leaf in zip(np.linspace(min_x, max_x, len(leaves)), leaves):
        leaf.xpos = xpos

####################################################

def add_node_x_positions(nodes):
    """
    Adds x positions in place
    """
    offset = len(nodes)-1
    nodes[-1].xpos = offset
    x_positions_traverse(nodes[-1], offset)

    shift_list = [0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.02]

    rescale_positions(nodes)

    [reduce_crossovers(nodes) for _ in range(5)]

    redistribute_leaves(nodes)
    [reduce_crossovers(nodes) for _ in range(5)]
    [reduce_total_branch_length(nodes, shift) for shift in shift_list]
    rescale_positions(nodes)
    [reduce_crossovers(nodes) for _ in range(5)]



    # rescale_positions(nodes)

    # shift_list = [0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.02]
    # # shift_list = [0.5, 0.4, 0.3, 0.2]

    # [reduce_crossovers(nodes) for _ in range(5)]

    # # [reduce_total_branch_length(nodes, shift) for shift in shift_list]
    # # [reduce_crossovers(nodes) for _ in range(4)]

    # redistribute_leaves(nodes)

    # [reduce_crossovers(nodes) for _ in range(5)]

    # # [reduce_total_branch_length(nodes, shift) for shift in shift_list]

    # [reduce_total_branch_length(nodes, shift) for shift in shift_list]

    # # [reduce_crossovers(nodes) for _ in range(5)]

    # # redistribute_leaves(nodes)
    # # [reduce_crossovers(nodes) for _ in range(5)]
    # # redistribute_leaves(nodes)
    # # [reduce_crossovers(nodes) for _ in range(5)]


    # rescale_positions(nodes)



#    print('====')

def get_arg_nodes(n=5, N=10000, r=1e-8, L=5e3, simulation="arg"):
    """
    Simulates an ARG
    """
    # because we use the sequence interval from 0 to 1
    r = r * L

#    print('\nSIM:')


    assert simulation in ["arg", "smcprime", "smc"]

    nodes = list()
    live = list() # live lineages

    for i in range(n):
        leaf = Leaf(i, height=0)
        nodes.append(leaf)
        lin = Lineage(lineageid=i, down=leaf, intervals=[(0, 1)])
        leaf.parent = lin
        live.append(lin)
        last_node = i # max node number used so far
        last_lineage = i
        
    while len(live) > 1:
        shuffle(live)

        coal_prob = (len(live) * len(live)-1) / 2 / (2*N) 
        tot_ancestral = sum(interval_sum(x.intervals) for x in live)
        rec_prob = r * tot_ancestral
        wating_time = exponential(1 / (coal_prob + rec_prob))
        height = wating_time + nodes[-1].height
        if random() < coal_prob / (coal_prob + rec_prob):
            # coalescence

            # if full arg:            
            a, b = 0, 1

            if simulation == "smcprime'":
                # intervals must overlap or be adjacent:
                while not (interval_intersect(live[a].intervals, live[b].intervals) \
                    or interval_any_shared_borders(live[a].intervals, live[b].intervals)):
                    shuffle(live)
            elif simulation == "smc":
                # intervals must overlap:
                while not interval_intersect(live[a].intervals, live[b].intervals):
                    shuffle(live)
    
            lin_a, lin_b = live.pop(0), live.pop(0)

            intervals = interval_union(lin_a.intervals, lin_b.intervals)

            # new node
            if lin_a.down is lin_b.down:
                # make sure diamond is symmetric
                children=[lin_a.down.left_parent, lin_b.down.right_parent]
            else:
                children=[lin_a, lin_b]  
            node_c = Coalescent(nodeid=last_node+1, height=height, 
                        children=children)
            
            last_node += 1
            nodes.append(node_c)

            # add node to top of coalescing lineages
            lin_a.up = node_c
            lin_b.up = node_c

            # new lineage
            lin_c = Lineage(lineageid=last_lineage+1, down=node_c, intervals=intervals) # fixme intervals
            last_lineage += 1
            live.append(lin_c)
            node_c.parent = lin_c
        else:
            # recombination

            # rec_lin = choices(live, weights=[interval_sum(x.intervals) for x in live], k=1)[0]
            rec_lin = choices(live, weights=[interval_span(x.intervals) for x in live], k=1)[0]
            live.remove(rec_lin)

            # total ancestral material
            total_anc = interval_sum(rec_lin.intervals)

            # recombination point in ancestral material
            recomb_point_anc = random() * total_anc

            # recombination point in full sequence
            cum = 0
            for s, e in rec_lin.intervals:
                cum += e - s
                if cum > recomb_point_anc:
                    recomb_point = e - (cum - recomb_point_anc)
                    break

            # recombination node
            rec_node = Recombination(nodeid=last_node+1, height=height, 
                child=rec_lin, recomb_point=recomb_point)
            last_node += 1
            nodes.append(rec_node)

            # two new lineages both refering back to recombination node
            intervals_a, intervals_b = interval_split(rec_lin.intervals, recomb_point)
            assert interval_sum(intervals_a) and interval_sum(intervals_b)
            assert sum(e-s for s, e in intervals_a) + sum(e-s for s, e in intervals_b) == sum(e-s for s, e in rec_lin.intervals)
            lin_a = Lineage(lineageid=last_lineage+1, down=rec_node, intervals=intervals_a)
            last_lineage += 1
            lin_b = Lineage(lineageid=last_lineage+1, down=rec_node, intervals=intervals_b)
            last_lineage += 1
            live.append(lin_a)
            live.append(lin_b)

            # add parents of node
            rec_node.left_parent = lin_a
            rec_node.right_parent = lin_b

            # add parent for child
            rec_lin.up = rec_node

    add_node_x_positions(nodes)

    return nodes

def _rescale_layout(pos, scale=1):
    # rescale to (0,pscale) in all axes

    # shift origin to (0,0)
    lim = 0  # max coordinate for all axes
    for i in range(pos.shape[1]):
        pos[:, i] -= pos[:, i].min()
        lim = max(pos[:, i].max(), lim)
    # rescale to (0,scale) in all directions, preserves aspect
    for i in range(pos.shape[1]):
        pos[:, i] *= scale / lim
    return pos

def get_positions(nodes):
    """
    Gets list of x,y positions for nodes
    """
    positions = list()
    for node in nodes:
        positions.append((node.xpos, node.height))
    return positions

def get_breakpoints(nodes):
    """
    Gets list of recombination break points
    """    
    # get list of break points
    breakpoints = list()
    for n in nodes:
        if type(n) is Recombination:
#        if type(n) is Recombination and all(p.up in nodes for p in n.parents):
            breakpoints.append(n.recomb_point)
    return sorted(breakpoints)

def get_parent_lineages(nodes, root=True):
    """
    Get list of parent lineages of nodes ordered by id
    """

    lineages = list()
    for node in nodes:
        if type(node) is Coalescent:
            lineages.append(node.parent)
        if type(node) is Recombination:
            lineages.append(node.left_parent)
            lineages.append(node.right_parent)
        if type(node) is Leaf:
            lineages.append(node.parent)

    # get unique and sort them so lineageid matches index
    lineages = list(set(lineages))
    lineages.sort(key=lambda x: x.lineageid)

    if not root:
        # remove dangling root lineage
        lineages = lineages[:-1]

    return lineages

def get_child_lineages(nodes):
    """
    Get list of child lineages of nodes ordered by id
    """

    lineages = list()
    for node in nodes:
        if type(node) is Coalescent:
            lineages.extend(node.children)
        if type(node) is Recombination:
            lineages.append(node.child)

    # get unique and sort them so lineageid matches index
    lineages = list(set(lineages))
    lineages.sort(key=lambda x: x.lineageid)

    return lineages

# def traverse_marginal(node, interval):
#     """
#     Recursive function for getting marginal tree/ARG
#     """    
#     node = deepcopy(node) # TODO: remove if input is a cloned arg.
#     tree_nodes = set()
#     if type(node) is Leaf:
#         tree_nodes.add(node)
#     if type(node) is Recombination:
#         if interval_intersect([interval], node.child.intervals):
#             tree_nodes.add(node)
#             tree_nodes.update(traverse_marginal(node.child.down, interval))
#     elif type(node) is Coalescent:
#         if node.parent is None or interval_intersect([interval], node.parent.intervals):
#             tree_nodes.add(node)
#         del_child = None
#         for i, child in enumerate(node.children):
#             if interval_intersect([interval], child.intervals):
#                 tree_nodes.update(traverse_marginal(child.down, interval))
#             else:
#                 del_child = i
#         if del_child is not None:
#             del node.children[del_child]
#     return tree_nodes

def _traverse_marginal(node, interval):
    """
    Recursive function for getting marginal tree/ARG
    """    
    tree_nodes = set()
    if type(node) is Leaf:
        tree_nodes.add(node)
    if type(node) is Recombination:
        if interval_intersect([interval], node.child.intervals):
            tree_nodes.add(node)
            tree_nodes.update(_traverse_marginal(node.child.down, interval))
    elif type(node) is Coalescent:
        if node.parent is None or interval_intersect([interval], node.parent.intervals):
            tree_nodes.add(node)
        del_child = None
        for i, child in enumerate(node.children):
            if interval_intersect([interval], child.intervals):
                tree_nodes.update(_traverse_marginal(child.down, interval))
            else:
                del_child = i
        if del_child is not None:
            del node.children[del_child]
    return tree_nodes       

def traverse_marginal(node, interval):
    clone = deepcopy(node) # TODO: remove if input is a cloned arg.
    return _traverse_marginal(clone, interval)


def remove_dangling_root(tree_nodes):
    """
    Remove the nodes of a marginal tree or ARG 
    from root to first coalescence
    """    
    for_del = list()
    for i in range(len(tree_nodes)-1, -1, -1):
        if type(tree_nodes[i]) is Coalescent and len(tree_nodes[i].children) == 2:
            break
        for_del.append(i)
    for i in for_del:
        del tree_nodes[i]

def marginal_arg(nodes, interval, strip_dangling_root=True):
    """
    Gets the marginal ARG given a sequene interval
    """

    # TODO: 
    # nodes = clone_arg(nodes)

    # get nodes for marginal
    marg_nodes = traverse_marginal(nodes[-1], list(interval))

    # set to list
    marg_nodes = list(marg_nodes)

    # sort on height
    marg_nodes.sort(key=lambda x: x.height)
    # prune top path above last coalescence
    if strip_dangling_root:
        remove_dangling_root(marg_nodes)

    return marg_nodes

def marginal_trees(nodes, interval, strip_dangling_root=False):
    """
    Gets list of marginal trees
    """
    tree_list = list()
    breakpoints = get_breakpoints(nodes)

    borders = [0] + breakpoints + [1]
    
    # for interval in zip(borders[:-1], borders[1:]):
    #     marg_nodes = marginal_arg(nodes, interval)
    #     tree_list.append(marg_nodes)
    # return tree_list

    interval_list = list()
    for interv in zip(borders[:-1], borders[1:]):
        marg_nodes = marginal_arg(nodes, interv, strip_dangling_root=strip_dangling_root)
        # TODO: make this correct:...
        if marg_nodes and interv[1] > interval[0] and interval[1] > interv[0]:
            tree_list.append(marg_nodes)
            interval_list.append(interv)
    return tree_list, interval_list

# def draw_graph(nodes):
#     """
#     Draws graph using matplotlib
#     """    
#     # make a graph and positions list
#     positions = list()
#     arg = nx.Graph()
#     for node in nodes:
#         arg.add_node(node.nodeid)
#         positions.append((node.xpos, node.height))
#         if isinstance(node, Recombination):
#             arg.add_edge(node.child.down.nodeid, node.nodeid)
#         elif isinstance(node, Coalescent):
#             for child in node.children:
#                 arg.add_edge(child.down.nodeid, node.nodeid)

#     positions = np.array(positions)
#     positions = dict(zip(arg.nodes(), positions))

#     #pos = nx.spring_layout(arg)
#     nx.draw(arg, positions, alpha=0.5, node_size=200, with_labels=True)
#     #with_labels=True, 
#     #connectionstyle='arc3, rad = 0.1', 
#     #arrowstyle='-')

#     # G = nx.DiGraph() #or G = nx.MultiDiGraph()
#     # G.add_node('A')
#     # G.add_node('B')
#     # G.add_edge('A', 'B', length = 2)
#     # G.add_edge('B', 'A', length = 3)
#     # pos = nx.spring_layout(G)
#     # nx.draw(G, pos, with_labels=True, connectionstyle='arc3, rad = 0.1')
#     # edge_labels=dict([((u,v,),d['length'])
#     #              for u,v,d in G.edges(data=True)])

#     plt.show()

def rescale_positions(nodes):
    """
    Rescales xpos and heights to 0,1 range
    """
    max_height = max(n.height for n in nodes)
    max_xpos = max(n.xpos for n in nodes)
    min_xpos = min(n.xpos for n in nodes)
    for node in nodes:
        node.xpos = (node.xpos - min_xpos) / (max_xpos-min_xpos)
        node.height = node.height / max_height

def arg2json(nodes):

    lineages = get_parent_lineages(nodes)

    json_data = dict(Coalescent=[], Recombination=[], Leaf=[], 
                        Lineage=[l.get_dict() for l in lineages])
    for node in nodes:
        json_data[node.__class__.__name__].append(node.get_dict())

    return json.dumps(json_data, indent=4)

def json2arg(json_str):

    nodes = list()

    data = json.loads(json_str)

    # make lineages (with indexes instead of node refs)
    lineages = [Lineage(**data) for data in data['Lineage']]

    for node_data in data['Coalescent']:            
        # make the node
        node = Coalescent(**node_data)
        # populate the parent and children with actual lineages
        node.parent = lineages[node.parent]
        node.children = [lineages[i] for i in node.children]
        nodes.append(node)

    for node_data in data['Recombination']:
        # make the node
        node = Recombination(**node_data)
        # populate the parents and child with actual lineages
        node.left_parent = lineages[node.left_parent]
        node.right_parent = lineages[node.right_parent]
        node.child = lineages[node.child]
        nodes.append(node)

    for node_data in data['Leaf']:
        # make the node
        node = Leaf(**node_data)
        # populate the parent with actual lineages
        node.parent = lineages[node.parent]
        nodes.append(node)

    nodes.sort(key=lambda x: x.nodeid)

    # populate up and down with the nodes
    for lineage in lineages:
        lineage.down = nodes[lineage.down]
        if lineage.up is not None:
            lineage.up = nodes[lineage.up]

    return nodes

if __name__ == '__main__':

    # get arg and add positions
    nodes = get_arg_nodes()

    print(nodes)

    print(deepcopy(nodes))

    sys.exit()

    json_str = arg2json(nodes)

    retrieved_nodes = json2arg(json_str)

    print(nodes)
    print(retrieved_nodes)
    print(nodes == retrieved_nodes)



    # get breakpoints
    breakpoints = get_breakpoints(nodes)
    print(breakpoints)

    # get marginal trees
    trees = marginal_trees(nodes)
    print(trees)

    # marginal arg for some consequtive intervals
    marg_arg = marginal_arg(nodes, [0, breakpoints[1]])


    # # draw graphs for testing
    # draw_graph(nodes)
    # draw_graph(marg_arg)

    # draw_graph(nodes)
    # for tree in marginal_trees(nodes):
    #     print([n.xpos for n in tree])
    #     draw_graph(tree)

