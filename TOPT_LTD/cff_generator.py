#!/usr/bin/env python3
import argparse
import itertools
from pprint import pprint, pformat

class cFFException(Exception):
    pass

class cFF_Analyser(object):

    def __init__(self, internal_edges, external_edges):
        self.ie = list(internal_edges)
        self.ee = list(external_edges)
        self.edges = internal_edges+external_edges
        self.iv = sorted(list(set(sum([[a,b] for a,b in self.ie],[]))))
        self.ev = {}
        for i_ee, (k, v) in enumerate(self.ee):
            if k in self.iv:
                self.ev[k] = (v,-1,i_ee)
                raise cFFException("Define all external edges as incoming.")
            elif v in self.iv:
                self.ev[v] = (k,1,i_ee)
            else:
                raise cFFException("External edges not connected to any internal vertex.")
        self.vertices = sorted(list(set(list(self.iv)+[v_infos[0] for v_infos in self.ev.values()])))

    def get_cluster_boundary(self, cluster, seed_nodes_visited):
        
        neighbours = []
        for v in cluster:
            for a, b in self.ie:
                if a in cluster:
                    if b not in cluster:
                        neighbours.append(b)
                else:
                    if b in cluster:
                        neighbours.append(a)
        return sorted([v for v in set(neighbours) if v not in seed_nodes_visited])

    def grow_cluster_list(self,cluster_list,seed_nodes_visited):
        """ Grow the list of clusters supplied to aquired one more neighbour."""

        grown_clusters = []
        for cluster in cluster_list:
            cluster_boundary = self.get_cluster_boundary(cluster, seed_nodes_visited)
            grown_clusters.extend( [cluster+[v,] for v in cluster_boundary] )

        return [ list(c) for c in set([tuple(sorted(c)) for c in grown_clusters]) ]


    def get_connected_clusters(self):
        """ Build all possible connected clusters of this graph."""

        seed_nodes_to_visit = list(self.iv)
        seed_nodes_to_visit.reverse()
        seed_nodes_visited = []

        all_clusters = []
        while len(seed_nodes_to_visit)>0:
            seed_nodes_visited.append(seed_nodes_to_visit.pop())
            grown_cluster_list = [ [seed_nodes_visited[-1],], ]
            while len(grown_cluster_list)>0:
                all_clusters.extend(grown_cluster_list)
                grown_cluster_list = self.grow_cluster_list(grown_cluster_list, seed_nodes_visited)
        
        return [ tuple(c) for c in all_clusters ]

    def enforce_connected_completement(self, connected_clusters):

        filtered_connected_clusters = []
        for cc in connected_clusters:
                
            complement = [v for v in self.iv if v not in cc]
            # Also remove the connected cluster corresponding to the full graph
            if len(complement) == 0:
                continue
            if tuple(complement) not in connected_clusters:
                continue
            filtered_connected_clusters.append(cc)

        return filtered_connected_clusters

    def is_completing_a_family(self, family, cluster):
        # TODO: Warning: Quite inefficient, find better approach
        
        # First test if the cluster is completed by the family
        possible_unions = []
        for n_clusters in range(1,len(family)+1):
            for comb in itertools.combinations(family,n_clusters):
                possible_unions.append(set().union(*comb))
                if possible_unions[-1] == cluster:
                    return True
        
        possible_unions = [ set(sorted_u) for sorted_u in set([tuple(sorted(list(u))) for u in possible_unions]) ] 

        #print(possible_unions)

        # Then check if helps completing any other set
        for cf in family:
            for union in possible_unions:
                #print(cf, len(cf) - len(union), len(cluster), cf.difference(union) , cluster)
                if len(cf) - len(union) != len(cluster):
                    continue
                if cf.difference(union) == cluster:
                    return True

        return False

    def is_violating_connectivity_count(self,family,cc,possible_unions_of_n_connected_clusters):

        all_internal_vertices = set(self.iv)

        for n_clusters in range(1,len(family)+1):
            for comb in itertools.combinations(family,n_clusters):
                if not all_internal_vertices.difference(set().union(*comb).union(cc)) in possible_unions_of_n_connected_clusters[n_clusters]:
                    return True
        
        return False

    def analyze(self, full_analysis=True):
        analysis = {
            'cross_free_family' : None,
            'connected_clusters' : None
        }

        # Create all connected clusters
        connected_clusters = self.get_connected_clusters()

        # Now filter clusters that leave a complement that is not connected
        filtered_connected_clusters = connected_clusters
        filtered_connected_clusters = self.enforce_connected_completement(connected_clusters)

        # Now take all combinations of v-1 number of connected clusters that have either full or no intersections with other ones.
        filtered_connected_clusters = [set(c) for c in filtered_connected_clusters]

        analysis['connected_clusters'] = filtered_connected_clusters
        if not full_analysis:
            return analysis

        # Prepare the data for all possible unions of N connected clusters
        # possible_unions_of_n_connected_clusters = {0:[set(),]}
        # for n_cc in range(1,len(self.iv)-1):
        #     possible_unions_of_n_connected_clusters[n_cc] = []
        #     for fcc in filtered_connected_clusters:
        #         for u in possible_unions_of_n_connected_clusters[n_cc-1]:
        #             possible_unions_of_n_connected_clusters[n_cc].append(u.union(fcc))
        #     possible_unions_of_n_connected_clusters[n_cc] = [ set(sorted_u) for sorted_u in 
        #         set([ tuple(sorted(list(u))) for u in possible_unions_of_n_connected_clusters[n_cc] ]) ]

        cross_free_families = [[],]
        n_clusters_combined = 0
        while n_clusters_combined < len(self.iv)-1:
            new_cross_free_families = []
            for i_cc, cc in enumerate(filtered_connected_clusters):
                print('n_clusters_combined = %d/%d, i_cc = %d/%d, n_cross_free_families = %d'%(
                    n_clusters_combined, len(self.iv)-1, i_cc, len(filtered_connected_clusters), len(cross_free_families)
                ), end='\r')
                for family in cross_free_families:
                    if cc in family:
                        continue
                    if any(cc.intersection(f) not in [cc,f,set([])] for f in family):
                        continue
                    # TODO: Warning: Quite inefficient, find better approach
                    if self.is_completing_a_family(family,cc):
                        continue
                    # TODO: Also potentially quite inefficient
                    # if self.is_violating_connectivity_count(family,cc,possible_unions_of_n_connected_clusters):
                    #     continue
                    new_cross_free_families.append(list(family)+[cc,])
            new_cross_free_families = set( tuple(sorted( (tuple(cc) for cc in ccf) )) for ccf in new_cross_free_families)
            cross_free_families = [ [ set(cc) for cc in cff ] for cff in new_cross_free_families ]
            n_clusters_combined += 1
        # Remove cross-free families completing the whole graph
        cross_free_families = [cff for cff in cross_free_families if set().union(*cff)!=set(self.iv)]

        analysis['cross_free_family'] = cross_free_families

        return analysis