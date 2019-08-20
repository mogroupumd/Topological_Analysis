from __future__ import division, unicode_literals

__author__ = "Xingfeng He, Yunsheng Liu"
__copyright__ = "Copyright 2019, UMD Mo. Group"
__version__ = "0.2"
__maintainer__ = "Yunsheng Liu"
__email__ = "yliu1240@umd.edu"
__date__ = "Apr 26th, 2019"

# System dependencies
from copy import deepcopy
from monty.json import MSONable
from collections import Counter
import time
import numpy as np

# Topological_Analysis dependencies
from Topological_Analysis.ZeoExtendFunctions import get_voronoi_percolate_nodes, get_voronoi_node_edge, get_percolated_node_edge

# Pymatgen dependencies
from pymatgen.io.zeopp import get_free_sphere_params
from pymatgen.analysis.bond_valence import calculate_bv_sum
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.periodic_table import Specie, Element
from pymatgen.core.structure import Structure
from pymatgen.io.cif import CifParser, CifWriter


class GetVoronoiNodes(MSONable):
    """
        A convenient class to get the framework structure and Voronoi node structure from the original structure.
        Target diffusion specie / element will be removed from original structure to get the framework structure.
        Voronoi node structure will only contain positions of all nodes (only positions are meaningful)
    """
    def __init__(self, structure, diff_specie, radii_dict):
        """
            Args:
                structure: the original structure;
                diff_specie: target diffusion specie / element;
                radii_dict: a dictionary which provides ionic radius for all species;
            
            ..Attribute::
                frame_structure: framework without target diffusion specie / element;
                final_structure: positions of all Voronoi nodes;
        """
        self.original_structure = structure.copy()
        self.sp = diff_specie
        self.rs = radii_dict
        
        self.frame_structure = None
        self.node_structure = None
        
        elems = [el for el in structure.composition.elements if el != self.sp]
        all_els = [el for el in structure.composition.elements]
        el_names = [str(el) for el in elems]
        ex_names = [n for n in self.rs.keys()]
        
        # check if the diffusion specie exists in the structure
        if not self.sp in all_els:
            raise Exception('Diffusion specie not found : {} / {}'.format(self.sp, structure.composition.reduced_formula))
            
        # check if the ionic radius is provided in dictionary for each specie
        if not (set(el_names) < set(ex_names)):
            # Note: we use < here because target diffusion specie is excluded from the element list
            print 'More ionic radius expected:'
            for i in (set(el_names) - set(ex_names)):
                print i
            raise Exception('Incomplete ionic radius dictionary.')
            
        frames = []
        for site in self.original_structure:
            if not str(self.sp) in site.species_string:
                frames.append(site)
        self.frame_structure = Structure.from_sites(frames, charge=None, validate_proximity=False, to_unit_cell=False)
        
        # run Voronoi analysis:
        error = True
        for attempts in range(5):
            try:
                node_struct = get_voronoi_node_edge(self.frame_structure.copy(), self.rs)
                error = False
                break
            except:
                time.sleep(5) # safety first
        
        if not error:
            self.node_structure = node_struct.copy()
    
    @property
    def final_structure(self):
        return self.node_structure

    
class OrderFrameworkFilter():
    """
        Since the Zeo++ cannot deal with disordered structures, we have to order the disordered structure first.
        Current implementation method is to replace the disordered sites with the specie whose ionic radius is smaller.
    """
    def __init__(self, structure, radii_dict, diff_specie):
        """
            Args:
                structure: a Structure object;
                diff_specie: target diffusion specie;
                radii_dict: a dictionary which provides all ionic radius of elements
            
            ..Attribute::
                framework: the original framework with disordered sites;
                virtual_framework: the framework where original disordered sites are replaced by smaller radius species;
                virtual_structure: the structure where original disordered non-diffusion sites are replaced.
        """
        self.original_structure = structure.copy()
        self.sp = diff_specie
        self.rd = radii_dict
        
        self.framework = self.frame.copy()
        if self.isordered(self.framework):
            self.virtual_structure = structure.copy()
            self.virtual_framework = self.frame.copy()
        else:
            self.virtual_framework = self.replace_frame.copy()
            new_structure = self.replace_frame.copy()
            for i in self.original_structure.copy():
                if str(self.sp) in i.species_string:
                    new_structure.append(i.species_and_occu, i.coords, coords_are_cartesian=True)
            new_structure.sort()
            self.virtual_structure = new_structure.copy()
        
    def isordered(self, structure):
        """
            Check if there's disordered sites in the structure
        """
        orders = True
        for i in structure.copy():
            if not str(self.sp) in i.species_string:
                if not i.is_ordered:
                    orders = False
                    break
        return orders
                    
    @property
    def frame(self):
        frame_list = []
        for i in self.original_structure.copy():
            if not str(self.sp) in i.species_string:
                frame_list.append(i)
        return Structure.from_sites(frame_list, charge=None, validate_proximity=False, to_unit_cell=False)
    
    @property
    def replace_frame(self):
        frame_list = []
        for i in self.original_structure.copy():
            if not str(self.sp) in i.species_string:
                if i.is_ordered:
                    frame_list.append(i)
                else:
                    radii_list = []
                    for sps in i.species_and_occu.keys():
                        if not sps in self.rd.keys():
                            print 'Warning! Ionic radius not found: {}'.format(sps)
                        else:
                            radii_list.append((sps, self.rd[sps]))
                    if len(radii_list) > 0:
                        radii_list = sorted(radii_list, key=lambda k: k[1]) # pick the specie with smallest radius
                    else:
                        radii_list.append((i.species_and_occu.keys()[0], 1)) # pick the first specie
                    new_i = PeriodicSite(radii_list[0][0], i.coords, i.lattice, to_unit_cell=False,
                                         coords_are_cartesian=True)
                    frame_list.append(new_i)
        return Structure.from_sites(frame_list, charge=None, validate_proximity=False, to_unit_cell=False)
    

class OxidationStateFilter():
    """
        To be convenient, the structure is required to be decorated with proper oxidation states to perform the analysis.
    """
    def __init__(self, structure):
        """
            Args:
                structure: a pymatgen Structure object;
            
            ..Attribute::
                final_structure: None if the structure doesn't have oxidation states;
                decorated: boolean to check whether the structure has oxidation states.
        """
        self.original_structure = structure.copy()
        self.final_structure = None
        if self.decorated:
            self.final_structure = structure.copy()
    
    @property
    def decorated(self):
        f = True
        for el, occu in self.original_structure.composition.items():
            if isinstance(el, Element):
                f = False
                break
        return f

    
class TAPercolateFilter(MSONable):
    def __init__(self, structure, radii_dict, diff_specie, percolate_r):
        """
            Args:
                structure: a Structure object with oxidation states. This means the 'elements' in structure are Specie objects;
                radii_dict: a dictionary which provides all ionic radius of elements;
                diff_specie: the target specie which intends to diffuse;
                percolate_r: the minimum percolation radius for target diffusion specie;
                
            ..Attribute::
                final_structure: None if the structure doesn't have enough room for percolation;
                analysis_results: structures and other results in Voronoi analysis.
                                  For more information, use attribute 'analysis_keys';
                analysis_keys: a list of usable results from Voronoi analysis.
        """
        self.rs = radii_dict
        self.sp = diff_specie
        self.percolate_r = percolate_r
        self.original_structure = structure.copy()
        self.final_structure = None
        self.voronoi_analysis = None
        
        all_els = [el for el in structure.composition.elements]
        elems = [el for el in structure.composition.elements if el != self.sp]
        el_names = [str(el) for el in elems]
        ex_names = [n for n in self.rs.keys()]
        
        # check if the diffusion specie exists in the structure
        if not self.sp in all_els:
            raise Exception('Diffusion specie not found : {} / {}'.format(self.sp, structure.composition.reduced_formula))
        
        # check if the ionic radius is provided in dictionary for each specie
        if not (set(el_names) < set(ex_names)):
            # Note: we use < here because target diffusion specie is excluded from the element list
            print 'More ionic radius expected:'
            for i in (set(el_names) - set(ex_names)):
                print i
            raise Exception('Incomplete ionic radius dictionary.')
            
        # check percolation
        voronoi_results = self.percolate()
        if voronoi_results:
            self.final_structure = structure.copy()
            self.voronoi_analysis = voronoi_results
    
    def percolate(self):
        """
            To decide whether the percolation radius in the structure is larger than the radius provided.
        """
        struct = self.original_structure.copy()
        # elems = [el for el in struct.composition.elements if el != self.sp]
        
        # create the framework structure without target diffusion specie, including disordered sites:
        framework = []
        for i in struct:
            if not str(self.sp) in i.species_string:
                framework.append(i)
        frame_structure = Structure.from_sites(framework, charge=None, validate_proximity=False, to_unit_cell=False)
        
        # free sphere parameters:
        free_sphere_params = get_free_sphere_params(frame_structure, rad_dict=self.rs, probe_rad=0.1)
        # update the sphere parameters to final analysis results
        results = free_sphere_params
        results['Framework'] = frame_structure.copy()
        
        if free_sphere_params['free_sph_max_dia'] < 2 * self.percolate_r:
            return False # no percolation path
        else:
            error = True
            # sometimes Voronoi analysis program may counter numerical errors, so better try running it for several times
            for attempts in range(5): # run for 5 times
                try:
                    node_struct = get_voronoi_node_edge(frame_structure, self.rs)
                    node_struct_access = get_percolated_node_edge(frame_structure, node_struct, self.rs, self.percolate_r)
                    error = False
                    break
                except:
                    time.sleep(5) # safety first when running heavy computation loads
                    
            if error: # can percolate, but cannot do Voronoi analysis
                f = open('./error.out', 'a')
                f.write('Error in Voronoi analysis. Attempt: {}\n'.format(attempts + 1))
                f.close()
                # update the voronoi analysis structures to None
                results['Voronoi_structure'] = None
                results['Voronoi_accessed_node_structure'] = None
                return results
            
            # can both percolate and do Voronoi analysis:
            """
                Note that in ZeoExtendFunctions, all node structures don't have specific site element. All are 'X'.
            """
            new_list = []
            for node in node_struct:
                properties = deepcopy(node.properties)
                new_node = PeriodicSite(str(self.sp), node.coords, node.lattice, to_unit_cell=False, coords_are_cartesian=True,
                                    properties=properties)
                new_list.append(new_node)
            results['Voronoi_structure'] = Structure.from_sites(new_list, charge=None, validate_proximity=False,
                                                                to_unit_cell=False).copy()
            
            new_list = []
            for node in node_struct_access:
                properties = deepcopy(node.properties)
                new_node = PeriodicSite(str(self.sp), node.coords, node.lattice, to_unit_cell=False, coords_are_cartesian=True,
                                    properties=properties)
                new_list.append(new_node)
            results['Voronoi_accessed_node_structure'] = Structure.from_sites(new_list, charge=None, validate_proximity=False,
                                                                              to_unit_cell=False).copy()
            return results
    
    @property
    def analysis_results(self):
        return self.voronoi_analysis
    
    @property
    def analysis_keys(self):
        return [n for n in self.voronoi_analysis.keys()]


class TACoulombReplusionFilter(MSONable):
    """
        A simple filter to eliminate Voronoi nodes which are too close to cations or anions.
        Physically, the strong Coulomb replusion will prevent atom being too close to ions with similar oxidation states.
        Here, we simply apply a minimum distance to imitate the replusion process.
    """
    def __init__(self, node_structure, frame_structure, prune='C', min_d_to_ion=0):
        """
            Args:
                node_structure: A Structure which contains all positions of Voronoi nodes.
                                Available from Voronoi analysis results in percolation filter or GetVoronoiNodes;
                frame_structure: the framework Structure which doesn't contain the target diffusion specie;
                prune: specify what kind of ions we're filtering. Can either be cations or anions.
                        Ususally, the prune type should be the same ion type to target diffusion specie;
                        e.g. both diffusion specie and prune type are cations.
                min_d_to_ion: the minimum distance permitted to prune type ions;
            
            ..Attribute::
                final_structure: A Structure which contains all positions of available nodes after Coulomb replusion filtering.
                                Similar to node_structure.
                
        """
        if not prune[0].lower() in ['c', 'a']:
            raise Exception('Can\'t recognize the provided prune type. Should either be cation or anion.')
        self.prune = prune
        self.min_d = min_d_to_ion
        self.original_nodes = node_structure.copy()
        self.framework = frame_structure.copy()
        
        # check the minimum distance to ions:
        if self.min_d >= 5:
            print 'Warning! Minimum distance more than 5 A...'
        
        # do replusion check:
        self.final_nodes = self.replusion_check()
        
    def replusion_check(self):
        nodes_struct = self.original_nodes.copy()
        frames = self.framework.copy()
        pruned_v_nodes = [] # our return structure
        
        for node_i, node in enumerate(nodes_struct):
            neighbor_list = []
            
            # if minimum distance is not provided, only use nearest neighbors.
            if self.min_d != 0: 
                sphere_neighbor_list = frames.get_sites_in_sphere(node.coords, self.min_d, include_index=False)
                for i in sphere_neighbor_list:
                    neighbor_list.append(i[0])
            else:
                sphere_neighbor_list = frames.get_sites_in_sphere(node.coords, 6, include_index=False)
                sphere_neighbor_list = sorted(sphere_neighbor_list, key=lambda k: k[1])
                neighbor_list.append(sphere_neighbor_list[0][0])
            
            good = True
            for site in neighbor_list:
                el, occu = site.species_and_occu.items()[0] # Note: el must be Specie, not Element
                if self.prune[0].lower() == 'c': # cations
                    if el.oxi_state > 0:
                        good = False
                        break
                elif self.prune[0].lower() == 'a': # anions
                    if el.oxi_state < 0:
                        good = False
                        break
                else:
                    raise Exception('Unexpected prune type. Must be either cation or anion...')
            if good:
                pruned_v_nodes.append(node)
        
        if pruned_v_nodes:
            tmp_structure = Structure.from_sites(pruned_v_nodes, charge=None, validate_proximity=False, to_unit_cell=False)
            tmp_structure.sort()
            return tmp_structure.copy()
        else:
            return None
    
    @property
    def final_structure(self):
        if self.final_nodes:
            return self.final_nodes.copy()
        else:
            return None

        
class TABvFilter(MSONable):
    """
        To restrict the bond valence range of sites in computed Voronoi nodes. The node structure and framework structure can
        be gathered from GetVoronoiNodes() class or TAPercolateFilter() class. The bond valence limitation here is a very good
        restriction to eliminate unreasonable sites.
    """
    def __init__(self, node_structure, frame_structure, bv_range, scale_factor=1):
        """
            Args:
                node_structure: A Structure which contains all positions of Voronoi nodes.
                                Available from Voronoi analysis results in percolation filter or GetVoronoiNodes;
                frame_structure: the framework Structure which doesn't contain the target diffusion specie;
                bv_range: Minimum and maximum of bond valences. e.g.[0.5, 1.0] means that minimum bond valence is 0.5 while
                          maximum is 1.0.
                scale_factor: A scaling factor for bond valences. 1.0 for experimental structures, 1.015 for GGA structures.
                              The scale factor is the same with the one used in Pymatgen calculate_bv_sum().
            
            ..Attribute::
                final_structure: A Structure which contains all positions of available nodes with appropriate bond valence
                                filtering. Similar to node_structure.
                                
        """
        self.min = bv_range[0]
        self.max = bv_range[1]
        self.scale = scale_factor
        self.original_nodes = node_structure.copy()
        self.framework = frame_structure.copy()
        
        self.final_nodes = self.bv_filter(self.original_nodes.copy(), self.framework, bvmin=self.min, bvmax=self.max)
        
    def bv_filter(self, nodes, frame, bvmin=0, bvmax=1, r=5):
        if nodes == None:
            return None
        
        good_node_list = []
        
        for site in nodes.copy():
            neighbor_list = frame.get_neighbors(site, r) # get the neighbor atom list for bond valence calculation
            bvsum = calculate_bv_sum(site, neighbor_list, scale_factor=self.scale)
            site_property = deepcopy(site.properties)
            site_property['valence_state'] = bvsum
            new_site = PeriodicSite(site.species_string, site.coords, site.lattice, to_unit_cell=False,
                                    coords_are_cartesian=True, properties=site_property)
            if (bvsum >= bvmin) and (bvsum <= bvmax):
                good_node_list.append(new_site)
        
        if len(good_node_list) != 0:
            new_struct = Structure.from_sites(good_node_list, charge=None, validate_proximity=False, to_unit_cell=False)
            new_struct.sort()
        else:
            new_struct = None
        
        return new_struct
    
    @property
    def final_structure(self):
        if self.final_nodes:
            return self.final_nodes.copy()
        else:
            return None
    
class TADenseNeighbor(MSONable):
    """
        Check neighbors and clustering neighboring nodes in Voronoi node structures.
    """
    def __init__(self, node_structure, close_criteria, big_node_radius, radius_range, use_radii_ratio=True,
                 use_pbc_dist=True):
        """
            Args:
                node_structure: Structure which contains all Voronoi nodes;
                close_criteria: A criteria to decide whether 2 nodes are neighbors.
                                
                                If using radius ratio algorithm, this variable will be a ratio of radius sum (see
                                'use_radii_ratio' below). It has no unit.
                                
                                If NOT using radius ratio algorithm, this variable will be a cut-off distance (see
                                'use_radii_ratio' below). Its unit will be A (Angstrom).
                                
                big_node_radius: A radius criteria to decide whether a node is a big node or not;
                radius_range: A shell radius range to decide number of neighbors. The parameter should have a minimum and a
                              maximum. Since there may still exist nodes that are very close to each other after clustering,
                              so a minimum radius is necessary here.
                              e.g.
                                                      A-B --- C
                                                      |
                                                      |
                                                      D
                                   Here, B is very close to A, so A and B should be one site, while C and D are A's neighbors.
                use_radii_ratio: Whether to use the radius ratio algorithm to decide whether 2 nodes are neighbors;
                                 For node A and node B:
                                 If True, the neighbor criteria will be defined as:
                                 
                                     Neighbor_Distance = ratio * (A_radius + B_radius)
                                     If AB < Neighbor_Distance:
                                         A, B are neighbors
                                     else:
                                         A, B are not neighbors
                                     Here, 'close_criteria' will act as 'ratio' variable;
                                     
                                 If False, the close neighbor criteria will be a distance cut-off:
                                 
                                     Neighbor_Distance = distance_cut_off
                                     If AB < Neighbor_Distance:
                                         A, B are neighbors
                                     else:
                                         A, B are not neighbors
                                     Here, 'close_criteria' will act as distance_cut_off.
                                     
                use_pbc_dist: Whether to use periodic boundary distance or not.
            
            ..Attribute::
                final_structure: A node structure which combines all neighbor nodes and stores node radii as in their site
                                 properties.
                averageNeighborNum: average neighbor number;
                clustering: make all nodes into a list of clusters;
                            e.g. cluster_list:
                                 [cluster 1, cluster 2, ......]
                                 cluster:
                                 [node 1, node 2, ......]
                close_to_cluster: decide whether a node belongs to the cluster.                
        """
        self.original_structure = node_structure.copy()
        self.cc = close_criteria
        self.bn_r = big_node_radius
        self.rr = radius_range
        self.use_ratio = use_radii_ratio
        self.use_pbc = use_pbc_dist
        
        # to prune unnecessary nodes
        node_struct = self.original_structure.copy()
        dense_struct = self.prune_neighbor_nodes(node_struct, self.cc, self.bn_r, self.use_ratio, self.use_pbc)
        # to get average neighbor atom number
        neighbor_num = 0
        shell_r = (self.rr[0] + self.rr[1]) / 2 # r = (max + min)/2
        shell_dr = abs(self.rr[0] - self.rr[1]) / 2 # dr = (max - min)/2
        for k in dense_struct:
            neighbor_num += len(dense_struct.get_neighbors_in_shell(k.coords, shell_r, shell_dr))
        self.avg_nn = neighbor_num / float(dense_struct.num_sites)
        if dense_struct:
            dense_struct.sort()
            self.final_structure = dense_struct.copy()
        else:
            self.final_structure = None
        
    def prune_neighbor_nodes(self, node_structure, close_criteria, big_radius, use_radii_ratio, use_pbc_dist):
        """
            To merge neighbor nodes and get a clean structure. Since an original node structure may contain many
            nodes at similar locations, this function aims to merge all these 'close' nodes by implementing a simple
            clustering algorithm.
            
            Args:
                node_structure: Original structure with all available nodes.
                close_criteria
                use_radii_ratio
                use_pbc_dist
            
            Return:
                A node structure which contains simplified nodes.
        """
        if node_structure == None:
            return None
        
        pruned_nodes = []
        cluster_list = self.clustering(node_structure.copy(), close_criteria, use_radii_ratio, use_pbc_dist)
        pruned_nodes_lists = [self.select_big_nodes(deepcopy(single_cluster), big_radius) for single_cluster in cluster_list]
        for i in pruned_nodes_lists:
            for j in i:
                pruned_nodes.append(j)
        pruned_structure = Structure.from_sites(pruned_nodes, charge=None, validate_proximity=False, to_unit_cell=False)
        pruned_structure.sort()
        return pruned_structure.copy()
    
    def clustering(self, node_structure, close_criteria, use_radii_ratio, use_pbc_dist):
        """
            An algorithm to merge all individual nodes into several clusters.
            The overall algorithm is very simple:
                1. If a node belongs to only one cluster, merge the node and the cluster;
                2. If a node belongs to more than one cluster, that means these clusters and nodes are connected.
                   Merge all clusters and nodes connected;
                3. If a node doesn't belong to any cluster, list an individual cluster for the node.
            
            Args:
                node_structure: A structure contains all available Voronoi nodes;
                close_criteria
                use_radii_ratio
                use_pbc_dist
            
            Return:
                A list of all clusters.
                e.g. [(node 1, node 2, node 3, ......), (node i, node j, node k, ......), ......]
        """
        clusters = []
        for node_index, node in enumerate(node_structure):
            """
                Create a list of indexes which represents all connected clusters;
            """
            cluster_index_list = []
            for c_index, cluster in enumerate(clusters):
                nn = self.close_to_cluster(node, cluster, close_criteria, use_radii_ratio, use_pbc_dist)
                if nn:
                    cluster_index_list.append(c_index)
            """
                Algorithm:
                1. If there's no connected cluster, the node is an individual cluster;
                2. If there's one cluster, the node belongs to that cluster;
                3. If there're multiple clusters, the node connects all these clusters and they should be merged.
                
                    e.g.        Cluster A -- node 1 -- Cluster B
                                                |
                                                |
                                            Cluster C
                        Then node 1 connects Cluster A, B and C, so all 3 clusters should be merged as one.
            """
            if len(cluster_index_list) == 0:
                clusters.append([node]) # make an individual cluster
            elif len(cluster_index_list) == 1:
                clusters[cluster_index_list[0]].append(node) # merge the specific cluster and the node
            else:
                new_clusters = [] # the new cluster list
                new_cluster = [] # the new cluster which contains all clusters connected by current node
                # for other unrelated clusters
                for i in range(0, len(clusters)):
                    if not i in cluster_index_list:
                        new_clusters.append(deepcopy(clusters[i]))
                # merge all connected clusters
                for i in cluster_index_list:
                    for j in clusters[i]:
                        new_cluster.append(j)
                new_cluster.append(node) # merge the cluster and the current node
                new_clusters.append(new_cluster)
                clusters = deepcopy(new_clusters) # update the cluster list
        return clusters
    
    def close_to_cluster(self, node, cluster, close_criteria, use_radii_ratio, use_pbc_dist):
        """
            To decide whether a node is too close to any group in the group list.
            
            Args:
                node: A Site object with voronoi_radius property.
                cluster: A group of nodes. All nodes in the same cluster are neighbors.
                         Nodes don't have to be neighbors to all the other nodes in the same cluster, but any member in the
                         cluster must be neighbor to at least one member in the same cluster. (it's a little bit hard
                         to understand at first, but please take a while).
                         e.g. Cluster:
                                         A -- B -- C
                                              |    |
                                              D -- E
                                               \  /
                                                F
                              Both A and F are members in the group, even though A is only the neighbor of B, while F is
                              a neighbor of D and E. A and F are not neighbors.
                                 
                              Of course, this algorithm can be further extended to higher restrictions (you may change
                              the code so that each member must be neighbors to at least 2 other members. In this case,
                              A will not be in the group while F will remain in the group).
                close_criteria: see 'close_criteria' reference in the __init__().
                use_radii_ratio
                use_pbc_dist
            
            Return:
                True / False -- whether the node is a neighbor of the group.
        """
        close_clusters = []
        n_r = node.properties['voronoi_radius'] # the Voronoi radius for the target node
        # Note: As long as the target node is a neighbor of any site in a group, we decide it's in the group.
        # This algorithm may be further changed into being a neighbor of any number of sites in a group (currently 1).
        for site_index, site in enumerate(cluster):
            """
                Whether use radius ratio algorithm or regular distance cut-off to decide neighbors.
                The neighbor criteria will be defined separately.
            """
            if use_radii_ratio:
                # Voronoi radius: a built-in property for each site from Voronoi calculation. All sites from
                # GetVoronoiNodes() or TAPercolateFilter()'s node structure have this property.
                site_radii = site.properties['voronoi_radius'] # the Voronoi radius for the site in group
                criteria = close_criteria * (n_r + site_radii)
            else:
                criteria = close_criteria
                
            """
                Whether use periodic boundary distance.
            """
            if use_pbc_dist:
                dist = site.distance(node)
            else:
                dist = site.distance(node, jimage=np.array([0, 0, 0]))
            if dist <= criteria:
                return True
        return False
    
    def select_big_nodes(self, cluster, big_criteria):
        """
            Select only big nodes.
            
            Args:
                cluster: a list of nodes;
                big_criteria: the radius cut-off for big nodes;
            
            Return:
                node_list: a list of all big nodes. If no node achieve the criteria of being a big node, return the largest
                           node only.
        """
        if len(cluster) == 0:
            return []
        
        node_list = []
        sorted_cluster = sorted(cluster, key=lambda k: k.properties['voronoi_radius'], reverse=True)
        if sorted_cluster[0].properties['voronoi_radius'] < big_criteria:
            return [sorted_cluster[0]]
        else:
            for node in sorted_cluster:
                if node.properties['voronoi_radius'] > big_criteria:
                    node_list.append(node)
            return node_list
    
    @property
    def averageNeighborNum(self):
        return self.avg_nn
            
class TALongFilter(MSONable):
    """
        Check whether the structure has long nodes. The length of nodes are key parameter here.
    """
    def __init__(self, node_structure, long_criteria, use_voro_radii=True):
        """
            Args:
                node_structure: Structure which contains all Voronoi nodes;
                long_criteria: A radius cut-off to decide whether a node is a long one or not;
                use_voro_radii: Whether to use Voronoi radius or not;
            
            Return:
                final_structure: The same node structure as the original one;
                longest_node_length: The maximum length of all node clusters;
                has_long_node: Boolean whether there're long nodes in the structure;
                clusters: a list of cluster coordinates and the length of each cluster.
                          e.g. [(coords, length), (coords, length), ......]
                get_cluster_length: get the length of a cluster;
                get_average_coords: get the central coordinates of a cluster;
        """
        self.l_c = long_criteria
        self.original_structure = node_structure.copy()
        self.use_voro_radii = use_voro_radii
        
        voro_dense = TADenseNeighbor(self.original_structure.copy(), close_criteria=1,
                                     big_node_radius=0, radius_range=[0, 0], use_radii_ratio=True)
        cluster_list = voro_dense.clustering(self.original_structure.copy(), 1, True, True)
        cluster_info = []
        for cluster in cluster_list:
            c_coords = self.get_average_coords(cluster)
            c_length = self.get_cluster_length(cluster, self.use_voro_radii)
            cluster_info.append((c_coords, c_length))
        self.cluster_info = sorted(cluster_info, key=lambda k: k[1], reverse=True)
        self.max_node = max([n[1] for n in cluster_info])
     
    def get_cluster_length(self, cluster, use_voro_radii):
        """
            An algorithm to calculate the length of a cluster of nodes.
            The length of a cluster is defined as the longest distance between 2 nodes in the cluster.
            If use Voronoi radius, it's the longest distance between 2 nodes in the cluster plus the Voronoi radius of 2 nodes.
            e.g.
                If we use Voronoi radius:
                
                                                - A --- B
                                                   \    |
                                                    \   |
                                                     \  |
                                                      \ |
                                                       C
                                                        \
                                                        
                                         Cluster_Length = AC + R_a + R_c
                                         where AC is the distance between A and C, R_a is the Voronoi radius of A, while
                                         R_c is the Voronoi radius of C.
                                         
                If we don't use Voronoi radius:
                                         
                                         Cluster_Length = AC
            
            Args:
                cluster: A list of nodes;
                use_voro_radii: Whether to use Voronoi radius or not;
            
            Return:
                Length of cluster, with unit A.
        """
        if len(cluster) == 0:
            return 0 # if there's no node, the length of cluster is 0
        elif len(cluster) == 1:
            # if there's only 1 node, the length of cluster is 2R when using Voronoi radius and 0 when not using radius.
            if use_voro_radii:
                return 2 * cluster[0].properties['voronoi_radius']
            else:
                return 0
        else:
            longest_d = 0
            for i in range(len(cluster)-1):
                for j in range(i+1, len(cluster)):
                    """
                        Note: the 'if / else' must be inside the loop because it can be that AB is very large while R_a and
                              R_b are very small, which leads to a small distance overall.
                    """
                    if use_voro_radii:
                        r1 = cluster[i].properties['voronoi_radius']
                        r2 = cluster[j].properties['voronoi_radius']
                        tmp_d = cluster[i].distance(cluster[j]) + r1 + r2
                    else:
                        tmp_d = cluster[i].distance(cluster[j])
                    if tmp_d >= longest_d:
                        longest_d = tmp_d
            return longest_d
    
    def get_average_coords(self, cluster):
        """
            A convenient method to get the center position of a cluster.
            Note: if a cluster of nodes is very long, or well connected by several sub-clusters, the center position may be
                  very weird.
            
            Args:
                cluster: A list of nodes;
            
            Return:
                Center coordinate of the cluster, in cartesian coordinate.
        """
        if len(cluster) == 0:
            return [0, 0, 0]
        elif len(cluster) == 1:
            return cluster[0].coords
        else:
            struct_temp = Structure.from_sites(cluster)
            longest_d = max([struct_temp[0].distance(site) for site in struct_temp[1:]])
            sites = struct_temp.get_sites_in_sphere(struct_temp[0].coords, longest_d + 0.01, include_index=True)
            if len(sites) != struct_temp.num_sites:
                print'There are periodic image of some sites exists, only use the nearest one'
                unique_sites_index = []
                unique_sites = []
                for i in sites:
                    if i[-1] not in unique_sites_index:
                        unique_sites_index.append(i[-1])
                        unique_sites.append(i)
                position = sum(i[0].coords for i in unique_sites) / (len(unique_sites))
            else:
                position = sum(i[0].coords for i in sites) / (len(sites))   
            return position
    
    @property
    def final_structure(self):
        return self.original_structure.copy()
    
    @property
    def longest_node_length(self):
        return self.max_node
    
    @property
    def has_long_node(self):
        return self.max_node >= self.l_c
    
    @property
    def clusters(self):
        return self.cluster_info

class TAOptimumSiteFilter(MSONable):
    """
        To optimize the site prediction, making clusters of sites into unique and separated sites.
        There will be several method to decide the priority of all available nodes. Bond valence method will keep the sites
        with best-fitting bond valence value. Radius method will keep the site with largest Voronoi radius first. A regular
        method will begin with random site and make sure that every 2 unique sites are separated by certain distance.
    """
    def __init__(self, structure, dist_criteria, diff_specie, sort_type='radius', use_exp_ordered_site=False):
        """
        Args:
            structure: The structure of material which contains both diffusion specie and framework ions. All sites in the
                       structure must be decorated with oxidation state.
            dist_criteria: A distance cut-off to decide whether 2 sites are far away enough from each other. If the distance
                           is greater than criteria, then these 2 sites are unique sites, otherwise they'll be merged as one;
            diff_specie: target diffusion specie. MUST be a Specie object if use bond valence sorting method;
            sort_type: Specify a method to decide the priority of sites. It currently supports Bond Valence Method (bv) and
                       Voronoi Radius Method (radius). Other values will activate regular random method;
            use_exp_ordered_site: Whether to use the original site from original structure. If True, all sites from original
                                  structure will be preserved.
        
        Return:
            final_structure: A structure which contains all possible sites for diffusion specie and framework ions;
            site_structure: A structure which contains all possible sites for diffusion specie;
            optimize_list: a list of nodes where all members are far away from each other;
            optimize_cluster / optimize_clusters
        """
        self.original_structure = structure.copy()
        self.d = dist_criteria
        self.sp = diff_specie
        self.bv = None
        self.sort = sort_type
        self.use_exp = use_exp_ordered_site
        
        if self.sort.lower() == 'bv':
            try:
                self.bv = self.sp.oxi_state
            except:
                raise Exception('Target diffusion specie without oxidation state!...')
        
        node_list = []
        frame_list = []
        for site in structure:
            if (str(self.sp) in site.species_string) and site.is_ordered:
                # Note that diffusion specie sites can be disordered
                properties = site.properties
                node_list.append(PeriodicSite(str(self.sp), site.coords, site.lattice, to_unit_cell=False,
                                              coords_are_cartesian=True, properties=properties))
            else:
                frame_list.append(site)
        self.framework = Structure.from_sites(frame_list, charge=None, validate_proximity=False, to_unit_cell=False)
        if self.use_exp:
            if self.sort.lower() != 'radius':
                self.site_structure = Structure.from_sites(node_list, charge=None, validate_proximity=False,
                                                           to_unit_cell=False)
            else:
                raise Exception('Experimental sites don\'t have Voronoi radius...')
        else:
            self.site_structure = []
        
        if self.use_exp:
            self.final_structure = self.original_structure.copy()
        else:
            self.final_structure = self.framework.copy()
    
    def optimize_list(self, cluster, nn_dist, sort_type):
        """
            To optimize a list of nodes / sites.
            All nodes / sites in the final list should be far away from each other (unique nodes / sites).
            
            Also depends on number of nodes in the cluster.
        """
        final_list = []
        if len(cluster) <= 1:
            return cluster
        else:
            # check whether the cluster is longer than length criteria
            longest_d = 0
            vr = 0
            for i in range(len(cluster)-1):
                for j in range(i+1, len(cluster)):
                    if cluster[i].distance(cluster[j]) >= longest_d:
                        longest_d = cluster[i].distance(cluster[j])
                        if sort_type.lower() == 'radius':
                            vr = (cluster[i].properties['voronoi_radius'] + \
                                  cluster[j].properties['voronoi_radius'])/2
            # if the node is shorter than length criteria, it has only 1 node available
            if longest_d <= nn_dist:
                properties = None
                if sort_type.lower() == 'radius':
                    properties = deepcopy(cluster[0].properties)
                    properties['voronoi_radius'] = vr
                final_list.append(PeriodicSite(cluster[0].species_string, self.get_average_coords(cluster),
                                               cluster[0].lattice, to_unit_cell=False, coords_are_cartesian=True,
                                               properties=properties))
                return final_list
        # the sorting process is done outside the function.
        for n in cluster:
            if len(final_list) == 0:
                final_list.append(n)
            else:
                nn = False
                for s in final_list: # to check whether the node / site is a neighbor to EVERY previous unique one.
                    if s.distance(n) < nn_dist:
                        nn = True
                        break
                if not nn:
                    final_list.append(n)
        return final_list
    
    def optimize_cluster(self, cluster, nn_dist, sort_type, target_bv=None):
        """
            Optimize a list of nodes / sites. Can apply Bond Valence Method (bv), Voronoi Radius Method (radius) and regular
            method.
            
            Args:
                cluster: A list of nodes / sites. Should have Voronoi radius as its property;
                nn_dist: the distance criteria to decide whether 2 sites are too close to each other;
                sort_type: Specify which sorting method to use;
                target_bv: appropriate bond valence value. Must be specified when using bond valence method.
            
            Return:
                final_list: A list of unique nodes / sites. All nodes / sites are quite far away from each other.
        """
        # create a temporary list of nodes and sort it according to specific method
        tmp_list = []
        # check the cluster situation:
        if len(cluster) <= 1:
            return self.optimize_list(cluster, nn_dist, sort_type)
        else:
            # check the node length
            longest_d = 0
            for i in range(len(cluster)-1):
                for j in range(i+1, len(cluster)):
                    if cluster[i].distance(cluster[j]) >= longest_d:
                        longest_d = cluster[i].distance(cluster[j])
            if longest_d <= nn_dist:
                return self.optimize_list(cluster, nn_dist, sort_type)
            
        if sort_type.lower() == 'bv':
            # calculate bond valence and sort all sites by how close they are to appropriate bond valence
            for node in cluster:
                nn_list = self.framework.get_neighbors(node, r_c)
                bv = calculate_bv_sum(node, nn_list, scale_factor=1)
                properties = node.properties
                properties.update({'valence_state': bv})
                tmp_list.append(PeriodicSite(str(self.sp), node.coords, node.lattice, to_unit_cell=False,
                                             coords_are_cartesian=True, properties=properties))
            tmp_list = sorted(tmp_list, key=lambda k: abs(target_bv - k.properties['valence_state']))
        elif sort_type.lower() == 'radius':
            # keep the largest node in the cluster
            tmp_list = sorted(cluster, key=lambda k: k.properties['voronoi_radius'], reverse=True)
        else:
            tmp_list = deepcopy(cluster)
        
        return self.optimize_list(tmp_list, nn_dist, sort_type)
    
    def optimize_clusters(self, cluster_list, nn_dist, sort_type, target_bv=None):
        """
            Optimize a list of clusters. Can apply Bond Valence Method (bv), Voronoi Radius Method (radius) and regular
            method. First, we will optimize each cluster, then we merge all unique nodes / sites into one list and optmize
            this final node / site list.
            
            Args:
                clusters: A list of clusters. Should have Voronoi radius as its property;
                nn_dist: the distance criteria to decide whether 2 sites are too close to each other;
                sort_type: Specify which sorting method to use;
                target_bv: appropriate bond valence value. Must be specified when using bond valence method.
            
            Return:
                final_list: A list of unique nodes / sites. All nodes / sites are quite far away from each other.
        """
        single_list = []
        for cluster in cluster_list:
            site_list = self.optimize_cluster(cluster, nn_dist, sort_type, target_bv)
            for node in site_list:
                single_list.append(node)
        
        if sort_type.lower() == 'bv':
            # calculate bond valence and sort all sites by how close they are to appropriate bond valence
            single_list = sorted(single_list, key=lambda k: abs(target_bv - k.properties['valence_state']))
        elif sort_type.lower() == 'radius':
            # keep the largest node in the cluster
            single_list = sorted(single_list, key=lambda k: k.properties['voronoi_radius'], reverse=True)
        
        return self.optimize_list(single_list, nn_dist, sort_type)
    
    def add_clusters(self, cluster_list):
        """
            An operation to add more clusters to site structure.
            Be aware that there's no sorting method in this part. Previous sites are always preferred to remain.
        """
        added_nodes = self.optimize_clusters(cluster_list, self.d, self.sort, self.bv)
        if self.site_structure:
            original_nodes = self.site_structure.copy()
        else:
            original_nodes = []
        
        single_list = []
        for node in original_nodes:
            single_list.append(node)
        for node in added_nodes:
            single_list.append(node)
        
        final_list = self.optimize_list(single_list, self.d, self.sort)
        # update the site structure
        self.site_structure = Structure.from_sites(final_list, charge=None, validate_proximity=False,
                                                   to_unit_cell=False).copy()
        # update the final structure
        new_frame = self.framework.copy()
        for node in self.site_structure:
            new_frame.append(node.species_string, node.coords, coords_are_cartesian=True)
        new_frame.sort()
        self.final_structure = new_frame.copy()
        
    def add_cluster(self, cluster):
        added_nodes = self.optimize_cluster(cluster, self.d, self.sort, self.bv)
        if self.site_structure:
            original_nodes = self.site_structure.copy()
        else:
            original_nodes = []
        
        single_list = []
        for node in original_nodes:
            single_list.append(node)
        for node in added_nodes:
            single_list.append(node)
        
        final_list = self.optimize_list(single_list, self.d, self.sort)
        # update the site structure
        self.site_structure = Structure.from_sites(final_list, charge=None, validate_proximity=False,
                                                   to_unit_cell=False).copy()
        # update the final structure
        new_frame = self.framework.copy()
        for node in self.site_structure:
            new_frame.append(node.species_string, node.coords, coords_are_cartesian=True)
        new_frame.sort()
        self.final_structure = new_frame.copy()
    
    def get_average_coords(self, cluster):
        if len(cluster) == 0:
            return [0, 0, 0]
        elif len(cluster) == 1:
            return cluster[0].coords
        else:
            struct_temp = Structure.from_sites(cluster)
            longest_d = max([struct_temp[0].distance(site) for site in struct_temp[1:]])
            sites = struct_temp.get_sites_in_sphere(struct_temp[0].coords, longest_d + 0.01, include_index=True)
            if len(sites) != struct_temp.num_sites:
                print'There are periodic image of some sites exists, only use the nearest one'
                unique_sites_index = []
                unique_sites = []
                for i in sites:
                    if i[-1] not in unique_sites_index:
                        unique_sites_index.append(i[-1])
                        unique_sites.append(i)
                position = sum(i[0].coords for i in unique_sites) / (len(unique_sites))
            else:
                position = sum(i[0].coords for i in sites) / (len(sites))   
            return position
