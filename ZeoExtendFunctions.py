from __future__ import division, unicode_literals

"""
This module defines more function for Zeo++
"""

__copyright__ = "Copyright 2019, UMD Mo. Group"
__author__ = "Xingfeng He"
__version__ = "0.1"
__maintainer__ = "Yunsheng Liu"
__email__ = "yliu1240@umd.edu"
__date__ = "Jul 9, 2015"

import subprocess
import copy
from monty.dev import requires
from monty.tempfile import ScratchDir
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.io.zeopp import ZeoCssr, ZeoVoronoiXYZ
from pymatgen.core.sites import PeriodicSite
from zeo.netstorage import AtomNetwork

process = subprocess.Popen('which network'.split(), stdout=subprocess.PIPE)
output = process.communicate()[0]
network_exe = (output != '')

@requires(network_exe, "get_voronoi_percolate_nodes requires Zeo++ c++ executable")
def get_voronoi_percolate_nodes(structure, rad_dict=None, probe_rad=0.1):
    """
    This function is used to get voronoi percolate nodes. Different from get_percolated_node_edge,
    the vor_accessible_node_struct returned by this one does not contain neighbor information. So if you need
    neighboring information of each accessible node, you may use get_percolated_node_edge function.
    Args:
        structure (Structure): Structure object for analysis
        rad_dict (dict): optional, dictionary of radii of elements in structures.
            If not given, Zeo++ default values are used.
            Note: Zeo++ uses atomic radii of elements.
            For ionic structures, pass rad_dict with ionic radii.
        probe_rad:

    Returns:
        vor_node_struct, vor_accessible_node_struct, vor_edgecenter_struct, vor_facecenter_struct (Structure):


    """
    with ScratchDir('.'):
        name = "temp_zeo2"
        zeo_inp_filename = name + ".cssr"
        ZeoCssr(structure).write_file(zeo_inp_filename)
        rad_file = None
        rad_flag = False
        if rad_dict:
            rad_file = name + ".rad"
            rad_flag = True
            with open(rad_file, 'w+') as fp:
                for el in rad_dict.keys():
                    fp.write("{} {}\n".format(el, rad_dict[el].real))

        atmnet = AtomNetwork.read_from_CSSR(
            zeo_inp_filename, rad_flag=rad_flag, rad_file=rad_file)
        vornet, vor_edge_centers, vor_face_centers = \
            atmnet.perform_voronoi_decomposition()
        vornet.analyze_writeto_XYZ(name, probe_rad, atmnet)
        voro_out_filename = name + '_voro.xyz'
        voro_node_mol = ZeoVoronoiXYZ.from_file(name + '_voro.xyz').molecule
        voro_accessible_node_mol = ZeoVoronoiXYZ.from_file(name + '_voro_accessible.xyz').molecule
    species = ["X"] * len(voro_node_mol.sites)
    coords = []
    prop = []
    for site in voro_node_mol.sites:
        coords.append(list(site.coords))
        prop.append(site.properties['voronoi_radius'])

    lattice = Lattice.from_lengths_and_angles(
        structure.lattice.abc, structure.lattice.angles)
    vor_node_struct = Structure(
        lattice, species, coords, coords_are_cartesian=True,
        to_unit_cell=False, site_properties={"voronoi_radius": prop})

    # percolate node struct
    species = ["X"] * len(voro_accessible_node_mol.sites)
    coords = []
    prop = []
    for site in voro_accessible_node_mol.sites:
        coords.append(list(site.coords))
        prop.append(site.properties['voronoi_radius'])

    lattice = Lattice.from_lengths_and_angles(
        structure.lattice.abc, structure.lattice.angles)
    vor_accessible_node_struct = Structure(
        lattice, species, coords, coords_are_cartesian=True,
        to_unit_cell=False, site_properties={"voronoi_radius": prop})

    # PMG-Zeo c<->a transformation for voronoi face centers
    rot_face_centers = [(center[1], center[2], center[0]) for center in
                        vor_face_centers]
    rot_edge_centers = [(center[1], center[2], center[0]) for center in
                        vor_edge_centers]

    species = ["X"] * len(rot_face_centers)
    prop = [0.0] * len(rot_face_centers)  # Vor radius not evaluated for fc
    vor_facecenter_struct_origin = Structure(
        lattice, species, rot_face_centers, coords_are_cartesian=True,
        to_unit_cell=False, site_properties={"voronoi_radius": prop})
    vor_facecenter_struct = Structure.from_sites(list(set([i for i in vor_facecenter_struct_origin])))
    species = ["X"] * len(rot_edge_centers)
    prop = [0.0] * len(rot_edge_centers)  # Vor radius not evaluated for fc
    vor_edgecenter_struct_origin = Structure(
        lattice, species, rot_edge_centers, coords_are_cartesian=True,
        to_unit_cell=False, site_properties={"voronoi_radius": prop})
    vor_edgecenter_struct = Structure.from_sites(list(set([i for i in vor_edgecenter_struct_origin])))
    return vor_node_struct, vor_accessible_node_struct, vor_edgecenter_struct, vor_facecenter_struct


@requires(network_exe, "get_voronoi_node_edge requires Zeo++ c++ executable")
def get_voronoi_node_edge(structure, rad_dict, write_nt2_file=False):
    """
    This function obtain the voronoi nodes and the edge information with 0 probe radius. It will contain neighboring
    information of each nodes.
    !!!Notice the nt2 file change the xyz to zxy. Thus the result need rotate. What a strange setting!!!
    Args:
        structure: orginal Structure
        rad_dict: {'Na':10,...}
    Returns:
        Structure: the voronoi node structure. For each site, there are voronoi_radius, neighbor_atoms, neighbor_nodes
            in the properties. 
            The neighbor_nodes in the format : [[node_No,channel_size,node_image]...] node image is list of int, 
                not array
            The neighbor_atoms in the format: [1,2,3,4] the number is the sequence number of sites in structure. So 
            take care of the structure, don't do any change on the sequence order.
    """
    with ScratchDir('.'):
        name = "temp_zeo3"
        cssr_file = name + ".cssr"
        rad_file = name + ".rad"
        mass_file = name + ".mass"
        out_file = name + ".nt2"
        cssr = ZeoCssr(structure)
        with open(cssr_file, 'w+') as fp:
            fp.write(str(cssr))
        # CifWriter(structure).write_file(cif_file) the xyz need to rotate to zxy, so just use existing ZeoCssr
        # rad_file
        with open(rad_file, 'w+') as fp:
            for el in rad_dict.keys():
                fp.write("{} {}\n".format(el, rad_dict[el].real))
        # mass_file
        struct_element_set = set(structure.composition.elements)
        with open(mass_file, 'w+') as fp:
            for el in struct_element_set:
                fp.write("{} {}\n".format(el, float(el.atomic_mass)))
        bashcmd = "network -r {} -mass {} -nt2 {} {}".format(rad_file, mass_file, out_file, cssr_file)
        process = subprocess.Popen(bashcmd.split(), stdout=subprocess.PIPE)
        zeo_plus_plus_output = process.communicate()[0]
        with open(out_file, 'r') as fp:
            nt2_file = fp.readlines()
    if write_nt2_file != False:
        with open(write_nt2_file, 'w+') as fp:
            for i in nt2_file:
                fp.write(i)
    node_info_string_list = nt2_file[1:nt2_file.index('\n')]
    edge_info_string_list = nt2_file[nt2_file.index('\n') + 2:]

    # node info string to Structure object
    species = ["X"] * len(node_info_string_list)
    coords = []
    voronoi_radius = []
    neighbor_atoms_num = []  # the four closest atoms, the number corresponding to the sequence number in structure
    neighbor_nodes = [list() for _ in range(len(node_info_string_list))]
    # the connected neighbor voronoi nodes, [node_No,channel_size,node_image], node image: [-1,0,0] 
    # means the actual connected neighbor will have a transition in -a direction.
    # final_frac_coords = node_frac_coords+node_image

    for site in node_info_string_list:
        coords_temp = map(float, site.split()[1:4])
        coords.append([coords_temp[1], coords_temp[2], coords_temp[0]])
        voronoi_radius.append(float(site.split()[4]))
        neighbor_atoms_num.append(map(int, site.split()[5:]))

    for neighbor in edge_info_string_list:
        valid_info = neighbor.split()
        valid_info.pop(1)  # there is a arrow in the output string
        image_temp = map(int, valid_info[3:6])
        valid_info = [int(valid_info[0]), int(valid_info[1]), float(valid_info[2]),
                      [image_temp[1], image_temp[2], image_temp[0]]]
        neighbor_nodes[valid_info[0]].append(copy.deepcopy(valid_info[1:]))
    lattice = Lattice.from_lengths_and_angles(
        structure.lattice.abc, structure.lattice.angles)
    vor_node_struct = Structure(
        lattice, species, coords, coords_are_cartesian=True,
        to_unit_cell=False, site_properties={"voronoi_radius": voronoi_radius,
                                             "neighbor_atoms": neighbor_atoms_num,
                                             "neighbor_nodes": neighbor_nodes})  # notice, to_unit_cell must be False,
    # or the image of neighbor_nodes need change
    return vor_node_struct


def get_percolated_node_edge(origin_struct, node_struct, rad_dict, probe_rad):
    """
    The function "get_voronoi_node_edge" gave all node information, but we need the percolated nodes. This function 
    conbine with function "get_voronoi_percolate_nodes" to delete the isolated nodes and correct the neighbor 
    information in the Structure from get_voronoi_node_edge.
    """
    _, vor_accessible_node_struct, _, _ = get_voronoi_percolate_nodes(origin_struct, rad_dict, probe_rad)
    # print vor_accessible_node_struct.lattice == node_struct.lattice
    # print vor_accessible_node_struct
    # node_struct_copy = copy.deepcopy(node_struct)
    unaccess_list = []
    access_list = []
    for i_index, i in enumerate(node_struct):
        match = 0
        for j_index, j in enumerate(vor_accessible_node_struct):
            if i.is_periodic_image(j, tolerance=1e-3, check_lattice=True):
                match = 1
                break
        if match == 0:
            unaccess_list.append(i_index)
        else:
            access_list.append(i_index)
    new_id = range(len(access_list))
    old_new_corr = {}
    for i_index, i in enumerate(access_list):
        old_new_corr[i] = new_id[i_index]
    access_node_site_list = [copy.deepcopy(site) for site_index, site in enumerate(node_struct) \
                             if site_index in access_list]
    final_access_node_site_list = []
    for i_index, i in enumerate(access_node_site_list):
        neighbor_nodes = [copy.deepcopy(neighbor) for neighbor in i.properties['neighbor_nodes'] \
                          if neighbor[0] in access_list]
        for j_index, j in enumerate(neighbor_nodes):
            new_id = old_new_corr[j[0]]
            neighbor_nodes[j_index] = copy.deepcopy([new_id] + j[1:])
        properties = i.properties
        properties['neighbor_nodes'] = neighbor_nodes
        node_site = PeriodicSite(i.species_string, i.coords, i.lattice, to_unit_cell=False,
                                 coords_are_cartesian=True, properties=properties)
        final_access_node_site_list.append(node_site)
    if final_access_node_site_list:
        new_node_edge = Structure.from_sites(final_access_node_site_list, validate_proximity=False, to_unit_cell=False)
    else:
        new_node_edge = None
    return new_node_edge
