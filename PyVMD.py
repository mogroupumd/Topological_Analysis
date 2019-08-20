import numpy as np

__author__ = "Xingfeng He, Yunsheng Liu"
__copyright__ = "Copyright 2019, UMD Mo. Group"
__version__ = "0.2"
__maintainer__ = "Yunsheng Liu"
__email__ = "yliu1240@umd.edu"
__date__ = "Jun 12th, 2019"

"""
    Several simple functions to generate VMD Tcl command lines.
    Not difficult any way you slice it......
"""

def cmd_by_radius(node_structure, r_cut, molecular_id=None, scale=0.5):
    """
        To generate Tcl command lines by Voronoi radius of each node.
        
        Args:
            node_structure: a list of nodes, whose properties must include 'voronoi_radius'.
            r_cut: the cut off radius. Nodes with Voronoi radius lower than this value will not be presented.
            molecular_id: the graphics id in VMD Tcl command line. This must exist to execute the commands.
            scale: use scaled radius to better visualize the site structure.
        
        Return:
            cmds: a list of string with VMD Tcl command lines.
    """
    cmds = []
    for site in node_structure:
        if site.properties['voronoi_radius'] < r_cut:
            continue
        rad = site.properties['voronoi_radius'] * scale
        coords = ' '.join(map(str, site.coords))
        if molecular_id:
            cmd = 'graphics {} sphere {{{}}} radius {} resolution 20\n'.format(molecular_id, coords, rad)
        else:
            cmd = 'draw sphere {{{}}} radius {} resolution 20\n'.format(coords, rad)
        cmds.append(cmd)
    return cmds

def cmd_by_bv(node_structure, molecular_id=None):
    """
        To generate Tcl command lines by bond valence of each site.
        
        Args:
            node_structure: a list of nodes, whose properties must include 'valence_state'.
            molecular_id: the graphics id in VMD Tcl command line. This must exist to execute the commands.
        
        Return:
            cmds: a list of string with VMD Tcl command lines.
    """
    cmds = []
    for site in node_structure:
        bv = round(site.properties['valence_state'], 2)
        coords = ' '.join(map(str, site.coords + np.array([0, 0, 1]))) # change the text coordinate a little bit
        if molecular_id:
            cmd = 'graphics {} text {{{}}} {}\n'.format(molecular_id, coords, bv)
        else:
            cmd = 'draw text {{{}}} {}\n'.format(coords, bv)
        cmds.append(cmd)
    return cmds

def cmd_by_edge(node_structure, r_cut, molecular_id=None, scale=0.5):
    """
        To generate Tcl command lines by nearest neighbor distance.
        
        Args:
            node_structure: a list of nodes, whose properties must include 'neighbor_nodes'.
            r_cut: the cut off radius. Nodes with neighbor distance lower than this value will not be presented.
            molecular_id: the graphics id in VMD Tcl command line. This must exist to execute the commands.
            scale: use scaled radius to better visualize the site structure.
        
        Return:
            cmds: a list of string with VMD Tcl command lines.
    """
    cmds = []
    for site_i, site in enumerate(node_structure):
        for nn in site.properties['neighbor_nodes']:
            if (nn[0] < site_i) and (nn[2] == [0, 0, 0]):
                continue # avoid repeating the same internal path
            elif nn[1] < r_cut:
                continue
            rad = nn[1] * scale
            start_coords = ' '.join(map(str, site.coords))
            end_coords_frac = node_structure[nn[0]].frac_coords + nn[2]
            end_coords_cart = np.dot(end_coords_frac, node_structure.lattice.matrix)
            if np.linalg.norm(end_coords_cart - site.coords) > 8: # the edge distance may be too large
                print 'Site Index {}; Neighbor Index {}'.format(site_i, nn[0])
                print 'Site Frac Coords {}; Neighbor Frac Coords {}'.format(site.frac_coords, node_structure[nn[0]].frac_coords)
                print 'Distance'.format(nn[2])
            end_coords = ' '.join(map(str, end_coords_cart))
            if molecular_id:
                cmd = 'graphics {} cylinder {{{}}} {{{}}} radius {} resolution 20\n'.format(molecular_id, start_coords, end_coords,
                                                                                       rad)
            else:
                cmd = 'draw cylinder {{{}}} {{{}}} radius {} resolution 20\n'.format(start_coords, end_coords, rad)
            cmds.append(cmd)
    return cmds
