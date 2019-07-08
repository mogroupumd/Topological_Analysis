# system & I/O dependencies
import os
import sys
import time
import argparse
import pandas as pds
from ruamel.yaml import YAML
from copy import deepcopy

# pymatgen dependencies
from pymatgen.io.cif import CifParser, CifWriter
from pymatgen.core.periodic_table import Specie, Element
from pymatgen.core.composition import Composition

# Topological Analyzer dependencies
from Topological_Analyzer.filters import *
from Topological_Analyzer.PyVMD import cmd_by_radius
# TAPercolateFilter, TABvFilter, TALongFilter, TAOptimumSiteFilter, TACoulombReplusionFilter, OxidationStateFilter,
# TADenseNeighbor

__author__ = "Xingfeng He, Yunsheng Liu"
__copyright__ = "Copyright 2019, UMD Mo. Group"
__version__ = "0.2"
__maintainer__ = "Yunsheng Liu"
__email__ = "yliu1240@umd.edu"
__date__ = "Jun 21st, 2019"

# Preparation
yaml = YAML()

# customized directory to specify Topological_Analyzer is if it's not in PYTHONPATH or other python directory
# base_dir = '/Users/yunshengliu'
# if not base_dir in sys.path:
#     sys.path.append(base_dir)

def Analyze_Voronoi_Nodes(args):
    """
    A standard process to apply all filters. Zeo++ finds all possible polyhedrons and corresponding sites while this class
    will screen bad sites and merge them. The program currently support CIF input files ONLY;
    
    Args:
        args.cif_file (str): Directory of CIF file
        args.input_file (yaml): Directory of input file which specify filter parameters. The input file must be a yaml:
        
                                Mandatory:
                                    1. SPECIE: a string of target diffusion specie, with oxidation state;
                                               e.g. Li+: Li specie with +1 oxidation state.
                                            
                                Optional: (Each parameter must be added according to filters specified)
                                    Overall:
                                    2. ANION: a string of potential anion type in the structure.
                                              This parameter will automatically specify parameters for further analysis:
                                                  BV_UP
                                                  BV_LW
                                                  R_CUT
                                              However, these parameters will be overwritten if they're explicitly assigned.
                                              
                                              e.g. S (sulfur) will have R_CUT: 1.5 A,
                                                   If input yaml file has another R_CUT to be 2 A, the final R_CUT will be 2 A.
                                              
                                              Currently support following anions:
                                                               |   BV_LW   |   BV_UP   |   R_CUT(A)   |
                                                  ------------------------------------------------------
                                                  S (sulfur)   |    0.4    |    1.1    |      2.5     |
                                                  O (oxygen)   |    0.5    |    1.2    |      2.3     |
                                              
                                    VoroPerco:
                                    3. PERCO_R: the percolation radius for diffusion specie;
                                    VoroBV:
                                    4. BV_UP: the maximum bond valence of a site which considered to be appropriate;
                                    5. BV_LW: the minimum bond valence of a site which considered to be appropriate;
                                    Coulomb:
                                    6. R_CUT: the minimum distance between target specie and nearest ion (either anion or
                                              cation);
                                    VoroLong:
                                    7. LONG: the criteria to decide whether a node is a long node or not. Unit: A;
                                    MergeSite:
                                    8. NEIGHBOR: the distance criteria to decide whether 2 sites / nodes are too close to 
                                                 each other. Unit: A.
                                    9. LONG
                                    
        args.filters (str): strings to specify which filter to use in analysis:
        
                            FILTER: filters applied. Currently support following filters:
                                    Ordered: OrderFrameworkFilter
                                    PropOxi: OxidationStateFilter
                                    VoroPerco: TAPercolateFilter
                                    Coulomb: TACoulombReplusionFilter
                                    VoroBV: TABvFilter
                                    VoroLong: TALongFilter
                                    MergeSite: OptimumSiteFilter
                                    VoroInfo: TALongFilter, but only output the center coordinates and length of each node
                                
    Output:
        CIF files after applying each filter. The predicted sites for target specie will be represented as sites with 50%
        partial occupancy.
        Note that some filters may be fundamental (decide whether they're good CIFs or not) and they may have
        no output structures.
        e.g. if applying OxidationStateFilter, TAPercolateFilter and TABvFilter,
             there will be 2 output CIF files:
                 1. CIF with all accessible sites;
                 2. CIF with all sites having good bond valence;
             OxidationStateFilter has no ouput structure.
    """
    import Topological_Analyzer
    
    # built-in radius for different species
    va_dir = os.path.dirname(VoronoiAnalyzer.__file__)
    radii_yaml_dir = os.path.join(va_dir, 'files/radii.yaml')
    with open(radii_yaml_dir, 'r') as f:
        radii = yaml.load(f)
    f.close()
    
    # read structure from CIF file
    name = args.cif_file[:-4] # the last 4 characters are '.cif'
    precif = CifParser(args.cif_file, occupancy_tolerance=2.0)
    structure = precif.get_structures(primitive = False)[0].copy()
    # for input parameter file
    with open(args.input_file, 'r') as f:
        input_parameters = yaml.load(f)
    f.close()
    
    # target specie
    sp = Specie.from_string(input_parameters['SPECIE'])
    
    # other possible parameters
    if 'ANION' in input_parameters.keys():
        if input_parameters['ANION'].lower() == 's':
            bv_range = (0.4, 1.1)
            rc = 2.5
        elif input_parameters['ANION'].lower() == 'o':
            bv_range = (0.5, 1.2)
            rc = 2.3
        else:
            print '##    Unsupported anion type: {}'.format(input_parameters['ANION'])
            bv_range = (0, 1.5)
            rc = 2.0
            
    if 'PERCO_R' in input_parameters.keys():
        pr = input_parameters['PERCO_R'] # percolation radius
    else:
        pr = None
        
    try:
        # these exist further bond valence limits to overwrite existing ones
        tmp = bv_range
        if 'BV_UP' in input_parameters.keys():
            bv_range = (bv_range[0], input_parameters['BV_UP'])
        if 'BV_LW' in input_parameters.keys():
            bv_range = (input_parameters['BV_LW'], bv_range[1])
    except:
        # these's no anion type to assign bond valence range
        if ('BV_UP' in input_parameters.keys()) and ('BV_LW' in input_parameters.keys()):
            bv_range = (input_parameters['BV_LW'], input_parameters['BV_UP'])
        else:
            bv_range = None
            
    if 'R_CUT' in input_parameters.keys():
        rc = input_parameters['R_CUT'] # cut-off distance of coulomb replusion
    else:
        try:
            tmp = rc # to check whether parameter exists, if it doesn't exist, set it to None.
                     # only necessary for bv_range and r_cut because these 2 may be set by ANION parameter
        except:
            rc = None
            
    if 'LONG' in input_parameters.keys():
        long = input_parameters['LONG'] # cut-off distance to decide whether a node is long or not
    else:
        long = None
    if 'NEIGHBOR' in input_parameters.keys():
        nn = input_parameters['NEIGHBOR'] # cut-off distance to decide whether 2 sites are neighbors
    else:
        nn = None
    
    # temporary parameters for filters applied
    frame_structure = None
    org_frame = None
    node_structure = None
    predicted_structure = None
    
    for f_index, f in enumerate(args.filters):
        print 'Step {}: {}'.format(f_index, f)
        
        if f.lower() == 'ordered':
            # Check whether the framework is ordered or not.
            print '#     Check framework disordering.'
            orderFrame = OrderFrameworkFilter(structure.copy(), radii, sp)
            org_structure = orderFrame.virtual_structure.copy()
            frame_structure = orderFrame.virtual_framework.copy()
            org_frame = orderFrame.framework.copy()
            print '#     Check finishes.'
            
        if f.lower() == 'propoxi':
            # Check oxidation states in structures. This is necessary for bond valence filter.
            print '#     Check oxidation states in structure.'
            PropOxi = OxidationStateFilter(org_structure.copy())
            if not PropOxi.decorated:
                print '##    Oxidation state check fails...'
                sys.exit()
            else:
                print '#     Check finishes.'
                
        elif f.lower() == 'voroperco':
            # Check whether there's enough space for percolation.
            print '#     Check Voronoi percolation raduis.'
            if pr:
                VoroPerco = TAPercolateFilter(org_structure.copy(), radii, sp, pr)
            else:
                print '##    No percolation radius provided...'
                sys.exit()
            
            if not VoroPerco.analysis_results:
                print '##    Cannot percolate...'
                sys.exit()
            else:
                """
                    The Voronoi analysis results include:
                        Voronoi_accessed_node_structure: A structure with all nodes (with Voronoi radius added to the property
                                                         of each node);
                        Voronoi_structure: A structure containing nodes whose Voronoi radius is greater than a certain value;
                        Framework: The framework structure with no target diffusion specie;
                        free_sph_max_dia: Maximum spherical diameter in the structure;
                        ......
                    To see other results, please use 'analysis_keys' attribute of the class.
                """
                results = deepcopy(VoroPerco.analysis_results)
                print '#     Percolation diameter (A): {}'.format(round(results['free_sph_max_dia'], 3))
                output_structure = org_frame.copy()
                if results['Voronoi_accessed_node_structure']:
                    node_structure = results['Voronoi_accessed_node_structure'].copy()
                    for nodes in node_structure.copy():
                        output_structure.append(str(sp), nodes.coords, coords_are_cartesian=True)
                    CifWriter(output_structure).write_file('{}_all_accessed_node.cif'.format(name))
                    print '#     Percolation check finishes.'
                else:
                    print '##    Errors in Voronoi analysis structure...'
                    
        elif f.lower() == 'coulomb':
            print '#     Check Coulomb replusion effects.'
            if (not frame_structure) or (not node_structure):
                print '##    No framework and node structure provided for Coulomb Replusion analysis...'
                sys.exit()
            elif not rc:
                print '##    No Coulomb replusion cut-off distance provided...'
                sys.exit()
            else:
                if sp.oxi_state < 0:
                    ion = 'anion'
                else:
                    ion = 'cation'
                print '#     Processing Coulomb replusion check.'
                print '#     {} effect detected, minimum distance to {}s is {} A.'.format(ion, ion, round(rc, 3))
                CoulRep = TACoulombReplusionFilter(node_structure.copy(), frame_structure.copy(), prune=ion, min_d_to_ion=rc)
                if CoulRep.final_structure:
                    node_structure = CoulRep.final_structure.copy()
                    output_structure = org_frame.copy()
                    for node in node_structure.copy():
                        output_structure.append(str(sp), node.coords, coords_are_cartesian=True)
                    CifWriter(output_structure).write_file('{}_coulomb_filtered.cif'.format(name))
                    print '#     Coulomb replusion check finishes.'
                else:
                    print '##    All available nodes will experience high Coulomb replusion...'
                    print '##    The structure is either unreasonable or the replusion radius cut-off is too large...'
                    sys.exit()
                    
        elif f.lower() == 'vorobv':
            print '#     Check bond valence limits.'
            if (not frame_structure) or (not node_structure):
                print '##    No framework and node structure provided for bond valence analysis...'
                sys.exit()
            elif not bv_range:
                print '##    No bond valence range provided...'
                sys.exit()
            else:
                print '#     Processing bond valence check.'
                print '#     Bond valence limitation: {} - {}'.format(bv_range[0], bv_range[1])

                VoroBv = TABvFilter(node_structure.copy(), frame_structure.copy(), bv_range)
                if VoroBv.final_structure:
                    node_structure = VoroBv.final_structure.copy()
                    output_structure = org_frame.copy() # output cif structure
                    output_doc = {} # output csv file
                    variables = ['Cartesian_Coords', 'Voronoi_R', 'Bond_Valence']
                    for i in variables:
                        output_doc[i] = []
                    for node in node_structure.copy():
                        output_structure.append(str(sp), node.coords, coords_are_cartesian=True)
                        tmp_coords = [round(n, 4) for n in node.coords]
                        output_doc['Cartesian_Coords'].append(tmp_coords)
                        output_doc['Voronoi_R'].append(round(node.properties['voronoi_radius'], 3))
                        output_doc['Bond_Valence'].append(round(node.properties['valence_state'], 2))
                        
                    CifWriter(output_structure).write_file('{}_bond_valence_filtered.cif'.format(name))
                    df = pds.DataFrame(data=output_doc).sort_values(by=['Voronoi_R'])
                    df = df.reindex(variables, axis=1)
                    df.to_csv('{}_bv_info.csv'.format(name))
                    
                    print '#     Bond valence check finishes.'
                else:
                    print '##    All available nodes are excluded for bad bond valences...'
                    print '##    The structure is either unreasonable or the bond valence range is bad...'
                    sys.exit()
        
        elif f.lower() == 'vorolong':
            print '#     Check long nodes in structure.'
            if not node_structure:
                print '##    No node structure provided for long Voronoi node analysis...'
                sys.exit()
            elif not long:
                print '##    No length provided to decide Voronoi node length...'
                sys.exit()
            else:
                print '#     Processing Voronoi length check.'
                print '#     Voronoi length limitation: {} A'.format(round(long, 3))
                VoroLong = TALongFilter(node_structure.copy(), long, use_voro_radii=True)
                print '#     Maximum node length detected: {} A'.format(round(VoroLong.longest_node_length, 3))
                output_doc = {}
                variables = ['Center_Coords', 'Node_Length']
                for i in variables:
                    output_doc[i] = []
                for i in VoroLong.clusters:
                    tmp_coords = [round(n, 4) for n in i[0]]
                    output_doc['Center_Coords'].append(tmp_coords)
                    output_doc['Node_Length'].append(round(i[1], 4))
                df = pds.DataFrame(data=output_doc).sort_values(by=['Node_Length'])
                df = df.reindex(variables, axis=1)
                df.to_csv('{}_node_length_info.csv'.format(name))
                print '#     Central node information written.'
                if VoroLong.has_long_node:
                    print '#     Long node check finishes.'
                else:
                    print '##    The structure has no long nodes or node length restriction is bad...'
                    print '##    Please check the node length CSV for more information...'
                    sys.exit()
                    
        elif f.lower() == 'voroinfo':
            print '#     Output the center coordinates and length of each node......'
            if not node_structure:
                print '##    No node structure provided for Voronoi information...'
                sys.exit()
            else:
                VoroLong = TALongFilter(node_structure.copy(), 0, use_voro_radii=True)
                print '#     Maximum node length detected: {} A'.format(round(VoroLong.longest_node_length, 3))
                output_doc = {}
                variables = ['Center_Coords', 'Node_Length']
                for i in variables:
                    output_doc[i] = []
                for i in VoroLong.clusters:
                    tmp_coords = [round(n, 4) for n in i[0]]
                    output_doc['Center_Coords'].append(tmp_coords)
                    output_doc['Node_Length'].append(round(i[1], 4))
                df = pds.DataFrame(data=output_doc).sort_values(by=['Node_Length'])
                df = df.reindex(variables, axis=1)
                df.to_csv('{}_node_length_info.csv'.format(name))
                print '#     Voronoi node information written.'
                    
        elif f.lower() == 'mergesite':
            # before we use TAOptimumSiteFilter, we need to have a list of different clusters,
            # thus must use TADenseNeighbor and TALongFilter. Also note that all clusters in the list must be. 
            if not node_structure:
                print '##    No node structure provided for optimizing sites...'
                sys.exit()
            if (not nn) or (not long):
                print '##    No neighbor distance cut-off and long node cut-off provided for site optimization...'
                sys.exit()
                
            voro_dense = TADenseNeighbor(node_structure.copy(), close_criteria=1,
                                         big_node_radius=0, radius_range=[0, 0], use_radii_ratio=True)
            voro_long = TALongFilter(node_structure.copy(), 0, use_voro_radii=True)
            cluster_list = voro_dense.clustering(node_structure.copy(), 1, True, True)
            
            long_list = []
            short_list = []
            for i in cluster_list:
                if voro_long.get_cluster_length(i, use_voro_radii=True) >= long:
                    long_list.append(i)
                else:
                    short_list.append(i)
            print '#     Processing site optimization: nearest neighbor cut-off {} A.'.format(round(nn, 3))
            OpSite = TAOptimumSiteFilter(org_structure.copy(), nn, sp, sort_type='None', use_exp_ordered_site=False)
            opt_long_list = []
            opt_short_list = []
            for i in long_list:
                tmp_list = OpSite.optimize_cluster(i, nn, sort_type='radius')
                for j in tmp_list:
                    opt_long_list.append(j)
            for i in short_list:
                tmp_list = OpSite.optimize_cluster(i, nn, sort_type='radius')
                for j in tmp_list:
                    opt_short_list.append(j)
            print '#     Long node number: {}'.format(len(opt_long_list))
            print '#     Short node number: {}'.format(len(opt_short_list))
            new_list = []
            for i in opt_long_list:
                new_list.append(i)
            for i in opt_short_list:
                new_list.append(i)
            OpSite.add_cluster(new_list)
            
            output_structure = OpSite.site_structure.copy()
            half_list = [] # it seems 50% occupancy sites are easier to see. You may directly use output_structure otherwise
            for i in output_structure:
                ppt = deepcopy(i.properties)
                new_i = PeriodicSite({str(sp): 0.5}, i.coords, i.lattice, to_unit_cell=False, coords_are_cartesian=True,
                                     properties=ppt)
                half_list.append(new_i)
            half_structure = Structure.from_sites(half_list, charge=None, validate_proximity=False, to_unit_cell=False)
            CifWriter(half_structure).write_file('{}_{}_optimized_sites.cif'.format(name, 'radius'))
            # CifWriter(output_structure).write_file('{}_{}_optimized_sites.cif'.format(name, 'radius'))
            
            # for predicted structure:
            tot_num = org_structure.composition[sp]
            current_num = OpSite.site_structure.composition.num_atoms
            ratio = tot_num / current_num
            if ratio > 1:
                print '##    Prediction error, please be cautious about the predicted results.'
                print '##    Please also double check whether the input parameters are reasonable...'
                ratio = 1
            prediction = org_frame.copy()
            for site in OpSite.site_structure.copy():
                prediction.append({str(sp): ratio}, site.coords, coords_are_cartesian=True)
            prediction.sort()
            predicted_structure = prediction.copy()
            print '#     Site optimization finishes.'
            
        else:
            print '##    Unsupported operation...'
    if predicted_structure:
        comp = org_structure.composition.reduced_formula
        CifWriter(predicted_structure).write_file('{}_{}_predicted.cif'.format(name, comp))
        cmds = cmd_by_radius(half_structure, 0.5)
        cmd_file = open('{}_cmd'.format(name), 'w')
        cmd_file.write('mol new\n')
        for lines in cmds:
            cmd_file.write(lines)
        cmd_file.close()
        
if __name__ == '__main__':
    start_time = time.time()
    parser = argparse.ArgumentParser(description='Voronoi analysis on structures')
    parser.add_argument("cif_file", type=str, help='CIF file directory')
    parser.add_argument("-i", "--input_file", type=str, help='Input yaml file to specify different parameters')
    parser.add_argument("-f", "--filters", nargs="+",
                        default=["Ordered", "PropOxi", "VoroPerco", "Coulomb", "VoroBV", "VoroInfo", "MergeSite"],
                        type=str, help="Default is PropOxi, VoroPerco, Coulomb, VoroBV, VoroInfo, MergeSite"
                                       "Ordered list of filters. Current only support 6 filters."
                                       "Please read README for further information")
    # default: "Vorolong"
    parser.set_defaults(func=Analyze_Voronoi_Nodes)
    args = parser.parse_args()
    args.func(args)
    print('Total used time: {}'.format(str(time.time() - start_time)))
