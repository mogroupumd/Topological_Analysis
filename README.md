# Topological_Analysis
Analyze material structures

# Introduction

Topological_Analyzer is a Python library to analyze interstitial sites in oxides and sulfides, based on Voronoi analysis. It specifies possible interstitial sites for a given specie, including its original sites in the structure. Topological Analyzer provides several Python filters for site selection and a simple executable script. The script produces several CIFs and CSV files which will show the Voronoi radius of interstitial sites and their positions. Bond valence and other information will also be included in CSV files. The Topological_Analyzer only supports Python 2 environment, recommending Python 2.7. It requires pre-installed Zeo++.

# Dependency

**The following softwares and Python libraries are required to use Topological_Analyzer:**
Zeo++ (http://pymatgen.org/_modules/pymatgen/io/zeopp.html#ZeoCssr)

Pymatgen(http://pymatgen.org/)

Numpy (http://www.numpy.org/)

os (https://docs.python.org/2/library/os.html)

argparse (https://docs.python.org/3/library/argparse.html)

PrettyTable (https://pypi.python.org/pypi/PrettyTable)

monty (https://materialsvirtuallab.github.io/monty/)

# Installation

1. Make sure all dependencies are installed and the executing environment is Python 2;
2. Add the library to your Python path: $export PYTHONPATH="/your/directory/Topological_Analysis/:$PYTHONPATH";
3. Add /Topological_Analyzer/scripts/ to your PATH: $export PATH="/your/directory/Topological_Analyzer/scripts/:$PATH";
4. Run 'analyze_Voronoi_nodes.py' for analysis or import filters in your Python programs.

# To Use the Script

The 'analyze_Voronoi_nodes.py' requires 2 input files. A CIF file for input structure and a Yaml file for analysis parameters. The Yaml file specifies percolation radius and other necessary information to use Topological_Analyzer. Please read the code comments for further details.

To run the script:
    
    $python analyze_Voronoi_nodes.py <cif_file> -i <yaml_file>

You may run -h to the script to read other documentations. An example Yaml parameter file (only list all usable parameters, example.yaml) and an executable input Yaml file (inputs.yaml for oxides) are included in script directory as well.

# Example

**Here, we will use LGPS (188887.cif) and Li2S (60432.cif) as examples. Input Yaml file (input)**
To analyze LGPS (188887.cif), simply execute:

    $cd /your/directory/Topological_Analyzer/examples/
    
    $python analyze_Voronoi_nodes.py 188887.cif -i inputs.yaml
    
Note: the input CIF file may contain disordered cation sites. Analysis codes in current version will replace the disordered cation sites with lowest radius cation ions. However, the output structure will preserve the original disordered sites. Please refer to the analysis comments for more information.

To use the VMD command lines, simply copy and paste the command line file contents to VMD terminal.

# Citing the Topological_Analyzer

If you use the Topological_Analyzer extensively, you may want to cite the original publication:
Xingfeng He1, Qiang Bai1, Yunsheng Liu, Adelaide M. Nolan, Chen Ling, Yifei Mo*, “Crystal Structural Framework of Lithium Super-Ionic Conductors”, Advanced Energy Materials 9, 1902078 (2019)  https://doi.org/10.1002/aenm.201902078
