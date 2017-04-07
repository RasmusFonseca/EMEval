
#from TEMPy.MapParser import readMRC
from TEMPy.MapParser import fromfile
from TEMPy.ScoringFunctions import ScoringFunctions
from TEMPy.StructureParser import PDBParser
import sys

if len(sys.argv)<3:
    print("Usage: python emeval.py <map-file> <structure-file> [<structure-file2> ...]")
    sys.exit(-1)

map_filename    = sys.argv[1:]
structure_files = sys.argv[2:]

scorer = ScoringFunctions()

target_map     = MapParser.readMRC(map_filename) #read target map
target_map_min = target_map.min()
target_map_max = target_map.max()

for structure_filename in structure_files:
    structure_instance=PDBParser.read_PDB_file(structure_filename, structure_file)

    #minimum density value based on protein molecular weight.
    min_thr=target_map.get_min_threshold(structure_instance.get_prot_mass_from_atoms(), target_map_min, target_map_max) 
    score_ENV=scorer.envelope_score(target_map, min_thr, structure_instance,norm=True)
    print("%s %s %.3f"%(map_filename, structure_filename, score_ENV))

