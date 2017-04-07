from TEMPy.MapParser import MapParser
from TEMPy.ScoringFunctions import ScoringFunctions
from TEMPy.StructureParser import PDBParser
from TEMPy.StructureBlurrer import StructureBlurrer
import os,sys
from TEMPy.class_arg import TempyParser
from traceback import print_exc

def print_usage():
  print 'Usage:',sys.argv[0],"[options]"
  print 'Where options can be'
  print '  --help for help'
  print '  -m <input map>           .. Mandatory'
  print '  -p <input pdb>           .. Mandatory'
  print '  -r <resolution>          .. Mandatory'
  print '  -t <density threshold>   .. Optional'
  print 'This program takes the input pdb and computes the goodness-of-fit'
  print 'to the input map using different methods.'
  sys.exit(-1)

if len(sys.argv)==1 or "--help" in sys.argv:
  print_usage()

tp = TempyParser()
tp.generate_args()
#COMMAND LINE OPTIONS
m = tp.args.inp_map
r = tp.args.res
c = tp.args.thr
p = tp.args.pdb


#calculate map contour
def map_contour(m,t=-1.):
  mName = os.path.basename(m).split('.')[0]
  #print 'reading map'
  emmap=MapParser.readMRC(m)
  c1 = None
  if t != -1.0:
    #print 'calculating contour'
    zeropeak,ave,sigma1 = emmap._peak_density()
    if not zeropeak is None: c1 = zeropeak+(t*sigma1)
    else:
      c1 = 0.0
  return mName,emmap,c1

#calculate model contour
def model_contour(p,res=4.0,emmap=False,t=-1.):
  pName,modelmap = blur_model(p,res,emmap)
  c1 = None
  if t != -1.0:
    #print 'calculating contour'
    c1 = t*emmap.std()#0.0
  return pName,modelmap,c1

def blur_model(p,res=4.0,emmap=False):
  pName = os.path.basename(p).split('.')[0]
  #print 'reading the model'
  structure_instance=PDBParser.read_PDB_file(pName,p,hetatm=False,water=False)
  #print 'filtering the model'
  blurrer = StructureBlurrer()
  if res is None: sys.exit('Map resolution required..')
  #emmap = blurrer.gaussian_blur(structure_instance, res,densMap=emmap_1,normalise=True)
  modelmap = blurrer.gaussian_blur_real_space(structure_instance, res,densMap=emmap,normalise=True) 
  return pName,modelmap


#Get input data for one map and model
#print 'reading map'
if c is None: Name1,emmap1,c1 = map_contour(m,t=1.5)
else:
    Name1 = os.path.basename(m).split('.')[0]
    emmap1=MapParser.readMRC(m)
if r is None: sys.exit('Input a map, a model, map resolution and contours (optional)')
if p is None: sys.exit('Input a map, a model, map resolution and contours (optional)')
#print 'reading model'
Name2,emmap2,c2 = model_contour(p,res=r,emmap=emmap1,t=0.5)



#print 'Scoring...'
if not None in [Name1,Name2]:
  scores = {}
  sc = ScoringFunctions()
  #OVR
  try:
    ccc_mask,ovr = sc.CCC_map(emmap1,emmap2,c1,c2,3)
    print 'Percent overlap:', ovr
    if ovr < 0.0: ovr = 0.0
  except:
    print 'Exception for lccc and overlap score'
    print_exc()
    ovr = 0.0
  scores['overlap'] = ovr
  if ovr < 0.02:
    sys.exit("Maps do not overlap.")
  #SCCC
  print 'Local correlation score: ', ccc_mask
  if ccc_mask < -1.0 or ccc_mask > 1.0:
    ccc_mask = 0.0
  scores['local_correlation'] = ccc_mask
  #LMI
  try:
    mi_mask = sc.MI(emmap1,emmap2,c1,c2,3)
    print 'Local Mutual information score: ', mi_mask
    if mi_mask < 0.0: mi_mask = 0.0
  except:
    print 'Exception for MI score'
    print_exc()
    mi_mask = 0.0
  scores['mutual_information'] = mi_mask

  #NMI
  try:
    nmi = sc.MI(emmap1,emmap2,c1,c2,1,None,None,True)
    print 'Normalized Mutual information score:', nmi
    if nmi < 0.0: nmi = 0.0
  except:
    print 'Exception for NMI score'
    print_exc()
    nmi = 0.0
  #CD
  try:
    chm = sc._surface_distance_score(emmap1,emmap2,c1,c2,'Minimum')
    if chm == 0.0 or chm is None:
      chm = sc._surface_distance_score(emmap1,emmap2,c1,c2,'Mean')
    print 'Surface distance score: ', chm
    if chm < 0.0: chm = 0.0
  except:
    print 'Exception for surface distance score'
    print_exc()
    chm = 0.0
  scores['surface_distance'] = chm
  
  try:
    nv = sc.normal_vector_score(emmap1,emmap2,float(c1),float(c1)+(emmap1.std()*0.05),'Minimum')
    if nv == 0.0 or nv is None: nv = sc.normal_vector_score(emmap1,emmap2,float(c1),float(c1)+(emmap1.std()*0.05))
    print 'Normal vector score: ', nv
    if nv < 0.0: 
      nv = 0.0
  except:
    print 'Exception for NV score'
    print_exc()
    nv = 0.0
  scores['normal_vector'] = nv
