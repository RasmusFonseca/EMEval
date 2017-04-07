# EMEval

Evaluate fit of PDB/mmCif to a Cryo-EM map using the [TEMPy](http://tempy.ismb.lon.ac.uk/overview.html) library. Make sure TEMPy is installed and simply run
```bash
$ python2.7 emeval.py map-file.mrc structure-file.pdb [structure-file2.pdb ...]
map-file.mrc structure-file.pdb 10.22
map-file.mrc structure-file2.pdb 10.22
...
```
Each stdout line indicates a map-file, structure file-name, and the corresponding cross-correlation goodness-of-fit.

