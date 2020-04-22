####### Libraries #######
from itertools import product
from pathlib import Path
import glob
import json
import os

####### Util functions #######
def expand_list(A,B):
    return([''.join(r) for r in product(A,B)])

def extractFilenames(fullnames,suffix):
    names = []
    for file in fullnames:
        names.append(os.path.basename(file).split(suffix)[0])
    return sorted(names)

def findLibraries(path,prefix,suffix):
    filenames_path = glob.glob(os.path.join(path,prefix) + "*" + suffix)
    names = []
    for file in filenames_path:
        library = os.path.basename(file).split(suffix)[0]
        if(library not in names):
            names.append(library)
    return sorted(names)

def loadGenome(ref):
    if(not ref.endswith(".json")):
        raise ValueError("expecting file with .json format for the reference genome")
    FA = None
    GTF = None
    with open(ref) as genome_data:
        data = json.load(genome_data)
        for i in data.keys():
            if(i.endswith(".fa.gz")):
                FA = Path(i).stem
            elif(i.endswith(".gtf.gz") or i.endswith(".gff3.gz")):
                GTF = Path(i).stem
    if((FA is None) or (GTF is None)):
        raise ValueError("reference genome file wrongly formatted")
    return FA, GTF

def which(file):
    for path in os.environ["PATH"].split(os.pathsep):
        if os.path.exists(os.path.join(path, file)):
            return path
    return None
