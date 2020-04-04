####### Libraries #######
from itertools import product
import glob
import os

####### Util functions #######
def expand_list(A,B):
    return([''.join(r) for r in product(A,B)])

def extractFilenames(fullnames,suffix):
    names = []
    for file in fullnames:
        names.append(os.path.basename(file).split(suffix)[0])
        return sorted(names)

def fastqc(dir,libs,format,ends=[],extra=[]):
    ends = expand_list(["_"],ends)
    suffix = expand_list(["_fastqc."],format)
    if(len(ends) > 0):
        suffix = expand_list(ends,suffix)
    if(len(extra) > 0):
        suffix = expand_list(extra,suffix)
    return(expand_list(dir,expand_list(libs,suffix)))

def findLibraries(path,prefix,suffix):
    filenames_path = glob.glob(os.path.join(path,prefix) + "*" + suffix)
    names = []
    for file in filenames_path:
        library = os.path.basename(file).split(suffix)[0]
        if(library not in names):
            names.append(library)
    return sorted(names)

def which(file):
    for path in os.environ["PATH"].split(os.pathsep):
        if os.path.exists(os.path.join(path, file)):
            return path
    return None
