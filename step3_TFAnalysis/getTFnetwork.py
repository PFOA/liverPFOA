
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from MulticoreTSNE import MulticoreTSNE as TSNE
import json
import base64
import zlib
import glob
import pickle
import pandas as pd
import numpy as np
from dask.diagnostics import ProgressBar
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2
from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
import  os,sys






TF=sys.argv[1].rstrip()
adjacencies = pd.read_csv("adj.csv", index_col=False, sep=',')
lf = lp.connect("Samplef_loom_path_scenic_output.loom", mode='r', validate=False)
meta = json.loads(zlib.decompress(base64.b64decode(lf.attrs.MetaData)))
exprMat = pd.DataFrame(lf[:, :], index=lf.ra.Gene, columns=lf.ca.CellID).T
modules = list(modules_from_adjacencies(adjacencies, exprMat))
tf = TF
tf_mods = [x for x in modules if x.transcription_factor == tf]

regulons = {}
for i,r in pd.DataFrame(lf.ra.Regulons,index=lf.ra.Gene).iteritems():
    regulons[i] =  list(r[r==1].index.values)
for i, mod in enumerate(tf_mods):
    print(f'{tf} module {str(i)}: {len(mod.genes)} genes')

# write these modules, and the regulon to files:
for i, mod in enumerate(tf_mods):
    with open(tf + '_module_' + str(i) + '.txt', 'w') as f:
        for item in mod.genes:
            f.write("%s\n" % item)

with open(tf + '_regulon.txt', 'w') as f:
    for item in regulons[tf + '(+)']:
        f.write("%s\n" % item)