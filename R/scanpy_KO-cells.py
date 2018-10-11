print('Import...')
import numpy as np
import pandas as pd
import scanpy.api as sc


print('Reading...')
inpath = "/ye/yelabstore2/dosageTF/tfko_140/combined/nsnp20.raw.sng.guide_sng.norm.h5ad"
adata_singlet = sc.read(inpath)
print(type(adata_singlet.obs))

file_out = pd.DataFrame(data = adata_singlet.obs)

file_out.to_csv('KO_cells.csv')
print('Done')
