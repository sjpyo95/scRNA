import scanpy as sc
import os, sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scvi

import warnings
warnings.simplefilter("ignore", FutureWarning)
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

import celltypist
from celltypist import models

def load_MPAL_ref(rds_file):
	import anndata2ri
	from rpy2.robjects import r
	from rpy2.robjects import pandas2ri
	import rpy2.robjects.packages as rpackages

	# Activate the conversion between AnnData and SingleCellExperiment
	anndata2ri.activate()

	# Activate pandas to R dataframe conversion
	pandas2ri.activate()

	# Import the necessary R packages
	seurat = rpackages.importr('Seurat')
	summarized_experiment = rpackages.importr('SummarizedExperiment')

	# Pass the file path to R and read the RDS file
	r(f'rse <- readRDS("{rds_file}")')

	# Create Seurat object
	r('seurat_object <- CreateSeuratObject(counts = assay(rse, "counts"), meta.data = as.data.frame(colData(rse)))')

	# Convert to SingleCellExperiment
	r('healthy <- as.SingleCellExperiment(seurat_object)')

	# Retrieve the SingleCellExperiment object as a pandas DataFrame (or AnnData object)
	healthy = r['healthy']
	
#	sc.pp.normalize_total(healthy, target_sum = 1e4) #Note this is only for cell annotation, recommended by authors but not best
#	sc.pp.log1p(healthy)
	
#	ref_model = celltypist.train(healthy, labels = 'BioClassification', n_jobs = 22,use_SGD = False, feature_selection = True, top_genes = 300)
#	ref_model.write('"/mnt/mone/PriFiles/sunyme95/singlecell/ref/ref.pkl')
	print('ref model loaded...\n')
	return healthy


## Load reference models from celltypist database (Immune_All_Low, Immune_All_High)  and created reference model (MPAL-Single-Cell Hematopoeisis cell)
model_low = models.Model.load(model="Immune_All_Low.pkl")
model_hi = models.Model.load(model="Immune_All_High.pkl")
ref_model = models.Model.load(model="/mnt/mone/PriFiles/sunyme95/singlecell/ref/ref.pkl")	#If you don't have it. Create it from load_MPAL_ref() function

def predict_cells(adata):
	sc.pp.filter_genes(adata, min_cells = 10)
	sc.pp.normalize_total(adata, target_sum=1e4) #not recommended for typical pp
	sc.pp.log1p(adata)

	adata.X = adata.X.toarray()
	
	predictions = celltypist.annotate(adata, model=model_low, majority_voting=False)
	predictions_adata = predictions.to_adata()
	adata.obs["low_label"] = predictions_adata.obs.loc[adata.obs.index, "predicted_labels"]
	adata.obs["low_score"] = predictions_adata.obs.loc[adata.obs.index, "conf_score"]
	
	predictions = celltypist.annotate(adata, model=model_hi, majority_voting=False)
	predictions_adata = predictions.to_adata()
	adata.obs["hi_label"] = predictions_adata.obs.loc[adata.obs.index, "predicted_labels"]
	adata.obs["hi_score"] = predictions_adata.obs.loc[adata.obs.index, "conf_score"]

	predictions = celltypist.annotate(adata, model=ref_model, majority_voting=False)
	predictions_adata = predictions.to_adata()
	adata.obs["ref_label"] = predictions_adata.obs.loc[adata.obs.index, "predicted_labels"]
	adata.obs["ref_score"] = predictions_adata.obs.loc[adata.obs.index, "conf_score"]
	
	#if more reference add

	return adata.obs



def main():
	inputdir = ag.inputdir
	outputdir = ag.outputdir
	os.makedirs(outputdir, exist_ok=True)
	rds_file = ag.rds

	healthy = load_MPAL_ref(rds_file)
	adatas = [sc.read_h5ad(inputdir + x) for x in os.listdir(inputdir) if x.endswith('.h5ad')]

	print('AnnData loaded...\n')

	predictions = [predict_cells(ad.copy()) for ad in adatas]
	predictions = pd.concat(predictions)[['low_label', 'low_score', 'hi_label', 'hi_score', 'ref_label', 'ref_score']]

	adata = sc.concat(adatas)

	predictions.to_csv(outputdir + '/PREDICTIONS.csv')
	
	print('PREDICTIONS.csv\n')
	##scVI label transfer

	adata.obs['CellType'] = 'Unknown'
	adata.obs['Batch'] = 'KOR'

	healthy.obs['CellType'] = healthy.obs['BioClassification']
	healthy.obs['Batch'] = 'ref'
	healthy.obs['Sample'] = healthy.obs.index.map(lambda x: x.split(':')[0])

	dater = sc.concat((adata, healthy))

	sc.pp.highly_variable_genes(dater, flavor = 'seurat_v3', n_top_genes=2000, batch_key='Batch', subset = True)

	print(dater)
	scvi.model.SCVI.setup_anndata(dater, batch_key='Batch', categorical_covariate_keys = ['Sample'])
	vae = scvi.model.SCVI(dater)
	vae.train()

	lvae = scvi.model.SCANVI.from_scvi_model(vae, adata = dater, unlabeled_category = 'Unknown',
											labels_key = 'CellType')
	
	lvae.train(max_epochs=20, n_samples_per_label=100)

	dater.obs['predicted'] = lvae.predict(dater)
	dater.obs['transfer_score'] = lvae.predict(soft = True).max(axis = 1)
	dater = dater[dater.obs.Batch == 'KOR']
	adata.obs = adata.obs.merge(right = dater.obs[['predicted', 'transfer_score']], left_index=True, right_index=True)

	adata.obs = adata.obs.merge(right = predictions, left_index=True, right_index=True)
	
	adata.write_h5ad(outputdir + '/unintigrated.h5ad')
	print('unintigrated.h5ad file generated')

if __name__ == '__main__':
	prs = argparse.ArgumentParser(description = 'Integration and automated annotation by using celltypist')
	prs.add_argument('-I', '--inputdir', help = 'Input directory of h5ad files', dest = 'inputdir')
	prs.add_argument('-O', '--outputdir', help = 'Output directory', default = '~/sc_anno/', dest = 'outputdir')
	prs.add_argument('-R', '--refernce', help = 'RDS file to add reference model', dest = 'rds')

	ag = prs.parse_args()

	main()
