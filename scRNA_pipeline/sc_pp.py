#QC and doublet removal for large scRNA datas
import scanpy as sc
import sys, os
import argparse
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import doubletdetection
from scipy.stats import median_abs_deviation as mad

import warnings
warnings.simplefilter("ignore", FutureWarning)
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

def load_it(adata, info_df):
	samp = adata.split('/')[-1].split('_')[0]
	age, sex, dx = info_df.loc[samp, ['Age', 'Sex', 'dx']]
	adata = sc.read_10x_h5(adata)
	adata.obs['Age'] = age
	adata.obs['Sex'] = sex
	adata.obs['DX'] = dx
	adata.obs['Sample'] = samp
	adata.obs.index = adata.obs.index + '-' + samp
	adata.var_names_make_unique()
	return adata
	

def qc(adata):
	
	#you could also use a whitelist of barcodes from the filtered barcodes for each sample
	sc.pp.filter_cells(adata, min_genes = 200) #if you use the whitelist, you can get rid of this
	adata.var["mt"] = adata.var_names.str.startswith("MT-")
	adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
	adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))
	sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)
	
	remove = ['total_counts_mt', 'log1p_total_counts_mt', 'total_counts_ribo', 
				'log1p_total_counts_ribo', 'total_counts_hb', 'log1p_total_counts_hb']
	
	adata.obs = adata.obs[[x for x in adata.obs.columns if x not in remove]]

	return adata

def qc_vis(df, value, outputdir):
	sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
	g = sns.FacetGrid(df, row="Sample", hue="Sample", aspect=15, height=0.5, palette="tab20")
	g.map(sns.kdeplot, value, clip_on=False, fill=True, alpha=1, linewidth=1.5)
	g.map(sns.kdeplot, value, clip_on=False, color="w", lw=2)
	g.map(plt.axhline, y=0, lw=2, clip_on=False)
	
	def label(x, color, label):
		ax = plt.gca()
		ax.text(0, .2, label, fontweight="bold", color=color, 
				ha="left", va="center", transform=ax.transAxes)
	g.map(label, value)
	g.figure.subplots_adjust(hspace=-.6)
	g.set_titles("")
	g.set(yticks=[], ylabel="")
	g.despine(bottom=True, left=True)
	
	for ax in g.axes.flat:
		ax.axvline(x=df[value].median(), color='r', linestyle='-')
	
	
	plt.savefig(outputdir + value + '.png')

def mad_outlier(adata, metric, nmads, upper_only = False):
	M = adata.obs[metric]

	if not upper_only:
		return (M < np.median(M) - nmads * mad(M)) | (M > np.median(M) + nmads * mad(M))
	
	return (M > np.median(M) + nmads * mad(M))

def pp(adata, clf):
	adata = adata[adata.obs.pct_counts_mt < 25]
	
	bool_vector = mad_outlier(adata, 'log1p_total_counts', 5) +\
		mad_outlier(adata, 'log1p_n_genes_by_counts', 5) +\
		mad_outlier(adata, 'pct_counts_in_top_20_genes', 5) +\
		mad_outlier(adata, 'pct_counts_mt', 3, upper_only = True)
	
	adata = adata[~bool_vector]
	
	adata.uns['cells_removed'] = sum(bool_vector)
	
	doublets = clf.fit(adata.X).predict(p_thresh=1e-3, voter_thresh=0.5)
	doublet_score = clf.doublet_score()
	
	adata.obs['doublet'] = doublets
	adata.obs['doublet_score'] = doublet_score
	
	adata.uns['doublets_removed'] = adata.obs.doublet.sum()
	adata = adata[adata.obs.doublet == 0]
	
	return adata

def make_unique(names):
	seen = {}
	result = []
	for name in names:
		if name in seen:
			seen[name] += 1
			result.append(f"{name}_{seen[name]}")
		else:
			seen[name] = 0
			result.append(name)
	return result

def main():

	inputdir = ag.inputdir
	outputdir = ag.outputdir
	os.makedirs(outputdir, exist_ok=True)
	info_df = pd.read_csv(ag.info, sep='\t')
	info_df.set_index('Sample_ID', inplace=True)
	## Load AnnData files
	adatas = [inputdir + '/' + x for x in os.listdir(inputdir) if x.endswith('filtered.h5')]
	adatas = [load_it(ad, info_df) for ad in adatas if ad.split('/')[-1].split('_')[0] in info_df.index]

	print("AnnData", len(adatas), "Samples Loaded...\n")
	adatas = [qc(ad) for ad in adatas]
	
	## Visualize Qualities of each samples
	df = pd.concat([x.obs for x in adatas])
	df.sort_values('Sample')
	
	soutputdir = outputdir + '/bfpp_qc/'
	os.makedirs(soutputdir, exist_ok=True)
	qc_vis(df, 'pct_counts_mt', soutputdir)
	qc_vis(df, 'n_genes', soutputdir)
	qc_vis(df, 'pct_counts_in_top_20_genes', soutputdir)
	qc_vis(df, 'log1p_total_counts', soutputdir)

	## Detect and remove doublets by using DoubletDetection
	clf = doubletdetection.BoostClassifier(
		n_iters=10, 
		clustering_algorithm="louvain",
		standard_scaling=True,
		pseudocount=0.1,
		n_jobs=-1)
	
	adatas = [pp(ad, clf) for ad in adatas]
	
	with open(outputdir + '/pp.stat', 'w') as f:
		f.write('Total\tCell_removed\tDoublet_removed')
		for adata in adatas:			
			f.write('\n'+'\t'.join(map(str, [len(adata), adata.uns['cells_removed'], adata.uns['doublets_removed']])))

	## Visualize Qualities after pp
	df2 = pd.concat([x.obs for x in adatas])
	df2.sort_values('Sample')
	
	soutputdir = outputdir + '/afpp_qc/'
	os.makedirs(soutputdir, exist_ok=True)
	qc_vis(df2, 'pct_counts_mt', soutputdir)
	qc_vis(df2, 'n_genes', soutputdir)
	qc_vis(df2, 'pct_counts_in_top_20_genes', soutputdir)
	qc_vis(df2, 'log1p_total_counts', soutputdir)

	# Check and resolve duplicates in variable names
	for i, adata in enumerate(adatas):
		if adata.var.index.duplicated().any():
			print(f"Resolving duplicates in variable names of adata{i+1}")
			adata.var.index = make_unique(adata.var.index)

	os.makedirs(outputdir+'/pp_adata/', exist_ok=True)
	for adata in adatas:
		samp = adata.obs['Sample'][0]
		adata.write_h5ad(outputdir+'/pp_adata/'+samp+'_pp.h5ad')

if __name__ == '__main__':
	prs = argparse.ArgumentParser(description = 'Preprocessing of scRNA-seq data by using Scanpy and DoubletDetection method')
	prs.add_argument('-I', '--inputdir', help = 'Input directory contains h5 file(s)', dest = 'inputdir')
	prs.add_argument('-O', '--outputdir', help = 'Output directory', default = '~/sc_pp/', dest = 'outputdir')
	prs.add_argument('-i', '--infodir', help = 'Directory contains Information tables of the data', dest = 'info')
	
	ag = prs.parse_args()
	main()
