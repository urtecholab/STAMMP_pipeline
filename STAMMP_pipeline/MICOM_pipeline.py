import pandas as pd
import numpy as np
import warnings
from micom.qiime_formats import load_qiime_medium
from micom.workflows import grow, build
from micom.interaction import interactions
from micom.workflows import load_results, GrowthResults
import argparse


import cobra
print("Default solver:", cobra.Configuration().solver)

import sys
import os

sys.path.append('..')
from utils import analysis

warnings.filterwarnings('ignore')


def main():
	# part 0 argument parsing
	parser = argparse.ArgumentParser()
	parser.add_argument(
        "-a", "--agora",
        required=True,
        help="Path to the AGORA CSV file"
    )
	parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path to where you want the model to be saved"
    )
	parser.add_argument(
        "-r", "--reference",
        required=True,
        help="Path to reference taxonomy files and genome sizes"
    )
	args = parser.parse_args()

	print("creating :: %s/micom_outputs/processed"%(args.output))
	print(os.path.dirname("%s/micom_outputs/models"%(args.output)))
	
	os.makedirs("%s/micom_outputs/processed"%(args.output), exist_ok=True)
	os.makedirs("%s/micom_outputs/models"%(args.output), exist_ok=True)
	# part 1 data loading
	reference_taxonomy = pd.read_csv('%s/reference_taxonomy.csv'%(args.reference)).drop(columns=['Unnamed: 0'])
	agora_genomes = pd.read_csv('%s/genome_size_ungapped.csv'%(args.reference)).drop(columns=['Unnamed: 0'])
	reference_taxonomy_counts = pd.merge(reference_taxonomy, agora_genomes, on='id', how='right')
	reference_taxonomy_counts['NCBI Taxonomy ID'] = reference_taxonomy_counts['NCBI Taxonomy ID'].astype(int).astype(str)
	abundance = pd.read_csv(args.agora) # output from agora pipeline is the abundance file
	abundance['taxon_id'] = abundance['taxon_id'].astype(str)

	cols_to_keep = ['taxon_id'] + [col for col in abundance.columns if col.startswith('D') or col.endswith('-0') or col.endswith('-6') ]
	#abundance = abundance[cols_to_keep]
	abundance = abundance.merge(reference_taxonomy_counts, left_on='taxon_id', right_on='NCBI Taxonomy ID', how='left')
	abundance.set_index('taxon_id', inplace=True)
	meta_cols = [
		'id', 'NCBI Taxonomy ID', 'kingdom', 'phylum', 'class', 'order',
		'family', 'genus', 'species', 'RefSeq', 'organism_name', 'genome_size_ungapped'
	]

	abundance['genome_size_ungapped'] = pd.to_numeric(abundance['genome_size_ungapped'], errors='coerce')
	data_cols = abundance.columns.difference(meta_cols)
	abundance[data_cols] = abundance[data_cols].div(abundance['genome_size_ungapped'], axis=0)
	abundance = abundance.drop(columns=['RefSeq', 'organism_name', 'genome_size_ungapped', 'NCBI Taxonomy ID']).reset_index()

	abundance.drop(columns=['kingdom', 'phylum', 'class', 
							'order', 'family', 
							'genus', 'species', 
							'id']).set_index('taxon_id').to_csv('%s/micom_outputs/processed/abundance_normalized_table.csv'%(args.output)) # edited

	abundance = abundance.drop(columns=['id']).melt(id_vars=['taxon_id', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'], var_name='sample_id', value_name='abundance').rename(columns={'taxon_id': 'id'})
	abundance['sample_id'] = abundance['sample_id'].astype(str)
	abundance.to_csv('%s/micom_outputs/processed/abundance_normalized.csv'%(args.output)) # edited # saving agora abundance file after normalization
	abundance

	# part 2 model building
	metadata = build(abundance, 
					model_db='%s/agora_species_103.qza'%(args.reference), # agora species reference
					out_folder='%s/micom_outputs/models'%(args.output),
					cutoff=0.0001,
					threads=12)
					
	# part 3 model growing and saving results
	western = load_qiime_medium('%s/western_diet_gut_agora.qza'%(args.reference)) # agora western diet reference
	res = grow(manifest=metadata,
			medium=western,
			model_folder='%s/micom_outputs/models'%(args.output),
			threads=12,
			tradeoff=0.7,
			# ---- stability knobs ----
			presolve=True,        # let the solver simplify the QP before solving
			atol=1e-8,            # absolute flux tolerance (bigger number = more permissive)
			rtol=1e-6             # relative tolerance   
			)
	res.save('%s/micom_outputs/models/res_western.zip'%(args.output))

if __name__ == "__main__":
    main()