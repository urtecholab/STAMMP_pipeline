import pandas as pd
import warnings
from micom import interaction
from micom.workflows import GrowthResults
from sklearn.metrics.pairwise import cosine_distances, euclidean_distances
from typing import List, Optional, Literal
from sklearn.metrics.pairwise import pairwise_distances

warnings.filterwarnings('ignore')

def custom_sort(taxa: List[str]) -> List[str]:
    """
    Custom sort function for sorting the taxa

    Acceptable formats: 

    - String Format: k__Bacteria, p__Firmicutes, c__Clostridia, o__Clostridiales, f__Lachnospiraceae, g__Blautia
    - String and Numeric Format: AZC3J_1, AZC3J_2, AZC3J_3, AZC3J_4, AZC3J_5, AZC3J_6

    If both formats are detected in the input list, raises a ValueError.
    If only string format is provided, the list is sorted lexicographically.
    If string and numeric format is provided, it sorts by the prefix (before '_') 
    and then by the numeric part.

    :param taxa: list of taxa
    """
    string_only = all('_' not in x for x in taxa)
    string_and_numeric = all('_' in x and x.split('_')[1].isdigit() for x in taxa)

    if not (string_only or string_and_numeric):
        raise ValueError('Taxa list contains mixed types or invalid labels')
    if string_only:
        return sorted(taxa)
    elif string_and_numeric:
        return sorted(taxa, key=lambda x: (x.split('_')[0], int(x.split('_')[1])))

def metabolic_interaction_distance(donor_df: pd.DataFrame,
                                   recipient_df: pd.DataFrame,
                                   taxa: Optional[List[str]] = None, 
                                   metabolites: Optional[List[str]] = None,
                                   method: Literal['euclidean', 'cosine', 'mahalanobis'] = 'euclidean') -> pd.DataFrame:

    if method not in ['euclidean', 'cosine', 'mahalanobis']:
        raise ValueError('Please provide a valid method for calculating the distance')

    if not all([col in donor_df.columns \
                for col in ['metabolite', 'focal', 'partner', 'sample_id', 'flux', 'class']]) and \
            not all([col in recipient_df.columns \
                     for col in ['metabolite', 'focal', 'partner', 'sample_id', 'flux', 'class']]):
        raise ValueError('Please provide a valid interaction dataframe with metabolite, focal, partner, sample_id, '
                         'flux, and class columns')

    int_donor = donor_df.groupby(['metabolite',
                                        'focal', 'partner', 'class'])['flux'].mean().to_frame().reset_index().rename(
        columns={'focal': 'taxon'})
    int_recipient = recipient_df.groupby(['metabolite',
                                        'focal', 'partner', 'class'])['flux'].mean().to_frame().reset_index().rename(
        columns={'focal': 'taxon'})

    if metabolites is None:
        metabolites = sorted(pd.concat([int_donor['metabolite'], int_recipient['metabolite']], axis=0).unique())

    if taxa is None:
        donor_taxa = int_donor['taxon'].unique()
        recipient_taxa = int_recipient['taxon'].unique()
    else:
        donor_taxa = taxa
        recipient_taxa = taxa

    int_donor = int_donor.pivot_table(index='taxon',
                                      columns='metabolite',
                                      values='flux',
                                      fill_value=0.0,
                                      aggfunc='sum')
    
    int_donor = int_donor.reindex(index=donor_taxa, columns=metabolites, fill_value=0.0)

    int_recipient = int_recipient.pivot_table(index='taxon',
                                              columns='metabolite',
                                              values='flux',
                                              fill_value=0.0,
                                              aggfunc='sum')
    
    int_recipient = int_recipient.reindex(index=recipient_taxa, columns=metabolites, fill_value=0.0)


    if method == 'euclidean':
        distance_matrix = euclidean_distances(int_donor.values, int_recipient.values)
    elif method == 'cosine':
        distance_matrix = cosine_distances(int_donor.values, int_recipient.values)
    elif method == 'mahalanobis':
        combined = pd.concat([int_donor, int_recipient], axis=0)
        from sklearn.covariance import LedoitWolf
        cov = LedoitWolf().fit(combined)
        S_inv = cov.precision_
        distance_matrix = pairwise_distances(int_donor.values,
                                             int_recipient.values,
                                             metric='mahalanobis', VI=S_inv)
    else:
        raise ValueError('Please provide a valid method for calculating the distance')

    dist_df = pd.DataFrame(distance_matrix, index=int_donor.index, columns=int_recipient.index)
    return dist_df

if __name__ == '__main__':
    pass
