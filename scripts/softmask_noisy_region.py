from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
from pybedtools import BedTool
from scipy.stats import zscore
import numpy as np
from sklearn.linear_model import LinearRegression
import warnings
import gzip

plt.style.use('seaborn-white')
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams["image.cmap"] = "Dark2"
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.cm.Dark2.colors)

# A small function to ensure proper header format for repeat file. 
def read_gzipped_tsv(file_path):

    # Define the header for the repeat file. 
    repeat_header = ['#bin', 'swScore', 'milliDiv', 'milliDel', 'milliIns', 'genoName',
                       'genoStart', 'genoEnd', 'genoLeft', 'strand', 'repName', 'repClass',
                       'repFamily', 'repStart', 'repEnd', 'repLeft', 'id']

    with gzip.open(file_path, 'rt') as f:
        first_line = f.readline().strip()

    # Check if the first line starts with the expected header's first column name
    if first_line.startswith(repeat_header[0]):
        # Header is present, load normally with compression specified
        df = pd.read_csv(file_path, compression='gzip', header=0, sep = "\t")
    else:
        # Header is not present, load without header and assign expected column names
        df = pd.read_csv(file_path, compression='gzip', header=None, names=repeat_header, sep = "\t")
    return df

def fit_regression(raw_counts, total_reads, pseudocount_nread = 10, rsquare_threshold = 0.8):
    ''' given raw counts, perform regression and return deviation z-score '''
    
    individual_coef_df = []
    y_dev_df = []
    pseudocount = total_reads.div(total_reads.sum())*pseudocount_nread
    
    for index, row in raw_counts.iterrows():
        
        frac_table = pd.concat([total_reads, 
                                (row+pseudocount)/total_reads],
                  axis = 1)
        
        frac_table.columns = ['#total reads', f'fraction']

        X=np.log(frac_table[['#total reads']])
        y=np.log(frac_table[[f'fraction']])
        reg = LinearRegression().fit(X,y)
        reg_score=reg.score(X,y)
        coef=reg.coef_
        intercept=reg.intercept_

        if reg_score > rsquare_threshold and coef < 0:
            # detect outlier
            y_hat = reg.predict(X)
            y_dev = y-y_hat
            frac_table[index]=y_dev

            y_dev_df.append(frac_table[index])

        individual_coef_df.append([index, reg_score, coef[0,0], intercept[0]])
    individual_coef_df = pd.DataFrame(individual_coef_df, columns = ['name', 'R_square', 'coef', 'intercept'])
    try:
        y_dev_df = pd.concat(y_dev_df, axis = 1).T
        
        # calculate deviation z-score
        deviation_zscore = np.reshape(zscore(y_dev_df.values.flatten()), y_dev_df.shape)
        deviation_zscore = pd.DataFrame(deviation_zscore, index = y_dev_df.index, columns = y_dev_df.columns)
    except:
        deviation_zscore = pd.DataFrame()
        warnings.warn("No family has good regression line")
    
    return deviation_zscore, individual_coef_df


def generate_repeat_family_mask(repeat_family_to_mask, repeat_annotation, total_reads,
                                raw_counts, 
                                pseudocount_nread = 10, rsquare_threshold = 0.8):
    '''
    Generate repeat family mask by checking all repeat families that has a high noise type.
    For each "repeat name", perform the regression, and save the deviation from the regression line, 
    as the standard to be masked or not
    '''
    repeat_name_to_mask = repeat_annotation.loc[repeat_annotation['repFamily'].isin(repeat_family_to_mask), 'repName'].unique()
    
    deviation_zscore, individual_coef_df = fit_regression(raw_counts.loc[raw_counts.index.isin(repeat_name_to_mask)]
                                                          , total_reads, pseudocount_nread = pseudocount_nread,
                                                         rsquare_threshold = rsquare_threshold)
    
    return deviation_zscore, individual_coef_df
    


    
def generate_genome_mask(genome_transcript_types_to_mask, genomic_annotation, total_reads, 
                                genome_counts,
                                pseudocount_nread = 10, rsquare_threshold = 0.8):
    '''
    Generate genome window mask by checking all transcript types that has a high noise type.
    For each "window", perform the regression, and save the deviation from the regression line, 
    as the standard to be masked or not
    '''
    windows_to_mask = genomic_annotation.loc[genomic_annotation['transcript_type_top'].isin(genome_transcript_types_to_mask), 'name']
    
    deviation_zscore, individual_coef_df = fit_regression(genome_counts.loc[windows_to_mask]
                                                          , total_reads, pseudocount_nread = pseudocount_nread,
                                                         rsquare_threshold = rsquare_threshold)
    return deviation_zscore, individual_coef_df

if __name__=='__main__':
    basedir=  Path(sys.argv[1])
    out_stem = Path(sys.argv[2])
    rsquare_threshold = 0.3 # how good the regression has to be
    zscore_cutoff = 1 # how much does the rbp-rna interaction have to deviate from the noise regression to not be masked
    genomic_annotation = pd.read_csv(Path(sys.argv[3]),
                                    sep = '\t')
    repeat_annotation = read_gzipped_tsv(Path(sys.argv[4]))
    
    ##### read in all the files #####
    # counts
    raw_counts = pd.read_csv(basedir/'counts'/'repeats'/'megatables'/'name'/f'{out_stem}.tsv.gz', sep = '\t', index_col = 0)
    genome_counts = pd.read_csv(basedir/'counts'/'genome'/'megatables'/f'{out_stem}.tsv.gz', sep = '\t')
    genome_counts.index = genome_counts.index+1

    total_reads = raw_counts.sum(axis=0)+genome_counts.sum(axis = 0)

    ##### fit for repeat families
    rep_family_counts = []
    for family in repeat_annotation['repFamily'].unique():
        cnt = raw_counts.loc[raw_counts.index.isin(repeat_annotation.loc[repeat_annotation['repFamily']==family, 'repName'])
                                        ].sum(axis = 0)
        cnt.name = family
        rep_family_counts.append(cnt)
    rep_family_counts = pd.concat(rep_family_counts, axis = 1).T
    _, rep_coef_df = fit_regression(rep_family_counts, total_reads, rsquare_threshold = rsquare_threshold)

    ##### fit for genomic transcript types

    genome_counts['transcript_type_top']=genomic_annotation.set_index('name')['transcript_type_top']
    transcript_type_counts = genome_counts.groupby(by = 'transcript_type_top').sum()
    genome_counts.drop('transcript_type_top', inplace = True, axis = 1)
    _, genome_coef_df = fit_regression(transcript_type_counts, total_reads, rsquare_threshold = rsquare_threshold)

    ##### plot & save results
    f, axes = plt.subplots(2,1, figsize = (3,6), sharex = True)
    genome_coef_df.loc[genome_coef_df['coef']<0].set_index('name')['R_square'].sort_values().iloc[-30:].plot.barh(
        ax = axes[0], color = 'lightgrey')
    rep_coef_df.loc[rep_coef_df['coef']<0].set_index('name')['R_square'].sort_values().iloc[-30:].plot.barh(
        ax = axes[1], color = 'lightgrey')
    axes[0].set_title('genomic transcript types')
    axes[1].set_title('repeat repFaily')
    _ = [ax.vlines(x=rsquare_threshold, ymin = 0, ymax = 30, color = 'grey', linestyle = '--') for ax in axes]
    plt.xlabel('regression R squared')
    sns.despine()
    plt.savefig(basedir / 'mask' / f'{out_stem}.category_coefficient.pdf')

    genome_coef_df['type']='genome'
    rep_coef_df['type']='repeat'
    pd.concat([genome_coef_df, rep_coef_df], axis = 0).to_csv(basedir / 'mask' / f'{out_stem}.category_coefficient.csv') # this saves the regression outputs


    repeat_family_to_mask = rep_coef_df.loc[(rep_coef_df['coef']<0)&(rep_coef_df['R_square'].ge(rsquare_threshold)),'name'].tolist()
    genome_transcript_types_to_mask = genome_coef_df.loc[(genome_coef_df['coef']<0)&(genome_coef_df['R_square'].ge(rsquare_threshold)),'name'].tolist()
    print(f'These family will be soft-masked: {repeat_family_to_mask} {genome_transcript_types_to_mask}')

    ##### flagged high noise family, test individual window
    # get individual repeat name coefficients and y_dev
    zscore_dev_df, individual_coef_df = generate_repeat_family_mask(repeat_family_to_mask, repeat_annotation, total_reads,
                                    raw_counts,  rsquare_threshold=rsquare_threshold)
    fraction_total_masked = zscore_dev_df.le(2).all(axis = 1).sum()/individual_coef_df.shape[0]
    repeats_total_masked = zscore_dev_df.loc[zscore_dev_df.le(zscore_cutoff).all(axis = 1)].index
    print(f'''
            =======================================
            fraction of repeat windows totally masked (no signal): {fraction_total_masked:.2f}.
            
            Here are the repName: {repeats_total_masked}
            '''
        )
    zscore_dev_df_genome, individual_coef_df_genome = generate_genome_mask(
        genome_transcript_types_to_mask, genomic_annotation, total_reads,
                                    genome_counts, rsquare_threshold=rsquare_threshold)
    fraction_total_masked = zscore_dev_df_genome.le(2).all(axis = 1).sum()/individual_coef_df_genome.shape[0]
    genomic_total_masked = genomic_annotation.loc[
        genomic_annotation['name'].isin(
        zscore_dev_df_genome.loc[zscore_dev_df_genome.le(zscore_cutoff).all(axis = 1)].index)
    , 'gene_name'].unique()
    print(f'''
            =======================================
            fraction of genomic windows totally masked (no signal): {fraction_total_masked:.2f}.
            
            Here are the gene names: {genomic_total_masked}
            '''
        )

    ##### Save
    zscore_dev_df_genome.to_csv(basedir / 'mask' / f'{out_stem}.genome_deviation_zscore.csv') # this saves the regression outputs
    zscore_dev_df_genome.ge(zscore_cutoff).to_csv(basedir / 'mask' / f'{out_stem}.genome_mask.csv') # this saves the regression outputs
    individual_coef_df_genome.to_csv(basedir / 'mask' / f'{out_stem}.genome_window_regcoef.csv') # this saves the regression outputs

    zscore_dev_df.to_csv(basedir / 'mask' / f'{out_stem}.repeat_deviation_zscore.csv') # this saves the regression outputs
    individual_coef_df.to_csv(basedir / 'mask' / f'{out_stem}.repeat_window_regcoef.csv') # this saves the regression outputs
    zscore_dev_df.ge(zscore_cutoff).to_csv(basedir / 'mask' / f'{out_stem}.repeat_mask.csv') # this saves the regression outputs