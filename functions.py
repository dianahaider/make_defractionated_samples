import pandas as pd
import numpy as np

def import_docs(comm):
    #import your ASV biom table
    biom = pd.read_csv('ASVs_'+comm+'.csv')

    #import your metadata
    md = pd.read_csv('Metadata_OG.csv') #take the habit of avoiding spaces in file names

    #melt your ASV table to attach dna metadata to ASVs
    biom.rename(columns={'#OTU ID': 'feature_id'}, inplace=True)
    md.rename(columns={'ID': 'sample_id'}, inplace=True)


    biomelted = biom.melt(id_vars=['feature_id'], var_name='sample_id', value_name='feature_frequency')
    
    return biomelted, md

def make_defract(comm, biomelted, md_SF):
    #remove rows where samples where not size fractionated, i.e. base/tray/top
    md_SF = md_SF[md_SF['SizeFraction2'].isin(['S', 'L'])]

    #make a new column of total [DNA] per sample that were size fractionated and need to be pooled
    md_SF['[DNAt]'] = md_SF.groupby(['ID3'])['Concentration'].transform('sum')

    #separate small and large size fraction
    sep_S = md_SF[md_SF.SizeFraction2 == 'S']
    sep_L = md_SF[md_SF.SizeFraction2 == 'L']

    #calculate DNA proportion per size fraction
    md_SF['DNApr'] = md_SF['Concentration']/md_SF['[DNAt]']

    #merge with separated on common columns to get corresponding rel. abundances
    #md_SF = md_SF[['sample_id', 'DNApr', '[DNAt]']].copy()
    merged = pd.merge(biomelted, md_SF, on=['sample_id'], how='left') #all_md is the metadata file

    #remove the ASVs with a null read count
    sepSLRA = merged[merged.feature_frequency != 0]

    #calculate corrected per sample ratio, and corrected feature frequency of de-fractionated samples
    sepSLRA['Newfeature_frequency'] = sepSLRA['feature_frequency'] * sepSLRA['DNApr']
    sepSLRA['Newff'] = sepSLRA.groupby(['feature_id', 'ID3'])['Newfeature_frequency'].transform('sum')

    #remove the rows where there was no size fractionation (base, tray, top..)
    sepSLRA = sepSLRA[sepSLRA['Newff'].notna()]

    #make a new id for the new combined samples
    sepSLRA['sampleid'] = sepSLRA['ID3'].astype(str) + "SL"

    #uncomment the line above if merging smallandlarge
    sepSLRA['SizeFraction'] = 'SL'

    #rename the columns
    sepSLRA.rename(columns={'feature_frequency':'old_feature_frequency'}, inplace=True)
    sepSLRA.rename(columns={'Newff':'feature_frequency'}, inplace=True)
    sepSLRA = sepSLRA.drop_duplicates()

    #recalculate ratios
    sepSLRA['Total'] = sepSLRA['feature_frequency'].groupby(sepSLRA['sampleid']).transform('sum')
    sepSLRA['ratio'] = sepSLRA['feature_frequency']/sepSLRA['Total'] #calculate the relative abundance of a feature (0-1 scale per sample)
    sepSLRA['nASVs'] = sepSLRA['feature_id'].groupby(sepSLRA['sampleid']).transform('nunique') #calculate the number of ASVs per sample

    sepSLRA = sepSLRA.drop_duplicates()

    #make a new biom table
    newbiom = sepSLRA[['sampleid', 'feature_id', 'feature_frequency']].copy()
    newbiom.drop_duplicates(inplace=True)
    newbiom = newbiom.pivot(index='feature_id', columns='sampleid', values='feature_frequency')
    newbiom = newbiom.fillna(0)
    
    #save outputs to csv
    newbiom.to_csv('newbiom_'+comm+'.csv')
    sepSLRA.to_csv('CombinedSL_metadata_'+comm+'.csv')
    
    return merged, sepSLRA, newbiom