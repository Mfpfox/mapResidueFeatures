
import sys
import os
sys.path.append("/Users/mariapalafox/Desktop/TOOLBOXPY")
from all_funx import *

# input: 
    # features_summary_QCd_sequences.csv
# output:
    # features_summary_QCd_02122022_proteinLevelColumns.csv
    # features_summary_QCd_NotIdenticalSequence_proteins.csv
    # features_summary_QCd_IdenticalSeq_proteins.csv

dirname = "Uniprot_Domains/20220212/"

df = pd.read_csv(dirname+"features_summary_QCd_sequences.csv")

# save protein level file
#  # protein -level merge later, 'Protein', 'Gene', 'Length','SubunitStructure', '3D', 'PDB', 'InteractsWith', 'PTM', 
proteinlevel = df[['UKBID', 
         'identicalToRefSequence',
           'query.minus.REF.difference', 
             'QC.posID.difference',
                 'Protein', 'Gene', 'Length','SubunitStructure',
               '3D', 'PDB', 'InteractsWith', 'PTM']].copy()

print('saving features_summary_QCd_proteinLevelColumns.csv in dated subfolder..')
proteinlevel.to_csv("features_summary_QCd_02122022_proteinLevelColumns.csv", index=False)

# simplifly df drop Sequence col
df = df[['UKBID', 
         'identicalToRefSequence',
           'query.minus.REF.difference', 
             'QC.posID.difference',
                 'Protein', 'Gene', 'Length','SubunitStructure',
               '3D', 'PDB', 'InteractsWith', 'PTM', 
             'SubstrateBindingSite',
             'ActiveSite', 
             'MetalBindSite', 
             'ModifiedResidue',
             'Lipidation',
             'Glycosylation', 
             'Mutagenesis', 
                'DisulfideBond',
    'NucleotideBindRegion', 'DNAbindRegion',
   'ZincFinger', 'Motif', 
         'DisorderedRegion', 'Domain', 
   'BetaStrand', 'Helix', 
         'Turn'
         ]].copy()

#qc
checkColumnValue(df, 'identicalToRefSequence')

# save dropped protein annotation rows 
falsedf =  df[~df['identicalToRefSequence']].copy()
print('saving features_summary_QCd_NotIdenticalSequence_proteins.csv - data shape:', falsedf.shape)
falsedf.to_csv("features_summary_QCd_NotIdenticalSequence_proteins.csv", index=False)
# saving falsedf , sequences dropped because they differed from the reference proteome

# drop rows with different sequences between query and ref for cpdaa
df =  df[df['identicalToRefSequence']].copy()
print('saving features_summary_QCd_IdenticalSeq_proteins.csv - data shape:', df.shape)
df.to_csv("features_summary_QCd_IdenticalSeq_proteins.csv", index=False)

