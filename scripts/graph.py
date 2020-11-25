import pandas as pd
import plotly.express as px
from IPython.display import display

# ======================================================================================================================
data = pd.read_csv('ITPA_scores.csv')

def mutate(mutation: str):
    # p.Ala33Orn --> A33O
    name3 = {'A': 'ALA',
          'C': 'CYS',
          'D': 'ASP',
          'E': 'GLU',
          'F': 'PHE',
          'G': 'GLY',
          'H': 'HIS',
          'I': 'ILE',
          'L': 'LEU',
          'K': 'LYS',
          'M': 'MET',
          'N': 'ASN',
          'P': 'PRO',
          'Q': 'GLN',
          'R': 'ARG',
          'S': 'SER',
          'T': 'THR',
          'V': 'VAL',
          'W': 'TRP',
          'Y': 'TYR',
          'O': 'ORN',
          'U': 'PYL',
          '*': 'TER'}
    try:
        name1=dict(zip(name3.values(), name3.keys()))
        if mutation[:2] == 'p.':
            mutation = mutation.replace('p.', '')
        if mutation[1].isdigit():
            return mutation
        else:
            return name1[mutation[:3].upper()] + mutation[3:-3] + name1[mutation[-3:].upper()]
    except Exception as error:
        print(f'Error {error.__class__.__name__}: {error}')
    return None
# ======================================================================================================================
gnomad = []
patho = []
# ======================================================================================================================
grouping = {**{mutate(k): 'gnomAD' for k in gnomad},
            **{mutate(k): 'pathogenic' for k in patho}}

ref = pd.read_csv('gnomAD_v3.1_ENSG00000125877_2020_11_25_11_15_10.csv')
freq = dict(zip(ref['Consequence'].apply(lambda k: mutate(k) if 'c.' not in k else '').values, ref['Allele Frequency'].values))

data = data.assign(group=data.mutation.apply(lambda x: grouping[x] if x in grouping else 'other'),
                   freq=data.mutation.apply(lambda x: freq[x] if x in freq else 0)
                  )
data = data.loc[~data.complex_mutant_dG.isna()]
data = data.drop_duplicates(subset='mutation', keep='last')
data['complex_ddG'] = data['complex_ddG'].astype(float)

display(data.loc[data.group=='pathogenic'])
# ======================================================================================================================
fig = px.scatter(data,
                 x="freq",
                 y="complex_ddG",
                 color="group",
                 #size='FA_RMSD',
                 hover_data=["mutation"]
                )
fig.update_layout(title='Difference in folding Gibbs free energy',
                  xaxis=dict(title='gnomAD allele frequency'),
                  yaxis=dict(title='∆∆G'))
fig.show(renderer="notebook")