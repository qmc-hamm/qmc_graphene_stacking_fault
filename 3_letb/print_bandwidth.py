import pandas as pd

d = pd.read_csv('processed/bandwidth.csv')[['twist_angle', 'potential', 'bandwidth']]
d['Potential'] = d['potential'].map({'qmc': 'KC-QMC', 'ouyang': 'KC-Ouyang', 'dft_d2': 'KC-DFT-D2', 'dft_d3': 'KC-DFT-D3'})
t = d.pivot(index='Potential', columns='twist_angle', values='bandwidth').reset_index(0)
t['sort'] = t['Potential'].map({'KC-QMC': 0, 'KC-Ouyang': 1, 'KC-DFT-D2': 2, 'KC-DFT-D3': 3})
t = t.sort_values(by='sort').drop('sort', axis=1).reset_index()
print(t.columns)
print(t.round(3).to_latex(index=False))
