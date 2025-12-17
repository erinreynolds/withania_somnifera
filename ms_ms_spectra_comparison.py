#adapted from code written by Trevor Moss in 2024

##### IMPORTS #####
import glob, os
import pandas as pd
import matplotlib.pyplot as plt

##### PARAMETERS #####
input_folder = '' #folder containing spectra in .csv form
output_folder = ''
xlsx_files = sorted(glob.glob(os.path.join(input_folder, '*.xlsx')))
colors = ['blue', 'red']
comparisonStyle = 'Reflected'
exteriorLegend = True
numLabeledPeaks = 8
verticalLineWidth = 0.5
xlabelSize = 7
annotationSize = 5

##### STYLING #####
plt.figure(figsize=(4, 2))
plt.axhline(y = 0, linewidth = verticalLineWidth, color = 'darkgray') 
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.xlabel('m/z', fontdict={'weight': 'normal', 'size': xlabelSize})

##### PLOTTING #####
for n, filepath in enumerate(xlsx_files):
    df = pd.read_excel(filepath, skiprows=1)
    df.loc[(df['m/z'] < 202.0779) & (df['m/z'] > 202.0775), ['Intensity']]=0 #specific to Weng lab orbitrap (internal calibration peak)
    df['Intensity'] = df['Intensity']/df.max().iloc[1]*100

    if n == 1 and comparisonStyle == 'Reflected':
        df['Intensity'] = df['Intensity'] * -1

    for x, y in zip(df['m/z'], df['Intensity']):
        g = plt.vlines(x, 0, y, color = colors[n], linestyle='-', linewidth=verticalLineWidth, alpha=1, label=filepath.split('/')[-1].split('-')[0])
            
    if n == 1 and comparisonStyle == 'Reflected':
        dfMaxes = df.sort_values(by='Intensity', ascending=True).head(numLabeledPeaks)
    else:
        dfMaxes = df.sort_values(by='Intensity', ascending=False).head(numLabeledPeaks)
    
    texts = [plt.text(x, y, round(x, 4), ha='center', fontdict={'weight': 'normal', 'size': annotationSize}) for x, y in zip(dfMaxes['m/z'], dfMaxes['Intensity'])]

    # Change xtick and ytick label size
    plt.xticks(fontsize=7)  # Change x-axis tick label size
    plt.yticks(fontsize=7)  # Change y-axis tick label size
    
    ax = plt.gca() #get current axis
    ax.spines['top'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['right'].set_linewidth(0.5)

# Save the plot
output_file_2 = os.path.join(output_folder, 'plotname.svg')
plt.savefig(output_file_2)
plt.show()