import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
from scipy.stats import mannwhitneyu
from matplotlib.ticker import ScalarFormatter

# User input:

# Input and output folders
input_file = 'withaferinA_peakarea.csv'
output_folder = 'wfa_graph'
    
# Import data and create a DataFrame
df = pd.read_csv(input_file)

#Calculate Mann U Whitney significance
ctrl_group = df.loc[df['condition'] == 'ctrl']
exp_group = df.loc[df['condition'] == 'silenced']

key_list = df['gene'].unique().tolist()

#for each key in key_list, calculate Mann Whitney U significance
for gene in key_list:
    ctrl = ctrl_group.loc[ctrl_group['gene'] == gene, 'peak area']
    exp = exp_group.loc[exp_group['gene'] == gene, 'peak area']

    # Perform Mann-Whitney U test
    stat, p = mannwhitneyu(ctrl, exp, alternative='two-sided')
    
    print (gene)
    print(f"U statistic = {stat}")
    print(f"p-value = {p}")

# Create color palette
palette = ["#029e73", "#d55e00"]
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['ytick.labelsize'] = 7

# Custom errorbar function: return (lower, upper)
def iqr_error(x):
    q1 = np.percentile(x, 25)
    q3 = np.percentile(x, 75)
    return q1, q3   # seaborn will use this as error bar bounds

# Create the bar plot
plt.figure(figsize=(3, 2), dpi =300)
sns.barplot(
    x='gene', 
    y='peak area', 
    hue='condition', 
    data=df,
    palette=palette,
    alpha = 0.4, 
    errwidth=0.7,
    capsize = 0.1,
    errcolor = 'black',
    dodge=True,
    linewidth = 0.5,
    estimator=np.median, #make the barplots use median
    errorbar=iqr_error #make the error bars interquartile range
    )

# Modify the legend labels
# Get current legend handles and labels
handles, labels = plt.gca().get_legend_handles_labels()

# Add individual data points
sns.stripplot(
    x='gene', 
    y='peak area', 
    hue='condition', 
    data=df, 
    color='black', 
    size=1, 
    alpha = 1, 
    linewidth=0.5, 
    dodge=True, 
    jitter=True
)

# Modify the legend labels, keeping the handles (colors)
new_labels = ['control', 'silenced']
plt.legend(handles=handles, labels=new_labels, title='', frameon=False, loc='upper right', fontsize=6)

ax = plt.gca()
ax.set_ylim(bottom=0)

class OneDecimalScalarFormatter(ScalarFormatter):
    def _set_format(self):
        self.format = "%.1f"  # force 1 decimal

ax = plt.gca()
formatter = OneDecimalScalarFormatter(useOffset=True, useMathText=True)
formatter.set_scientific(True)
ax.yaxis.set_major_formatter(formatter)

#format axes
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

# Get the offset text object
offset = ax.yaxis.get_offset_text()

# Move it to the top right corner of the axes
offset.set_ha('right')       # horizontal alignment
offset.set_va('bottom')      # vertical alignment
offset.set_x(0)            # place at far right of axis (1.0 = right edge in axes coords)
offset.set_y(1.3)            # place at top (1.0 = top in axes coords)

plt.gca().set_xticklabels([])
plt.gca().set_xticks([])
plt.yticks(fontsize=7)
plt.xlabel('')
plt.ylabel('Withaferin A peak area', fontsize=7)
sns.despine()

# Save the plot
output_file_2 = os.path.join(output_folder, 'WFA_medianbarplot.svg')
plt.savefig(output_file_2)
plt.show()
    

        

