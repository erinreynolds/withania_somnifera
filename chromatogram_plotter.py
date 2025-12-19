#Trevor Moss
#Weng Lab, Fall 2024
#Created: 20240916
#Last updated: 20251219

##### IMPORTS #####
import glob, os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

##### PARAMETERS #####
#input folder contains .csv files of chromatograms, named in the format "0-1, 1-1, 2-1..." where the first digit is a code for the yeast strain and second number indicates the replicate
input_folder = 'input_folder_file_path'


### The following customization variables have specific variable type requirements shown in each comment. ###
### If not using preference field, set variable to False. If using field, follow input format in example. ###

# Plot autosave #
output_folder = 'output_file_path'  # Specified where to save plots. 
plotFiletype = 'svg'  # Tested: 'png', 'svg', 'pdf', 'ps', 'tiff', 'webp'.

# Custom trace ordering #
# order is from bottom to top in the graph
custom_sort_order = ('1', '5', '4', '3', '2', '0')  # Specific order of files, ordered by filename before dash. Ex: Filenames: {0-etc, 1-etc, 2-etc} Input: (2, 0, 1) Sorts to {2-etc, 0-etc, 1-etc}. Type: tuple

# Axis parameters #
yAxisTitle = 397.3465  # Sets y-axis label to be: EIC (m/z = {input here}). Ex: 443.3156 results in EIC (m/z = 443.3156). Type: int/float
custom_x_limits = (6, 12)  # Zooms plot to specified x values. Tuple of two x coordinates. Ex: (5, 8). Type: tuple
rtDemarcation = [10.6]  # Vertical dashed lines to highlight regions along the x-axis. Ex: [6, 6.8841] or [4]. Type: list

# Misc. #
traceExtraOffsetFactor = 2  # Extra space between traces: highest peak in the entire plot, divided by the factor

##### STYLE SETTINGS & X-AXIS TRIM #####
sns.set_palette(sns.color_palette('colorblind'))

def graphingSetup():
    # Set font
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial']

    # EIC label on y-axis, or remove
    if type(yAxisTitle) == float or type(yAxisTitle) == int:
        g.set_ylabel(f'EIC (m/z = {yAxisTitle})', fontdict={'weight': 'normal'}, fontsize=7)
        g.spines['left'].set_linewidth(0.5)
    else:
        g.set_ylabel('')
        g.spines['left'].set_visible(False)
    g.yaxis.set_ticks([])

    # Customize x-axis: bold axis, ticks, axis label, and tick labels
    g.spines['bottom'].set_linewidth(0.5)
    g.xaxis.set_tick_params(width=0.5)
    plt.draw()
    g.set_xlabel('Retention time (min)', fontdict={'weight': 'normal'}, fontsize = 7)
    plt.xticks(fontsize=6)
   
    # Remove plot spines
    g.spines['top'].set_visible(False)
    g.spines['right'].set_visible(False)

def trim_dataframe_by_time(df, x_limits):
    # Filter the dataframe to keep rows where x is between start_time and end_time
    trimmed_df = df[(df['Time(min)'] >= x_limits[0]) & (df['Time(min)'] <= x_limits[1])]
    return trimmed_df

##### FILE IMPORT & OFFSET #####
if type(custom_sort_order) == tuple:  # CSV import, sorted according to custom order, or alphanumerically.
    order_lookup = {num: i for i, num in enumerate(custom_sort_order)}
    csv_files = sorted(glob.glob(os.path.join(input_folder, '*.csv')), key=lambda f: order_lookup[os.path.basename(f).split('-')[0]])
else:
    csv_files = sorted(glob.glob(os.path.join(input_folder, '*.csv')), reverse=True)

dfmaxes = []
for filepath in csv_files:  # Determines extra_offset to separate the traces. Without extra_offset, the top of one trace can touch the bottom of the trace above.
    df = pd.read_csv(filepath,skiprows=3)
    if type(custom_x_limits) == tuple:  # Accounting for custom x range. If not accounted for, offscreen peaks will determine the spacing.
        df = trim_dataframe_by_time(df, custom_x_limits)
    dfmaxes.append(df.max().iloc[1])
extra_offset = max(dfmaxes) / traceExtraOffsetFactor

##### PLOTTING #####
offset = 0
plt.figure(figsize=(3, 1.25)) 
for filepath in csv_files:  # Plots each trace, translated above the previous trace's highest onscreen peak.
    base_name = os.path.basename(filepath)
    strain_num = base_name.split('-')[0]
    df = pd.read_csv(filepath, skiprows=3)
    if type(custom_x_limits) == tuple:  # Accounting for custom x range. If not accounted for, offscreen peaks will determine the spacing.
        df = trim_dataframe_by_time(df, custom_x_limits)
    df['Relative Abundance'] = df['Relative Abundance'] + offset + extra_offset  # Adds highest peak from previous trace to all y values to separate the lines. Also adds extra offset to avoid line contact.
    offset = df.max().iloc[1]  # Sets the new offset for the next line/trace
    g = sns.lineplot(x='Time(min)', y='Relative Abundance', data=df, label=strain_num, linewidth = 0.75)  # Writes line to the plot

##### LEGEND #####
# remove legend
g.legend().remove()

##### X-AXIS ZOOM #####
if type(custom_x_limits) == tuple:
    g.set_xlim(custom_x_limits[0], custom_x_limits[1])
else:
    g.set_xlim(g.get_xticks()[1], g.get_xticks()[-2])  # The defalt xticks include buffer on either side, this sets the region to match the MS data window.

##### DASHED LINES TO HIGHLIGHT RTs #####
if type(rtDemarcation) == list:
    for rt in rtDemarcation:
        plt.axvline(rt, ymin=-0.00, color='grey', linewidth=0.5, linestyle='dashed', clip_on=False)  # Adds a vertical line at each x value in rtDemarcation, extending below the y-axis for clarity (ymin=-0.03).

##### CALL STYLE SETTINGS FUNCTION, SAVE AND SHOW PLOT #####
graphingSetup()  # Calls style function

if output_folder != False:  # Only saves if location is specified.
    os.makedirs(output_folder, exist_ok=True)
    output_file = os.path.join(output_folder, f'{input_folder.split("/")[-1]}.{plotFiletype}')  # Saves graph as specified filetype in folder specified by output_folder
    plt.savefig(output_file, format=plotFiletype)

plt.show()
