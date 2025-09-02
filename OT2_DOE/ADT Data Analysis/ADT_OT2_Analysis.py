import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.ticker import AutoMinorLocator, MultipleLocator


def best_fit(X,Y):
    
    xbar = sum(X)/len(X)
    ybar = sum(Y)/len(Y)
    n = len(X) # or len(Y)
    
    numer = sum([xi*yi for xi,yi in zip(X, Y)]) - n * xbar * ybar
    denum = sum([xi**2 for xi in X]) - n * xbar**2
    
    m = numer / denum
    b = ybar - m * xbar
    
    return m, b


def absorptivity(calibration_curve_csv):
    df = pd.read_csv(calibration_curve_csv)
    X = df['C_uM']
    Y = df['Abs']
    m, b = best_fit(X,Y)
    print('line of best fit:\ny = {:.4f}x + {:.3f}'.format(m,b))
    return(m, b)


def longform_record(filename, absorptivity, dilution):
    df = pd.read_csv(filename)
    df['C_sample'] = df['Abs']/absorptivity
    df['C_original'] = df['C_sample']*dilution
    df = df.round(decimals=3)
    return(df)


def shortform_C(filename, absorptivity, dilution):
    df = pd.read_csv(filename)
    ds = df.drop(['Time'], axis=1)
    data = pd.DataFrame(data=df['Time'])
    
    for column_name, column_data in ds.items():
        data[column_name] = ds[column_name]*dilution/absorptivity
        
    data = data.round(decimals=4)
    return data


def data_calculations(shortform_filename, absorptivity, dilution_factor, 
                     total_vol, area, H_cells, index_range):
    '''
    Performs all calculations for OT2 experimental data
    
    Parameters
    ----------
    shortform_filename: file name of the .csv for processing. Variable defined earlier.
    absorptivity = absorption coefficient (obtained from calibration curve)
    dilution_factor = dilution factor between original H-cells and UV-Vis samples
    
    total_vol: total solution volume in each h-cell in cubic microns
    area: exposed membrane area, in squared microns
    L_um: membrane thicknesses. Default is 120 microns (defined earlier), but this can 
        be specified for each membrane.
    H_cells: dictionary that correlates each H-cell position to its membrane sample
    index_range: the selection of data points to take the average diffusion coefficient from
    
    Returns
    ----------
    df: df with all original and calculated values
    D_aves: df with final D values and errors for each membrane sample
    '''
    shortform = shortform_C(shortform_filename, absorptivity, dilution_factor)
    df = pd.DataFrame(data=shortform)
    
    
    # diffusivities
    #H_cells is the dictionary
    columns = []
    
    for key, value in H_cells.items():
        for column_name, column_data in df.items():
            if column_name[:2] == key:
                if column_name[-2:] == 'C1':
                    H_cells[key]['cell 1'] = column_name
                elif column_name[-2:] == 'C2':
                    H_cells[key]['cell 2'] = column_name
                else:
                    pass
            else:
                pass
        C2 = H_cells[key]['cell 2']
        C1 = H_cells[key]['cell 1']
        C0 = df[C1][0]+df[C2][0]
        L_um = H_cells[key]['membrane_L']
        sample = H_cells[key]['sample']
        columns.append(sample)
        df[sample] = (-1)*(total_vol*L_um/(4*area*df['Time']*3600))*np.log(1-(2*df[C2]/C0))
        
    df = df.round(decimals = 2)
    
    # average diffusivities per individual sample
    D_aves = pd.DataFrame(columns=columns)
    D_aves.loc['mean'] = df.iloc[index_range].mean()
    D_aves.loc['std'] = df.iloc[index_range].std()
    D_aves = D_aves.round(decimals=1)
    return(df, D_aves)


def progress_plots_4(full_df, sample_names, max_time, max_D, color_set, font=18):
    '''
    Makes scatter plots of Concentration vs Time (left) and Diffusivity (right)
    
    Parameters
    ----------
    full_df: the df containing concentrations and permeability calculations
    plot_labels: [] containing labels for each membrane sample
    max_time: upper limit for x range
    max_D: upper limit for diffusivity range
    color_set: colors to be used for each sample.
    
    '''
    
    df = full_df
    
    D1 = sample_names['H1']['sample']
    D2 = sample_names['H2']['sample']
    D3 = sample_names['H3']['sample']
    D4 = sample_names['H4']['sample']
    
    fig, (ax1,ax2) = plt.subplots(1,2, figsize=(12,6), layout='constrained')
    
    donor = mlines.Line2D([], [], color='black', ls='', marker='v',
                          markersize=8, label='Donor Chamber (Cell 1)')
    receptor = mlines.Line2D([], [], color='black', ls='', marker='s',
                             markersize=8, label='Receptor Chamber (Cell 2)')
    diffusivity = mlines.Line2D([], [], color='black', ls='', marker='o',
                                markersize=8, label='Diffusion Coefficient')
    
    sample_1 = mpatches.Patch(color=color_set[0], label=D1)
    sample_2 = mpatches.Patch(color=color_set[1], label=D2)
    sample_3 = mpatches.Patch(color=color_set[2], label=D3)
    sample_4 = mpatches.Patch(color=color_set[3], label=D4)

    ax1.plot(df['Time'], df['H1_C1'], ms=8, marker='v', color=color_set[0], linewidth=0)
    ax1.plot(df['Time'], df['H2_C1'], ms=8, marker='v', color=color_set[1], linewidth=0)
    ax1.plot(df['Time'], df['H3_C1'], ms=8, marker='v', color=color_set[2], linewidth=0)
    ax1.plot(df['Time'], df['H4_C1'], ms=8, marker='v', color=color_set[3], linewidth=0)
    ax1.plot(df['Time'], df['H1_C2'], ms=8, marker='s', color=color_set[0], linewidth=0)
    ax1.plot(df['Time'], df['H2_C2'], ms=8, marker='s', color=color_set[1], linewidth=0)
    ax1.plot(df['Time'], df['H3_C2'], ms=8, marker='s', color=color_set[2], linewidth=0)
    ax1.plot(df['Time'], df['H4_C2'], ms=8, marker='s', color=color_set[3], linewidth=0)
    ax1.set_xlabel('Time (hr)', fontsize=font)
    ax1.set_ylabel(r'Concentration ($\mu M$)', fontsize=font)
    ax1.set_ylim(0,df['H1_C1'][0]+25)
    ax1.set_xlim(-2,max_time+5)
    ax1.xaxis.set_minor_locator(MultipleLocator(2))
    ax1.legend(handles=[donor, receptor], fontsize=font-4, edgecolor='inherit')
    ax1.tick_params(labelsize=font-2)

    ax2.plot(df['Time'][1:], df[D1][1:], ms=8, marker='o', color=color_set[0], linewidth=0)
    ax2.plot(df['Time'][1:], df[D2][1:], ms=8, marker='o', color=color_set[1], linewidth=0)
    ax2.plot(df['Time'][1:], df[D3][1:], ms=8, marker='o', color=color_set[2], linewidth=0)
    ax2.plot(df['Time'][1:], df[D4][1:], ms=8, marker='o', color=color_set[3], linewidth=0)
    ax2.set_xlabel('Time (hr)', fontsize=font)
    ax2.set_ylabel('Diffusivity ($\mu m^2/s$)', fontsize=font)
    ax2.set_ylim(0,max_D)
    ax2.set_xlim(0,max_time+5)
    ax2.xaxis.set_minor_locator(MultipleLocator(2))
    ax2.legend(handles=[diffusivity], fontsize=font-4, edgecolor='inherit')
    ax2.tick_params(labelsize=font-2)
    
    fig.legend(handles=[sample_1, sample_2, sample_3, sample_4], loc='outside upper center',
               ncols=4, fontsize=font-2, edgecolor='inherit')

    return(plt.show())



def progress_plots_8(full_df, plot_labels, sample_names, max_time, max_D, font=18):
    '''
    Makes scatter plots of Concentration vs Time (left) and Diffusivity (right)
    
    Parameters
    ----------
    full_df: the df containing concentrations and permeability calculations
    plot_labels: [] containing labels for each membrane sample
    max_time: upper limit for x range
    max_D: upper limit for diffusivity range
    color_set: colors to be used for each sample.
    
    '''
    
    df = full_df
    
    D1 = sample_names['H1']['sample']
    D2 = sample_names['H2']['sample']
    D3 = sample_names['H3']['sample']
    D4 = sample_names['H4']['sample']
    D5 = sample_names['H5']['sample']
    D6 = sample_names['H6']['sample']
    D7 = sample_names['H7']['sample']
    D8 = sample_names['H8']['sample']
    
    color_set = ['#332288','#117733','#44AA99','##88CCEE',
                 '#C5B044','#CC6677','#AA4499','#882255']
    
    fig, (ax1,ax2) = plt.subplots(1,2, figsize=(12,6), layout='constrained')
    
    donor = mlines.Line2D([], [], color='black', ls='', marker='v',
                          markersize=8, label='Donor Chamber (Cell 1)')
    receptor = mlines.Line2D([], [], color='black', ls='', marker='s',
                             markersize=8, label='Receptor Chamber (Cell 2)')
    diffusivity = mlines.Line2D([], [], color='black', ls='', marker='o',
                                markersize=8, label='Diffusion Coefficient')
    
    sample_1 = mpatches.Patch(color=color_set[0], label=D1)
    sample_2 = mpatches.Patch(color=color_set[1], label=D2)
    sample_3 = mpatches.Patch(color=color_set[2], label=D3)
    sample_4 = mpatches.Patch(color=color_set[3], label=D4)
    sample_5 = mpatches.Patch(color=color_set[4], label=D5)
    sample_6 = mpatches.Patch(color=color_set[5], label=D6)
    sample_7 = mpatches.Patch(color=color_set[6], label=D7)
    sample_8 = mpatches.Patch(color=color_set[7], label=D8)

    ax1.plot(df['Time'], df['H1_C1'], ms=8, marker='v', color=color_set[0], linewidth=0)
    ax1.plot(df['Time'], df['H2_C1'], ms=8, marker='v', color=color_set[1], linewidth=0)
    ax1.plot(df['Time'], df['H3_C1'], ms=8, marker='v', color=color_set[2], linewidth=0)
    ax1.plot(df['Time'], df['H4_C1'], ms=8, marker='v', color=color_set[3], linewidth=0)
    ax1.plot(df['Time'], df['H5_C1'], ms=8, marker='v', color=color_set[4], linewidth=0)
    ax1.plot(df['Time'], df['H6_C1'], ms=8, marker='v', color=color_set[5], linewidth=0)
    ax1.plot(df['Time'], df['H7_C1'], ms=8, marker='v', color=color_set[6], linewidth=0)
    ax1.plot(df['Time'], df['H8_C1'], ms=8, marker='v', color=color_set[7], linewidth=0)
    ax1.plot(df['Time'], df['H1_C2'], ms=8, marker='s', color=color_set[0], linewidth=0)
    ax1.plot(df['Time'], df['H2_C2'], ms=8, marker='s', color=color_set[1], linewidth=0)
    ax1.plot(df['Time'], df['H3_C2'], ms=8, marker='s', color=color_set[2], linewidth=0)
    ax1.plot(df['Time'], df['H4_C2'], ms=8, marker='s', color=color_set[3], linewidth=0)
    ax1.plot(df['Time'], df['H5_C2'], ms=8, marker='s', color=color_set[4], linewidth=0)
    ax1.plot(df['Time'], df['H6_C2'], ms=8, marker='s', color=color_set[5], linewidth=0)
    ax1.plot(df['Time'], df['H7_C2'], ms=8, marker='s', color=color_set[6], linewidth=0)
    ax1.plot(df['Time'], df['H8_C2'], ms=8, marker='s', color=color_set[7], linewidth=0)
    ax1.set_xlabel('Time (hr)', fontsize=font)
    ax1.set_ylabel(r'Concentration ($\mu M$)', fontsize=font)
    ax1.set_ylim(0,df['H1_C1'][0]+25)
    ax1.set_xlim(-2,max_time+5)
    ax1.xaxis.set_minor_locator(MultipleLocator(2))
    ax1.legend(handles=[donor, receptor], fontsize=font-4, edgecolor='inherit')
    ax1.tick_params(labelsize=font-2)

    ax2.plot(df['Time'][1:], df[D1][1:], ms=8, marker='o', color=color_set[0], linewidth=0)
    ax2.plot(df['Time'][1:], df[D2][1:], ms=8, marker='o', color=color_set[1], linewidth=0)
    ax2.plot(df['Time'][1:], df[D3][1:], ms=8, marker='o', color=color_set[2], linewidth=0)
    ax2.plot(df['Time'][1:], df[D4][1:], ms=8, marker='o', color=color_set[3], linewidth=0)
    ax2.plot(df['Time'][1:], df[D5][1:], ms=8, marker='o', color=color_set[4], linewidth=0)
    ax2.plot(df['Time'][1:], df[D6][1:], ms=8, marker='o', color=color_set[5], linewidth=0)
    ax2.plot(df['Time'][1:], df[D7][1:], ms=8, marker='o', color=color_set[6], linewidth=0)
    ax2.plot(df['Time'][1:], df[D8][1:], ms=8, marker='o', color=color_set[7], linewidth=0)
    ax2.set_xlabel('Time (hr)', fontsize=font)
    ax2.set_ylabel('Diffusivity ($\mu m^2/s$)', fontsize=font)
    ax2.set_ylim(0,max_D)
    ax2.set_xlim(0,max_time+5)
    ax2.xaxis.set_minor_locator(MultipleLocator(2))
    ax2.legend(handles=[diffusivity], fontsize=font-4, edgecolor='inherit')
    ax2.tick_params(labelsize=font-2)
    
    fig.legend(handles=[sample_1, sample_2, sample_3, sample_4,sample_5,sample_6,sample_7,sample_8],
               loc='outside upper center', ncols=4, fontsize=font-2, edgecolor='inherit')

    return(plt.show())


