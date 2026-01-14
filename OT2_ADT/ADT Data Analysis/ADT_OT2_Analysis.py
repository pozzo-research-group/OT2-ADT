import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.ticker import MultipleLocator

def best_fit(X,Y):   
    '''
    Compute the slope and intercept of the line of best fit using least squares.
    
    Parameters
    ----------
    X : array-like
        Independent variable values
    Y : array-like
        Dependent variable values
    
    Returns
    -------
    m : float
        Slope of the best-fit line
    b : float
        Intercept of the best-fit line
    '''
   
    xbar = sum(X)/len(X)
    ybar = sum(Y)/len(Y)
    n = len(X)  # or len(Y)
    
    numer = sum([xi*yi for xi,yi in zip(X, Y)]) - n * xbar * ybar
    denum = sum([xi**2 for xi in X]) - n * xbar**2
    
    m = numer / denum
    b = ybar - m * xbar
    
    return m, b

def absorptivity(calibration_curve_csv):
    '''
    Calculate the slope and intercept of the UV-Vis calibration curve

    Parameters
    ----------
    calibration_curve_csv : str
        Path to the CSV file containing the calibration curve data.
        Must contain 'C_uM' (concentration) and 'Abs' (Absorbance).

    Returns
    -------
    m : float
        Slope of the calibration curve (absorptivity coefficient)
    b : float
        Intercept of the calibration curve
    '''

    df = pd.read_csv(calibration_curve_csv)
    X = df['C_uM']
    Y = df['Abs']
    m, b = best_fit(X,Y)
    print('line of best fit:\ny = {:.4f}x + {:.3f}'.format(m,b))
    return(m,b)


def longform_record(filename, absorptivity, dilution):
    '''
    Generate a detailed dataframe of concentration calculations from
    UV-Vis data.

    Parameters
    ----------
    filename : str
        Path to the CSV file containing absorbance data.
    absorptivity : float
        Slope of the calibration curve.
    dilution : float
        Dilution factor applied to the sample aliquots.

    Returns
    -------
    df : pandas.DataFrame
        DataFrame containing calculated sample concentrations and
        original (undiluted) concentrations.
    '''

    df = pd.read_csv(filename)
    df['C_sample'] = df['Abs']/absorptivity
    df['C_original'] = df['C_sample']*dilution
    df = df.round(decimals=3)
    return(df)

def shortform_C(filename, absorptivity, dilution):
    '''
    Generate a simplified dataframe of concentration values from
    blank-subtracted Abs data.

    Parameters
    ----------
    filename : str
        Path to CSV file containing absorbance data.
    absorptivity : float
        Slope of the calibration curve/ absorptivity coefficient.
    dilution : float
        Dilution factor applied to the sample aliquots.
 
    Returns
    -------
    data : pandas.DataFrame
        DataFrame containing time points and calculated concentrations 
        for each sample.
    '''
    
    df = pd.read_csv(filename)
    ds = df.drop(['Time'], axis=1)
    data = pd.DataFrame(data=df['Time'])
    
    for column_name, column_data in ds.items():
        data[column_name] = ds[column_name]*dilution/absorptivity
        
    data = data.round(decimals=4)
    return data


def data_calculations(shortform_filename, absorptivity, dilution_factor,
                      H_cells, index_range):
    '''
    Calculates membrane permeability (1e-5 cm/s) and diffusivity (1e-7 cm^2/s).
    
    Parameters
    ----------
    shortform_filename : str
        Path to CSV file containing the shortform concentration data.
    absorptivity : float
        Absorption coefficient (obtained from calibration curve).
    dilution_factor : float
        Dilution factor between original H-cells and UV-Vis samples.
    H_cells : dict
        Dictionary mapping H-cell IDs to membrane properties.
        Must contain:
            - 'membrane_L' : float, membrane thickness in µm
            - 'sample' : str, label for the membrane sample
    index_range : list or slice
        Range of indices to use for averaging diffusion coefficients
    
    Returns
    -------
    df : pandas.DataFrame
        Full DataFrame containing original and calculated values
    D_aves : pandas.DataFrame
        DataFrame of averaged diffusion coefficients and standard deviations,
        based on the values selected in the index_range.
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
        vol_cc = H_cells[key]['total_volume']
        radius = H_cells[key]['radius']
        area_cm = np.pi*(radius)**2  # (cm^2)
        L_cm = H_cells[key]['membrane_L']/10000

        sample_D = H_cells[key]['sample'] + '_D'
        sample_P = H_cells[key]['sample'] + '_P'

        columns.append(sample_D)
        df[sample_D] = (-1)*(vol_cc*L_cm/(
            4*area_cm*df['Time']*3600))*np.log(1-(2*df[C2]/C0))*(1e7)

        columns.append(sample_P)
        df[sample_P] = (-1)*(vol_cc/(
            4*area_cm*df['Time']*3600))*np.log(1-(2*df[C2]/C0))*(1e5)

    df = df.round(decimals = 3)

    # average diffusivities per individual sample
    df_aves = pd.DataFrame(columns=columns)
    df_aves.loc['mean'] = df.iloc[index_range].mean()
    df_aves.loc['std'] = df.iloc[index_range].std()
    df_aves = df_aves.round(decimals=3)
    return(df, df_aves)


def progress_plots(full_df, variable, sample_names, 
                   color_set, y_subticks, markers=6, font=11):
    """
    Create side-by-side progress plots for any number of H-cells.

    Parameters
    ----------
    full_df : pd.DataFrame
        Data containing Time and all sample columns.
    variable : str
        Either 'P' or 'D', used to select the permeability or diffusivity columns.
    sample_names : dict
        Dictionary like H_cells with keys ('H1', 'H2', ...) and sample info.
    color_set : list
        List of colors for plotting each sample.
    y_subticks : float
        Minor tick interval for y-axis.
    markers : int, optional
        Marker size.
    font : int, optional
        Font size for labels.
    """
    df = full_df.copy()
    
    # Build a dictionary mapping each H# → "S#_variable"
    D = {key: val['sample'] + '_' + variable for key, val in sample_names.items()}
    H_keys = list(sample_names.keys())
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4), layout='constrained')

    # Legend markers
    donor = mlines.Line2D([], [], color='black', ls='', marker='v',
                          markersize=8, label='Donor Chamber')
    receptor = mlines.Line2D([], [], color='black', ls='', marker='s',
                             markersize=8, label='Receptor Chamber')
    permeability = mlines.Line2D([], [], color='black', ls='', marker='o',
                                 markersize=8, label='Permeability')

    # -------------------------
    # Plot concentration data
    # -------------------------
    for i, key in enumerate(H_keys):
        color = color_set[i % len(color_set)]
        ax1.plot(df['Time'], df[f'{key}_C1'], ms=markers, marker='v',
                 color=color, linewidth=0)
        ax1.plot(df['Time'], df[f'{key}_C2'], ms=markers, marker='s',
                 color=color, linewidth=0)
    
    ax1ymax = df[[f'{key}_C1' for key in H_keys]].max().max() / 0.9
    ax1.set_xlabel('Time (hr)', fontsize=font)
    ax1.set_ylabel(r'Concentration ($\mu M$)', fontsize=font)
    ax1.set_ylim(0, ax1ymax)
    ax1.set_xlim(-2, df['Time'].iloc[-1] + 5)
    ax1.xaxis.set_minor_locator(MultipleLocator(2))
    ax1.legend(handles=[donor, receptor], fontsize=font-1, edgecolor='inherit')
    ax1.tick_params(labelsize=font)

    # -------------------------
    # Plot permeability/diffusivity data
    # -------------------------
    for i, key in enumerate(H_keys):
        color = color_set[i % len(color_set)]
        ax2.plot(df['Time'][1:], df[D[key]][1:], ms=markers, marker='o',
                 color=color, linewidth=0)

    # Axis labels and limits
    ax2.set_xlabel('Time (hr)', fontsize=font)
    if variable == 'P':
        yaxisname = 'Permeability (1e-5 cm/s)'
    elif variable == 'D':
        yaxisname = 'Diffusivity (1e-7 $cm^2/s$)'
    else:
        yaxisname = 'Unknown'

    ymax = df[[D[k] for k in H_keys]][1:].max().max() / 0.75
    ax2.set_ylabel(yaxisname, fontsize=font)
    ax2.set_ylim(0, ymax)
    ax2.set_xlim(0, df['Time'].iloc[-1] + 5)
    ax2.xaxis.set_minor_locator(MultipleLocator(2))
    ax2.yaxis.set_minor_locator(MultipleLocator(y_subticks))
    ax2.legend(handles=[permeability], fontsize=font-1, edgecolor='inherit')
    ax2.tick_params(labelsize=font)

    # -------------------------
    # Unified sample legend
    # -------------------------
    sample_patches = [
        mpatches.Patch(color=color_set[i % len(color_set)], label=D[key])
        for i, key in enumerate(H_keys)
    ]
    fig.legend(handles=sample_patches,
               loc='outside upper center', ncols=len(H_keys),
               fontsize=font, edgecolor='inherit')

    plt.show()


def data_averages_1_4(df, H_cells, variable,
                      H1_index = [],
                      H2_index = [],
                      H3_index = [],
                      H4_index = []):
    '''
    Calculates average diffusivity or permeabilty for each sample using the
    assigned index ranges for each sample, specific to H-cells labeled 
    H1-H4.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing concentrations and calculated 
        diffusivity/permeability values.
    sample_names : dict
        Mapping of H-cell IDs to sample labels.
    variable : str
        Either "D" for diffusivity or "P" for permeability.
    H1_index = list or slice
        Range of indices to use for calculating D/P coefficients and 
        standard deviations for the membrane sample in position H1.
    H2_index = list or slice
        Range of indices to use for calculating D/P coefficients and 
        standard deviations for the membrane sample in position H2.
    H3_index = list or slice
        Range of indices to use for calculating D/P coefficients and 
        standard deviations for the membrane sample in position H3.
    H4_index = list or slice
        Range of indices to use for calculating D/P coefficients and 
        standard deviations for the membrane sample in position H4.
    
    Returns
    -------
    None
        Prints each calculated permeability/diffisivity per sample.
    '''
    
    H1 = H_cells['H1']['sample'] + '_' + variable
    H2 = H_cells['H2']['sample'] + '_' + variable
    H3 = H_cells['H3']['sample'] + '_' + variable
    H4 = H_cells['H4']['sample'] + '_' + variable
    
    H1_ave = df[H1][H1_index].mean()
    H2_ave = df[H2][H2_index].mean()
    H3_ave = df[H3][H3_index].mean()
    H4_ave = df[H4][H4_index].mean()
    
    H1_std = df[H1][H1_index].std()
    H2_std = df[H2][H2_index].std()
    H3_std = df[H3][H3_index].std()
    H4_std = df[H4][H4_index].std()
    
    print(str(H1), ' = ', H1_ave.round(3), '+/-', H1_std.round(3))
    print(str(H2), ' = ', H2_ave.round(3), '+/-', H2_std.round(3))
    print(str(H3), ' = ', H3_ave.round(3), '+/-', H3_std.round(3))
    print(str(H4), ' = ', H4_ave.round(3), '+/-', H4_std.round(3))
    return

def data_averages_5_8(df, H_cells, variable,
                      H5_index=[],
                      H6_index=[],
                      H7_index=[],
                      H8_index=[]):
    '''
    Calculates average diffusivity or permeabilty for each sample using the
    assigned index ranges for each sample, specific to H-cells labeled 
    H5-H8.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing concentrations and calculated 
        diffusivity/permeability values.
    sample_names : dict
        Mapping of H-cell IDs to sample labels.
    variable : str
        Either "D" for diffusivity or "P" for permeability.
    H1_index = list or slice
        Range of indices to use for calculating D/P coefficients and 
        standard deviations for the membrane sample in position H1.
    H2_index = list or slice
        Range of indices to use for calculating D/P coefficients and 
        standard deviations for the membrane sample in position H2.
    H3_index = list or slice
        Range of indices to use for calculating D/P coefficients and 
        standard deviations for the membrane sample in position H3.
    H4_index = list or slice
        Range of indices to use for calculating D/P coefficients and 
        standard deviations for the membrane sample in position H4.
    
    Returns
    -------
    None
        Prints each calculated permeability/diffisivity per sample.
    '''

    H5 = H_cells['H5']['sample'] + '_' + variable
    H6 = H_cells['H6']['sample'] + '_' + variable
    H7 = H_cells['H7']['sample'] + '_' + variable
    H8 = H_cells['H8']['sample'] + '_' + variable

    H5_ave = df[H5][H5_index].mean()
    H6_ave = df[H6][H6_index].mean()
    H7_ave = df[H7][H7_index].mean()
    H8_ave = df[H8][H8_index].mean()

    H5_std = df[H5][H5_index].std()
    H6_std = df[H6][H6_index].std()
    H7_std = df[H7][H7_index].std()
    H8_std = df[H8][H8_index].std()

    print(str(H5), ' = ', H5_ave.round(3), '+/-', H5_std.round(3))
    print(str(H6), ' = ', H6_ave.round(3), '+/-', H6_std.round(3))
    print(str(H7), ' = ', H7_ave.round(3), '+/-', H7_std.round(3))
    print(str(H8), ' = ', H8_ave.round(3), '+/-', H8_std.round(3))
    return

def data_averages_8(df, H_cells, variable,
                    H1_index = [],
                    H2_index = [],
                    H3_index = [],
                    H4_index = [],
                    H5_index = [],
                    H6_index = [],
                    H7_index = [],
                    H8_index = []):
    '''
    Calculates average diffusivity (1e-7 cm^2/s) or permeabilty (1e-5 cm/s) for 
    each sample using the assigned index ranges for each sample.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing concentrations and calculated 
        diffusivity/permeability values.
    sample_names : dict
        Mapping of H-cell IDs to sample labels.
    variable : str
        Either "D" for diffusivity or "P" for permeability.
    H1_index = list or slice
        Range of indices to use for calculating D/P coefficients and 
        standard deviations for the membrane sample in position H1.
    H2_index = list or slice
        Range of indices to use for calculating D/P coefficients and 
        standard deviations for the membrane sample in position H2.
    H3_index = list or slice
        Range of indices to use for calculating D/P coefficients and 
        standard deviations for the membrane sample in position H3.
    H4_index = list or slice
        Range of indices to use for calculating D/P coefficients and 
        standard deviations for the membrane sample in position H4.
    H5_index = list or slice
        Range of indices to use for calculating D/P coefficients and 
        standard deviations for the membrane sample in position H5.
    H6_index = list or slice
        Range of indices to use for calculating D/P coefficients and 
        standard deviations for the membrane sample in position H6.
    H7_index = list or slice
        Range of indices to use for calculating D/P coefficients and 
        standard deviations for the membrane sample in position H7.
    H8_index = list or slice
        Range of indices to use for calculating D/P coefficients and 
        standard deviations for the membrane sample in position H8.
    
    Returns
    -------
    None
        Prints each calculated permeability/diffisivity per sample.
    '''
    
    H1 = H_cells['H1']['sample'] + '_' + variable
    H2 = H_cells['H2']['sample'] + '_' + variable
    H3 = H_cells['H3']['sample'] + '_' + variable
    H4 = H_cells['H4']['sample'] + '_' + variable
    H5 = H_cells['H5']['sample'] + '_' + variable
    H6 = H_cells['H6']['sample'] + '_' + variable
    H7 = H_cells['H7']['sample'] + '_' + variable
    H8 = H_cells['H8']['sample'] + '_' + variable
    
    H1_ave = df[H1][H1_index].mean()
    H2_ave = df[H2][H2_index].mean()
    H3_ave = df[H3][H3_index].mean()
    H4_ave = df[H4][H4_index].mean()
    H5_ave = df[H5][H5_index].mean()
    H6_ave = df[H6][H6_index].mean()
    H7_ave = df[H7][H7_index].mean()
    H8_ave = df[H8][H8_index].mean()
    
    H1_std = df[H1][H1_index].std()
    H2_std = df[H2][H2_index].std()
    H3_std = df[H3][H3_index].std()
    H4_std = df[H4][H4_index].std()
    H5_std = df[H5][H5_index].std()
    H6_std = df[H6][H6_index].std()
    H7_std = df[H7][H7_index].std()
    H8_std = df[H8][H8_index].std()
    
    print(str(H1), ' = ', H1_ave.round(3), '+/-', H1_std.round(3))
    print(str(H2), ' = ', H2_ave.round(3), '+/-', H2_std.round(3))
    print(str(H3), ' = ', H3_ave.round(3), '+/-', H3_std.round(3))
    print(str(H4), ' = ', H4_ave.round(3), '+/-', H4_std.round(3))
    print(str(H5), ' = ', H5_ave.round(3), '+/-', H5_std.round(3))
    print(str(H6), ' = ', H6_ave.round(3), '+/-', H6_std.round(3))
    print(str(H7), ' = ', H7_ave.round(3), '+/-', H7_std.round(3))
    print(str(H8), ' = ', H8_ave.round(3), '+/-', H8_std.round(3))
    return