import time
import string
import numpy as np
import pandas as pd

# -------------------------------------------------------------------
# Core Utility Functions
# -------------------------------------------------------------------

def get_Hcell_labware(loaded_dict):
    """
    Extract and organize H-cell wells and associated sample plate
    wells from the loaded labware.
    
    This function separates all destination wells into either "H-cell wells"
    or standard 96-well plate "sample wells". All H-cell wells are assorted 
    column-wise to keep H-cell donor and receptor chambers paired together. 
    The sample wells are sorted by group (row A, row B, etc) so that each 
    H-cell chamber has a designated row in the sample plate.

    Parameters
    ----------
    loaded_dict : dict
        Dictionary containing labware and well mappings loaded by the OT2.
        Must contain a key 'Destination Wells', mapping to the list of 
        H-cell and sample plate hardware.
        
    Returns
    -------
    h_cell_wells : list
        Ordered list of H-cell wells arranged column-wise.
        E.g., ['A1', 'B1', 'A2', 'B2']
    wellplates: list of list
        Nested list of sample wells grouped row-wise corresponding to each
        H-cell chamber. Each sublist corresponds to the row dedicated to one
        H-cell well.
        
    Notes
    -----
    - H-cell chambers are arranged column-wise, pairing donor and receptor chambers.
    - Each H-cell chamber has a dedicated row on the sample plate
      (e.g., H1 donor → row A, H1 receptor → row B). 

    """
    
    dest_lab = loaded_dict['Destination Wells']
    
    h_cell_wells = []
    wellplates = []
    for row in dest_lab:
        if 'hcell' in str(row):
            h_cell_wells.append(row)
        else:
            wellplates.append(row)
    
    # Organize the wells in the H-cell list to be column-wise
    h_cell_wells= rearrange_hcell_wells(h_cell_wells)
    
    # Group sample wells corresponding to each H-cell
    # Each H-cell chamber has a dedicated row (i.e. H1 donor goes to row A)
    number_of_hcells = len(h_cell_wells)//2
    
    wellplates = get_Hcell_sample_wells(wellplates, number_of_hcells)
    
    return h_cell_wells, wellplates


def rearrange_hcell_wells(h_cell_wells):
    """
    Arrange H-cell wells in a consistent column-wise order.

    Given a flat list of H-cell wells, this function organizes them so that
    each donor/receptor pair is grouped together. Within each plate of four H-cell
    chambers, wells are sorted according to the expected pattern:
    A1 → B1 → A2 → B2.

    Parameters
    ----------
    h_cell_wells : list
        Unordered list of H-cell wells (e.g., ['A1', 'A2', 'B1', 'B2']).

    Returns
    -------
    new_h_cell_wells : list
        Ordered list of H-cell wells sorted column-wise.
        Example:
            ['A1', 'B1', 'A2', 'B2']

    Notes
    -----
    - Each H-cell plate contains four chambers: two donor and two receptor.
    - Ensures donor and receptor chambers for the same H-cell are adjacent.
    """
    
    # H-cell units are arranged column-wise
    # This pairs each donor and receptor chamber
    order = ['A1', 'B1', 'A2', 'B2']
    
    h_cell_plates = len(h_cell_wells)//4
    
    new_h_cell_wells = []
    for i in range(h_cell_plates):
        for well in order:
            for chamber in h_cell_wells[i * 4: (i + 1) * 4]:
                if well in str(chamber):
                    new_h_cell_wells.append(chamber)
    
    return new_h_cell_wells


def get_Hcell_sample_wells(wellplates, number_of_hcells):
    """
    Group sample plate wells row-wise for each H-cell chamber.

    This function assigns the correct rows from 96-well plates to each H-cell.
    Each H-cell chamber (donor or receptor) gets a dedicated row of sample wells.
    If multiple 96-well plates are present, wells are grouped by plate and row.

    Parameters
    ----------
    wellplates : list
        Flat list of wells from all sample plates (excluding H-cell wells).
    number_of_hcells : int
        Number of H-cell chambers. Used to determine well-to-chamber mapping.

    Returns
    -------
    sample_wells : list of list
        Nested list where each sublist corresponds to one row of a 96-well plate
        mapped to a specific H-cell chamber.
        Example:
            [
                ['A1', 'A2', 'A3', ...],  # H1 donor wells
                ['B1', 'B2', 'B3', ...],  # H1 receptor wells
                ...
            ]

    Notes
    -----
    - For four H-cells (8 chambers total), the mapping allows two sample plates
    to be used, for up to 24 samples per chamber. In this case:
      H1 donor → Row A, H1 receptor → Row B, H2 donor → Row C, etc. for both plates.
    - Supports setups with multiple 96-well plates (either max 8 H-cells, or <=4 with
    extended sample collection).
    """
    
    # assigning wells in the 96-well sample plate
    n_plates = len(wellplates)//96
    row_letter = list(string.ascii_uppercase[:8])
    
    sample_wells = []
    
    if number_of_hcells == 4:
        for letter in row_letter:
            wells = []
            for row in wellplates:
                if letter in str(row).split(' ')[0]:
                    wells.append(row)
            sample_wells.append(wells)
    else:
        for i in range(n_plates):
            for letter in row_letter:
                wells = []
                for row in wellplates[i * 96 : (i+1) * 96]:
                    if letter in str(row).split(' ')[0]:
                        wells.append(row)
                sample_wells.append(wells)
                
    return sample_wells
    


def get_time_schedule(total_time, interval_time, interval_frequency):
    """
    Generates cumulative sampling timepoints and delays.

    Parameters
    ----------
    total_time : float
        Total duration of the experiment in hours.
    interval_time : list of float
        Duration of each sampling phase in hours.
    interval_frequency : list of float
        Sampling frequency per phase in hours.

    Returns
    -------
    schedule : list
        Cumulative list of sampling time points (hours).
    delays : list
        Delay durations between sampling points (hours).
    """
    assert sum(interval_time) == total_time, \
        "Mismatch between total experimental time and provided intervals."
    assert len(interval_time) == len(interval_frequency), \
        "Mismatch between interval counts and frequencies."

    schedule = []
    delays = []
    
    schedule.append(0) # sample aliquot at time zero
    delays.append(0) #sample aliquot at time zero
    hr = 0
    
    for i, time in enumerate(interval_time):
        instances = time//interval_frequency[i]
        for n in range(1, instances+1):
            hr = hr + interval_frequency[i]
            schedule.append(hr)
            delays.append(interval_frequency[i])

    return schedule, delays # in hours


def get_hcell_name(h_cell_wells):
    """
    Generates H-cell chamber names and mapping.

    Parameters
    ----------
    h_cell_wells : list
        Ordered list of H-cell wells (chambers), from get_Hcell_labware().

    Returns
    -------
    h_cell_name : list
        List of H-cell chamber names, e.g. ['H1_1', 'H1_2', 'H2_1', ...].
    h_cell_dict : dict
        Mapping of H-cell chamber names to corresponding wells.
    """
    
    h_cell_dict = {}
    n_hcell = len(h_cell_wells) // 2
    h_cell_name = []
    for i in range(n_hcell):
        h_cell_name.append('H{}_1'.format(i+1))
        h_cell_name.append('H{}_2'.format(i+1))
    
    for i, cell_well in enumerate(h_cell_wells):
        h_cell_dict[h_cell_name[i]] = cell_well
        
    return h_cell_name, h_cell_dict


def get_sample_schedule(h_cell_name, sample_wells, time_schedule, delay_time):
    """
    Creates a sampling schedule DataFrame for aliquot collection.

    Parameters
    ----------
    h_cell_name : list
        Names of each H-cell chamber.
    sample_wells : list
        Nested list of sample wells grouped by H-cell.
    time_schedule : list of float
        List of scheduled sampling points (hours).
    delay_time : list of float
        List of delays between scheduled time points (hours).

    Returns
    -------
    final_schedule : pandas.DataFrame
        DataFrame containing well positions, time stamps, and delays.
    """
    data = {}
    for i, name in enumerate(h_cell_name):
        data[name] = sample_wells[i]
    sample_schedule = pd.DataFrame(data)
    final_schedule = sample_schedule.truncate(after=len(time_schedule)-1)
    final_schedule['Time_Stamp_hr'] = list(map(float, time_schedule))
    final_schedule['Delay_hr'] = list(map(float, delay_time))
    final_schedule['Delay_min'] = [t * 60. for t in delay_time]
    
    return final_schedule


# -------------------------------------------------------------------
# Main Protocol Function
# -------------------------------------------------------------------

def H_cell_protocol(protocol, loaded_labware_dict,
                    sample_schedule, h_cell_info, aliquot_volume=20,
                    sample_dilution_volume=180, max_stock_vol=17,
                    log_file_info={}, log_filename='test',
                    dye_well_num=None, dye_volume=None):
    """
    Runs the OT-2 H-cell protocol including aliquot sampling and optional dye injection.

    Parameters
    ----------
    protocol : opentrons.protocol_api.ProtocolContext
        The OT-2 protocol context.
    loaded_labware_dict : dict
        Dictionary containing loaded labware and pipettes.
    sample_schedule : pandas.DataFrame
        Sampling schedule with well mappings and delays.
    h_cell_info : dict
        Mapping of H-cell compartment names to their wells.
    aliquot_volume : float, optional
        Volume (µL) for each aliquot. Default is 20.
    sample_dilution_volume : float, optional
        Volume (µL) of diluent per sample. Default is 180.
    max_stock_vol : float, optional
        Maximum volume (mL) per stock well before switching. Default is 17.
    log_file_info : dict, optional
        Metadata for logging.
    log_filename : str, optional
        Name of the log file. Default is 'test'.
    dye_well_num : int, optional
        Well number containing dye stock in the 12-well water stock plate.
    dye_volume : float, optional
        Volume of dye injected (µL) on the first iteration to each donor chamber.

    Notes
    -----
    - The dye injection step occurs **only during the first iteration**.
    """
    
    # create a log file that keeps track of protocol steps and saves to directory
    log_file = open(log_filename+ '.txt', 'a')
    log_file.write('Author: {}\n'.format(log_file_info['author']))
    log_file.write('Sample: {}\n'.format(log_file_info['sample_name']))
    log_file.write('Membrane: {}\n'.format(log_file_info['hcell_membrane']))
    log_file.write('Solvent: {}\n'.format(log_file_info['solvent']))
    log_file.write('Restock: {}\n'.format(log_file_info['restock_solvent']))
    log_file.write('---\n')

    protocol.home()
    start = time.time()

    start_message = 'experiment started: {}'.format(time.ctime())
    print(start_message)

    log_file.write(start_message)
    log_file.write('\n')

    log_file.close()
        
    # Retrieve pipettes and plates
    small_pipette = loaded_labware_dict['Small Pipette']
    large_pipette = loaded_labware_dict['Large Pipette']
    small_tiprack = loaded_labware_dict['Small Tiprack']
    large_tiprack = loaded_labware_dict['Large Tiprack']
    destination_plate = loaded_labware_dict['Destination Wells']
    
    # Checking if the OT2 has tips attached prior
    if small_pipette.has_tip:
        small_pipette.drop_tip()
    if large_pipette.has_tip:
        large_pipette.drop_tip()

    # Rename top/bottom halves to donor/receptor compartments
    h_cell_names = list(h_cell_info.keys())
    n_hcells = len(h_cell_names) // 2
    
    # variables to keep track of overall water volume dispensed
    # to avoid use of an empty stock well, once water_vol = max_stock_vol
    # pipette will switch to secondary water supply. ww is the well #
    water_vol = 0
    ww = 0
    
    # execution delay: accounts for the time required to complete a single
    # iteration and subtracts from overall delay time between iterations
    h_cell_col = [ col for col in sample_schedule.columns if col.startswith('H') ]
    elapsed_time= False
    time_stamp_schedule = np.zeros((len(sample_schedule), len(h_cell_col)))
    

    # Main protocol loop
    for i in range(len(sample_schedule)):
        
        iteration_message = 'ITERATION {}'.format(i)
        log_file = open(log_filename + '.txt', 'a')
        log_file.write(iteration_message)
        log_file.write('\n')
        log_file.close()
        
        water_well = loaded_labware_dict['Stock Wells'][ww]
        sample_wells = sample_schedule.iloc[i]
        
        tiprack = large_tiprack
        pipette = large_pipette



        # Perform aliquot sampling for each H-cell
        
        # apply delay to sample collection
        if elapsed_time != False:
            new_time= (sample_wells['Delay_min']- elapsed_time)
            print('Corrected Delay {:.2f}'.format(new_time)) 
            protocol.delay(minutes=(sample_wells['Delay_min']- elapsed_time))
        else:
            protocol.delay(minutes=sample_wells['Delay_min'])
            
        # Distribute dilution water to the sample wells first 
        s = time.time() # start time of iteration
        
        pipette.pick_up_tip(tiprack[i])
        
        pipette.transfer(sample_dilution_volume,
                               water_well,
                               [sw.top(-3) for sw in sample_wells[h_cell_col].to_list()],
                               new_tip='never',
                               blow_out=True,
                               blowout_location='destination well')
        large_pipette.return_tip()
        
        water_message = 'Water distribution completed.'
        log_file = open(log_filename + '.txt', 'a')
        log_file.write(water_message)
        log_file.write('\n')
        log_file.close()
        print(water_message)
            
        for m, n in enumerate(h_cell_col):
            
            # Inject dye during first iteration only
            if i == 0 and dye_volume and dye_well_num:
                dye_well = loaded_labware_dict['Stock Wells'][dye_well_num]

                if n.endswith('_1'):
                    tiprack = large_tiprack
                    pipette = large_pipette

                    pipette.transfer(dye_volume,
                                     dye_well.bottom(10),
                                     h_cell_info[n].top(-2),
                                     new_tip = 'always',
                                     blow_out = True,
                                     blowout_location = 'destination well')
                    dye_message = f"Injected {dye_volume}µL dye into {n}"
                    log_file = open(log_filename + '.txt', 'a')
                    log_file.write(dye_message)
                    log_file.write('\n')
                    log_file.close()
                    print(dye_message)
                else:
                    pass
            else:
                pass                    
            
            # collect aliquot sample from H-cell chamber
            pipette = small_pipette
            tiprack = small_tiprack
            
            if pipette == large_pipette and pipette.has_tip is True:
                pipette.return_tip()
            
            pipette.transfer(aliquot_volume,
                             h_cell_info[n].top(-60),
                             sample_wells[n].top(-7),
                             mix_after = (8,20),
                             blow_out = True,
                             blowout_location = 'destination well')
            aliquot_message = 'Aliquot # {} from {} was sampled at {}'.format(i, n, time.ctime())
            time_stamp_schedule[i,m] = time.time()
            
            log_file = open(log_filename + '.txt', 'a')
            log_file.write(aliquot_message)
            log_file.write('\n')
            log_file.close()
            print(aliquot_message)
            
            water_vol += aliquot_volume
        
        tiprack = large_tiprack
        pipette = large_pipette
        
        # replenish H-cells with water using the same aliquot volume
        pipette.pick_up_tip(tiprack[i])
        
        pipette.transfer(aliquot_volume,
                         water_well.bottom(5),
                         [hw.top(-30) for hw in list(h_cell_info.values())],
                         new_tip = 'never',
                         blow_out = True,
                         blowout_location = 'destination well')
        
        e = time.time() #end time for current iteration
        elapsed_time = round(e-s)/60 #time elapsed during current iteration
        
        water_vol += (sample_dilution_volume + aliquot_volume)*len(h_cell_names)
        
        # once the first water stock is empty, moves to second water well (+1 well position)
        if water_vol >= max_stock_vol:
            ww += 1
            water_vol = 0
        else:
            pass
        
        log_file = open(log_filename+ '.txt', 'a')
        water_message= 'The H-cells were replenished with water at {}'.format(
            time.ctime())
        print(water_message)

        log_file.write(water_message)
        log_file.write('\n')

        log_file.close()
        
        if small_pipette.has_tip is True:
            small_pipette.drop_tip()
        if large_pipette.has_tip is True:
            large_pipette.drop_tip()

    #for line in protocol.commands():
    #    print(line)
        
    # Keeping track of execution time. Will print total run time in minutes
    end = time.time()
    time_consumed = end-start

    print("This protocol took \033[1m{}\033[0m minutes to execute".format(
        np.round(time_consumed/60, 3)))
    
    