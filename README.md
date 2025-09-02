# OT2-ADT
A fork of the original OT2-DOE repo, which contains a group of python modules and notebooks made for high throughput measurement and analysis of samples made through a liquid handling robot (Opentrons or OT2). The ADT version includes specific modules, notebooks, and designs to implement an automated diffusion testing (ADT) setup on an OT2 platform.

    Note: As of 02/01/22 this repo is still being developed with new functions, bugs and documentation updated constantly. 
    Feel free to address issues in the issues tab or edit them yourself so long as you document and justify the change. 
    Currently the framework is still not entirely up and running as testing allows for adaptation for the groups use. 

This workflow has been designed to facilitate the testing and characterization of separations-based membranes using dual-compartment H-cells. 

![Depiction of using ]("assets\automated_diffusion_testing.png")


## Introduction
Automatic handling robots (ALH) are one of many high throughput tools that allow for increases discovery of solution-based chemistry. They can also allow for the automation of repetitive tasks, such as sample collection of multi-hour permeation testing. We have developed a workflow that allows for up to 8 H-cell membrane tests to be conducted in parallel. The robot (an OT2) extracts aliquots from each H-cell chamber at pre-designated frequencies and prepares them in a 96-well sample plate for UV-Vis characterization. The robot then replenishes the lost aliquot volume in each H-cell with restock solution (in this case, water).


## To get started
You can install this package using the GitHub repository in following steps: 
* In your terminal, run git clone https://github.com/pozzo-research-group/OT2-ADT.git
* Change the directory to `OT2_DOE` root directory, by running `cd OT2_DOE`
* (Recommended)- Create an environment using the provided `environment.yml` file. To do so, run the following lines:

	`conda env create --file environment.yml`
	
	`conda activate OT2_DOE`
* Install the package by running `python setup.py install` in your terminal

## Directory Organization
This repo contains all OT2-DOE content from the original pozzo-research-group/OT2-DOE. The relevant content for the ADT workflow is as follows:
- `OT2-DOE/`: directory with all files and modules
  - `ADT Data Analysis/`: directory with notebooks and scripts for analyzing UV-Vis data collected from the automated testing workflow. A walkthrough notebook with explanations has been provided, along with sample data in the resepctive sub-folders. A blank template of the noteobook is also available.
  - `Custom Labware/`: directory with customized hardware created for use with the OT2. The following .json files are the labware definitions created using the Opentrons custom labware tool. Each hcell file is unique to a specific dual H-cell assembly, as the custom builds and glassware all have slightly different measurements.
    - `20mlscintillationeven_12_wellplate_18000ul.json`
    - `adt_hcell_1_2.json` (for the 22 mL total volume H-cells)
    - `adt_hcell_3_4.json` (for the 22 mL total volume H-cells)
    - `adt_hcell_5_6.json` (for the 30 mL total volume H-cells)
    - `adt_hcell_7_8.json` (for the 30 mL total volume H-cells)
    - `adt_hcell_9_10.json` (for the 30 mL total volume H-cells)
  - `Labware Designs/`:
    - contains the .stl and .step files for all 3D-printed hardware used in this workflow.
  - `Plan/`: 
    - `CreateSamples.py`
  - `Prepare/`:
    - `OT2Commands.py`
    - `OT2_Hcell_commands.py`
  - `Chemical Database.csv`
  - `H_cell_protocol_4hcell.csv`
  - `H_cell_protocol_8hcell.csv`
  - `ADT_Hcell_Testing.ipynb`

### Important!

- The OT2 layout protocols (`H_cell_protocol_4hcell.csv`, `H_cell_protocol_8hcell.csv`) do not contain the robot-specific calibration offsets. These must be added before implementing the workflow.
- the custom labware definitions for the H-cells (e.g., `adt_hcell_1_2.json`) ARE NOT equivalent to the wellplate designs in `Labware Designs/`. The opentrons defintion consists of the **entire dual H-cell assembly**. The labware definition simplifies the 4-component assembly into one very tall wellplate with round, flat-bottom wells. Therefore, all well spacing and depth measurements must be taken from this full assembly.

![Figure depiction of a "dual H-cell assembly", which will be defined as a well plate in the Opentrons custom labware library.]("assets\dual_hcell_assembly.png")

- The H-cell-to-sample plate mapping has been designed so that each H-cell chamber has a designated row in the sample plate, and each iteration (round of sampling) has a designated column.

|                | Iteration | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H-cell chamber |           | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10| 11 | 12 |
| H1_C1          |    A      |   |   |   |   |   |   |   |   |   |   |    |    |
| H1_C2          |    B      |   |   |   |   |   |   |   |   |   |   |    |    |
| H2_C1          |    C      |   |   |   |   |   |   |   |   |   |   |    |    |
| H2_C2          |    D      |   |   |   |   |   |   |   |   |   |   |    |    |
| H3_C1          |    E      |   |   |   |   |   |   |   |   |   |   |    |    |
| H3_C2          |    F      |   |   |   |   |   |   |   |   |   |   |    |    |
| H4_C1          |    G      |   |   |   |   |   |   |   |   |   |   |    |    |
| H4_C2          |    H      |   |   |   |   |   |   |   |   |   |   |    |    |


## Using this workflow
### In advance: hardware printing and labware definitions

The 3D designs for the customized hardware can be found in `Labware Designs/`. The H-cell well plates are designed based on the H-cell specifications that can be found in the `assets/` folder. Notably, the 22 mL H-cells do not sit flush to a surface when assembled with the clamp, hence the added material below each chamber in the well plate design. 

Dual stir plate construction:

The dual stir plates were constructed by 3D printing the Dual_Stir_Plate_base and Dual_Stir_Plate_lid, then outfitting each unit with the components of two commerically available stir plates. The base and lid of the unit were designed as separate pieces so that the motors are connected to the lid and can separated from the base. The support columns on the lid were designed specifically for the stir plate design from [INTLLAB](https://www.amazon.com/Magnetic-stirrer-magnetic-Stirring-Capacity/dp/B072K24X5P). If using different hardware, the lid can be edited using the .step file. Images of the interiors of the stir plates can be found in the `assets/` folder.

Once the dual H-cell assembly is created, a custom labware definition must be created for each unique assembly and added to the opentrons library. You can either create a new file using the labware definition creator tool from Opentrons, or edit the existing files in `Custom Labware/`. Calibration and offsets must also be completed.

### Before starting the workflow

Membranes should be assembled in the H-cells. H-cells should be placed in the H-cell holders and outfitted with stir bars, then placed on top of the dual stir plates. The water/solutions should be added immediately before initiating the protocol. If utilizing the dye injection step in the automated protocol, then make sure to account for the added volume in each donor cell when that step is executed. 

- e.g. for a 22 mL H-cell in which 80 uL of the stock dye will be used:
  - with NO injection step: add 10.920 mL water to donor chamber (C1) and 11 mL water to receptor chamber (C2). Immediately before executing the protocol, add the 80 uL of stock dye to the donor chamber.
  - WITH the injection step: add 10.920 mL and 11 mL water to the donor and receptor chambers, respectively. Arrange the H-cells in the pre-determined layout in the OT2, and initiate the protocol. The robot will distribute the 80 uL of stock dye to each donor chamber (C1).
  
The `ADT_Hcell_Testing.ipynb` notebook explains each step of the execution. There are several experiment-specific parameters that should be defined in advance:
- `path`: path used to obtain the protocol. The protocol could be `H_cell_protocol_4hcell.csv`, `H_cell_protocol_8hcell.csv`, or your own protocol with the same dictionary definitions and stucture as the provided protocols.
- `log_metadata`: a dictionary containing experiment-specific metadata, useful for record keeping purposes. The information provided here will be saved to the output log file created during the experiment.
- `aliquot_volume`: volume to be extracted from each H-cell chamber (in uL)
- `sample_dilution_volume`: volume of water to be distributed to each sample well in the given iteration, to which the sample aliquot will be mixed and diluted with. Aliquot volume + sample dilution volume must be less than the max volume of the sample well (here, 300 uL).
- `time_schedule, delay_time`: need inputs for total experimental time (hrs), how long each sampling freuquency will last (hrs), and which sampling frequencies to use (hrs). The total number of samples collected (total time / sampling frequency) cannot exceed available spots in the well plate. If using only 4 H-cells, up to 24 samples can be collected from each H-cell chamber. If using more than 4 H-cells (up to 8), there is a maximum of 12 samples per chamber.
- `HC.H_cell_protocol()`: the only inputs in this function that need to be edited are the `log_filename`, `dye_well_num` (location of the dye stock, if using the dye injection option), and `dye_volume` (amount of dye stock to be distributed to each donor chamber). Defaults for the last two inputs are `None`.
  
### Implementing the protocol

One the full OT2 layout is assembled and all parameters are determined in the testing notebook, run all of the cells in the notebook. The robot will begin immediately with iteration 0. It will distribute the dilution water to each sample well, then switch pipettes to extract each aliquot, distribute to the assigned sample well, and mix the diluted aliquot. After all aliquots have been collected, the robot will switch pipette tips again and distribute the equivalent volume in water to the aliquot volume, to each H-cell chamber.

Upon completion of the iteration, the robot will calculate the amount of time remaining until the next round of sampling, based on the pre-defined sample frequency. During these "down times," the sample well plates can be removed from the plaform for UV-Vis characterization, then returned to the platform.

The `ADT Testing Analysis/` folder has been provided as an optional guidance for processing and analyzing the data collected from this workflow. 
