from opentrons import protocol_api

metadata = {"protocolName": "ADT_4_Hcell_Offset_Calibration",
            "author": "Pozzo Research Group",
            "description": "Protocol to use for establishing the individual labware offsets for the 4-Hcell diffusion testing layout.",
            "apiLevel": '2.8'
}


def test_pipette(pipette, tipracks):
    '''check the calibration for the pipette and associated tip racks
    Parameters. Can be modified to pick up and reeturn the tip in the H12
    position, too.
    -----------
    pipette : pipette to be tested
    tipracks : array of tipracks for the given pipette
    
    Returns
    ----------
    Using the given pipette, picks up the tip in A1 of each tip rack and
    returns it to the same spot.
    '''
    for rack in tipracks:
        pipette.pick_up_tip(rack['A1'])
        pipette.return_tip()
        #
        # pipette.pick_up_tip(rack['H12'])
        # pipette.return_tip()
    return


def run(protocol: protocol_api.ProtocolContext):
    ## LEFT PIPETTE ##

    # define all tipracks for the left pipette
    left_tiprack1 = protocol.load_labware('opentrons_96_tiprack_20ul', 8)
    left_tiprack2 = protocol.load_labware('opentrons_96_tiprack_20ul', 11)
    # define left pipette
    left_single = protocol.load_instrument(
        'p20_single_gen2', 'left' ,
         tip_racks = [left_tiprack1,left_tiprack2])
    # Pick up pipette tip in LTRs (left tiprack) A1 and H12

    test_pipette(left_single, [left_tiprack1, left_tiprack2])

    ## RIGHT PIPETTE ##

    # define all tipracks for the right pipette
    right_tiprack1 = protocol.load_labware('opentrons_96_tiprack_300ul', 5)

    # define left pipette
    right_single = protocol.load_instrument(
        'p300_single', 'right' ,
         tip_racks = [right_tiprack1])
    # Pick up pipette tip in LTRs (left tiprack) A1 and H12
    test_pipette(right_single, [right_tiprack1])

    ## STOCK LABWARES ##
    stock_1 = protocol.load_labware('20mlscintillationeven_12_wellplate_18000ul', 2)

    ## DESTINATION LABWARE ##
    dest_1 = protocol.load_labware('adt_hcell_5_6', 1)
    dest_2 = protocol.load_labware('adt_hcell_7_8', 4)
    dest_3 = protocol.load_labware('corning_96_wellplate_360ul_flat', 6)
    dest_4 = protocol.load_labware('corning_96_wellplate_360ul_flat', 9)
