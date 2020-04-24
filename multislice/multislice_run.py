from multislice import Multislice


def Si100_mulslic():
    multi=Multislice('../decks/Silicon/Si110',keV=200,
        repeat=[6,6,40],NxNy=[512,512],slice_thick=1.3575,
        mulslice=False,v='D')
    #multi.run()         # creates the deck and run it
    # multi.beam_vs_thickness()
    # multi.structure_factor()
    # multi.diffraction_pattern()
    # multi.diffraction_section()
    return multi


multi=Si100_autoslic()
