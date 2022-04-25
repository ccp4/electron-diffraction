# Electron diffraction tools

Some python wrappers for common electron diffraction simulation codes and utilities for continuous rotation electron diffraction.


## Some figures produced
Bloch wave simulations performed for diamond.

N denotes the configuration :

- N=2 for 2 strongly excited beams,
- N=4+2w for 4 strongly exited beams and 2 weak contributing beams,
r1 for a random orientation.


N   | a  | b  | c
--- | -- | -- | --
2   | [![](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_2beam_Sw.svg)](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_2beam_Sw.svg)     | [![](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_2beam_Iz.svg)](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_2beam_Iz.svg) | [![](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_2beam_beams0.svg)](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_2beam_beams0.svg)
3   | [![](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_3beam_Sw.svg)](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_3beam_Sw.svg)     | [![](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_3beam_Iz.svg)](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_3beam_Iz.svg) | [![](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_3beam_beams0.svg)](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_3beam_beams0.svg)
3+1w| [![](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_3_1beam_Sw.svg)](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_3_1beam_Sw.svg) | [![](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_3_1beam_Iz.svg)](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_3_1beam_Iz.svg) | [![](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_3_1beam_beams0.svg)](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_3_1beam_beams0.svg)
4   | [![](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_4beam_Sw.svg)](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_4beam_Sw.svg)     | [![](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_4beam_Iz.svg)](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_4beam_Iz.svg) | [![](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_4beam_beams0.svg)](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_4beam_beams0.svg)
4+2w| [![](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_4_2beam_Sw.svg)](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_4_2beam_Sw.svg) | [![](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_4_2beam_Iz.svg)](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_4_2beam_Iz.svg) | [![](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_4_2beam_beams0.svg)](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_4_2beam_beams0.svg)
r1  | [![](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_r1_Sw.svg)](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_r1_Sw.svg) | [![](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_r1_Iz.svg)](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_r1_Iz.svg) | [![](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_r1_beams0.svg)](https://www.ccp4.ac.uk/ccp4-ed/figures/bloch/diamond_r1_beams0.svg)

a) Setup displaying excitation errors (orange map) and structure factor (diamond) at this beam orientation.
b) Intensity I_z as function of thickness for selected beams.
c) Rocking curve for selected beam at various thicknesses.
