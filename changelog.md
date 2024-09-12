# Changelog Summaries
## 1.6
- installation bug fixes 

## 1.3.1

## 1.2
- EDutils 1.2:
  - pets :
  - dials :
    - dev : fine orientation import and `get_RUB`
    - ref : more robust import of reflection.txt
  - rotate_exp :
    - dev : added `kin` flag to `plot_rocking`
    - dev : added filtering for full rocking curve and some info into `beams`
    - dev : added `show_excitation_map`
    - dev : added frame support
    - dev : `params` and `vals` (such as `hkl`) can be added to the continuous sweep
  - utilities :
    - dev : added `get_uvw` about an orientation    
    - dev : added export function `to_shelx`
    - dev : added `axis_and_angle_from_R`
- blochwave 1.0:
  - bug fix : fixed kinematic intensities
  - dev : added frame parameter and custom excitation error calculation
  - dev : added perturbation strength to Vg
  - dev : added `solve_strong_beam` to iteratively remove unnecessary beams
  - dev : added `show_df_G`
  - ref : - added structure factor dataframe `pd.read_pickle(b0.get_Fhkl_pkl())` as loadable structure factors
  - ref : changes to `get_lattice`
## 1.1.0

### EDutils
- xds importer
- dials importer
- EDimporter (base class for xds and dials)
- pets:
  - added negative sign already in `pets.uvw`
  - added dyngo functionalities
- rotate_exp
  - added dataframe based `integrate` function
- utilities :
  - added rodrigues rotation matrix
  - added convert `to_shelx`
### blochwave
- added dyngo support
- changes to frame_generation
- changes to load bloch

##1.0.7
### blochwave
- felix wrapper from Bloch._solve_Felix
- bloch_util.strong_beam condition selector
- bug fix in Bloch.convert_tiff  
