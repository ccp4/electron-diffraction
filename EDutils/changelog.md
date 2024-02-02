# Changelog
## 1.2dev
### Fri  2 Feb 14:52:12 GMT 2024
- rotate_exp :
  - added filtering for full rocking curve and some info into `beams`
  - added `show_excitation_map` 
- utilities :
  - added `get_uvw` about an orientation
  - added `import_fcf`
  - in `to_shelx` rescale max to 9999.99
### Fri 12 Jan 16:40:48 GMT 2024
- utilities :
  - added axis_and_angle_from_R
  - shelx misc tweak
  - compute_B from parameters (Dials convention)
- importED : get_UB and get_RUB base definition
- dials_utils :
  - fine orientation import and get_RUB

### Thu  7 Dec 08:18:49 GMT 2023
- added frame support for rotate_exp
- dials_utils : more robust import of reflection.txt

### Thu 30 Nov 15:21:35 GMT 2023
- removed minor changes in rotate_exp

### Wed 15 Nov 14:32:42 GMT 2023
- pets :
  - more stable import pets info
  - fix hkl S line in
  - nxy from tiff
- dials
  - added import capabilities for multi-panel detectors
  - import refl.txt from integrated.expt
  - get image_size
  - Sw excitation error
- rotate_exp
  - change_path fix
  - frame parameter `sweep_var`
### Mon 25 Sep 17:05:13 BST 2023
- fix 'img' reader bug so it does not automatically fails if TYPE not in header
### Tue 18 Jul 14:05:19 BST 2023
- added path option to import pets
- added manual integration to dials

### Fri  9 Jun 17:51:51 BST 2023
- quick update to the smv reader

### Fri  9 Jun 15:08:57 BST 2023
- added image readers in readers.py
- ported read_space_group
- added test for image_readers
- added smv/img format reader

### Thu  8 Jun 10:38:03 BST 2023
- added change path functionality
