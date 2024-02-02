### 1.0
###Fri  2 Feb 15:18:40 GMT 2024
- add quantities in df_G
- add `show_df_G`
### Fri 26 Jan 13:11:38 GMT 2024
- change structure_factor dataframe callback
### Wed 24 Jan 10:26:49 GMT 2024
- added structure factor dataframe `pd.read_pickle(b0.get_Fhkl_pkl())``
### Fri 12 Jan 16:40:48 GMT 2024
- get_lattice is not loaded any more but computed with ut.get_lattice to fix bug when updating Nmax
### Thu  7 Dec 08:17:10 GMT 2023
- added frame parameter and custom excitation error calculation
### Thu 30 Nov 15:21:35 GMT 2023
- initial changelog for blochwave library
- structure factor storage for continuous rotation


##TODO :
- change structure factor to dataframe in _solve
