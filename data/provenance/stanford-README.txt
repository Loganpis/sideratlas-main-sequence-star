Stanford Digital Repository (https://purl.stanford.edu/dv883vb9686) entry for
 
"A comprehensive three-dimensional radiative magnetohydrodynamic simulation of a solar flare", by
M. C. M. Cheung, M. Rempel, G. Chintzoglou, F. Chen, P. Testa, J. Martinez-Sykora, A. Sainz Dalda, M. L. DeRosa, 
A. Malanushenko, V. Hansteen, B. De Pontieu, M. Carlsson, B. Gudiksen & S. W. McIntosh
published in Nature Astronomy on Nov 26 2018
https://www.nature.com/articles/s41550-018-0629-3
https://doi.org/10.1038/s41550-018-0629-3

Corresponding author: Mark Cheung (mark.cheung@gmail.com)

(Reminder for the HGCR team: The following data set is from the run HGCR_Flare_Rerun_Jan_24_2017.)

Instructions for loading simulation data into IDL:

3D MHD data cubes
-----------------
See preamble of hgcr_muram_cube.pro for instructions on how to load 3D snapshots into IDL. 

2D slices of MHD quantities sampled at constant geometrical height or optical depth 
-----------------------------------------------------------------------------------
Make sure files such as xz_slice_0116* and tau_slice_0.100* are in a subdirectory named "2D" under the working directory.
See preamble of hgcr_tauslice.pro for instructions on how to load 2D slices of MHD quantities sampled at tau=0.1 (in the photosphere). 
See preamble of hgcr_hslice.pro for instructions on how to load 2D slices of MHD quantities sampled at constant height. For photosphere (height where <tau> = 1), set input parameter k = 116. 

2D slices of vertical emergent grey intensity
---------------------------------------------
Make sure files such as I_out.* are in a subdirectory named "2D" under the working directory. 
See preamble of hgcr_iout.pro for instructions on how to load intensity images.

Emission Measure Cubes
----------------------
Make sure files such as corona_emission_adj_dem_y* and corona_emission_adj_dem_z* are in a subdirectory named "2D" under the working directory.
See preamble of hgcr_muram2em.pro for instructions on how to load these cubes into IDL. 

Timestamps of snapshots
-----------------------
See preamble of hgcr_tauslice.pro or hgcr_hslice.pro to see how to get timestamps of each snapshot. 

Python
------
Use the module muram_io.py provided as part of the repository.
