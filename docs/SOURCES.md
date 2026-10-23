# Sources, rights and exact reuse

## Calculated plasma models

- **Gudiksen et al. (2011)**, *The stellar atmosphere simulation code Bifrost*, A&A 531, A154 — https://doi.org/10.1051/0004-6361/201116520
- **Hayek et al. (2010)**, radiative transfer with scattering — https://doi.org/10.1051/0004-6361/201014210
- **Carlsson et al. (2016)**, *A publicly available simulation of an enhanced network region of the Sun*, A&A 585, A4 — https://doi.org/10.1051/0004-6361/201527226; numeric source https://sdc.uio.no/vol/simulations/en024048_hion/atmos/
- **Cheung et al. (2019, online 2018)**, *A comprehensive three-dimensional radiative magnetohydrodynamic simulation of a solar flare*, Nature Astronomy 3, 160–166 — https://doi.org/10.1038/s41550-018-0629-3; numerical data https://purl.stanford.edu/dv883vb9686, **CC BY 3.0**. Only the photospheric numerical outgoing radiation and magnetic slices enter this movie; no source movie is used
- **Rempel (2017)**, *Extension of the MURaM Radiative MHD Code for Coronal Simulations* — https://doi.org/10.3847/1538-4357/834/1/10
- **Vögler et al. (2005)**, MURaM magnetoconvection — https://doi.org/10.1051/0004-6361:20041507
- **Martínez-Sykora, De Pontieu, Hansteen, Rouppe van der Voort, Carlsson & Pereira (2017)**, *On the generation of solar spicules and Alfvénic waves*, Science 356, 1269–1272 — https://doi.org/10.1126/science.aah5412
- **Martínez-Sykora, De Pontieu, Carlsson, Hansteen, Nóbrega-Siverio & Gudiksen (2017)**, *Two-dimensional Radiative Magnetohydrodynamic Simulations of Partial Ionization in the Chromosphere. II. Dynamics and Energetics of the Low Solar Atmosphere* — https://doi.org/10.3847/1538-4357/aa8866; full methods and public-data description https://arxiv.org/html/1708.06781; actual numeric source https://sdc.uio.no/vol/simulations/en096014_gol/atmos/
- **Martínez-Sykora et al. (2018)**, *Impact of Type II Spicules in the Corona: Simulations and Synthetic Observables* — https://arxiv.org/abs/1805.06475. Establishes the source model's multi-thermal behavior and discusses approximate equilibrium optically thin observables. It does not validate this project's finite-depth reconstruction

**Required MURaM attribution:** M. C. M. Cheung, M. Rempel, G. Chintzoglou, F. Chen, P. Testa, J. Martínez-Sykora, A. Sainz Dalda, M. L. DeRosa, A. Malanushenko, V. Hansteen, B. De Pontieu, M. Carlsson, B. Gudiksen and S. W. McIntosh. **Provided by the Stanford University Libraries.** Changes: radiative/color approximation, local-source sampling, spherical reconstruction, magnification, enhanced contrast, prescribed viewing geometry, optical imaging and video export. Original numerical data retain CC BY 3.0; this project's MIT code license does not replace it.

Bifrost data are publicly provided by ITA/University of Oslo and the IRIS modeling project for scientific reuse with citation of the code and experiment papers. No blanket MIT license is assigned to these data. The large source FITS files are retrieved by the scripts, not redistributed in the package. Exact variable, snapshot, byte range, header and payload hashes are retained in the source manifest.

## Radiation, color and imaging

- **CHIANTI atomic database**, freely available with requested scientific acknowledgement — https://www.chiantidatabase.org/ and https://www.chiantidatabase.org/chianti_download.html. Atomic files used are version **11.0.2**, from the public SolarSoft mirror https://sohoftp.nascom.nasa.gov/solarsoft/packages/chianti/dbase/. The exact C IV `elvlc`, `scups`, `wgfa` and `chianti.ioneq` files and checksums are retained. CHIANTI is a collaborative project involving George Mason University, University of Michigan, University of Cambridge and NASA Goddard
- **Dere et al. (1997)**, initial CHIANTI database, A&AS 125, 149; **Dere et al. (2023)**, version 10.1 revised ionization/recombination rates — https://arxiv.org/abs/2305.15221
- **Dufresne et al. (2024)**, CHIANTI XVIII, version 11 advanced ionization models — https://doi.org/10.3847/1538-4357/ad6765. **Limit:** this film uses the density-independent `chianti.ioneq` coronal approximation, not the advanced charge-transfer/density model; the paper explains material differences this can cause
- **Asplund et al. (2009)**, solar abundance reference — https://doi.org/10.1146/annurev.astro.46.060407.145222. Adopted carbon abundance is 8.43
- **Coddington et al. (2023)**, TSIS-1 HSRS v2 — https://doi.org/10.1029/2022EA002637; source spectral file and hash in the optical manifest
- **CIE 1931 standard observer**, CC BY-SA 4.0 data — https://cie.co.at/datatable/cie-1931-colour-matching-functions-2-degree-observer; provider metadata retained
- H⁻ coefficient tables follow Geltman (1962), Stilley & Callaway (1970) and the **MIT-licensed Lightweaver** implementation — https://github.com/Goobley/Lightweaver/blob/master/Source/Background.cpp. The original MIT notice is in `licenses/Lightweaver-MIT.txt`
- **Latin Modern** fonts retain their GUST Font License in `fonts/`

The source papers establish methods and model behavior. They do not independently validate this visualization's atmosphere assembly, display grading, interpolated frames, absolute C IV photometry or new numerical convergence claims. No scientists were contacted and no observed solar movie, NASA animation or generative-image asset is used as a texture.
