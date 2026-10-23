# Third-party material

- `continuum-coefficients.json`: numerical tables adapted from the MIT-licensed Lightweaver `Source/Background.cpp`, copyright Christopher M. J. Osborne (2019). The MIT notice is retained in `licenses/Lightweaver-MIT.txt`. Scientific origins: Geltman (1962), Stilley & Callaway (1970), Mihalas (1978)
- `CIE_xyz_1931_2deg.csv`: unmodified CIE 1931 2° color-matching functions, DOI [10.25039/CIE.DS.xvudnb9b](https://doi.org/10.25039/CIE.DS.xvudnb9b), International Commission on Illumination, 2019. CC BY-SA 4.0 as specified in the accompanying official metadata. Numerical integration and display conversion are this project's calculations
- Latin Modern font files: GUST Font License, retained under `fonts/`
- Bifrost atmosphere subsets: retrieved from the University of Oslo's public archive; citation/public-use statement described in `docs/SOURCES.md`. Raw source arrays are not included in the distributable package
- TSIS-1 HSRS v2: LASP reference spectrum, Coddington et al. (2023); retrieved by `fetch_auxiliary.py`; the exact reference file is included for reproducibility. Institutional provenance, official URL and SHA-256 are retained in `source-manifest.json`

The MIT license for this project's original code does not relicense any third-party scientific data.

- Active-region MURaM data: Cheung et al. (2019; online 2018), Stanford Digital Repository, [CC BY 3.0](https://creativecommons.org/licenses/by/3.0/), [source](https://purl.stanford.edu/dv883vb9686). Attribution, author credit, and modifications are detailed in `docs/ACTIVE-REGION-ATTRIBUTION.md`. Numeric gray radiation and magnetic fields are downloaded from this source; the documentation figure is derived from those data. The source authors’ movie footage is not used

CHIANTI atomic files and the public GOL model retain their original provider rights; see docs/SOURCES.md.
