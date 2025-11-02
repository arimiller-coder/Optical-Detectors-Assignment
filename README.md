# HR Diagram Pipeline for NGC 1261

This project creates a Hertzsprung-Russell (HR) diagram from Hubble Space Telescope WFPC2 observations of the globular cluster NGC 1261 using F336W and F555W filters.

## Overview

The pipeline performs the following steps:
1. **Cosmic Ray Removal**: Median-combines multiple exposures to remove cosmic rays
2. **Star Finding**: Detects stellar sources using DAOStarFinder
3. **Quality Filtering**: Applies quality checks (flux, roundness, sharpness, edge masking)
4. **Source Cataloging**: Creates catalogs with object IDs and coordinates
5. **Cross-Matching**: Matches sources between the two filters
6. **Aperture Photometry**: Measures stellar magnitudes in both filters
7. **HR Diagram**: Creates color-magnitude diagram

### Output Files

The pipeline generates:
- `F336W_catalog.txt`, `F555W_catalog.txt` - Individual filter catalogs
- `matched_catalog.txt` - Cross-matched sources
- `final_photometry_catalog.txt` - Final catalog with magnitudes in both filters
- `hr_diagram.png` - HR diagram plot


## References

- Stetson, P. B. (1987). DAOPHOT: A Computer Program for Crowded-Field Stellar Photometry. PASP, 99, 191. https://doi.org/10.1086/131977
- Holtzman, J. A., et al. (1995). The Photometric Performance and Calibration of WFPC2. PASP, 107, 1065. https://doi.org/10.1086/133664
- DAOStarFinder documentation: https://photutils.readthedocs.io/en/stable/api/photutils.detection.DAOStarFinder.html
- sigma_clipped_stats documentation: https://docs.astropy.org/en/stable/api/astropy.stats.sigma_clipped_stats.html


