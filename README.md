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


## Functions

### `median_combine_fits(directory, file_pattern='*.fits')`
Median-combines FITS files to remove cosmic rays.

### `find_stars(image, fwhm, threshold_sigma)`
Detects stellar sources using DAOStarFinder algorithm.

### `filter_sources(sources, image_shape, ...)`
Applies quality checks to remove spurious detections.

### `create_catalog(sources)`
Creates clean catalog with object IDs and coordinates.

### `save_catalog(catalog, filename)`
Saves catalog to text file.

### `match_catalogs(catalog1, catalog2, max_separation, ...)`
Matches sources between two catalogs by proximity.

### `perform_photometry(image, catalog, exptime, zp, ...)`
Performs aperture photometry with local background subtraction.

### `create_photometry_catalog(matched_catalog, phot_f336w, phot_f555w, ...)`
Combines photometry from both filters into final catalog.

### `create_hr_diagram(catalog, output_filename)`
Creates and saves HR diagram plot.

## References

- DAOStarFinder: Stetson (1987), PASP, 99, 191
- WFPC2 Photometry: Holtzman et al. (1995), PASP, 107, 1065
