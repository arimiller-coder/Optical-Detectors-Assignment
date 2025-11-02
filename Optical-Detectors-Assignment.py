import numpy as np
from astropy.io import fits
import glob
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
from astropy.visualization import ZScaleInterval
from astropy.table import Table
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry
def median_combine_fits(directory, file_pattern='*.fits'):
    """
    Median-combine all FITS files in a specified directory to remove cosmic rays.
    """

    # Find all FITS files in the directory
    search_pattern = directory + '/' + file_pattern
    fits_files = sorted(glob.glob(search_pattern))

    # Load all images into a list
    image_stack = []
    header = None

    for i, fits_file in enumerate(fits_files):
        with fits.open(fits_file) as hdul:
            data = hdul[1].data  # Data is in extension 1, not 0

            # Store the header from the first valid file
            if header is None:
                header = hdul[1].header.copy()  # Get header from extension 1

            image_stack.append(data)

    # Convert list to numpy array
    image_stack = np.array(image_stack)

    # Perform median combination along the first axis
    combined_data = np.median(image_stack, axis=0)

    return combined_data, header

def find_stars(image, fwhm, threshold_sigma):
    """
    Find stellar sources in an image using DAOStarFinder.

    """
    # Step 1: Estimate the background using sigma-clipped statistics
    # This removes bright stars from background calculation
    mean, median, std = sigma_clipped_stats(image)

    # Step 2: Set up DAOStarFinder
    # threshold = median + (threshold_sigma * std)
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold_sigma * std)

    # Step 3: Find sources (DAOStarFinder automatically subtracts the background)
    sources = daofind(image - median)  # Subtract median background

    return sources, mean, median, std

def filter_sources(sources, image_shape, flux_min=0.0, ellipticity_max=0.5,
                   sharpness_min=0.2, sharpness_max=1, left_buffer=40, right_buffer=0, 
                   top_buffer=0, bottom_buffer=50):
    """
    Apply quality checks to filter detected sources and remove spurious detections.

    """

    # Create boolean mask for sources that pass all checks
    mask = np.ones(len(sources), dtype=bool)

    # Check 1: Positive flux only
    flux_check = sources['flux'] > flux_min
    mask &= flux_check

    # Check 2: Ellipticity/Roundness
    # DAOStarFinder provides 'roundness1' and 'roundness2' 
    # Both should be close to 0 for circular sources
    roundness1_check = np.abs(sources['roundness1']) <= ellipticity_max
    roundness2_check = np.abs(sources['roundness2']) <= ellipticity_max
    roundness_check = roundness1_check & roundness2_check
    mask &= roundness_check

    # Check 3: Sharpness (filters extended sources and cosmic rays)
    # Sharpness should be around 0.2-1.0 for stars
    sharpness_check = (sources['sharpness'] >= sharpness_min) & (sources['sharpness'] <= sharpness_max)
    mask &= sharpness_check

    # Check 5: Edge masking: reject sources too close to edges
    height, width = image_shape
    x = sources['xcentroid']
    y = sources['ycentroid']
    edge_check = ((x > left_buffer) & (x < width - right_buffer) & 
                  (y > bottom_buffer) & (y < height - top_buffer))
    mask &= edge_check

    # Apply the combined mask
    filtered_sources = sources[mask]

    return filtered_sources

def create_catalog(sources):
    """
    Create a catalog of sources with unique IDs and coordinates.

    """

    # Create object IDs starting from 1
    object_ids = np.arange(1, len(sources) + 1)

    # Extract x and y centroids
    x_centers = sources['xcentroid']
    y_centers = sources['ycentroid']

    # Create new catalog table with column names
    catalog = Table()
    catalog['Object_ID'] = object_ids
    catalog['x_center'] = x_centers
    catalog['y_center'] = y_centers

    return catalog

def save_catalog(catalog, filename):
    """
    Save a source catalog to a text file.

    """

    # Save as space-delimited text file with header
    catalog.write(filename, format='ascii.fixed_width', overwrite=True)

def match_catalogs(catalog1, catalog2, max_separation=2.5, 
                   filter1_name='F336W', filter2_name='F555W'):
    """
    Match sources between two catalogs based on coordinate proximity.

    """

    # Extract coordinates
    x1 = np.array(catalog1['x_center'])
    y1 = np.array(catalog1['y_center'])
    x2 = np.array(catalog2['x_center'])
    y2 = np.array(catalog2['y_center'])

    # Track which sources have been matched
    matched_indices1 = []
    matched_indices2 = []
    separations = []

    # For each source in catalog1, find closest source in catalog2
    for i in range(len(catalog1)):
        # Calculate distances to all sources in catalog2
        dx = x2 - x1[i]
        dy = y2 - y1[i]
        distances = np.sqrt(dx**2 + dy**2)

        # Find closest source
        min_dist_idx = np.argmin(distances)
        min_dist = distances[min_dist_idx]

        # Check if within matching radius and not already matched
        if min_dist <= max_separation and min_dist_idx not in matched_indices2:
            matched_indices1.append(i)
            matched_indices2.append(min_dist_idx)
            separations.append(min_dist)

    # Create matched catalog
    matched_catalog = Table()
    matched_catalog['Object_ID'] = np.arange(1, len(matched_indices1) + 1)

    # Average the positions from both filters
    matched_catalog['x_center'] = (x1[matched_indices1] + x2[matched_indices2]) / 2.0
    matched_catalog['y_center'] = (y1[matched_indices1] + y2[matched_indices2]) / 2.0

    # Store original IDs for reference
    matched_catalog['f1_id'] = catalog1['Object_ID'][matched_indices1]
    matched_catalog['f2_id'] = catalog2['Object_ID'][matched_indices2]
    matched_catalog['separation'] = separations

    return matched_catalog

def perform_photometry(image, catalog,exptime,zp=21.1, aperture_radius=4.5, 
                       annulus_inner=6.0, annulus_outer=8.0):
    """
    Perform aperture photometry on sources in a catalog.
    """

    # Create apertures at source positions
    positions = np.transpose(np.array([catalog['x_center'], catalog['y_center']]))
    apertures = CircularAperture(positions, r=aperture_radius)
    annulus = CircularAnnulus(positions, r_in=annulus_inner, r_out=annulus_outer)

    # Perform aperture photometry
    phot_table = aperture_photometry(image, apertures)
    bkg_table = aperture_photometry(image, annulus)

    # Calculate mean background per pixel in the annulus
    annulus_area = annulus.area
    aperture_area = apertures.area
    background_mean = bkg_table['aperture_sum'] / annulus_area

    # Calculate total background in the source aperture
    background_sum = background_mean * aperture_area

    # Subtract background from aperture sum to get source flux
    flux = phot_table['aperture_sum'] - background_sum

    # Convert flux to magnitude: m = -2.5 × log10(flux)
    # Handle negative/zero fluxes by setting them to NaN
    mag = np.where(flux > 0, -2.5 * np.log10(flux/exptime)+zp, np.nan)

    # Add results to table
    phot_table['background_mean'] = background_mean
    phot_table['background_sum'] = background_sum
    phot_table['flux'] = flux
    phot_table['mag'] = mag

    # Remove the default id column and add catalog info
    phot_table.remove_column('id')
    phot_table['Object_ID'] = catalog['Object_ID']
    phot_table['x_center'] = catalog['x_center']
    phot_table['y_center'] = catalog['y_center']

    # Count valid magnitudes
    n_valid = np.sum(~np.isnan(mag))

    return phot_table

def create_photometry_catalog(matched_catalog, phot_f336w, phot_f555w, 
                              aperture_radius=4.5):
    """
    Create final photometry catalog combining measurements from both filters.
    """

    # Create new table with matched catalog info
    final_catalog = Table()
    final_catalog['Object_ID'] = matched_catalog['Object_ID']
    final_catalog['x_center'] = matched_catalog['x_center']
    final_catalog['y_center'] = matched_catalog['y_center']
    final_catalog['aperture_radius'] = aperture_radius

    # Add magnitudes from both filters
    final_catalog['mag_F336W'] = phot_f336w['mag']
    final_catalog['mag_F555W'] = phot_f555w['mag']

    # Remove sources with invalid magnitudes in either filter
    valid_mask = ~np.isnan(final_catalog['mag_F336W']) & ~np.isnan(final_catalog['mag_F555W'])
    final_catalog = final_catalog[valid_mask]

    return final_catalog


def create_hr_diagram(catalog, output_filename=None):
    """
    Create a Hertzsprung-Russell (HR) Diagram from photometry catalog.

    """

    # Calculate color (difference between filters)
    color = catalog['mag_F336W'] - catalog['mag_F555W']

    # Use F336W magnitude for brightness
    magnitude = catalog['mag_F336W']

    # Create the plot
    plt.figure(figsize=(10, 8))

    # Plot points
    plt.scatter(color, magnitude, s=5, alpha=0.6, c='black')

    # Invert y-axis (brighter stars at top)
    plt.gca().invert_yaxis()

    # Labels
    plt.xlabel('Color (F336W - F555W)', fontsize=14)
    plt.ylabel('Magnitude (F336W)', fontsize=14)
    plt.title('Hertzsprung-Russell Diagram', fontsize=16)

    # Grid for readability
    plt.grid(True, alpha=0.3)

    # Save if filename provided
    if output_filename:
        plt.savefig(output_filename, dpi=300, bbox_inches='tight')
        print(f"Saved HR diagram to: {output_filename}")

    plt.show()

    return color, magnitude

if __name__ == "__main__":
    # Step 1: Cosmic Ray Removal
    combined_f336w, header_f336w = median_combine_fits('data/F336W')
    combined_f555w, header_f555w = median_combine_fits('data/F555W')

    # Step 2: Star Finding
    sources_f336w, mean, median, std = find_stars(combined_f336w, fwhm=2.5, threshold_sigma=3)
    sources_f555w, mean, median, std = find_stars(combined_f555w, fwhm=2.5, threshold_sigma=3)

    # Step 3: Quality Filtering
    filtered_f336w = filter_sources(sources_f336w, combined_f336w.shape)
    filtered_f555w = filter_sources(sources_f555w, combined_f555w.shape)

    # Step 4: Create Individual Catalogs
    catalog_f336w = create_catalog(filtered_f336w)
    catalog_f555w = create_catalog(filtered_f555w)
    save_catalog(catalog_f336w, 'F336W_catalog.txt')
    save_catalog(catalog_f555w, 'F555W_catalog.txt')

    # Step 5: Match Sources Between Filters
    matched_catalog = match_catalogs(catalog_f336w, catalog_f555w, max_separation=1.5)
    save_catalog(matched_catalog, 'matched_catalog.txt')

    # Step 6: Photometry
    phot_f336w = perform_photometry(combined_f336w, matched_catalog,exptime=700, aperture_radius=4.5)
    phot_f555w = perform_photometry(combined_f555w, matched_catalog,exptime=70, aperture_radius=4.5)

    # Step 7: Create Final Photometry Catalog
    final_catalog = create_photometry_catalog(matched_catalog, phot_f336w, phot_f555w, aperture_radius=4.5)
    save_catalog(final_catalog, 'final_photometry_catalog.txt')

    # Step 8: Create HR Diagram
    color, magnitude = create_hr_diagram(final_catalog, output_filename='hr_diagram.png')



