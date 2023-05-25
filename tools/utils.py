"""
Module for utilities used for CHMAP Package from PSI
"""

import sys
import numpy as np
import sunpy 
import matplotlib.pyplot as plt
import matplotlib as mpl
import tools.chmap.datatypes as psi_d_types
from tools.chmap.info import DTypes
import tools.chmap.ezseg.ezsegwrapper as ezsegwrapper
from tools.chmap.map_manip import combine_cr_maps, combine_maps

def carrington_rotation_number_relative(time, lon):
    """
    A function that returns the decimal carrington rotation number for a spacecraft position
    that may not be at the same place at earth. In this case you know the carrington longitude
    of the spacecraft, and want to convert that to a decimal carrington number that is within
    +0.5 and -0.5 of the decimal rotation for the earth-based longitude.

    :param time: an astropy Time object indicating the time the position is known.
    :param lon: the carrington longitude of the spacecraft position.
    :return: the decimal_carrington number.
    """
    # get the decimal carrington number for Earth at this time
    cr_earth = sunpy.coordinates.sun.carrington_rotation_number(time)

    # convert that to the earth longitude (this should match sunpy.coordinates.sun.L0(time))
    cr0 = np.floor(cr_earth)
    lon_earth = np.mod((1 - (cr_earth - cr0)*360), 360)

    # compute the angular difference and the modulus
    diff = lon_earth - lon
    mod = np.mod(diff, 360.)

    # compute the fractional rotation offset, which depends on where the periodic boundary is.
    offset = 0.0
    if lon_earth < 180 and mod < 180 and diff < 0:
        offset = +1.0
    if lon_earth >= 180 and mod >= 180 and diff >= 0:
        offset = -1.0
    cr_now = cr0 + np.mod(1.0 - lon/360., 360.) + offset

    debug = False
    if debug:
        print('{: 7.3f} {: 7.3f} {: 7.3f} {: 7.3f} {: 7.3f} {: 7.3f}'.format(lon, diff, mod, cr_now, cr_earth,
                                                                             cr_now - cr_earth))
        print(cr_earth, cr0, lon_earth, sunpy.coordinates.sun.L0(time).value, lon, cr_now)

    return cr_now


def get_metadata(map):
    """
    This function gets the metadata we need and then creates a dictionary at the end.
    - If we want to enforce specific types for each tag we "may" want to define a class
      that defines the specific metadata tags and corresponding types a priori and then
      this class is instantiatied and then populated in a subroutine like this.
    - however, this would have to be compatible with how the record type is defined in
      SQL and might be somewhat of a pain? i'm not sure what the best solution is
    """
    # Observation time is saved as a Time object
    time_object = map.date
    # For SQL, we want this in the Python native 'datetime' format
    time_datetime = map.date.datetime

    # Get the time as a string
    time_string = map.date.isot

    # get the time as a floating point julian date
    time_float = time_object.jd

    # get the wavelength as an integer (map.wavelength is an astropy quantity)
    # here I am converting the astropy distance quantity to angstrom and then a float to be sure
    wavelength = int(map.wavelength.to("angstrom").value)

    # make a string that gives a unique observatory/instrument combo [remove whitespace]
    # o_str = map.observatory.replace(" ","")
    # d_str = map.detector.replace(" ","")
    # instrument = o_str+'_'+d_str

    # or just use the sunpy nickname, which is also unique (i think i like this more...)
    instrument = map.nickname

    # get the distance of the observer (in km) from "observer_coordinate" (a SkyCoord object)
    # here I am converting the astropy distance quantity to km and then a float
    d_km = map.observer_coordinate.radius.to("km").value

    # get the carringtion longitude and latitude in degrees
    cr_lon = map.carrington_longitude.to("degree").value
    cr_lat = map.carrington_latitude.to("degree").value

    # get the decimal carrington rotation number (for this central longitude, not earth).
    cr_rot = carrington_rotation_number_relative(time_object, cr_lon)

    # now build a dictionary with the information
    # the idea here is to formalize the metadata components as a dictionary, which can
    # be used to create nice sliceable dataframes later with pandas
    metadata = dict()

    metadata['date_string'] = time_string
    metadata['datetime'] = time_datetime
    metadata['jd'] = time_float
    metadata['wavelength'] = wavelength
    metadata['instrument'] = instrument
    metadata['distance'] = d_km
    metadata['cr_lon'] = cr_lon
    metadata['cr_lat'] = cr_lat
    metadata['cr_rot'] = cr_rot

    return metadata


def PlotMap(map_plot, title=None, save_map=False, save_dir='maps/synoptic/'):
    """
    Super simple plotting routine for PsiMap objects.
    imshow() should be replaced with pcolormesh() for generalizing to non-uniform rectilinear grids
    OR use Ron's Plot2D from PSI's 'tools'
    """
    # set color palette and normalization (improve by using Ron's colormap setup)
    norm = mpl.colors.LogNorm(vmin=np.nanpercentile(map_plot.data.flatten(),5), vmax=np.nanpercentile(map_plot.data.flatten(),99.9))
    im_cmap = plt.get_cmap('sdoaia193')
    plot_mat = map_plot.data

    #  convert map x-extents to degrees
    x_range = [180 * map_plot.x.min() / np.pi, 180 * map_plot.x.max() / np.pi]

    # setup xticks
    xticks = np.arange(x_range[0], x_range[1] + 1, 30)

    # labels and title
    plt.xlabel("Carrington Longitude")
    plt.ylabel("Sine Latitude")
    plt.xticks(xticks)

    # plot the EUV data
    plt.imshow(plot_mat, extent=[x_range[0], x_range[1], map_plot.y.min(), map_plot.y.max()],
                    origin="lower", cmap=im_cmap, aspect=90.0, norm=norm)

    # plot the CHD data
    if map_plot.chd is not None:
        x_extent = np.linspace(x_range[0], x_range[1], len(map_plot.x))
        plt.contour(x_extent, map_plot.y, map_plot.chd, origin="lower", colors='white',
                        extent=[x_range[0], x_range[1], map_plot.y.min(), map_plot.y.max()], linewidths=0.5)

    if title is not None:
        plt.title(title)

    plt.show(block=False)

    # save the figure
    if save_map:
        plt.savefig(save_dir + str(map_plot.data_info.date_obs[0]))

    return None


def combine_map(euv_map, chd_map, euv_combined, chd_combined, mu_cutoff=0.0, mu_merge_cutoff=None,
           del_mu=None):
    
    # create map lists
    euv_maps = [euv_map, ]
    chd_maps = [chd_map, ]
    if euv_combined is not None:
        euv_maps.append(euv_combined)
    if chd_combined is not None:
        chd_maps.append(chd_combined)
     
    # combine maps with minimum intensity merge
    if del_mu is not None:
        euv_combined = combine_maps(euv_maps, del_mu=del_mu, mu_cutoff=mu_cutoff)
    else:
        euv_combined = combine_maps(euv_maps, mu_cutoff=mu_cutoff, mu_merge_cutoff=mu_merge_cutoff)
    

    return euv_combined


def cr_map(euv_map, euv_combined, n_images, mu_cutoff=0.0, mu_merge_cutoff=None,
           del_mu=None):
    
    # create map lists
    euv_maps = [euv_map, ]
    if euv_combined is not None:
        euv_maps.append(euv_combined)
     
    # combine maps with minimum intensity merge
    if del_mu is not None:
        euv_combined = combine_cr_maps(n_images, euv_maps, del_mu=del_mu, mu_cutoff=mu_cutoff)
    else:
        euv_combined = combine_cr_maps(n_images, euv_maps, mu_cutoff=mu_cutoff, mu_merge_cutoff=mu_merge_cutoff)
    

    return euv_combined


def chd(los_image, map_x, map_y, thresh1, thresh2, nc=3, iters=100, R0=1.01):
    # define chd parameters
    los_image.get_coordinates()
    mu_array = los_image.mu
    data = los_image.data

    # get mu index locations where mu is valid.
    use_indices = np.logical_and(mu_array > 0., data > 2.)

   # define chd parameters
    use_chd = use_indices.astype(int)
    use_chd = np.where(use_chd == 1, use_chd, los_image.no_data_val)
    nx = los_image.x.size
    ny = los_image.y.size

    # calculate new threshold parameters based off reference (AIA) instrument
    t1 = thresh1 # * ref_alpha + ref_x
    t2 = thresh2 # * ref_alpha + ref_x

    # fortran chd algorithm
    np.seterr(divide='ignore')
    ezseg_output, iters_used = ezsegwrapper.ezseg(np.log10(data), use_chd, nt=nx, np=ny, thresh1=t1, thresh2=t2, nc=nc, iters=iters)
    chd_result = np.logical_and(ezseg_output == 0, use_chd == 1)
    chd_result = chd_result.astype(int)

    # create CHD image
    chd_image = psi_d_types.create_chd_image(los_image, chd_result)
    chd_image.get_coordinates()

    # create EUV and CHD map
    euv_map = los_image.interp_to_map(R0=R0, map_x=map_x, map_y=map_y)
    chd_map = chd_image.interp_to_map(R0=R0, map_x=map_x, map_y=map_y)
    euv_map.chd = chd_map.data

    return euv_map