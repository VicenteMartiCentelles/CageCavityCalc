from scipy.optimize import basinhopping, minimize
from scipy.spatial.distance import cdist
import numpy as np
from scipy.spatial import distance_matrix

from CageCavityCalc.data import vdw_radii
from CageCavityCalc.log import logger
#This is literally copied from cgbind

def get_max_escape_sphere(cage_coords, atom_names, basinh=False, max_dist_from_metals=10):
    """
    Get the maximum radius of a sphere that can escape from the centroid of
    the cage – will iterate through all theta/phi

    :param basinh: (bool) Find the true maximum escape sphere by basin
                   hopping on the surface
    :param max_dist_from_metals: (float) Distance in Å on top of the
                                 average M-M distance that will be used for
                                the search for the maximum escape sphere
    :return: (float) Volume of the maximum escape sphere in Å^3
    """
    # logger.info('Getting the volume of the largest sphere that can escape '
    #            'from the cavity')

    max_sphere_escape_r = 99999999999999.9

    # we approximate maximum distance by taking "at random" (first atom acctually) and checking the futherest atom, and its furthest atom
    furthest_atom = np.argmax(distance_matrix([cage_coords[0]], cage_coords)[0])

    avg_m_m_dist = np.max(distance_matrix([cage_coords[furthest_atom]], cage_coords)[0])

    centroid = np.mean(cage_coords, axis=0)  # centroid

    cage_coords = np.array([coord - centroid for coord in cage_coords])

    # For a distance from the origin (the cage centroid) calculate the
    # largest sphere possible without hitting atoms
    opt_theta_phi, opt_r = np.zeros(2), 0.0
    for r in np.linspace(0.0, avg_m_m_dist + max_dist_from_metals, 30):
        if basinh:
            opt = basinhopping(get_max_sphere_negative_radius,
                               x0=opt_theta_phi, stepsize=1.0, niter=5,
                               minimizer_kwargs={'args': (r, cage_coords),
                                                 'method': 'BFGS'})
        else:
            opt = minimize(get_max_sphere_negative_radius,
                           x0=opt_theta_phi,
                           args=(r, cage_coords),
                           method='BFGS')

        opt_theta_phi = opt.x

        # This is the correct way round because we want the largest sphere
        #  that CAN escape
        if -opt.fun < max_sphere_escape_r:
            max_sphere_escape_r = -opt.fun
            opt_r = r

    # Get the atom id that the max escape sphere hits into
    sphere_point = spherical_to_cart(r=opt_r, theta=opt_theta_phi[0], phi=opt_theta_phi[1])
    atom_id = np.argmin([np.linalg.norm(coord - sphere_point) for coord in cage_coords])

    radius = max_sphere_escape_r - vdw_radii[atom_names[atom_id]]
    # logger.info(f'Radius of largest sphere that can escape from the '
    #            f'cavity = {radius}')

    return radius


def get_max_sphere_negative_radius(theta_and_phi, r, cage_coords):
    """
    Get the maximum sphere radius that is possible at a point defined by the
    spherical polar coordinates theta, phi and r. This amounts to finding the
    minimum pairwise distance between the point and the rest of the cage. The
    negative radius is returned as it will be fed into scipy.optmize.minimise

    :param theta_and_phi: (list(float))
    :param r: (float)
    :param cage_coords: (np.ndarray) n_atoms x 3
    :return: (float)
    """

    theta, phi = theta_and_phi
    # Convert the point in spherical polars to Cartesian so the distances to
    # the rest of the cage can be calculated
    # needs to be a 1 x 3 matrix to use cdist
    point = np.array([spherical_to_cart(r=r, theta=theta, phi=phi)])

    return -np.min(cdist(point, cage_coords))


def spherical_to_cart(r, theta, phi):
    return np.array([r * np.cos(theta) * np.sin(phi),
                     r * np.sin(theta) * np.sin(phi),
                     r * np.cos(phi)])
