import numpy


def compute_focus_squared(au, pu, av, pv):
    """Computes the square of the focus of the ellipse from the
    amplitide and phase of x and y component of velocity. The focus f
    is the distance between the centre and each of the foci of the
    ellipse. If a and b are the major and minor radii, then this returns:
        f**2 = a**2 - b**2"""
    return numpy.sqrt(au**4+av**4+2*au**2*av**2*numpy.cos(2*(pu-pv)))


def tidal_ellipse_parameters(au, pu, av, pv):
    """Computes the major and minor radii, direction of, and phase wrt
    principle axis from the amplitude and phase of the u and v components of velocity.
    Accepts numpy arrays to compute parameters for multiple consituents
    at once."""
    f2 = compute_focus_squared(au, pu, av, pv)
    a = numpy.sqrt(0.5*(au**2+av**2+f2))
    b = numpy.sqrt(0.5*(au**2+av**2-f2))

    # from Pugh (A3:10)
    gc = numpy.arctan2(au*numpy.sin(pu)+av*numpy.cos(pv), au*numpy.cos(pu)-av*numpy.sin(pv))
    gac = numpy.arctan2(-au*numpy.sin(pu)+av*numpy.cos(pv), au*numpy.cos(pu)+av*numpy.sin(pv))
    theta = 0.5 * (gc+gac)
    g = 0.5 * (gc-gac) % (2*numpy.pi)
    return a, b, theta, g
