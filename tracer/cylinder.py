"""
Geometry managers based on a cylinder along the Z axis.

References:
.. [1] http://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter4.htm
"""

from .quadric import QuadricGM
import numpy as N

class InfiniteCylinder(QuadricGM):
    """
    A cylindrical surface infinitely long on the Z axis.
    """
    def __init__(self, diameter):
        self._R = diameter/2.
        QuadricGM.__init__(self)
    
    def get_ABC(self, bundle):
        """
        Finds the coefficients of the quadratic equation for the intersection
        of a ray and the cylinder. See [1]_.
        
        Arguments:
        bundle - a RayBundle instance with rays for which to get the
            coefficients.
        
        Returns:
        A, B, C - satisfying A*t**2 + B*t + C = 0 such that the intersection
            points are at distance t from each ray's vertex.
        """
        # Transform the the direction and position of the rays temporarily into the
        # frame of the paraboloid for calculations
        d = N.dot(self._working_frame[:3,:3].T, bundle.get_directions())
        proj = N.round(N.linalg.inv(self._working_frame), decimals=9)
        v = N.dot(proj, N.vstack((bundle.get_vertices(), N.ones(d.shape[1]))))[:3]
        
        A = N.sum(d[:2]**2, axis=0)
        B = 2*N.sum(d[:2]*v[:2], axis=0)
        C = N.sum(v[:2]**2, axis=0) - self._R**2
        
        return A, B, C
    
    def _normals(self, verts, dirs):
        # Move to local coordinates
        proj = N.round(N.linalg.inv(self._working_frame), decimals=9)
        hit = N.dot(proj, N.vstack((verts.T, N.ones(verts.shape[0]))))
        dir_loc = N.dot(self._working_frame[:3,:3].T, dirs.T)
        
        # The local normal is made from the X,Y components of the vertex:
        local_norm = N.vstack((hit[:2], N.zeros(hit.shape[1])))
        local_norm /= N.sqrt(N.sum(hit[:2]**2, axis=0))
        
        # Choose whether the normal is inside or outside:
        local_norm[:, N.sum(local_norm[:2] * dir_loc[:2], axis=0) > 0] *= -1
        
        # Back to global coordinates:
        return N.dot(self._working_frame[:3,:3], local_norm)


class FiniteCylinder(InfiniteCylinder):
    """
    This geometry manager represents a cylinder with a given height, centered
    on its origin, (so that the top Z is at height/2).
    """
    def __init__(self, diameter, height):
        self._half_h = height/2.
        InfiniteCylinder.__init__(self, diameter)
    
    def _select_coords(self, coords, prm):
        """
        Choose between two intersection points on a quadric surface.
        This implementation extends QuadricGM's behaviour by not choosing
        intersections higher or lower than half the cylinder height.

        Arguments:
        coords - a 2 by 3 by n array whose each column is the global coordinates
            of one intersection point of a ray with the sphere.
        prm - the corresponding parametric location on the ray where the
            intersection occurs.

        Returns:
        The index of the selected intersection, or None if neither will do.
        """
        select = N.empty(prm.shape[1])
        select.fill(N.nan)
        
        proj = N.round(N.linalg.inv(self._working_frame), decimals=9)
        height = N.sum(proj[None,2,:,None] * \
            N.concatenate((coords, N.ones((2,1,coords.shape[-1]))), axis=1), axis=1)
        inside = (abs(height) <= self._half_h)
        positive = prm > 1e-10
        
        hitting = inside & positive
        select[N.logical_and(*hitting)] = 1
        one_hitting = N.logical_xor(*hitting)
        select[one_hitting] = N.nonzero(hitting.T[one_hitting,:])[1]

        return select
    
    def mesh(self, resolution = None):
        """
        Represent the surface as a mesh in local coordinates. Uses cylindrical
        bins, i.e. the points are equally distributed by angle and axial 
        location, not by x,y.

        Arguments:
        resolution - in points per unit length (so the number of points 
            returned is O(height*pi*diameter*resolution**2))

        Returns:
        x, y, z - each a 2D array holding in its (i,j) cell the x, y, and z
            coordinate (respectively) of point (i,j) in the mesh.
        """
        if resolution is None:
            angres = 2*N.pi / 40
        else:
            angres = 2*N.pi * (resolution / (2*N.pi*self._R))

        # note: current mesh has no detail in the axial direction, just start/end points -- JP
        h = N.r_[[-self._half_h, self._half_h]]
        ang_end = 2*N.pi
        angs = N.r_[0:ang_end+angres:angres]

        x = N.tile(self._R*N.cos(angs), (len(h), 1))
        y = N.tile(self._R*N.sin(angs), (len(h), 1))
        z = N.tile(h[:,None], (1, len(angs)))

        return x, y, z

# vim: et:ts=4
