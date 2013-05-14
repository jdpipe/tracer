# Implements a circular conical surface

import numpy as N
from quadric import QuadricGM
import pivy.coin as coin

class InfiniteCone(QuadricGM):
    # TODO reimplement with support for z0, location of apex
    """
    Implements the geometry of an infinite circular conical surface. That
    means the sloping side-walls of the cone, doesn't include the base.
    """
    def __init__(self, c = 1.):
        """          
        Arguments: 
        c - constant in x**2 + y**2 = (c*z)**2, c = r/h (base radius/height)
        
        Private attributes:                                                                  
        c - constant in x**2 + y**2 = (c*z)**2
        """ 
        QuadricGM.__init__(self)
        self.c = c

    def _normals(self, hits, directs):
        """
        Finds the normal to the cone in a bunch of intersection points, by
        taking the derivative and rotating it. Used internally by quadric.
        
        Arguments:
        hits - the coordinates of intersections, as an n by 3 array.
        directs - directions of the corresponding rays, n by 3 array.
        """
        hit = N.dot(N.linalg.inv(self._working_frame), 
            N.vstack((hits.T, N.ones(hits.shape[0]))))
        dir_loc = N.dot(self._working_frame[:3,:3].T, directs.T)
        partial_x = 2*hit[0]
        partial_y = 2*hit[1]
        partial_z = -2*self.c**2*hit[2]

        #FIXME check what happens with inside and outside reflections?

        local_normal = N.vstack((partial_x, partial_y, partial_z))
        local_unit = local_normal/N.sqrt(N.sum(local_normal**2, axis=0))
        down = N.sum(dir_loc * local_unit, axis=0) > 0
        local_unit[:,down] *= -1
        normals = N.dot(self._working_frame[:3,:3], local_unit)
        
        return normals  
    
    def get_ABC(self, ray_bundle):
        # FIXME check this with a test case.
        # FIXME TODO XXX the calculation of intersections goes wrong if A,B,C are all negated
        """
        TODO Determines the variables forming the relevant quadric equation, [1]
        """
        # Transform the the direction and position of the rays temporarily into the
        # frame of the paraboloid for calculations
        d = N.dot(self._working_frame[:3,:3].T, ray_bundle.get_directions())
        v = N.dot(N.linalg.inv(self._working_frame), 
            N.vstack((ray_bundle.get_vertices(), N.ones(d.shape[1]))))[:3]
        
        A = d[0]**2 + d[1]**2 - (self.c*d[2])**2
        B = 2*(v[0]*d[0] + v[1]*d[1] - self.c**2*v[2]*d[2])
        C = v[0]**2 + v[1]**2 - (self.c*v[2])**2

        mult = 1. - 2.*(A<0)
        A *= mult
        B *= mult
        C *= mult

        print "mult",mult

        return A, B, C


class Cone(InfiniteCone):
    """
    Implements a finite cone. Parameters are r (base radius) and h (cone
    height). The cone is aligned with the (positive) z axis, and the apex of the
    cone is at the origin.
    """
    def __init__(self, r, h):
        if h < 0 or r < 0:
            raise AttributeError
        self.h = h
        self.r = r
        InfiniteCone.__init__(self, r/h)
    
    def _select_coords(self, coords, prm):
        """
        Refinement of QuadricGM._select_coords; we want to choose the correct
        intersection point for a set of rays and our *truncated* quadric cone
        surface.

        Arguments:
        coords - a 2 by 3 by n array whose each column is the global coordinates
            of one intersection point of a ray with the sphere.
        prm - a 2 by n array (CHECK THIS) giving the parametric location on the
            ray where the intersection occurs.

        Returns:
        The index of the selected intersection, or None if neither will do.
        """
        select = N.empty(prm.shape[1])
        select.fill(N.nan)

        print "prm",prm        
        print "A**-1 =",N.linalg.inv(self._working_frame) 
        print "A trim",N.linalg.inv(self._working_frame)[None,2,:,None]       
        print "coords",coords
        print "concat",N.concatenate((coords, N.ones((2,1,coords.shape[-1]))), axis=1)

        # FIXME CHECK THIS...
        height = N.sum(N.linalg.inv(self._working_frame)[None,2,:,None] * \
            N.concatenate((coords, N.ones((2,1,coords.shape[-1]))), axis=1), axis=1)

        print "height=\n",height
        #print "self.h =",self.h

        inside = (height >= 0) & (height <= self.h)

        print "inside=\n",inside

        positive = prm > 0

        # Assumption here seems to be that the first-given of each 'hit' is the nearer one -- JP
        hitting = inside & positive
        select[N.logical_and(*hitting)] = 0
        print "*hitting",N.logical_xor(*hitting)
        one_hitting = N.logical_xor(*hitting)
        print "one_hitting=\n",one_hitting
        select[one_hitting] = N.nonzero(hitting.T[one_hitting,:])[1]
        print "select",select

        return select
    
    def get_scene_graph(self, resolution = 1.):
        n = coin.SoSeparator()
        ro = coin.SoRotationXYZ()
        ro.axis = coin.SoRotationXYZ.X
        ro.angle = -N.pi / 2.
        n.addChild(ro)
        tr = coin.SoTranslation()
        tr.translation = (0,-self.h/2.,0)
        n.addChild(tr)
        co = coin.SoCone()
        co.bottomRadius = self.r
        # NOTE FOR LATER: to render the inside of the cone, we can set co.height < 0
        co.height = self.h
        co.parts = co.SIDES
        n.addChild(co)
        return n


class ConicalFrustum(InfiniteCone):
    """
    Implements a conical frustum, with base radius r, (virtual) apex at the
    original, small face at z=h1, large face at z=h2, axis on the z-axis.
    """
    def __init__(self, r, h1, h2):
        if r <= 0 or h1 <= 0 or h2 <= 0:
            raise AttributeError
        if h2 <= h1:
            raise AttributeError
        self.h1 = h1
        self.h2 = h2
        c = r / h2
        InfiniteCone.__init__(self, c)
    
    def _select_coords(self, coords, prm):
        """
        Refinement of QuadricGM._select_coords; we want to choose the correct
        intersection point for a set of rays and our *truncated* quadric cone
        surface.

        Arguments:
        coords - a 2 by 3 by n array whose each column is the global coordinates
            of one intersection point of a ray with the sphere.
        prm - a 2 by n array (CHECK THIS) giving the parametric location on the
            ray where the intersection occurs.

        Returns:
        The index of the selected intersection, or None if neither will do.
        """
        select = N.empty(prm.shape[1])
        select.fill(N.nan)

        print "prm",prm        
        print "A**-1 =",N.linalg.inv(self._working_frame) 
        print "A trim",N.linalg.inv(self._working_frame)[None,2,:,None]       
        print "coords",coords
        print "concat",N.concatenate((coords, N.ones((2,1,coords.shape[-1]))), axis=1)

        # FIXME CHECK THIS...
        height = N.sum(N.linalg.inv(self._working_frame)[None,2,:,None] * \
            N.concatenate((coords, N.ones((2,1,coords.shape[-1]))), axis=1), axis=1)

        print "height=\n",height
        #print "self.h =",self.h

        inside = (self.h1 <= height) & (height <= self.h2)

        print "inside=\n",inside

        positive = prm > 0

        # Assumption here seems to be that the first-given of each 'hit' is the nearer one -- JP
        hitting = inside & positive
        select[N.logical_and(*hitting)] = 0
        print "*hitting",N.logical_xor(*hitting)
        one_hitting = N.logical_xor(*hitting)
        print "one_hitting=\n",one_hitting
        select[one_hitting] = N.nonzero(hitting.T[one_hitting,:])[1]
        print "select",select

        return select

    def mesh(self, resolution):
        """
        Represent the surface as a mesh in local coordinates. Uses polar
        bins, i.e. the points are equally distributed by angle and radius,
        not by x,y.
        
        Arguments:
        resolution - in points per unit length (so the number of points 
            returned is O(A*resolution**2) for area A)
        
        Returns:
        x, y, z - each a 2D array holding in its (i,j) cell the x, y, and z
            coordinate (respectively) of point (i,j) in the mesh.
        """

        # Generate a circular-edge mesh using polar coordinates.    
        r1 = self.c * self.h1
        r2 = self.c * self.h2
        rs = N.r_[r1,r2]
        # Make the circumferential points at the requested resolution.
        ang_end = 2*N.pi + 1./(r2*resolution)
        angs = N.r_[0:ang_end:1./(r2*resolution)]

        x = N.outer(rs, N.cos(angs))
        y = N.outer(rs, N.sin(angs))
        z = 1/self.c * N.sqrt(x**2 + y**2)
        
        return x, y, z
    
    #use get_scene_graph from GeometryManager, since we've implemented mesh.

# vim: et:ts=4
