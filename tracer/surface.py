# -*- coding: utf-8 -*-
# Define some basic surfaces for use with the ray tracer. From this minimal 
# hierarchy other surfaces should be derived to implement actual geometric
# operations.
#
# References:
# [1] John J. Craig, Introduction to Robotics, 3rd ed., 2005. 

import numpy as N
from has_frame import HasFrame
import pivy.coin as coin

class Surface(HasFrame):
    """
    Defines the base of surfaces that interact with rays.
    """
    def __init__(self, geometry, optics, location=None, rotation=None):
        """
        Arguments:
        geometry - a GeometryManager object responsible for finding ray 
            intersections with the surface.
        optics - a callable that gets the geometry manageri, bundle and
            selector, and returns the outgoing ray bundle generated by the
            geometry and bundle.
        location, rotation - passed directly to the HasFrame constructor.
        """
        HasFrame.__init__(self, location, rotation)
        self._geom = geometry
        self._opt = optics
        
    def get_optics_manager(self):
        """
        Returns the optics-manager callable. May be useful for introspection.
        Note that it is a read-only attribute.
        """
        return self._opt
    
    def get_geometry_manager(self):
        """
        Returns the geometry-manager instance. May be useful for introspection.
        Note that it is a read-only attribute.
        """
        return self._geom
    
    def register_incoming(self, ray_bundle):
        """
        Records the incoming ray bundle, and uses the geometry manager to
        return the parametric positions of intersection with the surface along
        the ray.
        
        Arguments:
        ray_bundle - a RayBundle object with at-least its vertices and
            directions specified.
        
        Returns
        A 1D array with the parametric position of intersection along each of
            the rays. Rays that missed the surface return +infinity.
        """
        self._current_bundle = ray_bundle
        return self._geom.find_intersections(self._temp_frame, ray_bundle)
    
    def select_rays(self, idxs):
        """
        Informs the geometry manager that only the specified rays are to be
        used henceforth.
        
        Arguments:
        idxs - an array with indices into the last registered ray bundle,
            marking rays that will be used.
        """
        self._selected = idxs
        self._geom.select_rays(idxs)
    
    def get_outgoing(self):
        """
        Generates a new ray bundle, which is the reflections/refractions of the
        user-selected rays out of the incoming ray-bundle that was previously
        registered.
        
        Returns: 
        a RayBundle object with the new bundle, with vertices on the surface
            and directions according to optics laws.
        """
        return self._opt(self._geom, self._current_bundle, self._selected)
    
    def done(self):
        """
        When this is called, the surface will no longer be queried on the
        results of the latest trace iteration, so it can discard internal
        data to relieve memory pressure.
        """
        if hasattr(self, '_current_bundle'):
            del self._current_bundle
        self._geom.done()
        
    def global_to_local(self, points):
        """
        Transform a set of points in the global coordinates back into the frame
        used during tracing.
        
        Arguments:
        points - a 3 x n array for n 3D points
        
        returns:
        local - a 3 x n array with the respective points in local coordinates.
        """
        return N.dot(N.linalg.inv(self._temp_frame), 
            N.vstack((points, N.ones(points.shape[1]))))
    
    def mesh(self, resolution):
        """
        Represent the surface as a mesh in global coordinates.
        
        Arguments:
        resolution - in points per unit length (so the number of points 
            returned is O(A*resolution**2) for area A)
        
        Returns:
        x, y, z - each a 2D array holding in its (i,j) cell the x, y, and z
            coordinate (respectively) of point (i,j) in the mesh.
        """
        # The geometry manager has the local-coordinates mesh.
        x, y, z = self._geom.mesh(resolution)
        local = N.array((x, y, z, N.ones_like(x)))
        glob = N.tensordot(self._temp_frame, local, axes=([1], [0]))
        
        return glob[:3]

    def get_scene_graph(self,resolution=None):
        """
        Any object that provides a nice QuadMesh from the previous code should be able to render in Coin3D with with the following...
        """
    
        n = self.get_scene_graph_transform()

        o = self.get_optics_manager()
        if hasattr(o,'get_all_hits'):
            e, h = o.get_all_hits()
            # plot the histogram into the scenegraph
            """
                How to do this?
                * is the surface has a mesh, we could use IndexedFaceSet with
                  colouring set for each face.
                * if it's a flat surface, we can put a texture on the face
                * but we want to do receivers with concave surfaces
                * how to 'unwrap' the receiver won't always be obvious either
                * we could use axisymmetric assumptions to simplify this
                * maybe start with the flat surface case???
            """
            n.addChild(self._geom.get_scene_graph(resolution))
        else:         
            n.addChild(self._geom.get_scene_graph(resolution))

        return n


# vim: et:ts=4
