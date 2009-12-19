# An assembly modeling the parabolic dish for building-integrated installations,
# as developped in Tel Aviv University's Faculty of Engineering.
#
# References:
# [1] Kribus A., et al, A miniature concentrating photovoltaic and thermal 
#     system, Energy Conversion and Management, Volume 47, Issue 20, December
#     2006, Pages 3582-3590, DOI: 10.1016/j.enconman.2006.01.013.

import numpy as N

from .. import spatial_geometry as sp
from .. import optics_callables as opt
from ..assembly import Assembly
from ..surface import Surface
from ..paraboloid import ParabolicDishGM
from ..object import AssembledObject

from .one_sided_mirror import one_sided_receiver
from .homogenizer import rect_homogenizer

class MiniDish(Assembly):
    def __init__(self, diameter, focal_length, dish_opt_eff,\
        receiver_pos, receiver_side, homogenizer_depth, homog_opt_eff):
        """
        Arguments:
        diameter, focal_length - of the parabolic dish
        dish_opt_eff - the optical efficiency of the dish
        receiver_pos - the distance along the optical axis from the dish to the
            receiver's end surface - the PV panel (should be about the focal 
            length)
        receiver_side - the receiver is square, with this side length.
        homogenizer_depth - the homogenizer has base dimensions to fit the PV
            square, and this height.
        homog_opt_eff - the optical efficiency of each mirror in the homogenizer
        """
        self._side = receiver_side
        self._rec, rec_obj = one_sided_receiver(self._side, self._side)
        receiver_frame = N.dot(sp.translate(0, 0, receiver_pos), sp.rotx(N.pi))
        rec_obj.set_transform(receiver_frame)
        
        homogenizer = rect_homogenizer(self._side, self._side, \
            homogenizer_depth, homog_opt_eff)
        homogenizer.set_transform(receiver_frame)
        
        dish_surf = Surface(ParabolicDishGM(diameter, focal_length), 
            opt.gen_reflective(1 - dish_opt_eff))
        dish = AssembledObject(surfs=[dish_surf])
        
        Assembly.__init__(self, objects=[rec_obj, dish], subassemblies=[homogenizer])
    
    def get_receiver_surf(self):
        """for anyone wishing to directly access the receiver"""
        return self._rec
    
    def histogram_hits(self, bins=50):
        """
        Generates a 2D histogram of energy absorbed at the receiver surface,
        assuming a trace has been run using this assembly.
        
        Arguments:
        bins - hom many bins per axis to use (default 50)
        
        Returns:
        H - a 2D array with the energy falling on each bean, x axis along
            the first dimension, y along second
        xbins, ybins - the edges of the bins (so one point more than the number
            of beans for each axis)
        
        See Also:
        numpy.histogram2D()
        """
        energy, pts = self._rec.get_optics_manager().get_all_hits()
        x, y = self._rec.get_geometry_manager().global_to_local(pts)[:2]
        rng = self._side/2.
        
        H, xbins, ybins = N.histogram2d(x, y, bins, \
            range=([-rng,rng], [-rng,rng]), weights=energy)
        return H, xbins, ybins