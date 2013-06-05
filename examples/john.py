
import numpy as np
import pylab
np.set_printoptions(linewidth=140)

from tracer.surface import *
from tracer.cone import *
from tracer.cylinder import *
from tracer.sphere_surface import *
from tracer.paraboloid import *
from tracer.flat_surface import *
from tracer.assembly import *
from tracer.optics_callables import *
from tracer.object import *
from tracer.spatial_geometry import *
from tracer.sources import *
from tracer.tracer_engine import *
from tracer.models.one_sided_mirror import *
import types

import pivy.coin as coin
SOGUI_BINDING="SoQt"
from pivy.sogui import *


A = Assembly()

# Paraboloidal dish...
alpha = 0 # dish absorptivity
d = 22 # dish diameter
f = 13.4 # dish focal length
tr0 = translate(z=-f)
P = AssembledObject(surfs=[Surface(ParabolicDishGM(d, f), Reflective(alpha))], transform=tr0)
A.add_object(P)

# Cylinder
#alpha3 = 0.5
#tr = N.dot(rotx(N.pi), translate(z=-0.3))
#CY1 = AssembledObject(surfs=[Surface(FiniteCylinder(0.8, 0.8), Reflective(alpha3))], transform=tr)
#A.add_object(CY1)

# A beautiful thick, double frustum with 2 different sides. Still leaks.

alpha2 = 0.5
tr = N.dot(rotx(N.pi), translate(z=-0.3))
width = 1e-10 # receiver thickness at frustii junction

CO1front = Surface(ConicalFrustum(z1=-0.5,r1=0.01,z2=0,r2=0.7), ReflectiveReceiver(alpha2))
CO1back = Surface(ConicalFrustum(z1=-0.5,r1=0.01,z2=0,r2=0.7+width), Reflective(alpha2))
CO1 = AssembledObject(surfs=[CO1front,CO1back], transform=tr)

CO1.surfaces_for_next_iteration = types.MethodType(surfaces_for_next_iteration, CO1, CO1.__class__)

CO2front = Surface(ConicalFrustum(z1=0,r1=0.7,z2=0.3,r2=0.4), ReflectiveReceiver(alpha2))
CO2back = Surface(ConicalFrustum(z1=0,r1=0.7+width,z2=0.3,r2=0.4), ReflectiveReceiver(alpha2))
CO2 = AssembledObject(surfs=[CO2front,CO2back], transform=tr)

CO2.surfaces_for_next_iteration = types.MethodType(surfaces_for_next_iteration, CO2, CO2.__class__)

A.add_object(CO1)
A.add_object(CO2)

# A target surface
#rw = 1.
#rh = 1.
#R = AssembledObject(surfs=[Surface(RectPlateGM(rw,rh), LambertianReceiver(alpha))], transform=translate(0,0,f))
#A.add_object(R)

# do a raytrace
cr = np.array([[0,0,2*f]]).T
dr = np.array([0,0,-1])
ar = 5e-3 # radians, sun rays angular range (what's the correct value?)
G = 1000. # W/m2 solar flux
nrays = 10000 # number of rays escaping the source
#TODO code in the Buie sunshape instead of a pillbox
src = solar_disk_bundle(nrays, cr, dr, d*1., ar, G)

engine = TracerEngine(A)
engine.ray_tracer(src, 100, 0.001)


def show_rays(engine, escaping_len=20., highlight_level=None):
    """
    Function to draw the rays to a Coin3D scenegraph.
    Adapted from ../tracer/mayavi/scene_view.py
    """
    tree = engine.tree
    no = coin.SoSeparator()
    print tree.num_bunds()
    
    # loop through the reflection sequences?
    co = [] # regular lines
    co_h = [] # highlighted lines
    hist = {} # ray histories, for highlighted rays
    for level in xrange(tree.num_bunds()):
        print "bundle",level
        start_rays = tree[level]
        print "start_rays",start_rays.get_num_rays()
        sv = start_rays.get_vertices()
        sd = start_rays.get_directions()
        se = start_rays.get_energy()
        
        if level == tree.num_bunds() - 1:
            parents = []
        else:
            end_rays = tree[level + 1]
            ev = end_rays.get_vertices()
            parents = end_rays.get_parents()

        # loop through individual rays in this bundle
        for ray in xrange(start_rays.get_num_rays()):
            if se[ray] == 0:
                # ignore rays with no starting energy
                continue
            
            if ray in parents:
                # Has a hit on another surface
                first_child = N.where(ray == parents)[0][0]
                c1 = sv[:,ray]
                c2 = ev[:,first_child]
                #endpoints = N.c_[sv[:,ray], ev[:,first_child]]
            else:
                l = escaping_len
                if level == 0:
                    l = 0.1
                # Escaping ray.
                c1 = sv[:,ray]
                c2 = sv[:,ray] + sd[:,ray]*l
            if level == highlight_level:
                hist[ray] = tree.ray_history(ray,level)
                co_h += [(c1[0],c1[1],c1[2]), (c2[0],c2[1],c2[2])]
            else:
                co += [(c1[0],c1[1],c1[2]), (c2[0],c2[1],c2[2])]

    #print "num of vertices",len(co)
    #print "co_h",len(co_h)
    #print "co",len(co)
    #print "hist",hist

    #print "ray history:"
    #for i in hist:
    #    detailed_ray_history(engine,hist[i])

    # normal rays
    def plot_rays_color(co, color=(1,1,0.5)):
        """
        Add ray set of line color `color` to scenegraph `node`. Rays `co`
        should be stored as sequences of 3-vector pairs in a list, eg
        [(x1,y1,z1),(x2,y2,z2),...]
        """
        no1 = coin.SoSeparator()

        ma1 = coin.SoMaterial()
        ma1.diffuseColor = color
        no1.addChild(ma1)

        ds = coin.SoDrawStyle()
        ds.style = ds.LINES
        ds.lineWidth = 2
        no1.addChild(ds)

        coor = coin.SoCoordinate3()
        coor.point.setValues(0, len(co), co)
        no1.addChild(coor)

        ls = coin.SoLineSet()
        ind = [2] * (len(co)/2)
        ls.numVertices.setValues(0, len(ind), ind)
        no1.addChild(ls)

        return no1;
    
    no.addChild(plot_rays_color(co))
    no.addChild(plot_rays_color(co_h, color=(1,1,1)))

    return no

def detailed_ray_history(engine,seq):
    """
    Print a detailed ray history for a particular ray, being a tuple returned 
    from TraceTree.ray_history.
    """
    tree = engine.tree
    n = len(seq)
    print seq
    for i in range(n):
        bund = tree[i]
        ray = seq[(n-1)-i]
        print "...bund",i,"ray",ray
        print "...from",bund.get_vertices()[:,ray]
        print "...direction",bund.get_directions()[:,ray]

def axis_labels(length=10):
    """
    Create axis arrows/labels for addition to a Coin3D scenegraph.
    """
    r = coin.SoSeparator()
    st = coin.SoDrawStyle()
    r.addChild(st)
    st.lineWidth=3
    data = {'x':(1,0,0), 'y':(0,1,0), 'z':(0,0,1)}
    for k in data:
        vx,vy,vz = data[k]
        vec = (length*vx, length*vy, length*vz)
        s1 = coin.SoSeparator()
        tr1 = coin.SoTranslation()
        tr1.translation = vec
        s1.addChild(tr1)
        la = coin.SoLabel()
        la.label = k
        s1.addChild(la)
        r.addChild(s1)
        ma = coin.SoMaterial()
        ma.diffuseColor = data[k]
        r.addChild(ma)
        co = coin.SoCoordinate3()
        co.point.setValues(0,2,[(0,0,0),vec])
        r.addChild(co)
        ls = coin.SoLineSet()
        ls.numVertices.setValues(0,1,[2])
        r.addChild(ls)
    return r


# render the scene with Pivy
r = coin.SoSeparator()
r.addChild(axis_labels())
r.addChild(show_rays(engine, highlight_level=2))
r.addChild(A.get_scene_graph())

win = SoGui.init("hello")
viewer = SoGuiExaminerViewer(win)
viewer.setSceneGraph(r)
viewer.setTitle("Examiner Viewer")
viewer.viewAll()
viewer.show()

SoGui.show(win)
SoGui.mainLoop()


# vim: et:ts=4
