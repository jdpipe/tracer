
import numpy as np
import pylab

from tracer.surface import *
from tracer.cone import *
from tracer.sphere_surface import *
from tracer.paraboloid import *
from tracer.flat_surface import *
from tracer.assembly import *
from tracer.optics_callables import *
from tracer.object import *
from tracer.spatial_geometry import *
from tracer.sources import *
from tracer.tracer_engine import *

import pivy.coin as coin
SOGUI_BINDING="SoQt"
from pivy.sogui import *

A = Assembly()

# Paraboloidal dish...
alpha = 0.04 # dish absorptivity
d = 22 # dish diameter
f = 13.4 # dish focal length
P = AssembledObject(surfs=[Surface(ParabolicDishGM(d, f), Reflective(alpha))], transform=rotx(0*N.pi/6))
#A.add_object(P)


# a beautiful cone...
CO = AssembledObject(surfs=[Surface(ConicalFrustum(11.,8.,11.), Reflective(alpha))], transform=rotx(N.pi/3))
A.add_object(CO)


r = 0.1
for z in range(14):
	tr = translate(0,0,z)
	#print "translation",y,"=",tr
	S = AssembledObject(surfs=[Surface(SphericalGM(r), Reflective(alpha))], transform=tr)
	#A.add_object(S)

# A target surface
rw = 1.
rh = 1.
R = AssembledObject(surfs=[Surface(RectPlateGM(rw,rh), LambertianReceiver(alpha))], transform=translate(0,0,f))
#A.add_object(R)

# do a raytrace
cr = np.array([[0,0,2*f]]).T
dr = np.array([0,0,-1])
ar = 0*5e-3 # radians, sun rays angular range (what's the correct value?)
G = 1000. # W/m2 solar flux
#TODO code in the Buie sunshape instead of a pillbox
src = solar_disk_bundle(5000, cr, dr, d*1., ar, G)

engine = TracerEngine(A)
engine.ray_tracer(src, 100, 0.001)

# draw the rays... adapted from ../tracer/mayavi/scene_view.py
def show_rays(engine, escaping_len=20.):
    tree = engine.tree
    no = coin.SoSeparator()
    print tree.num_bunds()
    
    # loop through the reflection sequences?
    co = []
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
                co.append((c1[0],c1[1],c1[2]))
                co.append((c2[0],c2[1],c2[2]))
                #endpoints = N.c_[sv[:,ray], ev[:,first_child]]
            else:
                l = escaping_len
                if level == 0:
                    l = 0.1
                # Escaping ray.
                c1 = sv[:,ray]
                c2 = sv[:,ray] + sd[:,ray]*l
                co.append((c1[0],c1[1],c1[2]))
                co.append((c2[0],c2[1],c2[2]))

            #lines.append(scene.mlab.plot3d(*endpoints, tube_radius=None))
    #print "num of vertices",len(co)
    #print co

    ma = coin.SoMaterial()
    ma.diffuseColor = (1,1,0.5)
    no.addChild(ma)

    ds = coin.SoDrawStyle()
    ds.style = ds.LINES
    ds.lineWidth = 2
    no.addChild(ds)

    coor = coin.SoCoordinate3()
    coor.point.setValues(0, len(co), co)
    #print [x.getValue() for x in coor.point.getValues()]
    no.addChild(coor)

    ls = coin.SoLineSet()
    ind = [2] * (len(co)/2)
    #print "ind =",ind
    ls.numVertices.setValues(0, len(ind), ind)
    no.addChild(ls)
    
    return no

# render the scene with Pivy
r = coin.SoSeparator()

# create axes
def axis_labels(length=10):
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

r.addChild(axis_labels())

r.addChild(show_rays(engine))

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
