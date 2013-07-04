'''
__________________________________________________________________________________________________
************ 4 Parameters cavity receiver reflective losses calculation program *************
__________________________________________________________________________________________________
'''
import numpy as N
import pylab
import time
N.set_printoptions(linewidth=140)

from tracer.surface import *
from tracer.quadric import *
from tracer.paraboloid import *
from tracer.cone import *
from tracer.assembly import *
from tracer.optics_callables import *
from tracer.object import *
from tracer.spatial_geometry import *
from tracer.sources import *
from tracer.tracer_engine import *

from emissive_losses.emissive_losses import *

import pivy.coin as coin
SOGUI_BINDING="SoQt"
from pivy.sogui import *

'''
__________________________________________________________________________________________________
Study definition
__________________________________________________________________________________________________
'''
# timer start
t0 = time.clock()

#-------------------------Simulation parameters--------------------
# Receiver: 4 parameters cavity geometry ---------------------------------
apertureDiameter = 0.3 # (m)
depth = 0.5 # (m)
frustumAngle = 50.*N.pi/180 # (rad) frustum half-angle
coneAngle = 120.*N.pi/180 # (rad) cone half angle
# Receiver: optics ------------------------------------------------
absReceiver = 0.87 # Receiver absorptivity
emsReceiver = absReceiver # Receiver emissivity
# Receiver: Temperature boundary condition -------------------------------
tempAmbiant = 20+273.15 # Ambiant temperature in K
tempReceiver = 650+273.15 # Receiver internal wall temperature in K
# Receiver number of elements (for emissive losses) -------------------------
n = 20
# Collector: paraboloidal dish geometry ----------------------------------
dishDiameter = 22. # (m) dish diameter
dishFocus = 13.4 # (m) dish focal length
# Collector: optics -----------------------------------------------
absDish = 0.05 # mirror absorptivity
sigma = 4e-3 # (rad) local x and y axis gaussian shape error distribution parameter

# Surface definition parameters
d1 = (depth*N.tan(coneAngle)-apertureDiameter/2.)/(N.tan(coneAngle)+N.tan(frustumAngle))
d2 = depth-d1
rc = d2*N.tan(coneAngle)

'''
_____________________________________________________________________________________________
Assembly construction
_____________________________________________________________________________________________
'''
Sim = Assembly()
Receiver = Assembly()

# Paraboloid dish
# - Positioning
trd = translate(z=-dishFocus)
# - Construction
DISH = AssembledObject(surfs=[Surface(ParabolicDishGM(dishDiameter, dishFocus), RealReflectiveReceiver(absDish, sigma))], transform=trd)
Sim.add_object(DISH)
# 4 Parameters cavity receiver
# - Positioning of the receiver assembly
trr = None
# - Construction of the receiver surfaces
# -- Frustum surface
# --- Positioning of the frustum surface in the receiver frame
trf = None
# --- Construction of the frustum surface in the receiver frame
FRU = AssembledObject(surfs=[Surface(ConicalFrustum(z1=0., r1=apertureDiameter/2., z2=d1, r2=rc), LambertianReceiver(absReceiver))], transform=trf)
Receiver.add_object(FRU)
# --- Construction of the frustum external reflective surface in the receiver frame
FRUe = AssembledObject(surfs=[Surface(ConicalFrustum(z1=0., r1=apertureDiameter/2., z2=d1, r2=rc+1e-6), ReflectiveReceiver(0))], transform=trf)
Receiver.add_object(FRUe)
# -- Cone surface
# --- Positioning of the Cone surface in the receiver frame
trc = N.dot(rotx(N.pi), translate(z=-depth))
# --- Construction of the Cone surface in the receiver frame
CON = AssembledObject(surfs=[Surface(FiniteCone(r=rc, h=d2), LambertianReceiver(absReceiver))], transform=trc)
Receiver.add_object(CON)
# --- Construction of the Cone external surface in the receiver frame
CONe = AssembledObject(surfs=[Surface(FiniteCone(r=rc+1e-6, h=d2), ReflectiveReceiver(0))], transform=trc)
Receiver.add_object(CONe)
# - Completion of the global assembly and positioning of the receiver.
Sim.add_assembly(Receiver, transform=trr)

'''
_____________________________________________________________________________________________
Source declaration
_____________________________________________________________________________________________
'''
sourceCenter = N.array([[0,0,2.*dishFocus]]).T # Source center position
sourceDirection = N.array([0,0,-1.]) # Source normal direction
sourceRadius = 0.6*dishDiameter # m, Source radius
sourceAngle = 4.5e-3 # radians, sun rays angular range
G = 1000. # W/m2, solar flux
nrays = 100 # number of rays escaping the source
SOURCE = solar_disk_bundle(nrays, sourceCenter, sourceDirection, sourceRadius, sourceAngle, G)

'''
_____________________________________________________________________________________________
Raytrace!
_____________________________________________________________________________________________
'''
engine = TracerEngine(Sim)
itmax = 100 # stop iteration after this many ray bundles were generated (i.e. 
	        # after the original rays intersected some surface this many times).
minener = 0.0001 # minimum energy threshold
engine.ray_tracer(SOURCE, itmax, minener, tree=True)

# Raytrace runtime calculation
t1 = time.clock()-t0
print 'Raytrace calculation time: ', t1

'''
_____________________________________________________________________________________________
Re-emissive losses calculations
Functions called:
FourParameterCavity - Calculates the view factor of the four parameter cavity using the viewax program.
radiosity - Calculates the radiosity equation to get the emissive losses.
_____________________________________________________________________________________________
'''
# Calculate geometry view factor
VF = FourParameterCavity(n,coneAngle,frustumAngle,apertureDiameter,depth)
# Solve radiosity problem
AA,bb,J,Eb,T,q = radiosity(VF,emsReceiver, tempAmbiant, tempReceiver)

# Emissive losses calculation time
t2 = time.clock()-t1
print 'Re-emissive losses calculation time:', t2

'''
_____________________________________________________________________________________________
Results interpretation
_____________________________________________________________________________________________
'''
# REFLECTIVE LOSSES
Cone_acc = AbsorptionAccountant.get_all_hits(CON.get_surfaces()[0].get_optics_manager())
Cone_abs, Cone_hits = Cone_acc[0], Cone_acc[1]
Frustum_acc = AbsorptionAccountant.get_all_hits(FRU.get_surfaces()[0].get_optics_manager())
Frustum_abs, Frustum_hits = Frustum_acc[0], Frustum_acc[1]
Dish_acc = AbsorptionAccountant.get_all_hits(DISH.get_surfaces()[0].get_optics_manager())
Dish_abs, Dish_hits = Dish_acc[0], Dish_acc[1]

# Energy balance over ray tracing (as a check)
# Source emitted rays
S = N.sum(SOURCE.get_energy())

# Total absorbed in the scene
A = N.sum(N.concatenate((Cone_abs,Frustum_abs,Dish_abs)))
ReceiverA = N.sum(N.concatenate((Cone_abs,Frustum_abs)))
print 'Total flux absorbed by the receiver', ReceiverA,'W'
# Total escaping from the scene
E = N.sum(N.concatenate(engine.lost_energy))
# Total reflective losses
Ref_losses_tot = N.sum(N.concatenate(engine.lost_energy[1:]))
print 'Total reflective losses :', Ref_losses_tot,'W'
# Total raytrace energy balance check.
B = S-A-E
print 'Raytrace balance check: ', B,'W','Fraction lost: ', B/S*100,'%. Should be close to 0 to certify no energy is lost.'

# RE-EMISSIVE LOSSES
print 'Emissive losses through the aperture :', q[0], 'W'
# Net absorbed flux = receiver absorbed flux - re emitted flux.
RealReceiverA = ReceiverA+q[0]
print 'Receiver net absorbed flux: ', RealReceiverA, 'W'

# System efficiency = (flux hitting after leaving source - reflective losses - emissive losses)/flux hitting after leaving source
flux_input = N.sum(engine.tree._bunds[1].get_energy())
System_efficiency = (flux_input-Ref_losses_tot+q[0])/flux_input
print 'System efficiency :', System_efficiency*100, '%'


'''
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
__________________________________________________________________________________________________
Rendering:

Renders the scene. Offers the option to highlight specific rays according to the number of times they have been 
reflected.
__________________________________________________________________________________________________
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''

def show_rays(engine, escaping_len=5., highlight_level=None):
    """
    Function to draw the rays to a Coin3D scenegraph.
    """
    tree = engine.tree
    no = coin.SoSeparator()
    #print tree.num_bunds()
    
    # loop through the reflection sequences?
    co = [] # regular lines
    co_h = [] # highlighted lines
    pos = [] # 2D level text position
    text = [] # 2D level text
    hist = {} # ray histories, for highlighted rays

    for level in xrange(tree.num_bunds()):
        start_rays = tree[level]
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
            if se[ray] <= minener:
                # ignore rays with starting energy smaller than energy cutoff
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
                # Highlight rays that have the highlight level set.
                hist[ray] = tree.ray_history(ray,level)
                co_h += [(c1[0],c1[1],c1[2]), (c2[0],c2[1],c2[2])]
            else:
                co += [(c1[0],c1[1],c1[2]), (c2[0],c2[1],c2[2])]
            # Position and text of the level 2D text. ratio is the parameter of the text position on the ray.
            # eg. 1/5 is equivalent to a text position at one fifth of the total ray length in the scene.
            ratio = 1./5
            c3 = ratio*(((1./ratio)-1.)*c1+c2)
            pos += [(c3[0],c3[1],c3[2])]        
            text.append(str(level))

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

        return no1

    def plot_level_number(text_pos, level_number):
        """
        Shows the number of reflections a ray has had as a 2D text on the scene.
        Arguments:
        text_pos - the position of the 2D text over the ray
        level_number - the number of reflections already encountered by the ray according to ray history.
        """
        no2 = coin.SoSeparator()
             
        tr = coin.SoTransform()
        tr.translation.setValue(text_pos)
        no2.addChild(tr)

        fo = coin.SoFont()
        fo.name.setValue("Arial-Bold")
        fo.size.setValue(15)
        no2.addChild(fo)   

        ma2 = coin.SoMaterial()
        ma2.diffuseColor.setValue(1,0,1)
        no2.addChild(ma2) 

        tx = coin.SoText2()      
        tx.string = level_number        
        no2.addChild(tx)
                
        return no2       
    
    no.addChild(plot_rays_color(co))
    no.addChild(plot_rays_color(co_h, color=(1,1,1)))
    num = len(text)    
    #for ref in range(num):
    #    no.addChild(plot_level_number(pos[ref], text[ref]))
    print "Number of reflections", num
    return no

def detailed_ray_history(engine,seq):
    """
    Print a detailed ray history for a particular ray, being a tuple returned 
    from TraceTree.ray_history.
    """
    tree = engine.tree
    n = len(seq)
    for i in range(n):
        bund = tree[i]
        ray = seq[(n-1)-i]
        print "...bund",i,"ray",ray
        print "...from",bund.get_vertices()[:,ray]
        print "...direction",bund.get_directions()[:,ray]

def axis_labels(length=1):
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

        la = coin.SoLabel()
        la.label = k
        s1.addChild(la)

        tr1 = coin.SoTranslation()
        tr1.translation = vec
        s1.addChild(tr1)

        r.addChild(s1)
    
        s2 = coin.SoSeparator()

        tr2 = coin.SoTransform()
        tr2.translation.setValue(data[k])
        s2.addChild(tr2)

        matxt = coin.SoMaterial()
        matxt.diffuseColor = data[k]
        s2.addChild(matxt)
       
        txaxis = coin.SoText2()      
        txaxis.string = k       
        s2.addChild(txaxis)

        r.addChild(s2)

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

# Render the scene with Pivy
r = coin.SoSeparator()

r.addChild(axis_labels())
r.addChild(show_rays(engine, highlight_level=None))
r.addChild(Sim.get_scene_graph())

win = SoGui.init("hello")
viewer = SoGuiExaminerViewer(win)
viewer.setSceneGraph(r)
viewer.setTitle("Examiner Viewer")
viewer.viewAll()
viewer.show()

SoGui.show(win)
# Rendering timer
t3 = time.clock() - t2
print '3D rendering calculation time:', t3

SoGui.mainLoop()


