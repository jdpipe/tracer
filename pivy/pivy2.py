
import pivy.coin as coin
SOGUI_BINDING="SoQt"
from pivy.sogui import *

r = coin.SoSeparator()

r.addChild(coin.SoCube())

import numpy as np
# Note: pivy doesn't speak numpy yet... wouldn't that be great though?
m = ((2,0,0.4,0),(0,3,0,0),(-0.4,0,4,0),(1,5,6,1))
# How do we convert a matrix to the list structure above in a fast/efficient way?

tr = coin.SoTransform()
m1 = coin.SbMatrix(m)
print "MESH\n\n",m1.getValue()

# Note: pivy doesn't have a good string representation of SbVec3f or SbRotation.
mt,mr,ms,msa = m1.getTransform()
print "TRANSL = ",mt.getValue()
print "ROT = ",mr.getValue()
print "SCALE = ",ms.getValue()
print "MSA = ",msa.getValue()

tr.setMatrix(m1)
r.addChild(tr)

r.addChild(coin.SoCube())

win = SoGui.init("hello")

viewer = SoGuiExaminerViewer(win)
viewer.setSceneGraph(r)
viewer.setTitle("Examiner Viewer")
viewer.viewAll()
viewer.show()

SoGui.show(win)
SoGui.mainLoop()


