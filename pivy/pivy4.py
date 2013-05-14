
import pivy.coin as coin
SOGUI_BINDING="SoQt"
from pivy.sogui import *

r = coin.SoSeparator()

ds = coin.SoDrawStyle()
ds.style = coin.SoDrawStyle.LINES
ds.lineWidth = 2
#ds.pointSize = 10
r.addChild(ds)
print "line width =",ds.lineWidth.getValue()

# NUMPY BIT...
A = [(0,0,0), (1,0,0), (0,5,0), (0,0,0)]

coor = coin.SoCoordinate3()
coor.point.setValues(0,4,A)
r.addChild(coor)

#B = numpy.array([2,2,2,2]).flatten()
B = (2,2)
ls = coin.SoLineSet()
ls.numVertices.setValues(0,2,B)
r.addChild(ls)

win = SoGui.init("hello")

viewer = SoGuiExaminerViewer(win)
viewer.setSceneGraph(r)
viewer.setTitle("Examiner Viewer")
viewer.viewAll()
viewer.show()

SoGui.show(win)
SoGui.mainLoop()


