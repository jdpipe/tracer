
import pivy.coin as coin
SOGUI_BINDING="SoQt"
from pivy.sogui import *

r = coin.SoSeparator()

co = coin.SoMaterial()
co.diffuseColor = (0.5,0.5,0.5)
co.specularColor = (1,0.5,0.5)
r.addChild(co)

# NUMPY BIT...
import numpy as np
nr = 10
nc = 10
x = np.linspace(5,-5,nc)
y = np.linspace(-5,5,nr)
X,Y = np.meshgrid(x,y)
Z = -np.sqrt(10**2 - X**2 - Y**2)

# from numpy back to standard python; pivy doesn't speak numpy yet
A = [(X.flat[i],Y.flat[i],Z.flat[i]) for i in range(len(X.flat))]

coor = coin.SoCoordinate3()
coor.point.setValues(0, len(A), A)

r.addChild(coor)

sh = coin.SoShapeHints()
sh.shapeType = coin.SoShapeHintsElement.UNKNOWN_SHAPE_TYPE
sh.vertexOrdering = coin.SoShapeHintsElement.COUNTERCLOCKWISE
r.addChild(sh)

qm = coin.SoQuadMesh()
qm.verticesPerRow = nc
qm.verticesPerColumn = nr

r.addChild(qm)

# note that we can manipulate nodes already added!
co.specularColor = (0.5,0.5,1)

win = SoGui.init("hello")

viewer = SoGuiExaminerViewer(win)
viewer.setSceneGraph(r)
viewer.setTitle("Examiner Viewer")
viewer.viewAll()
viewer.show()

SoGui.show(win)
SoGui.mainLoop()


