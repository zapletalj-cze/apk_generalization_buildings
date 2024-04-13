from PyQt6.QtCore import *
from PyQt6.QtGui import *
from PyQt6.QtGui import QMouseEvent, QPaintEvent
from PyQt6.QtWidgets import *
import shapefile


class Draw(QWidget):
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        #definice našeho bodu a jeho velikosti

        self.building = QPolygonF()
        self.ch = QPolygonF()
        self.mbr = QPolygonF()
        self.pol = QPolygonF()
        self.polygons = []
        self.features = None
        self.min_max = [0,0,10,10]
        self.q = QPointF(-100, -100)
        self.pol = QPolygonF()
        self.result = []

    # function for loading input data
    def loadData(self, filename):   
        # load objects from shapefile
        filename, _ = QFileDialog.getOpenFileName(self, "Open file", "", "Shapefile (*.shp)")
        if filename:
            shp = shapefile.Reader(filename)
            self.features = shp.shapes()

            # calculate minimum and maximum coordinates
            min_x, min_y, max_x, max_y = float('inf'), float('inf'), -float('inf'), -float('inf')
            for feature in self.features:
                for point in feature.points:
                    min_x = min(min_x, point[0])
                    min_y = min(min_y, point[1])
                    max_x = max(max_x, point[0])
                    max_y = max(max_y, point[1])
            # store minimum and maximum coordinates
            self.min_max = [min_x, min_y, max_x, max_y]
        
    # function for rescaling data to fit Canvas  
    def resizeData(self):
        # get width and height of the widget
        width = self.frameGeometry().width()
        height = self.frameGeometry().height()
        
        # initialize list for storing polygons
        self.polygons = [None] * len(self.features)

        # rescale data and create polygons
        for k, feature in enumerate(self.features):
            self.polygons[k] = QPolygonF()
            for point in feature.points:
                x = int(((point[0] - self.min_max[0]) / (self.min_max[2] - self.min_max[0]) * width))
                y = int((height - (point[1] - self.min_max[1]) / (self.min_max[3] - self.min_max[1]) * height))
                p = QPointF(x, y)
                self.polygons[k].append(p)        

    def mousePressEvent(self, e: QMouseEvent):
        #get cursor position
        x = e.position().x()
        y = e.position().y()
        
        # #Add new vertex
        # p = QPointF(x, y)
            
        # #Add p to polygon
        # self.building.append(p)
            
        #Repaint screen
        self.repaint()
        
    def paintEvent(self, e: QPaintEvent):
        #draw situation
        
        #create new object
        qp = QPainter(self)
        #start drawing
        qp.begin(self)
        # #set graphical attributes
        qp.setPen(Qt.GlobalColor.black)
        qp.setBrush(Qt.GlobalColor.cyan)
        
        #draw polygon
        qp.drawPolygon(self.building)
        
        for index, polygon in enumerate(self.polygons):
            # set graphical attributes
            qp.setPen(QPen(Qt.GlobalColor.black))
            qp.setBrush(Qt.GlobalColor.lightGray)
            qp.drawPolygon(polygon)
        
        #draw convex hull
        qp.setPen(Qt.GlobalColor.red)
        qp.setBrush(Qt.GlobalColor.transparent)
        
        #draw minimum bounding box
        qp.drawPolygon(self.mbr)
        
        
        qp.end()
        
    def getBuilding(self):
        # Return building
        return self.building
    
    def setMBR(self, mbr):
        #Set MBR
        self.mbr = mbr
        
    
    
    def clearData(self):
        #Clear building
        self.building.clear()
        
        #Clear CH
        self.ch.clear()
        
        #Clear MBR
        self.mbr.clear()
                
        #Repaint screen
        self.repaint()
        
