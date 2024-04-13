from PyQt6.QtCore import *
from PyQt6.QtGui import *
from PyQt6.QtWidgets import *


from scipy.linalg import *
from numpy import *
from math import *

#processing data
class Algorithms:
    
    def __init__(self):
        pass
    
    def computeAngle(self, p1: QPointF, p2: QPointF, p3: QPointF, p4: QPointF):
        # compute angle between two lines
        # p1, p2 - first line
        # p3, p4 - second line
        # return angle in degrees
        ux = p2.x() - p1.x()
        uy = p2.y() - p1.y()
        vx = p4.x() - p3.x()
        vy = p4.y() - p3.y()
        
        #Dot product
        dot = ux * vx + uy * vy
        
        #Vector norms
        norm_u = sqrt(ux*ux + uy*uy)
        norm_v = sqrt(vx*vx + vy*vy)
        
        # Angle
        angle = dot/(norm_u*norm_v)
        if angle > 1:
            angle = 1
        elif angle < -1: #### přidat argmin argmax
            angle = -1
        
        #Return angle
        return acos(angle)
        
    # create convex hull with Jarvis algorithm
    def jarvisScan(self, pol: QPolygonF):
        
        # inicialize convex hull
        convexHull = QPolygonF()
        
        # find pivot
        q = min(pol, key=lambda k: k.y())
        # lambda k: k.y() is used as the key function, 
        # which extracts the y() attribute from each element k in the list pol
        s = min(pol, key=lambda k: k.x())
        
        # inicialize last two points
        qj = q
        qj1 = QPointF(s.x(), q.y())
        
        # add pivot to convex hull
        convexHull.append(qj)
        
        # process all points
        while True:
            # inicilalize maximum angle and index
            maxAngle = 0
            maxIndex = -1

            for i in range(len(pol)):
                # find point with maximum angle
                if qj != pol[i]:
                # compute angle
                    angle = self.computeAngle(qj, qj1, qj, pol[i])
                    # find maximum angle
                    if angle > maxAngle:
                        maxAngle = angle
                        maxIndex = i

            # add point to convex hull
            convexHull.append(pol[maxIndex])            
            # stop condition if the last point is the pivot
            if pol[maxIndex] == q:
                break
            
            # update last two points (last added edge)
            qj1 = qj
            qj = pol[maxIndex]
            
        # return convex hull
        return convexHull
    
    # create minimum bounding rectangle
    def minMaxBox(self, pol: QPolygonF):
        
        px_min = min(pol, key=lambda k: k.x()).x()
        px_max = max(pol, key=lambda k: k.x()).x()
        py_min = min(pol, key=lambda k: k.y()).y()
        py_max = max(pol, key=lambda k: k.y()).y()
                
        v1 = QPointF(px_min, py_min)
        v2 = QPointF(px_max, py_min)
        v3 = QPointF(px_max, py_max)
        v4 = QPointF(px_min, py_max)
        
        # create minimum bounding rectangle
        mbr = QPolygonF([v1, v2, v3, v4])
        
        return mbr
    
    def rotatePolygon(self, pol: QPolygonF, angle: float):
        # rotate polygon by angle
        # pol - polygon
        # angle - angle in degrees
        # return rotated polygon
        pol_r = QPolygonF()
        
        #Process all points
        for p in pol:
            #Rotate point
            x_r = p.x() * cos(angle) - p.y() * sin(angle)
            y_r = p.x() * sin(angle) + p.y() * cos(angle)
    
            #Ceate rotated point
            p_r = QPointF(x_r, y_r)
            
            #Add to polygon
            pol_r.append(p_r)

        
        return pol_r
    
    def getArea(self, pol: QPolygonF):
        # compute area of polygon
        # pol - polygon
        # return area
        area = 0
        # polygon length
        n = len(pol)
        for i in range(n):
            area += pol[i].x() * (pol[(i+1) % n].y() - pol[(i-1+n) % n].y())
        return abs(area)/2
    
    def resizeRectangle(self, rect: QPolygonF, build: QPolygonF):
        # resize rectangle to fit the area of the building
        #input: building and rectangle (area)
        
        
        #compute are of the building and rectangle
        ab = self.getArea(build)
        a = self.getArea(rect)
        
        # compute ratio
        k = ab/a
        
        # center of the rectangle
        tx = (rect[0].x() + rect[1].x() + rect[2].x() + rect[3].x())/4
        ty = (rect[0].y() + rect[1].y() + rect[2].y() + rect[3].y())/4
        
        # vectors
        u1x = rect[0].x() - tx
        u1y = rect[0].y() - ty
        u2x = rect[1].x() - tx
        u2y = rect[1].y() - ty
        u3x = rect[2].x() - tx
        u3y = rect[2].y() - ty
        u4x = rect[3].x() - tx
        u4y = rect[3].y() - ty
        
        #new vertices
        v1x = tx + sqrt(k) * u1x
        v1y = ty + sqrt(k) * u1y
        v2x = tx + sqrt(k) * u2x
        v2y = ty + sqrt(k) * u2y
        v3x = tx + sqrt(k) * u3x
        v3y = ty + sqrt(k) * u3y
        v4x = tx + sqrt(k) * u4x
        v4y = ty + sqrt(k) * u4y
        
        v1 = QPointF(v1x, v1y)
        v2 = QPointF(v2x, v2y)
        v3 = QPointF(v3x, v3y)
        v4 = QPointF(v4x, v4y)
        
        # add vertices to the resized rectangle
        resizedRect = QPolygonF([v1, v2, v3, v4])
        return resizedRect
    
    def mbr(self, pol: QPolygonF):
        # minimum area enclosing rectangle

        # convex hull
        ch = self.jarvisScan(pol)

        # number of points
        n = len(ch)

        # minimum bounding box
        mmb_min = self.minMaxBox(ch)
        # area of the minimum bounding box
        area_mmb_min = self.getArea(mmb_min)

        for i in range(n):
            # coordinate differences
            dx = ch[(i + 1) % n].x() - ch[i].x()
            dy = ch[(i + 1) % n].y() - ch[i].y()

            # direction
            sigma = atan2(dy, dx)

            # rotate convex hull
            ch_rot = self.rotatePolygon(ch, -sigma)

            # find minimum bounding box
            mmb_rot = self.minMaxBox(ch_rot)
            area_rot = self.getArea(mmb_rot)

            # compare areas
            if area_rot < area_mmb_min:
                mmb_min = mmb_rot
                area_mmb_min = area_rot
                sigma_min = sigma

        mmb_rot_back = self.rotatePolygon(mmb_min, sigma_min)
        mmb_resized = self.resizeRectangle(mmb_rot_back, pol)
        return mmb_resized
    
    def createERPCA(self, pol: QPolygonF):
        #Create enclosing area using PCA 

        #lists of coordinates
        x = []
        y = []

        #add coordinates to lists
        for p in pol:
            x.append(p.x())
            y.append(p.y())

        #Create array
        P = array([x, y])

        #compute covariation matrix
        C = cov(P)

        #singular decompozition
        [U, S, V] = svd(C)

        #Copute sigma
        sigma = atan2(U[0][1], U[0][0])

        #rotate polygon by angle -sigma
        pol_unrot = self.rotatePolygon(pol, -sigma)

        #find min max box
        mmb = self.minMaxBox(pol_unrot)

        #ROTATE MIN MAX BOR
        er = self.rotatePolygon(mmb, sigma)

        #resize min max box (er_r - resized enclosing rectangle)
        er_r = self.resizeRectangle(er, pol)

        return er_r
    

    
    def longestEdge(self, pol: QPolygonF):
        # find the longest edge of the polygon
        # pol - polygon
        # return the longest edge
        n = len(pol)
        longest_edge = 0
        for i in range(n):
            # edge length
            edge = sqrt((pol[(i+1) % n].x() - pol[i].x())**2 + (pol[(i+1) % n].y() - pol[i].y())**2)
            #find longest edge
            if edge > longest_edge:
                longest_edge = edge
                #find slope of the longest edge
                angle = atan2(pol[(i+1) % n].y() - pol[i].y(), pol[(i+1) % n].x() - pol[i].x())
        #rotate building
        pol_r = self.rotatePolygon(pol, -angle)
        
        #find min max box
        mmb, area = self.minMaxBox(pol_r)
        #rotate min max box
        er = self.rotatePolygon(mmb, angle)
        #resize building
        er_r = self.resizeRectangle(er, pol)
        return er_r
    
    def wallAverage(self, pol: QPolygonF):
        # find the average wall length of the building
        # pol - polygon
        # return the average wall length
        # compute slope 
        dx = pol[1].x() - pol[0].x()
        dy = pol[1].y() - pol[0].y()
        angle = atan2(dy, dx)
        
        #process all edges
        n = len(pol)
        
        #inicilize remainder of the angle
        wall_sum= 0
        
        for i in range(1,n):
            # compute angle for current edge
            dx_i = pol[(i+1) % n].x() - pol[i].x()
            dy_i = pol[(i+1) % n].y() - pol[i].y()
            angle_i = atan2(dy_i, dx_i)
            
            # compute angle difference
            angle_diff = abs(angle_i - angle)
            
            # compute fraction of pi
            ki = 2*angle_diff/pi
            
            #compute remainder
            rem_i = angle_diff - ki*pi/2
            
            # update wall sum
            wall_sum = wall_sum + rem_i
            
        #compute average direction
        angle_ave = angle + wall_sum/n
        
        #rotate building
        pol_r = self.rotatePolygon(pol, -angle_ave)
        #find min max box
        mmb, area = self.minMaxBox(pol_r)
        #rotate min max box
        er = self.rotatePolygon(mmb, angle_ave)
        #resize building
        er_r = self.resizeRectangle(er, pol)
        return er_r
            
        
        

