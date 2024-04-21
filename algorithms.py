from PyQt6.QtCore import *
from PyQt6.QtGui import *
from PyQt6.QtWidgets import *
from math import *
from numpy import *
from scipy.linalg import *
import sys


# Processing data
class Algorithms:

    def __init__(self):
        pass

    def analyzePointPolygonPosition(self, q: QPointF, pol: QPolygonF):

        # Inicialize amount of intersections
        k = 0

        # Amount of vertices
        n = len(pol)

        # Process all segments
        for i in range(n):
            # Reduce coordinates
            xir = pol[i].x() - q.x()
            yir = pol[i].y() - q.y()

            xi1r = pol[(i + 1) % n].x() - q.x()
            yi1r = pol[(i + 1) % n].y() - q.y()

            # Suitable segment?
            if ((yi1r > 0) and (yir <= 0)) or ((yir > 0) and (yi1r <= 0)):

                # Compute intersection
                xm = (xi1r * yir - xir * yi1r) / (yi1r - yir)

                # Right half plane
                if xm > 0:
                    k += 1

        # Point q inside polygon?
        if k % 2 == 1:
            return 1

        # Point q outside polygon
        return 0

    def get2LineAngle(self, p1: QPointF, p2: QPointF, p3: QPointF, p4: QPointF):
        # Get 2 line angle
        ux = p2.x() - p1.x()
        uy = p2.y() - p1.y()

        vx = p4.x() - p3.x()
        vy = p4.y() - p3.y()

        # Dot product
        dot = ux * vx + uy * vy

        # Vector norms
        nu = (ux ** 2 + uy ** 2) ** (1 / 2)
        nv = (vx ** 2 + vy ** 2) ** (1 / 2)
        # Correct interval
        arg = dot / (nu * nv)
        return acos(max(min(arg, 1), -1))

    def distance_points(self, point1, point2):
        x1, y1 = point1
        x2, y2 = point2
        return sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

    def cHull(self, pol: QPolygonF):
        # CH construction using Jarvis scan algorithm
        ch = QPolygonF()

        # Find pivot 1
        q = min(pol, key=lambda k: k.y())

        # Find pivot 2
        s = min(pol, key=lambda k: k.x())

        # Initialize last 2 points of CH
        qj = q
        qj1 = QPointF(s.x(), q.y())

        # Add pivot to CH
        ch.append(q)

        # Find all points of CH
        while True:
            # Maximum and its index
            omega_max = 0
            index_max = -1

            # Process all point
            for i in range(len(pol)):

                # Compute angle
                if qj != pol[i] and qj1 != qj:
                    omega = self.get2LineAngle(qj, qj1, qj, pol[i])

                    # Update maximum
                    if omega > omega_max:
                        omega_max = omega
                        index_max = i

            # Append point to CH
            ch.append(pol[index_max])

            # We found pivot again
            if pol[index_max] == q:
                break

                # Update last segment of CH
            qj1 = qj
            qj = pol[index_max]

        return ch

    def createMMB(self, pol: QPolygonF):
        # Compute points with extreme coordinates
        px_min = min(pol, key=lambda k: k.x())
        px_max = max(pol, key=lambda k: k.x())
        py_min = min(pol, key=lambda k: k.y())
        py_max = max(pol, key=lambda k: k.y())

        # Compute min-max box points
        v1 = QPointF(px_min.x(), py_min.y())
        v2 = QPointF(px_max.x(), py_min.y())
        v3 = QPointF(px_max.x(), py_max.y())
        v4 = QPointF(px_min.x(), py_max.y())

        # Create min max box
        box = QPolygonF([v1, v2, v3, v4])

        return box

    def rotate(self, pol: QPolygonF, sig: float):
        # Rotate polygon by given angle

        pol_r = QPolygonF()

        # Process all points
        for p in pol:
            # Rotate point
            x_r = p.x() * cos(sig) - p.y() * sin(sig)
            y_r = p.x() * sin(sig) + p.y() * cos(sig)

            # Ceate rotated point
            p_r = QPointF(x_r, y_r)

            # Add to polygon
            pol_r.append(p_r)

        return pol_r

    def getArea(self, pol: QPolygonF):
        # Return polygon area
        area = 0
        n = len(pol)

        # Proccesing of vertexes
        for i in range(n):
            area = area + pol[i].x() * (pol[(i + 1) % n].y() - pol[(i - 1 + n) % n].y())

        return abs(area) / 2

    def resizeRectangle(self, rect: QPolygonF, build: QPolygonF):
        # Resize rectangle to fit area of the building

        # Compute areas
        Ab = self.getArea(build)
        A = self.getArea(rect)
        if A > 0:
            # Compute ratio
            k = Ab / A

            # Center of mass
            tx = (rect[0].x() + rect[1].x() + rect[2].x() + rect[3].x()) / 4
            ty = (rect[0].y() + rect[1].y() + rect[2].y() + rect[3].y()) / 4

            # Vectors
            u1x = rect[0].x() - tx
            u1y = rect[0].y() - ty
            u2x = rect[1].x() - tx
            u2y = rect[1].y() - ty
            u3x = rect[2].x() - tx
            u3y = rect[2].y() - ty
            u4x = rect[3].x() - tx
            u4y = rect[3].y() - ty

            # New vertices
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

            # Add vertices to polygon
            rectR = QPolygonF([v1, v2, v3, v4])
        else:
            rectR = None

        return rectR

    def createMBR(self, pol: QPolygonF):
        # Create minimum bounding rectangle

        # Compute convex hull
        ch = self.cHull(pol)

        # Initialization
        mmb_min = self.createMMB(ch)
        area_min = self.getArea(mmb_min)
        sigma_min = 0

        # Process all segments of CH
        n = len(ch)
        for i in range(n):
            # Coordinate differences
            dx = ch[(i + 1) % n].x() - ch[i].x()
            dy = ch[(i + 1) % n].y() - ch[i].y()

            # Direction
            sigma = atan2(dy, dx)

            # Rotate convex hull by -sigma
            ch_rot = self.rotate(ch, -sigma)

            # Find mmb and its area
            mmb_rot = self.createMMB(ch_rot)
            area_rot = self.getArea(mmb_rot)

            # Is it a better approximation?
            if area_rot < area_min:
                mmb_min = mmb_rot
                area_min = area_rot
                sigma_min = sigma

        # Back rotation
        mmb_unrot = self.rotate(mmb_min, sigma_min)

        # Resize rectangle
        mmb_res = self.resizeRectangle(mmb_unrot, pol)

        return mmb_res

    def createERPCA(self, pol: QPolygonF):
        # Create enclosing rectangle using PCA
        x = []
        y = []

        # Add coordinates to lists
        for p in pol:
            x.append(p.x())
            y.append(p.y())

        # Create array
        P = array([x, y])

        # Compute covariation matrix
        C = cov(P)

        # Singular value decomposition
        [U, S, V] = svd(C)

        # Compute sigma
        sigma = atan2(V[0][1], V[0][0])

        # Rotate polygon by minus sigma
        pol_unrot = self.rotate(pol, -sigma)

        # Find min-max box
        mmb = self.createMMB(pol_unrot)

        # Rotate min-max box (create enclosing rectangle)
        er = self.rotate(mmb, sigma)

        # Resize enclosing rectangle
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
        pol_r = self.rotate(pol, -angle)
        
        #find min max box
        mmb = self.createMMB(pol_r)
        #rotate min max box
        er = self.rotate(mmb, angle)
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
        pol_r = self.rotate(pol, -angle_ave)
        #find min max box
        mmb = self.createMMB(pol_r)
        #rotate min max box
        er = self.rotate(mmb, angle_ave)
        #resize building
        er_r = self.resizeRectangle(er, pol)
        return er_r

    def weighted_bisector(self, pol: QPolygonF):
        dist = {}
        points = [(point.x(), point.y()) for point in pol]
        # find two longest diagonals and its lengths
        for i in range(len(points)):
            for j in range(i + 1, len(points)):
                point1 = points[i]
                point2 = points[j]
                key = (point1, point2)
                dist[key] = self.distance_points(point1, point2)
        dist_sort = dict(sorted(dist.items(), key=lambda item: item[1], reverse=True))
        distances = []
        angles = []
        i = 0
        # calculate angle between x-axis and diagonals
        for key, value in dist_sort.items():
            points = key
            point1, point2 = [QPointF(*point) for point in points]
            delta_y = point2.y() - point1.y()
            delta_x = point2.x() - point1.x()
            angle_x = atan2(delta_y, delta_x) - pi/2
            angles.append(angle_x)
            distances.append(value)
            i += 1
            if i > 1:
                break
        # Calculate rotation using a weighted average, where the weight is the length of the diagonals.
        angle_weighted = (distances[0] * angles[0] + distances[1] * angles[1])/(distances[0]+distances[1])
        # rotate building
        pol_r = self.rotate(pol, -angle_weighted)
        # find min max box
        mmb = self.createMMB(pol_r)
        # rotate min max box
        er = self.rotate(mmb, angle_weighted)
        # resize building
        er_r = self.resizeRectangle(er, pol)
        return er_r

    # def calculate_accuracy(self, pol: QPolygonF, er: QPolygonF):
    #
    #     # calculate main direction of the enclosing rectangle based on the longest edge
    #     k = len(er)
    #     longest_edge = 0
    #     point1 = QPointF()
    #     point2 = QPointF()
    #     angle
    #
    #     for i in range(k):
    #         dx = er[(i + 1) % k].x() - er[i].x()
    #         dy = er[(i + 1) % k].y() - er[i].y()
    #         edge_length = sqrt(dx ** 2 + dy ** 2)
    #         if edge_length > longest_edge:
    #             longest_edge = edge_length
    #             point1 = er[i]
    #             point2 = er[(i + 1) % k]
    #
    #     # compute angle of the longest edge
    #     main_direction = atan2(point2.y() - point1.y(), point2.x() - point1.x())
    #
    #     # compute angle of the building
    #     # process all edges
    #     n = len(pol)
    #
    #     # inicilize remainder of the angle
    #     wall_sum = 0
    #
    #     for i in range(1, n):
    #         # compute angle for current edge
    #         dx_i = pol[(i + 1) % n].x() - pol[i].x()
    #         dy_i = pol[(i + 1) % n].y() - pol[i].y()
    #         sigma_i = atan2(dy_i, dx_i)
    #
    #         # compute angle difference
    #         # angle_diff = abs(angle_i - main_direction)
    #
    #         # compute fraction of pi
    #         ki = 2 * sigma_i / pi
    #         ki_floor = floor(ki)
    #
    #
    #         # compute remainder
    #         ri = (ki - ki_floor) * pi / 2
    #         angle_diff = abs(ri-main_direction)
    #         wall_sum += angle_diff
    #         angle_ave = pi/2*n * wall_sum
    #         angle_ave_deg = angle_ave*180/pi
    #     print(angle_ave_deg)
    #     # return angle_ave


    def calculate_accuracy(self, pol: QPolygonF, er: QPolygonF):
        # calculate main direction of the enclosing rectangle based on the longest edge's direction
        n_er = len(er)
        longest_edge = 0
        point1 = QPointF()
        point2 = QPointF()

        for i in range(n_er):
            dx = er[(i + 1) % n_er].x() - er[i].x()
            dy = er[(i + 1) % n_er].y() - er[i].y()
            edge_length = sqrt(dx ** 2 + dy ** 2)
            if edge_length > longest_edge:
                longest_edge = edge_length
                point1 = er[i]
                point2 = er[(i + 1) % n_er]

        # compute angle of the longest edge (consider quadrant handling)
        main_direction = atan2(point2.y() - point1.y(), point2.x() - point1.x())

        # process all edges
        n_pol = len(pol)

        # initialize remainder of the angle
        angular_deviations = []

        for i in range(n_pol):
            # compute current edge vector
            dx_i = pol[(i + 1) % n_pol].x() - pol[i].x()
            dy_i = pol[(i + 1) % n_pol].y() - pol[i].y()
            angular_deviation = atan2(dy_i, dx_i) - main_direction
            if angular_deviation > pi:
                angular_deviation -= 2 * pi

            angular_deviations.append(angular_deviation)

        mean_angular_deviation = sum(angular_deviations) / n_pol
        mean_angular_deviation = mean_angular_deviation * 180/pi

        # check if the mean angular deviation is greater than 10 degrees and return signal to canvas
        if mean_angular_deviation > 10:
            return 0
        else:
            return 1



            


                
        
        

