def point_within_triangle(pt1, pt2, pt3):
   import random
   """
   Random point within the triangle with vertices pt1, pt2 and pt3.
   """
   s, t = sorted([random.random(), random.random()])
   return (s * pt1[0] + (t-s)*pt2[0] + (1-t)*pt3[0],
         s * pt1[1] + (t-s)*pt2[1] + (1-t)*pt3[1])