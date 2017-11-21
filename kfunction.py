import sys
sys.path.append('../geom')
sys.path.append('../indexing')

from geom.point import *
from indexing.extent import *
from indexing.kdtree1 import *
from indexing.kdtree2b import *

import random
from math import sqrt, pi
import numpy as np

def kfunc(tree, p, d, density):
    """
    Input
      tree:    a k-D tree
      p:       a point where the K-function is computed
      d:       radius of the circle around p
      density: density of points in the area

    Return
      n: count of points in the circle
      ld: L(d) value
    """
    neighbors = []
    range_query_circular(tree, p, d, neighbors)
    n = len(neighbors)
    kd = float(n)/density
    ld = sqrt(kd/pi)
    return n, ld

def kfunc_monte_carlo(n, area, radii, density, rounds=100):
    """
    Input
      n:            number of points
      area:         Extent object defining the area
      radii:        list containing a set of radii of circles
      density:      density of point events in the area
      rounds:       number of simulations
    Return
      percentiles:  a list of 2.5th and 97.5th percentiles
                    for each d in radii
    """
    alllds = []
    for test in range(rounds):
        points = [ Point(random.uniform(area.xmin, area.xmax),
                         random.uniform(area.ymin, area.ymax))
                   for i in range(n) ]
        t = kdtree2(points)
        lds = [0 for d in radii]
        for i, d in enumerate(radii):
            for p in points:
                ld = kfunc(t, p, d, density)[1]
                lds[i] += ld
        lds = [ld/n for ld in lds]
        alllds.append(lds)
    alllds = np.array(alllds)
    percentiles = []
    for i in range(len(radii)):
        percentiles.append([np.percentile(alllds[:,i], 2.5),
                            np.percentile(alllds[:,i], 97.5)])
    return percentiles

def test(points, area):
    """
    Input
      points: a list of Point objects
      area: an Extent object

    Return
      ds: list of radii
      percentiles: Monte Carlo result of using radii in ds
      lds: L(d) values for each radius in ds
    """
    n = len(points)
    density = n/area.area()
    t = kdtree2(points)
    d = min([area.xmax-area.xmin,area.ymax-area.ymin])*2/3/10
    ds = [ d*(i+1) for i in range(10)]
    lds = [0 for d in ds]
    for i, d in enumerate(ds):
        for p in points:
            ld = kfunc(t, p, d, density)[1]
            lds[i] += ld
    lds = [ld/n for ld in lds]
    percentiles = kfunc_monte_carlo(n, area, ds, density)
    return ds, percentiles, lds

if __name__ == "__main__":
    n = 100
    # random patter
    points = [ Point(random.uniform(20, 30),
                     random.uniform(30, 40))
               for i in range(n) ]
    xmin = min([p.x for p in points])
    ymin = min([p.y for p in points])
    xmax = max([p.x for p in points])
    ymax = max([p.y for p in points])
    area = Extent(xmin=xmin, xmax=xmax,
                  ymin=ymin, ymax=ymax)
    ds, percentiles, lds = test(points, area)
    print("Random")
    for i, v in enumerate(percentiles):
        print(ds[i], v[0], lds[i], v[1])

    # three blocks of points
    points1 = [ Point(random.uniform(10, 20),
                      random.uniform(10, 20))
                for i in range(n/3) ]
    points2 = [ Point(random.uniform(30, 40),
                      random.uniform(10, 20))
                for i in range(n/3) ]
    points3 = [ Point(random.uniform(20, 30),
                      random.uniform(30, 40))
                for i in range(n/3) ]
    points = points1 + points2 + points3
    xmin = min([p.x for p in points])
    ymin = min([p.y for p in points])
    xmax = max([p.x for p in points])
    ymax = max([p.y for p in points])
    area = Extent(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)
    ds, percentiles, lds = test(points, area)
    print("Three blocks")
    for i, v in enumerate(percentiles):
        print(ds[i], v[0], lds[i], v[1])
