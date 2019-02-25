
import numpy as np
from ..tools.rotation import a_on_b

def populate_circle(radius, distance):
    # generate 3d-(x,y)-grid with hexagonal shape

    x0 = np.mgrid[0:radius+distance:distance]
    x0 = np.append(-x0[-1:0:-1], x0)

    y0 = np.mgrid[0:radius+distance:distance * 2*np.cos(np.deg2rad(30))]
    y0 = np.append(-y0[-1:0:-1], y0)

    x0, y0 = np.meshgrid(x0, y0)

    dy = np.cos(np.deg2rad(30))*distance
    dx = np.sin(np.deg2rad(30))*distance
    x = np.concatenate((x0.flatten(), x0.flatten()+dx)).tolist()
    y = np.concatenate((y0.flatten(), y0.flatten()+dy)).tolist()

    # delete all points outside the circle
    r2 = radius**2
    for i in range(len(x)-1, -1, -1):
        if x[i]**2+y[i]**2 > r2:
            del x[i]
            del y[i]

    return np.array([np.array(x), np.array(y), np.zeros(len(x))]).T

def populate_object(object, population, radius_fibers):
    '''
    object: 2d array of np.array([x,y,z,r]) coordinates
    points: 2d array of np.array([x,y,z]) coordinates
    '''

    if population.shape[1] == 2:
        points = np.append(population, np.zeros((population.shape[0], 1)), axis=1)
    elif population.shape[1] == 3:
        points = population.copy()
    if population.shape[1] != 3:
        raise TypeError('wrong input format')

    if object.shape[1] != 4:
        raise TypeError('wrong input format')

    if object.shape[0] < 2:
        raise TypeError('object size is to small')


    fibers = np.empty([points.shape[0], object.shape[0], 4])
    tangent_old = np.array([0, 0, 1], float)
    radius_old = object[0, -1]

    for i in range(0, object.shape[0]):
        if i == 0:
            tangent_new = 0.5*(object[i+1, :-1]-object[i, :-1])
        elif i == object.shape[0]-1:
            tangent_new = 0.5*(object[i, :-1]-object[i-1, :-1])
        else:
            tangent_new = 0.5*(object[i+1, :-1]-object[i-1, :-1])

        if np.all(tangent_new == 0):
            raise TypeError('same point:', i)

        tangent_new = tangent_new / np.linalg.norm(tangent_new)
        R = a_on_b(tangent_old, tangent_new)

        for j in range(0, points.shape[0]):
            points[j,:] = np.dot(R, points[j,:]) / radius_old * object[i, -1]
            fibers[j, i,:] = np.append(points[j,:] + object[i, :-1], radius_fibers)

        tangent_old = tangent_new.copy()
        radius_old = object[i, -1]

    return fibers
