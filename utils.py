import numpy as np
import numpy.random as rd
from scipy import integrate

class rect_vector():
    def __init__(self, x, y):
        self.x = x
        self.y = y
        return

    def __add__(self, other):
        return rect_vector(self.x + other.x, self.y + other.y)

    def __mul__(self, factor):
        return rect_vector(self.x*factor, self.y*factor)

    @property
    def get_length(self):
        return (self.x**2 + self.y**2)**0.5

    @property
    def polar_vec(self):
        return polar_vector(self.get_length, np.arctan2(self.y, self.x))

class polar_vector():
    def __init__(self, rho, theta):
        self.rho = rho
        self.theta = theta
        return

    @property
    def rect_vec(self):
        return rect_vector(self.rho*np.cos(self.theta), self.rho*np.sin(self.theta))

class find_bin():
    def __init__(self):
        r'''
        This class find location of a variable in a array-like sequence, or a location of a pixel in a matrix-like object.
        '''
        pass

    @staticmethod
    def find_bin1d(x, seq):
        if x <= seq[0]:
            return 0
        else:
            return (x>seq).sum() - 1

    @staticmethod
    def find_bin2d(x, y, seqx, seqy):
        bx = find_bin.find_bin1d(x, seqx)
        by = find_bin.find_bin1d(y, seqy)
        return [bx, by]

class rdg():
    def __init__(self, function, xlim, ylim):
        self.function = function
        self.xlim = xlim
        self.ylim = ylim
        return
    
    def get_random(self):
        while True:
            x = rd.uniform(self.xlim[0], self.xlim[1])
            f = rd.uniform(self.xlim[0], self.xlim[1])
            fx = self.function(x)
            if f <= fx:
                return x
            
class nucleon_generator():
    def __init__(self, woods_saxon):
        self.A = woods_saxon.A
        self.r = []
        ymax = woods_saxon.function(0)
        ws_rdg = rdg(woods_saxon.function, [0, 15], [0, ymax])
        for _ in range(self.A):
            self.r.append(ws_rdg.get_random())
        self.r = np.array(self.r)
        self.theta = rd.uniform(0, np.pi, self.A)
        self.phi = rd.uniform(0, 2*np.pi, self.A)
        self.s = self.r*np.sin(self.theta)
        self.x = self.s*np.cos(self.phi)
        self.y = self.s*np.sin(self.phi)
        return

    def get_particle(self, mode='polar'):
        r'''
        Return a 2d array: rho and phi of nuclei matter of this particle.
        Or x and y.
        '''
        if mode == 'polar':
            res = np.concatenate((self.s.reshape(-1, 1), self.phi.reshape(-1, 1)), 1)
        elif mode == 'rect':
            res = np.concatenate((self.x.reshape(-1, 1), self.y.reshape(-1, 1)), 1)
        else:
            raise Exception('Invalid mode for nucleon_generator.get_particle.')
        return res

class bsq():
    def __init__(self, bmax):
        self.bmax = bmax
        self.a = 1.0 / bmax**3 * 3
        return

    def function(self, x):
        return self.a*x**2

class distance():
    @staticmethod
    def Euclidean(projectile, target):
        r'''
        projectile: a [1, 2] array, for (x, y) of projectile
        target:     a [n, 2] array, for a set of (x, y) of target
        Return:     a [n,] array, for Euclidean distance
        '''
        delta = projectile - target
        return np.sqrt(delta[:, 0]**2 + delta[:, 1]**2)