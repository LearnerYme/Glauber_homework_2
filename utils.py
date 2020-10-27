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
    
    def get_random(self, *args):
        while True:
            x = rd.uniform(self.xlim[0], self.xlim[1])
            f = rd.uniform(self.xlim[0], self.xlim[1])
            fx = self.function(x, *args)
            if f <= fx:
                return x
            
class mc_sample():
    @staticmethod
    def sample_rho(r, woods_saxon):
        return r*r*woods_saxon.function(r)
    
    @staticmethod
    def sample_theta(theta):
        return np.sin(theta)

    @staticmethod
    def sample_phi(low, high):
        r'''
        Not recommand for using this, please just sample phi by rd.uniform().
        '''
        return 1.0/(high - low)
        

class nucleon_generator():
    def __init__(self, woods_saxon, ws_max=None):
        r'''
        Note that, the variables can be sampled respected to following distributions:
        rho:    rho^2*Woods_Saxon
        theta:  sin(theta)
        phi:    uniform
        If you know the maximum of Woods Saxon value, tell the function to save computing resources.
        '''
        self.A = woods_saxon.A
        self.r = []
        self.theta = []
        self.phi = []
        #find the max boundary of r*r*ws
        x = np.linspace(0, 15, 150)
        if ws_max is None:
            ws_max = mc_sample.sample_rho(x, woods_saxon)
        rho_rdg = rdg(mc_sample.sample_rho, [0, 15], [0, ws_max])
        theta_rdg = rdg(mc_sample.sample_theta, [0, np.pi], [0, 1])
        #phi_rdg = rdg(mc_sample.sample_phi, [0, 2*np.pi], [0, 1])
        for _ in range(self.A):
            self.r.append(rho_rdg.get_random(woods_saxon))
            self.theta.append(theta_rdg.get_random())
        self.r = np.array(self.r)
        self.theta = np.array(self.theta)
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