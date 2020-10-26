from utils import rect_vector, polar_vector, find_bin, rdg, integrate, nucleon_generator, bsq, distance
import numpy as np
import matplotlib.pyplot as plt

r'''
If you wanna show figures instead of saving them, 
set args of corresponding custom function 'save' as False, 
and comment the following line.
'''
plt.switch_backend('agg')

#Optical part
class custom_function():
    r'''
    Custom functions have those basic methods:
    __init__:   set default arguments
    function:   the function of this class, e.g. Woods Saxon, thickness
    get_value:  assign self.x and self.y via apply function of this class
    save_data:  save the function value (discreted points)
    plot_func:  draw the plot pf this class, need to set self.args
    '''
    def __init__(self):
        self.args = {
            'title':'title',
            'save':False,
            'path':'./res.png',
            'data_path':'./data_path.csv'
        }
        self.x = None
        self.y = None
        self.valued = False
        return

    def function(self, *args, **kwargs):
        pass

    def get_value(self, *args, **kwargs):
        pass

    def save_data(self, *args, **kwargs):
        res = np.concatenate((self.x.reshape(-1, 1), self.y.reshape(-1, 1)), 1)
        np.savetxt(self.args['data_path'], res, fmt='%.4f', delimiter=',')
        return

    def plot_func(self, *args, **kwargs):
        pass

class Woods_Saxon(custom_function):
    def __init__(self, nuclei:dict):
        super(Woods_Saxon, self).__init__()
        self.radius = nuclei['radius']
        self.d = nuclei['d']
        self.A = nuclei['A']
        self.rho_0 = 0.17
        self.args['x'] = np.linspace(0, 15, 150)
        self.args['label'] = 'Woods Saxon'
        self.args['data_path'] = './Woods_Saxon_data.csv'
        return

    def function(self, r):
        r'''
        Here, r is a scalar.
        Input the length of it if it's a vector. 
        '''
        return self.rho_0/(1+np.exp((r - self.radius)/self.d))

    def get_value(self):
        if self.valued is False:
            args = self.args
            self.x = args['x']
            self.y = self.function(self.x)
        self.valued = True
        return

    def plot_func(self):
        r'''
        Note: if your computer has not install SciencePlots, 
        replace the list following plt.style.context by ['classic'].
        '''
        args = self.args
        self.get_value()
        self.save_data()
        with plt.style.context(['ieee', 'science', 'no-latex', 'grid']):
            fig, ax = plt.subplots()
            ax.plot(self.x, self.y, label=args['label'])
            ax.set_xlabel(r'$r\ \mathrm{fm}$')
            ax.set_ylabel(r'$\rho_A(r)$')
            ax.legend()
            ax.set_title('%s Radius=%.4f, d=%.4f'%(args['title'], self.radius, self.d))
            if args['save']:
                fig.savefig(args['path'])
            else:
                fig.show()
        return

class thickness(custom_function):
    def __init__(self, func:Woods_Saxon):
        super(thickness, self).__init__()
        self.ws = func
        self.args['smin'] = 0
        self.args['smax'] = 10
        self.args['sbins'] = 1000
        self.args['zmin'] = 0
        self.args['zmax'] = 15
        self.args['label'] = 'Thickness Function'
        self.args['data_path'] = './thickness_data.csv'
        self.A = self.ws.A
        return

    def integ_func(self, z, s):
        r = (z**2 + s**2)**0.5
        return self.ws.function(r)*2

    def function(self, s):
        r'''
        This method returns a 1d array res.
        Here, argument s needs to be self.x
        '''
        res = np.zeros_like(s)
        for idx, item in enumerate(s, 0):
            res[idx] =  integrate.quad( self.integ_func, 
                                            self.args['zmin'], 
                                            self.args['zmax'],
                                            args=(item,))[0]
        return res
    
    def get_value(self):
        if self.valued is False:
            args = self.args
            self.x = np.linspace(self.args['smin'], self.args['smax'], self.args['sbins'])
            self.y = self.function(self.x)
        self.valued = True
        return

    def plot_func(self):
        r'''
        Note: if your computer has not install SciencePlots, 
        replace the list following plt.style.context by ['classic'].
        '''
        args = self.args
        self.get_value()
        self.save_data()
        with plt.style.context(['ieee', 'science', 'no-latex', 'grid']):
            fig, ax = plt.subplots()
            ax.plot(self.x, self.y, label=args['label'])
            ax.set_xlabel(r'$\vec{r}_\perp\ \mathrm{fm}$')
            ax.set_ylabel(r'$T_A(\vec{r}_\perp)\ \mathrm{fm}^{-2}$')
            ax.legend()
            ax.set_title('%s'%args['title'])
            if args['save']:
                fig.savefig(args['path'])
            else:
                fig.show()
        return

class overlap(custom_function):
    def __init__(self, func1:thickness, func2:thickness):
        super(overlap, self).__init__()
        self.tk1 = func1
        self.tk1.get_value()
        self.tk1x = self.tk1.x
        self.tk1y = self.tk1.y
        self.tk2 = func2
        self.tk2.get_value()
        self.tk2x = self.tk2.x
        self.tk2y = self.tk2.y
        self.args['bmin'] = 0
        self.args['bmax'] = 20
        self.args['bbins'] = 100
        self.args['smin'] = 0
        self.args['smax'] = 15
        self.args['label'] = 'Overlap Function'
        self.args['data_path'] = './overlap_data.csv'
        self.A1 = self.tk1.A
        self.A2 = self.tk2.A
        return

    def integ_func(self, rho, theta, b):
        x1 = self.tk1x
        y1 = self.tk1y
        x2 = self.tk2x
        y2 = self.tk2y
        s = rho
        spb = (polar_vector(rho, theta).rect_vec + rect_vector(b, 0)).get_length
        sbin1 = find_bin.find_bin1d(s, x1)
        T1 = y1[sbin1]
        sbin2 = find_bin.find_bin1d(spb, x2)
        T2 = y2[sbin2]
        return T1*T2

    def function(self, b):
        r'''
        This method returns a 1d array res.
        Here, argument b needs to be self.x
        '''
        res = np.zeros_like(b)
        for idx, item in enumerate(b, 0):
            res[idx] =  integrate.dblquad(  self.integ_func, 
                                            self.args['smin'], self.args['smax'],
                                            lambda x:0, lambda x:2*np.pi,
                                            args=(item,), epsabs=1e-1)[0]
        return res

    def get_value(self):
        if self.valued is False:
            args = self.args
            self.x = np.linspace(self.args['bmin'], self.args['bmax'], self.args['bbins'])
            self.y = self.function(self.x)
        self.valued = True
        return

    def plot_func(self):
        r'''
        Note: if your computer has not install SciencePlots, 
        replace the list following plt.style.context by ['classic'].
        '''
        args = self.args
        self.get_value()
        self.save_data()
        with plt.style.context(['ieee', 'science', 'no-latex', 'grid']):
            fig, ax = plt.subplots()
            ax.plot(self.x, self.y, label=args['label'])
            ax.set_xlabel(r'$b\ \mathrm{fm}$')
            ax.set_ylabel(r'$T_{AA}(\vec{r}_\perp)\ \mathrm{fm}^{-2}$')
            ax.legend()
            ax.set_title('%s'%args['title'])
            if args['save']:
                fig.savefig(args['path'])
            else:
                fig.show()
        return

class Ncoll(custom_function):
    def __init__(self, func:overlap, sigma):
        super(Ncoll, self).__init__()
        self.ol = func
        self.ol.get_value()
        self.olx = self.ol.x
        self.oly = self.ol.y
        self.args['label'] = r'$<N_{\mathrm{coll}}>$'
        self.args['sigma'] = sigma
        self.args['data_path'] = './Ncoll.csv'
        return
    
    def function(self):
        pass

    def get_value(self):
        if self.valued is False:
            self.x = self.olx
            self.y = self.oly * self.args['sigma']
        self.valued = True
        return

    def plot_func(self):
        r'''
        Note: if your computer has not install SciencePlots, 
        replace the list following plt.style.context by ['classic'].
        '''
        args = self.args
        self.get_value()
        self.save_data()
        with plt.style.context(['ieee', 'science', 'no-latex', 'grid']):
            fig, ax = plt.subplots()
            ax.plot(self.x, self.y, label=args['label'])
            ax.set_xlabel(r'$b\ \mathrm{fm}$')
            ax.set_ylabel(r'$<N_{\mathrm{coll}}>$')
            ax.legend()
            ax.set_title('%s'%args['title'])
            if args['save']:
                fig.savefig(args['path'])
            else:
                fig.show()
        return

class Npart(custom_function):
    def __init__(self, func1:thickness, func2:thickness, sigma):
        super(Npart, self).__init__()
        self.tk1 = func1
        self.tk1.get_value()
        self.tk1x = self.tk1.x
        self.tk1y = self.tk1.y
        self.tk2 = func2
        self.tk2.get_value()
        self.tk2x = self.tk2.x
        self.tk2y = self.tk2.y
        self.args['bmin'] = 0
        self.args['bmax'] = 20
        self.args['bbins'] = 100
        self.args['smin'] = 0
        self.args['smax'] = 15
        self.args['label'] = r'$<N_{\mathrm{part}}>$'
        self.args['sigma'] = sigma
        self.args['data_path'] = './Npart_data.csv'
        self.A1 = self.tk1.A
        self.A2 = self.tk2.A
        return

    def integ_func1(self, rho, theta, b):
        x1 = self.tk1x
        y1 = self.tk1y
        x2 = self.tk2x
        y2 = self.tk2y
        s = rho
        spb = (polar_vector(rho, theta).rect_vec + rect_vector(+b, 0)).get_length
        sbin1 = find_bin.find_bin1d(s, x1)
        T1 = y1[sbin1]
        sbin2 = find_bin.find_bin1d(spb, x2)
        T2 = y2[sbin2]
        expr = T1*(1 - (1 - T2/self.A2*self.args['sigma'])**self.A2)
        return expr

    def integ_func2(self, rho, theta, b):
        x1 = self.tk1x
        y1 = self.tk1y
        x2 = self.tk2x
        y2 = self.tk2y
        s = rho
        spb = (polar_vector(rho, theta).rect_vec + rect_vector(b, 0)).get_length
        sbin1 = find_bin.find_bin1d(s, x2)
        T1 = y2[sbin1]
        sbin2 = find_bin.find_bin1d(spb, x1)
        T2 = y1[sbin2]
        expr = T1*(1 - (1 - T2/self.A1*self.args['sigma'])**self.A1)
        return expr

    def function(self, b):
        r'''
        This method returns a 1d array res.
        Here, argument b needs to be self.x
        '''
        res = np.zeros_like(b)
        for idx, item in enumerate(b, 0):
            res[idx] =  integrate.dblquad(  self.integ_func1, 
                                            self.args['smin'], self.args['smax'],
                                            lambda x:0, lambda x:2*np.pi,
                                            args=(item,), epsabs=1e-1)[0] + \
                        integrate.dblquad(  self.integ_func2, 
                                            self.args['smin'], self.args['smax'],
                                            lambda x:0, lambda x:2*np.pi,
                                            args=(item,), epsabs=1e-1)[0]           
        return res        

    def get_value(self):
        if self.valued is False:
            args = self.args
            self.x = np.linspace(self.args['bmin'], self.args['bmax'], self.args['bbins'])
            self.y = self.function(self.x)
        self.valued = True
        return

    def plot_func(self):
        r'''
        Note: if your computer has not install SciencePlots, 
        replace the list following plt.style.context by ['classic'].
        '''
        args = self.args
        self.get_value()
        self.save_data()
        with plt.style.context(['ieee', 'science', 'no-latex', 'grid']):
            fig, ax = plt.subplots()
            ax.plot(self.x, self.y, label=args['label'])
            ax.set_xlabel(r'$b\ \mathrm{fm}$')
            ax.set_ylabel(r'$<N_{\mathrm{part}}>$')
            ax.legend()
            ax.set_title('%s'%args['title'])
            if args['save']:
                fig.savefig(args['path'])
            else:
                fig.show()
        return

#Monte-Carlo Glauber part
class Nbp_mc(custom_function):
    def __init__(self, nuclei1, nuclei2, sigma, entries=1000):
        r'''
        Nbp_mc would generate a set of nucleon matters for two nuclei respectively, 
        however, in this homework they are all Au.
        '''
        super(Nbp_mc, self).__init__()
        self.nuclei1 = nuclei1
        self.nuclei2 = nuclei2
        self.nuleon1 = None
        self.nuleon2 = None
        self.A1 = nuclei1['A']
        self.A2 = nuclei2['A']
        self.sigma = sigma
        self.entries = entries
        self.args['bmin'] = 0
        self.args['bmax'] = 14
        self.args['bbins'] = 15
        self.args['label1'] = r'$<N_{\mathrm{coll}}>$'
        self.args['label2'] = r'$<N_{\mathrm{part}}>$'
        self.args['data_path'] = './Nbp_data.csv'
        return

    def simulation(self, b):
        self.nuleon1 = nucleon_generator(Woods_Saxon(self.nuclei1)).get_particle('rect')
        self.nuleon2 = nucleon_generator(Woods_Saxon(self.nuclei2)).get_particle('rect')
        self.nuleon2[:, 0] += b#shift target nucleon's x, and y will not change
        Ncoll = 0
        Npart = 0
        for nid_1 in range(self.A1):
            projectile = self.nuleon1[nid_1, :].reshape(1, -1)
            dist = distance.Euclidean(projectile, self.nuleon2)
            dist = np.where(np.pi*dist**2<=self.sigma, 1, 0)
            num = dist.sum()
            Ncoll += num
            if num > 0:
                Npart += 1
        for nid_2 in range(self.A2):
            projectile = self.nuleon2[nid_2, :].reshape(1, -1)
            dist = distance.Euclidean(projectile, self.nuleon1)
            dist = np.where(dist**2*np.pi<=self.sigma, 1, 0)
            num = dist.sum()
            if num > 0:
                Npart += 1
        return Ncoll, Npart
    
    def get_value(self):
        if self.valued is False:
            self.x = np.linspace(self.args['bmin'], self.args['bmax'], self.args['bbins'])
            self.Ncoll = []
            self.Npart = []
            for item in self.x:
                print('\nb = %.3f'%item)
                Ncoll_list = []
                Npart_list = []
                for entry in range(self.entries):
                    print('\r %d of %d'%(entry+1, self.entries), end='')
                    res = self.simulation(item)
                    Ncoll_list.append(res[0])
                    Npart_list.append(res[1])
                self.Ncoll.append(np.mean(Ncoll_list))
                self.Npart.append(np.mean(Npart_list))
            self.Ncoll = np.array(self.Ncoll)
            self.Npart = np.array(self.Npart)
        self.valued = True
        return

    def save_data(self):
        res = np.concatenate((self.x.reshape(-1, 1), self.Ncoll.reshape(-1, 1), self.Npart.reshape(-1, 1)), 1)
        np.savetxt(self.args['data_path'], res, fmt='%.4f', delimiter=',')
        return

    def plot_func(self):
        r'''
        Note: if your computer has not install SciencePlots, 
        replace the list following plt.style.context by ['classic'].
        '''
        args = self.args
        self.get_value()
        self.save_data()
        with plt.style.context(['ieee', 'science', 'no-latex', 'grid']):
            fig, ax = plt.subplots()
            ax.plot(self.x, self.Ncoll, 'r:', label=args['label1'])
            ax.plot(self.x, self.Npart, 'b--', label=args['label2'])
            ax.set_xlabel(r'$b\ \mathrm{fm}$')
            ax.set_ylabel(r'$<N_{\mathrm{coll}}>\ \mathrm{and}\ <N_{\mathrm{part}}>$')
            ax.legend()
            ax.set_title('%s'%args['title'])
            if args['save']:
                fig.savefig(args['path'])
            else:
                fig.show()
        return

class prob(custom_function):
    def __init__(self, bmax, load_path, entries=100000):
        super(prob, self).__init__()
        data = np.loadtxt(load_path, delimiter=',')
        self.b = data[:, 0]
        self.Ncoll = data[:, 1]
        self.Npart = data[:, 2]
        self.bsq = bsq(bmax)
        self.rdg = rdg(self.bsq.function, [0, bmax], [0, self.bsq.function(bmax)])
        self.entries = entries
        self.args['bbins'] = 14
        self.args['Ncollbins'] = 14
        self.args['Npartbins'] = 14
        self.args['label1'] = 'b'
        self.args['label2'] = r'$<N_{\mathrm{coll}}>$'
        self.args['label3'] = r'$<N_{\mathrm{part}}>$'
        self.args['data_path'] = './prob_data.csv'        
        return

    def simulation(self):
        r'''
        This method returns a index,
        use this index to derive b, Ncoll and Npart of this simulation.
        '''
        return find_bin.find_bin1d(self.rdg.get_random(), self.b)

    def get_value(self):
        if self.valued is False:
            b = []
            Ncoll = []
            Npart = []
            for entry in range(self.entries):
                print('\rEntry: %d of %d...'%(entry+1, self.entries), end='')
                idx = self.simulation()
                b.append(self.b[idx])
                Ncoll.append(self.Ncoll[idx])
                Npart.append(self.Npart[idx])
            self.b_res = np.array(b)
            self.Ncoll_res = np.array(Ncoll)
            self.Npart_res = np.array(Npart)
        self.valued = True
        return

    def save_data(self):
        res = np.concatenate((self.b_res.reshape(-1, 1), self.Ncoll_res.reshape(-1, 1), self.Npart_res.reshape(-1, 1)), 1)
        np.savetxt(self.args['data_path'], res, fmt='%.4f', delimiter=',')
        return

    def plot_func(self):
        r'''
        Note: if your computer has not install SciencePlots, 
        replace the list following plt.style.context by ['classic'].
        '''
        args = self.args
        self.get_value()
        self.save_data()
        with plt.style.context(['ieee', 'science', 'no-latex', 'grid']):
            fig, ax = plt.subplots(2, 2)
            fig.set_size_inches(8, 8)
            ax[0, 0].hist(self.b_res, histtype='step', density=True, label=args['label1'])
            ax[0, 1].hist(self.Ncoll_res, histtype='step', density=True, label=args['label2'])
            ax[1, 1].hist(self.Npart_res, histtype='step', density=True, label=args['label3'])
            ax[0, 0].set_xlabel(r'$b\ \mathrm{fm}$')
            ax[0, 0].set_ylabel(r'$P(b)$')
            ax[0, 1].set_xlabel(r'$<N_{\mathrm{coll}}>$')
            ax[0, 1].set_ylabel(r'$P(<N_{\mathrm{coll}}>)$')
            ax[1, 1].set_xlabel(r'$<N_{\mathrm{part}}>$')
            ax[1, 1].set_ylabel(r'$P(<N_{\mathrm{part}}>$)')
            ax[0, 0].legend()
            ax[0, 1].legend()
            ax[1, 1].legend()
            ax[1, 0].text(0.5, 0.5, 'made by yghuang', horizontalalignment='center', verticalalignment='center', transform=ax[1, 0].transAxes)
            fig.tight_layout()
            if args['save']:
                fig.savefig(args['path'])
            else:
                fig.show()
        return
            

