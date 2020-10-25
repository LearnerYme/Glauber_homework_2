from functions import Woods_Saxon, thickness, overlap, Ncoll, Npart

#initialize nuclei dictionary
Au = {'name':'Au', 'radius':6.38, 'd': 0.535, 'A':197}
Cu = {'name':'Cu', 'radius':4.20641, 'd':0.5977, 'A':'Unknwon'}
Pb = {'name':'Pb', 'radius':6.62, 'd':0.546, 'A':'Unknwon'}

#draw Woods Saxon and thickness
def job1(nuclei):
    r'''
    Return the instance of thickness for job2 (only for Au)
    '''
    #problem 1
    #Woods Saxon
    inst_woods_saxon = Woods_Saxon(nuclei)
    inst_woods_saxon.args['title'] = nuclei['name']
    inst_woods_saxon.args['save'] = True
    inst_woods_saxon.args['path'] = './%s_Woods_Saxon.png'%nuclei['name']
    inst_woods_saxon.args['data_path'] = './%s_Woods_Saxon_data.csv'%nuclei['name']
    inst_woods_saxon.plot_func()
    #problem 2
    #thickness
    inst_thickness = thickness(inst_woods_saxon)
    inst_thickness.args['title'] = nuclei['name']
    inst_thickness.args['save'] = True
    inst_thickness.args['path'] = './%s_thickness.png'%nuclei['name']
    inst_thickness.args['data_path'] = './%s_thickness_data.csv'%nuclei['name']
    inst_thickness.plot_func()
    return inst_thickness

#draw overlap, Ncoll and Npart, only for Au
def job2(nuclei):
    #problem 1 and 2
    inst_thickness = job1(nuclei)
    #problem 3
    #overlap
    inst_overlap = overlap(inst_thickness, inst_thickness)
    inst_overlap.args['title'] = 'Au + Au at 200GeV'
    inst_overlap.args['save'] = True
    inst_overlap.args['path'] = './AuAu_overlap.png'
    inst_overlap.args['data_path'] = './AuAu_overlap_data.csv'
    inst_overlap.plot_func()
    #N collision
    inst_Ncoll = Ncoll(inst_overlap, 4.2)
    inst_Ncoll.args['title'] = 'Au + Au at 200 GeV'
    inst_Ncoll.args['save'] = True
    inst_Ncoll.args['path'] = './AuAu_Ncoll.png'
    inst_Ncoll.args['data_path'] = './AuAu_Ncoll_data.csv'
    inst_Ncoll.plot_func()
    #problem 4
    #N part
    inst_Npart = Npart(inst_thickness, inst_thickness, 4.2)
    inst_Npart.args['title'] = 'Au + Au at 200GeV'
    inst_Npart.args['save'] = True
    inst_Npart.args['path'] = './AuAu_Npart.png'
    inst_Npart.args['data_path'] = './AuAu_Npart_data.csv'
    inst_Npart.plot_func()
    return

for nuclei in [Cu, Pb]:
    job1(nuclei)
job2(Au)
