from functions import Nbp_mc

#initialize nuclei dictionary
Au = {'name':'Au', 'radius':6.38, 'd': 0.535, 'A':197}

inst_Nbp = Nbp_mc(Au, Au, 4.2)
inst_Nbp.args['title'] = 'Au + Au at 200GeV (MC)'
inst_Nbp.args['save'] = True
inst_Nbp.args['path'] = './AuAuNbp_mc.png'
inst_Nbp.args['data_path'] = './AuAuNbp_mc_data.csv'
inst_Nbp.plot_func()