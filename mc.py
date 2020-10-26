from functions import Nbp_mc, prob

#initialize nuclei dictionary
Au = {'name':'Au', 'radius':6.38, 'd': 0.535, 'A':197}

#problem 5
#N coll and N part
inst_Nbp = Nbp_mc(Au, Au, 4.2)
inst_Nbp.args['title'] = 'Au + Au at 200GeV (MC)'
inst_Nbp.args['save'] = True
inst_Nbp.args['path'] = './AuAuNbp_mc.png'
inst_Nbp.args['data_path'] = './AuAuNbp_mc_data.csv'
inst_Nbp.plot_func()

#problem 6
inst_prob = prob(14, './AuAuNbp_mc_data.csv')
inst_prob.args['save'] = True
inst_prob.args['path'] = './AuAuprob.png'
inst_prob.args['data_path'] = './AuAuprob_data.csv'
inst_prob.plot_func()