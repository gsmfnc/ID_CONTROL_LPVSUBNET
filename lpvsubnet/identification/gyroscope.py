import deepSI
import numpy as np
from scipy.io import loadmat
from matplotlib import pyplot as plt

from modified_lpvsubnet import LPV_single_encoder_mod
from modified_lpvsubnet import LPV_single_encoder_mod_double_input

out = loadmat('../../data/ML_estim.mat')
u = np.concatenate([out['exp_u'], \
    np.transpose([np.array(out['exp_p'][:, 0])])], axis = 1)
y = out['exp_y']
train = deepSI.System_data(u = u, y = y)

out = loadmat('../../data/ML_valid.mat')
u = np.concatenate([out['exp_u'], \
    np.transpose([np.array(out['exp_p'][:, 0])])], axis = 1)
y = out['exp_y']
test = deepSI.System_data(u = u, y = y)

# Self-scheduling, input=[i_2]

sys = LPV_single_encoder_mod(nx=5, Np=2, na=10, nb=10, feedthrough=True, \
            include_u_in_p=True, f_net_kwargs=dict(F=10), \
            e_net_kwargs=dict(n_nodes_per_layer = 64, n_hidden_layers = 2), \
            p_net_kwargs=dict(n_nodes_per_layer = 64, n_hidden_layers = 2))
sys.init_model(sys_data=train)
sys.unique_code = "D5fiSA"
sys.fit(train, train[-25000:], epochs = 1000, loss_kwargs=dict(nf=80, \
    loss_nf_cutoff=1), print_full_time_profile = True)
sys.checkpoint_save_system('_best', 'results')

sys.checkpoint_load_system("_best", "results")

print("LPV single encoder, nx=5, Np=2, na=nb=10")
res = sys.apply_experiment(test)
print("Test NRMS: " + str(res.NRMS(test) * 100) + "%")

# External-scheduling, input=[i_2]

from modified_lpvsubnet import LPV_multi_encoder_mod

sys_m = LPV_multi_encoder_mod(nx=5, Np=2, na=10, nb=10, feedthrough=True, \
            include_u_in_p=True, f_net_kwargs=dict(F=10), \
            e_net_kwargs=dict(n_nodes_per_layer = 64, n_hidden_layers = 2), \
            p_net_kwargs=dict(n_nodes_per_layer = 64, n_hidden_layers = 2))
sys_m.init_model(sys_data=train)
sys_m.unique_code = "wds1ew"
sys_m.fit(train, train[-25000:], epochs = 1000, loss_kwargs=dict(nf=80, \
    loss_nf_cutoff=1), print_full_time_profile = True)
sys_m.checkpoint_save_system('_best', 'results')

sys_m.checkpoint_load_system("_best", "results")

print("LPV multi encoder, nx=5, Np=2, na=nb=10")
res = sys_m.apply_experiment(test)
print("Test NRMS: " + str(res.NRMS(test) * 100) + "%")

# Self-scheduling, input=[i_2, i_1]

out = loadmat('../../data/DBL_estim.mat')
u = np.concatenate([out['exp_u_est'], \
    np.transpose([np.array(out['exp_p_est'][:, 0])])], axis = 1)
y = out['exp_y_est']
train = deepSI.System_data(u = u, y = y)

out = loadmat('../../data/DBL_valid.mat')
u = np.concatenate([out['exp_u_val'], \
    np.transpose([np.array(out['exp_p_val'][:, 0])])], axis = 1)
y = out['exp_y_val']
test = deepSI.System_data(u = u, y = y)

sys = LPV_single_encoder_mod_double_input(nx=5, Np=2, na=10, nb=10, \
            feedthrough=True, \
            include_u_in_p=True, f_net_kwargs=dict(F=10), \
            e_net_kwargs=dict(n_nodes_per_layer = 64, n_hidden_layers = 2), \
            p_net_kwargs=dict(n_nodes_per_layer = 64, n_hidden_layers = 2))
sys.init_model(sys_data=train)
sys.unique_code = "n31BRg"
sys.fit(train, train[-25000:], epochs = 1000, loss_kwargs=dict(nf=80, \
    loss_nf_cutoff=1), print_full_time_profile = True)
sys.checkpoint_save_system('_best', 'results')

sys.checkpoint_load_system("_best", "results")

print("LPV single encoder, nx=5, Np=2, na=nb=10")
res = sys.apply_experiment(test)
print("Test NRMS: " + str(res.NRMS(test) * 100) + "%")
