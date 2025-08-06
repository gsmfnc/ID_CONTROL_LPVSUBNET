import deepSI
import numpy as np
from scipy.io import loadmat, savemat
from modified_lpvsubnet import LPV_single_encoder_mod
from modified_lpvsubnet import LPV_multi_encoder_mod

#sys = LPV_single_encoder_mod(nx=5, Np=2, na=10, nb=10, feedthrough=True, \
#            include_u_in_p=True, f_net_kwargs=dict(F=10), \
#            e_net_kwargs=dict(n_nodes_per_layer = 64, n_hidden_layers = 2), \
#            p_net_kwargs=dict(n_nodes_per_layer = 64, n_hidden_layers = 2))
#sys.unique_code = "D5fiSA"
#sys.checkpoint_load_system("_best", "results")

sys = LPV_multi_encoder_mod(nx=5, Np=2, na=10, nb=10, feedthrough=True, \
            include_u_in_p=True, f_net_kwargs=dict(F=10), \
            e_net_kwargs=dict(n_nodes_per_layer = 64, n_hidden_layers = 2), \
            p_net_kwargs=dict(n_nodes_per_layer = 64, n_hidden_layers = 2))
sys.unique_code = "wds1ew"
sys.checkpoint_load_system("_best", "results")

A = sys.fn.A
As = sys.fn.As
B = sys.fn.B
Bs = sys.fn.Bs
C = sys.hn.A
Cs = sys.hn.As
D = sys.hn.B
Ds = sys.hn.Bs

def save_to_mat(M, varname, filename):
    n = len(M)
    m = len(M[0])
    string = varname + " = [\n"
    for i in range(0, n):
        for j in range(0, m):
            string = string + str(M[i][j].item())
            if j < m - 1:
                string = string + ", "
        string = string + "\n"
    string = string + "];\n"
    with open(filename, "a") as f:
        f.write(string)

def save_to_arr(M, varname, filename):
    n = len(M)
    string = varname + " = [\n"
    for i in range(0, n):
        string = string + str(M[i].item()) + "\n"
    string = string + "];\n"
    with open(filename, "a") as f:
        f.write(string)

import os
os.system('rm matrices.m')

save_to_mat(A, "A", "matrices.m")
for i in range(0, len(As)):
    save_to_mat(As[i], "A" + str(i + 1), "matrices.m")

save_to_arr(B[:, 0], "B", "matrices.m")
for i in range(0, len(Bs)):
    save_to_arr(Bs[i][:, 0], "B" + str(i + 1), "matrices.m")

save_to_mat(C, "C", "matrices.m")
for i in range(0, len(Cs)):
    save_to_mat(Cs[i], "C" + str(i + 1), "matrices.m")

save_to_arr(D[:, 0], "D", "matrices.m")
for i in range(0, len(Ds)):
    save_to_arr(Ds[i][:, 0], "D" + str(i + 1), "matrices.m")
            
#pnet
lin_part = sys.pnet.get_submodule('net_lin')
save_to_mat(lin_part.weight, "V", "matrices.m")
save_to_arr(lin_part.bias, "c", "matrices.m")

nonlin_part = sys.pnet.get_submodule('net_non_lin').get_submodule('net')
linear0 = nonlin_part.get_submodule('0')
tanh1 = nonlin_part.get_submodule('1')
linear2 = nonlin_part.get_submodule('2')
tanh3 = nonlin_part.get_submodule('3')
linear4 = nonlin_part.get_submodule('4')

save_to_mat(linear0.weight, "W0", "matrices.m")
save_to_mat(linear2.weight, "W1", "matrices.m")
save_to_mat(linear4.weight, "W2", "matrices.m")

save_to_arr(linear0.bias, "b0", "matrices.m")
save_to_arr(linear2.bias, "b1", "matrices.m")
save_to_arr(linear4.bias, "b2", "matrices.m")

#encoder
enc = sys.encoder.get_submodule('net')

lin_part = enc.get_submodule('net_lin')
save_to_mat(lin_part.weight, "Ve", "matrices.m")
save_to_arr(lin_part.bias, "ce", "matrices.m")

nonlin_part = enc.get_submodule('net_non_lin').get_submodule('net')
linear0 = nonlin_part.get_submodule('0')
tanh1 = nonlin_part.get_submodule('1')
linear2 = nonlin_part.get_submodule('2')
tanh3 = nonlin_part.get_submodule('3')
linear4 = nonlin_part.get_submodule('4')

save_to_mat(linear0.weight, "W0e", "matrices.m")
save_to_mat(linear2.weight, "W1e", "matrices.m")
save_to_mat(linear4.weight, "W2e", "matrices.m")

save_to_arr(linear0.bias, "b0e", "matrices.m")
save_to_arr(linear2.bias, "b1e", "matrices.m")
save_to_arr(linear4.bias, "b2e", "matrices.m")

out = loadmat('../../data/ML_estim.mat')
u = np.concatenate([out['exp_u'], \
    np.transpose([np.array(out['exp_p'][:, 0])])], axis = 1)
y = out['exp_y']
train = deepSI.System_data(u = u, y = y)
trainn = sys.norm.transform(train)

save_to_mat(trainn.u, "train_u", "matrices.m")
save_to_arr(trainn.y, "train_y", "matrices.m")

out = loadmat('../../data/ML_valid.mat')
u = np.concatenate([out['exp_u'], \
    np.transpose([np.array(out['exp_p'][:, 0])])], axis = 1)
y = out['exp_y']
test = deepSI.System_data(u = u, y = y)
testt = sys.norm.transform(test)

save_to_mat(testt.u, "valid_u", "matrices.m")
save_to_arr(testt.y, "valid_y", "matrices.m")

string = "As = {A, A1, A2};\n"
string = string + "Bs = {B, B1, B2};\n"
string = string + "Cs = {C, C1, C2};\n"
string = string + "Ds = {D, D1, D2};\n"
with open("matrices.m", "a") as f:
    f.write(string)
