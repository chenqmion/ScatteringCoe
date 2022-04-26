import numpy as np
import matplotlib.pyplot as plt
import Toolkit_Simulation as tk

L = 5E-3
C_1 = 1e-14
load = 'lamb4'
w_list = np.linspace(6.637,6.683,10001) * 1E9 * (2*np.pi)

# L = 10E-3
# C_1 = 2e-14
# geometry = 'hanger'
# load = 'lamb2'
# w_list = np.linspace(6.637,6.683,1001) * 1E9 * (2*np.pi)

Num_S = []
Ana_S = []
for w in w_list:
    #%% Numerical
    Zr = tk.Waveguide_Z(w,L, load=load)
    Z = 1/(1j*w*C_1) + Zr
    M = tk.TM_Vert(Z)
    Num_S.append(tk.Num_SC(M))
    
    #%% analytical calculations
    params = tk.Ana_Params_Hanger(L,C_1, load=load)
    # params[1] = params[1] + 1E6*(np.pi)
    Ana_S.append(tk.Ana_SC_Hanger(w, params))

Num_S = np.array(Num_S)
Ana_S = np.array(Ana_S)

# np.save('Num_S_hanger.npy', Num_S)
# np.save('Ana_S_hanger.npy', Ana_S)

#%% plots
Num_S = Num_S.reshape(len(w_list),4)
Ana_S = Ana_S.reshape(len(w_list),4)
Labels = [r'$S_{11}$', r'$S_{12}$', r'$S_{21}$', r'$S_{22}$']

w_min = float("{:.2f}".format(w_list[0]/(2*np.pi)*1E-9))
w_max = float("{:.2f}".format(w_list[-1]/(2*np.pi)*1E-9))
w_range = w_max - w_min

fig = plt.figure( figsize=(7,3) )
grid = plt.GridSpec(2,4, wspace=0.05, hspace=0.05)
for num_1 in range(0,4):
    Num_data = Num_S[:,num_1]
    Num_Am = np.abs(Num_data)
    Num_Ph = np.angle(Num_data)
    
    Ana_data = Ana_S[:,num_1]
    Ana_Am = np.abs(Ana_data)
    Ana_Ph = np.angle(Ana_data)
    
    ax1 = plt.subplot(grid[0,num_1])
    ax1.plot(w_list/(2*np.pi)*1E-9, Num_Am, linewidth=2, color='red', label='Numerical')
    ax1.plot(w_list/(2*np.pi)*1E-9, Ana_Am, linewidth=2, linestyle='--', color='black', label='Analytical')
    ax1.set_xlim(w_min-w_range/7, w_max+w_range/7)
    ax1.set_ylim(-0.1, 1.1)
    
    ax1.set_xticks([])
    ax1.set_title(Labels[num_1])
        
    ax2 = plt.subplot(grid[1,num_1])
    ax2.plot(w_list/(2*np.pi)*1E-9, Num_Ph, linewidth=2, color='green', label='Numerical')
    ax2.plot(w_list/(2*np.pi)*1E-9, Ana_Ph, linewidth=2, linestyle='--', color='black', label='Analytical')
    ax2.set_xlim(w_min-w_range/7, w_max+w_range/7)
    ax2.set_ylim(-3.5, 3.5)
    
    ax2.set_xticks([w_min, w_max])
    ax2.set_xticklabels([str(w_min), str(w_max)])

    if (num_1 == 0):
        ax1.set_yticks([0, 0.5, 1])
        ax1.set_yticklabels(['0.0','0.5', '1.0'])
        ax1.set_ylabel('Amplitude $|S|$')
        
        ax2.set_yticks([-np.pi, 0, np.pi])
        ax2.set_yticklabels([r'$-\pi$', r'$0$', r'$\pi$'])
        ax2.set_ylabel(r'Phase $\angle S$ (rad)')
    else:
        ax1.set_yticks([])
        ax2.set_yticks([])
    
fig.text(0.5, -0.02, r'Frequency $\omega$ ($2\pi\times GHz$)', ha='center')    
plt.savefig('Fig_hanger.pdf', bbox_inches='tight')
plt.show()
        

