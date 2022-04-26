import numpy as np
import matplotlib.pyplot as plt
import Toolkit_Simulation as tk

L = 10E-3
C_1 = 1e-14
load = 'lamb2'
w_list = np.linspace(6.57,6.75,1001) * 1E9 * (2*np.pi)

N_wg = 4
Z0 = 50

Num_S = []
Ana_S = []
for w in w_list:
    #%% Numerical
    M_c = tk.TM_Hori(1/(1j*w*C_1))
    Mr = tk.TM_Waveguide(w,L)
    
    M_chain =  M_c @ Mr
    for num_r in range(0,N_wg-1):
        M_chain =  M_chain @ M_c @ Mr
    M_chain = M_chain @ M_c
    Num_S.append(tk.Num_SC(M_chain))
    
    #%% analytical calculations
    params = tk.Ana_Params_Necklace(L,C_1,C_1, load=load)
    w_0, w_r, Qi, Qe1, Qe2, Ql = params
    gamm_a = w_r/Qi
    gamm = w_r/Qe1
    g = (w_r**2)*Z0*C_1/np.pi

    sum_1 = 0
    sum_2 = 0
    sum_12 = 0
    for num_k in range(0,N_wg):
        Delta_k = (w_r-w) - 2*g*np.cos((num_k+1)*np.pi/(N_wg+1))
        gamm_k = (2*gamm/(N_wg+1)) * (np.sin((num_k+1)*np.pi/(N_wg+1)))**2
        
        sum_1 +=  gamm_k/(1j*Delta_k + gamm_a/2)
        sum_2 +=  gamm_k/(1j*Delta_k + gamm_a/2)
        sum_12 +=  (-1)**(num_k+1) * gamm_k/(1j*Delta_k + gamm_a/2)
        
    Val_down = (1+sum_1/2)*(1+sum_2/2) - (sum_12/2)**2
    Val_S11 = (1+sum_2/2) * sum_1 - (1/2) * (sum_12)**2
    Val_S21 = sum_12
    Val_S12 = Val_S21
    Val_S22 = (1+sum_1/2) * sum_2 - (1/2) * (sum_12)**2
    
    S_11 = 1 - Val_S11/Val_down 
    S_12 = -Val_S12/Val_down 
    S_21 = -Val_S21/Val_down 
    S_22 = 1 - Val_S22/Val_down     
    
    Ana_S.append(np.conjugate([[S_11,S_12],[S_21,S_22]]))

Num_S = np.array(Num_S)
Ana_S = np.array(Ana_S)

# np.save('Num_S_necklace.npy', Num_S)
# np.save('Ana_S_necklace.npy', Ana_S)

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
plt.savefig('Fig_necklace_chain.pdf', bbox_inches='tight')
plt.show()
        
