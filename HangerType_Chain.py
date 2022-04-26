import numpy as np
import matplotlib.pyplot as plt
import Toolkit_Simulation as tk

L = 5E-3
C_1 = 1e-14
load = 'lamb4'
w_list = np.linspace(6.637,6.683,10001) * 1E9 * (2*np.pi)

thet_wg = np.pi/3
N_wg = 4

Num_S = []
Ana_S = []
for w in w_list:
    #%% Numerical
    Zr = tk.Waveguide_Z(w,L, load=load)
    Z = 1/(1j*w*C_1) + Zr
    Mr = tk.TM_Vert(Z)
    
    M_wg = tk.TM_Waveguide(w,(thet_wg/2/np.pi)*(4*L))
    M_chain = Mr
    for num_r in range(0,N_wg-1):
        M_chain = M_chain @ M_wg @ Mr
    
    Num_S.append(tk.Num_SC(M_chain))
    
    #%% analytical calculations
    params = tk.Ana_Params_Hanger(L,C_1, load=load)
    w_0, w_r, Qi, Qe1, Ql = params
    gamm_a = w_r/Qi
    gamm = w_r/(2*Qe1)
    
    # form 1
    A = np.eye(N_wg) * (1j*(w_r-w) + gamm_a/2)
    b_lin = np.zeros(N_wg, dtype=complex)
    b_rin = np.zeros(N_wg, dtype=complex)
    for num1 in range(1,N_wg+1):
        b_lin[num1-1] = -np.sqrt(gamm) * np.exp(-1j*(num1-N_wg)*thet_wg)
        b_rin[num1-1] = -np.sqrt(gamm) * np.exp(1j*(num1-1)*thet_wg)
        for num2 in range(1,N_wg+1):
            A[num1-1,num2-1] += gamm * np.exp(1j*np.abs(num1-num2)*thet_wg)
            
    a_lin = np.linalg.solve(A, b_lin)
    a_rin = np.linalg.solve(A, b_rin)
    
    S_11 = 0+0j
    S_21 = np.exp(1j*(N_wg-1)*thet_wg)
    S_12 = np.exp(1j*(N_wg-1)*thet_wg)
    S_22 = 0+0j
    for num1 in range(1,N_wg+1):
        S_11 += np.sqrt(gamm)*np.exp(1j*(num1-1)*thet_wg)*a_rin[num1-1]
        S_12 += np.sqrt(gamm)*np.exp(1j*(num1-1)*thet_wg)*a_lin[num1-1]
        S_21 += np.sqrt(gamm)*np.exp(1j*(N_wg-num1)*thet_wg)*a_rin[num1-1]
        S_22 += np.sqrt(gamm)*np.exp(1j*(N_wg-num1)*thet_wg)*a_lin[num1-1]
        
    # # form 2
    # x = gamm/(1j*(w_r-w) - gamm + gamm_a/2)
    # A = np.zeros((N_wg,N_wg), dtype=complex)
    # b_lin = np.zeros(N_wg, dtype=complex)
    # b_rin = np.zeros(N_wg, dtype=complex)
    # for num1 in range(1,N_wg+1):
    #     if (num1 == 1):
    #         A[0,0] = -(1+2*x)/x
    #         A[0,1] = (1+x)/x
    #         b_lin[0] = np.exp(1j*N_wg*thet_wg)
    #         b_rin[0] = np.exp(1j*thet_wg)
    #     elif (num1 == N_wg):
    #         A[N_wg-1,N_wg-2] = (1+x)*np.exp(2j*thet_wg)/x
    #         A[N_wg-1,N_wg-1] = -(1+2*x+np.exp(2j*thet_wg))/x
    #         b_lin[N_wg-1] = -np.exp(1j*N_wg*thet_wg)*(np.exp(2j*thet_wg)-1)
    #     else:
    #         A[num1-1,num1-2] = (1+x)*np.exp(2j*thet_wg)/x
    #         A[num1-1,num1-1] = -(1+2*x+np.exp(2j*thet_wg))/x
    #         A[num1-1,num1] = (1+x)/x
    #         b_lin[num1-1] = -np.exp(1j*N_wg*thet_wg)*(np.exp(2j*thet_wg)-1)
            
    # c_lin = np.linalg.solve(A, b_lin)
    # c_rin = np.linalg.solve(A, b_rin)
    
    # A = np.zeros((N_wg,N_wg), dtype=complex)
    # b_lin = np.zeros(N_wg, dtype=complex)
    # b_rin = np.zeros(N_wg, dtype=complex)
    # for num1 in range(1,N_wg+1):
    #     if (num1 == 1):
    #         A[0,0] = -(1+2*x+np.exp(2j*thet_wg))/x
    #         A[0,1] = (1+x)*np.exp(2j*thet_wg)/x
    #         b_rin[0] = -(np.exp(1j*thet_wg)-np.exp(-1j*thet_wg))
    #     elif (num1 == N_wg):
    #         A[N_wg-1,N_wg-2] = (1+x)/x
    #         A[N_wg-1,N_wg-1] = -(1+2*x)/x
    #         b_lin[N_wg-1] = np.exp(-1j*N_wg*thet_wg)
    #         b_rin[N_wg-1] = np.exp(-1j*thet_wg)
    #     else:
    #         A[num1-1,num1-2] = (1+x)/x
    #         A[num1-1,num1-1] = -(1+2*x+np.exp(2j*thet_wg))/x
    #         A[num1-1,num1] = (1+x)*np.exp(2j*thet_wg)/x
    #         b_rin[num1-1] = -(np.exp(1j*thet_wg)-np.exp(-1j*thet_wg))
            
    # d_lin = np.linalg.solve(A, b_lin)
    # d_rin = np.linalg.solve(A, b_rin)
    
    # S_11 = np.exp(-1j*thet_wg)*c_rin[0]
    # S_21 = np.exp(1j*(N_wg-1)*thet_wg) + np.exp(-1j*thet_wg)*c_lin[0]
    # S_12 = np.exp(1j*(N_wg-1)*thet_wg) + np.exp(1j*N_wg*thet_wg)*d_rin[N_wg-1]
    # S_22 = np.exp(1j*N_wg*thet_wg)*d_lin[N_wg-1]

    Ana_S.append(np.conjugate([[S_11,S_12],[S_21,S_22]]))

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
plt.savefig('Fig_hanger_chain.pdf', bbox_inches='tight')
plt.show()
        

