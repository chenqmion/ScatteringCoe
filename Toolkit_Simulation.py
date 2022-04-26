import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const

alph= 5E-3
v_ph = 1.35E8
Z0 = 50

def Waveguide_Z(w,L, alph=alph,v_ph=v_ph,Z0=Z0,load='lamb4'):
    beta = w/v_ph
    gamm = alph + 1j*beta
    
    if (load=='lamb4'):
        Zr = Z0*np.tanh(gamm*L)
    elif (load=='lamb2'):
        Zr = Z0/np.tanh(gamm*L)
    else:
        print('geometry error')
    return Zr

#%% ABCD matrices
def TM_Hori(Z):
    return np.array([[1.0,Z],[0.0,1.0]])

def TM_Vert(Z):
    return np.array([[1.0,0.0],[1/Z,1.0]])

def TM_Waveguide(w,L, alph=alph,v_ph=v_ph,Z0=Z0):
    beta = w/v_ph
    gamm = alph + 1j*beta
    
    A = np.cosh(gamm*L)
    B = Z0*np.sinh(gamm*L)
    C = (1/Z0)*np.sinh(gamm*L)
    D = np.cosh(gamm*L)
    return np.array([[A,B],[C,D]])

def Num_SC(M, Z0=Z0):
    A, B = M[0]
    C, D = M[1]
    
    S_11 = (A + B/Z0 - C*Z0 - D)/(A + B/Z0 + C*Z0 + D)
    S_12 = 2.0*(A*D-B*C)/(A + B/Z0 + C*Z0 + D)
    S_21 = 2.0/(A + B/Z0 + C*Z0 + D)
    S_22 = (-A + B/Z0 - C*Z0 + D)/(A + B/Z0 + C*Z0 + D)
    return np.array([[S_11,S_12],[S_21,S_22]])

#%% Hanger
def Ana_Params_Hanger(L,C_1, alph=alph,v_ph=v_ph,Z0=Z0, load='lamb4'):
    if (load=='lamb4'):
        w_0 = v_ph/(4*L) * (2*np.pi)
        w_r = w_0 - 2*(w_0**2)*Z0*C_1/np.pi
        # w_r = w_0 - 2*(w_0*w_r)*Z0*C_1/np.pi
        
        Qi = np.pi/(4*alph*L)
        Qe1 = np.pi/(2*(w_r**2)*(Z0**2)*(C_1**2))
        Ql = Qe1*Qi/(Qe1+Qi)
    elif (load=='lamb2'):
        w_0 = v_ph/(2*L) * (2*np.pi)
        w_r = w_0 - (w_0**2)*Z0*C_1/np.pi
        # w_r = w_0 - (w_0*w_r)*Z0*C_1/np.pi
        
        Qi = np.pi/(2*alph*L)
        Qe1 = np.pi/((w_r**2)*(Z0**2)*(C_1**2))
        Ql = Qe1*Qi/(Qe1+Qi)
    else:
        print('geometry error')
    return np.array([w_0, w_r, Qi, Qe1, Ql])

def Ana_SC_Hanger(w, params):
    w_0, w_r, Qi, Qe1, Ql = params
    delt = w/w_r - 1 
    
    S_11 = -Ql/Qe1/(1+2j*Ql*delt)
    S_12 = 1-Ql/Qe1/(1+2j*Ql*delt)
    S_21 = 1-Ql/Qe1/(1+2j*Ql*delt)
    S_22 = -Ql/Qe1/(1+2j*Ql*delt)
    return np.array([[S_11,S_12],[S_21,S_22]])

#%% Necklace
def Ana_Params_Necklace(L,C_1,C_2, alph=alph,v_ph=v_ph,Z0=Z0, load='lamb4'):
    if (load=='lamb4'):
        w_0 = v_ph/(4*L) * (2*np.pi)
        w_r = w_0 - 2*(w_0**2)*Z0*C_1/np.pi
        # w_r = w_0 - 2*(w_0*w_r)*Z0*C_1/np.pi
        
        Qi = np.pi/(4*alph*L)
        Qe1 = np.pi/(4*(w_r**2)*(Z0**2)*(C_1**2))
        Qe2 = 0
        Ql = Qe1*Qi/(Qe1+Qi)
    elif (load=='lamb2'):
        w_0 = v_ph/(2*L) * (2*np.pi)
        w_r = w_0 - (w_0**2)*Z0*(C_1+C_2)/np.pi
        # w_r = w_0 - (w_0*w_r)*Z0*(C_1+C_2)/np.pi
        
        Qi = np.pi/(2*alph*L)
        Qe1 = np.pi/(2*(w_r**2)*(Z0**2)*(C_1**2))
        Qe2 = np.pi/(2*(w_r**2)*(Z0**2)*(C_2**2))
        Qe = Qe1*Qe2/(Qe1+Qe2)
        Ql = Qe*Qi/(Qe+Qi)
    else:
        print('geometry error')
        
    return np.array([w_0, w_r, Qi, Qe1, Qe2, Ql])

def Ana_SC_Necklace(w, params):
    w_0, w_r, Qi, Qe1, Qe2, Ql = params
    delt = w/w_r - 1 
    
    S_11 = 1-2*Ql/Qe1/(1+2j*Ql*delt)
    S_12 = 2*Ql/np.sqrt(Qe1*Qe2)/(1+2j*Ql*delt)
    S_21 = 2*Ql/np.sqrt(Qe1*Qe2)/(1+2j*Ql*delt)
    S_22 = 1-2*Ql/Qe2/(1+2j*Ql*delt)
    return np.array([[S_11,S_12],[S_21,S_22]])

#%% Cross
def Ana_Params_Cross(L,C_1,C_2, alph=alph,v_ph=v_ph,Z0=Z0, load='lamb4'):
    if (load=='lamb4'):
        w_0 = v_ph/(4*L) * (2*np.pi)
        w_r = w_0 - 2*(w_0**2)*Z0*(C_1+C_2)/np.pi
        # w_r = w_0 - 2*(w_0*w_r)*Z0*(C_1+C_2)/np.pi
        
        Qi = np.pi/(4*alph*L)
        Qe1 = np.pi/(4*(w_r**2)*(Z0**2)*(C_1**2))
        Qe2 = np.pi/(4*(w_r**2)*(Z0**2)*(C_2**2))
        Qe = Qe1*Qe2/(Qe1+Qe2)
        Ql = Qe*Qi/(Qe+Qi)
    elif (load=='lamb2'):
        w_0 = v_ph/(2*L) * (2*np.pi)
        w_r = w_0 - (w_0**2)*Z0*(C_1+C_2)/np.pi
        # w_r = w_0 - (w_0*w_r)*Z0*(C_1+C_2)/np.pi
        
        Qi = np.pi/(2*alph*L)
        Qe1 = np.pi/(2*(w_r**2)*(Z0**2)*(C_1**2))
        Qe2 = np.pi/(2*(w_r**2)*(Z0**2)*(C_2**2))
        Qe = Qe1*Qe2/(Qe1+Qe2)
        Ql = Qe*Qi/(Qe+Qi)
    else:
        print('geometry error')
        
    return np.array([w_0, w_r, Qi, Qe1, Qe2, Ql])

def Ana_SC_Cross(w, params):
    w_0, w_r, Qi, Qe1, Qe2, Ql = params
    delt = w/w_r - 1 
    
    S_11 = 1-2*Ql/Qe1/(1+2j*Ql*delt)
    S_12 = -2*Ql/np.sqrt(Qe1*Qe2)/(1+2j*Ql*delt)
    S_21 = -2*Ql/np.sqrt(Qe1*Qe2)/(1+2j*Ql*delt)
    S_22 = 1-2*Ql/Qe2/(1+2j*Ql*delt)
    return np.array([[S_11,S_12],[S_21,S_22]])
