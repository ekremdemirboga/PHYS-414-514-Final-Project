import numpy as np 
from scipy.integrate import solve_ivp

def TOV_solver(rho_c,k=50):
    
    
    def EOS(p, K = k, Gamma=2):
        rho = (p/K)**(1/Gamma)
        return rho
 
    def TOV(r, y):
        m = y[0]
        v = y[1]
        p = y[2]
        baryonic_m = y[3]
        dmdr = 4*np.pi*(r**2)*EOS(p)
        dvdr = (2*(m + 4*np.pi*r**3*p))/(r*(r-2*m))
        dpdr = -0.5*(p+EOS(p))*dvdr
        dbmdr = 4*np.pi*(1 - 2*m/r)**(-1/2)*r**2*EOS(p)
        return [dmdr, dvdr, dpdr, dbmdr]
    
    def stop_condition(r, y):
        return y[2]
    stop_condition.terminal = True
    
    dt = 1e-10
    p_c = rho_c**2*k
    tspan = np.linspace(0+dt,100,1000)  
    yinit = [0,0,p_c,0]    
    sol = solve_ivp(lambda r, y: TOV(r, y),[tspan[0], tspan[-1]], yinit,method='RK45',events=stop_condition)
    
    
    return sol
