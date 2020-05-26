# Set of right-hand-side equations
def SIR_equations(t, y, ts, R0s, gamma):
    
   # Unpack y variable
   S = y[0]
   I = y[1]
   R = y[2]
   
   N = S + I + R
   
   if (t < ts[0]):
      R0 = R0s[0]
   elif (t < ts[1]):
      R0 = R0s[1]
   else:
      R0 = R0s[2]
   
   # Equations (note that gamma factorises out; it has units of 1/time)
   Sdot = -gamma*R0*I*S/N
   Idot = gamma*R0*I*S/N - gamma*I
   Rdot = gamma*I
   
   return [Sdot, Idot, Rdot]

# Routine to solve the SIR equations
def solve_SIR(S_initial, I_initial, R_initial, ts, R0s, gamma):

   import numpy as np
   from scipy.integrate import solve_ivp
    
   # Time span
   ti = 0.
   tf = 30.*gamma
   nt = 1024
   t = np.linspace(ti, tf, nt)
   
   # Solve ODE system
   solution = solve_ivp(SIR_equations, (ti, tf), (S_initial, I_initial, R_initial), 
                        method='RK45', 
                        t_eval=t,
                        args=(ts, R0s, gamma) 
                     )

   # Break solution down into SIR
   S = solution.y[0]
   I = solution.y[1]
   R = solution.y[2]
   
   return t, S, I, R

# Takes S(n), I(n), and R(n) and updates to S(n+1), I(n+1), and R(n+1)
def iterate_discrete_SIR_system(S, I, R, R0):
   N = S+I+R
   In = int(round(R0*I*S/N))
   Rn = I+R
   Sn = N-In-Rn
   if(Sn < 0):
      In = In+Sn
      Sn = 0
   return Sn, In, Rn

# Solves the discrete SIR model
# Si - Initial number of susceptible people
# Ii - Initial number of infected people
# Ri - Initial number of recovered (and therefore immune) people
# n - Number of time steps
def solve_discrete_SIR(Si, Ii, Ri, R0, n):

   import numpy as np
    
   # Empty lists for solution
   S = np.zeros(n+1, dtype=int)
   I = np.zeros(n+1, dtype=int)
   R = np.zeros(n+1, dtype=int)
   
   # Empty lists for solution
   S[0] = Si
   I[0] = Ii
   R[0] = Ri
   
   # Loop over steps and update
   for i in range(n):
      S[i+1], I[i+1], R[i+1] = iterate_discrete_SIR_system(S[i], I[i], R[i], R0)
            
   # Return lists of values at each time step
   return S, I, R