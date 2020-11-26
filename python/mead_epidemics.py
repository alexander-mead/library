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

# Solves the discrete SIR model
# Si - Initial number of susceptible people
# Ii - Initial number of infected people
# Ri - Initial number of recovered (and therefore immune) people
# R0 - Average number of people an infected person infects
# Ti - Duration of infection [days]
# n - Number of time steps [days]
def solve_discrete_SIR(Si, Ii, Ri, R0, Ti, n):

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

      N = S[i]+I[i]+R[i]
      
      new_infections = int(round(R0*I[i]*S[i]/(Ti*N)))

      if (i-Ti < 0):
         new_recoveries = 0
      else:
         new_recoveries = int(round(R0*I[i-Ti]*S[i-Ti]/(Ti*N)))

      S[i+1] = S[i]-new_infections
      I[i+1] = I[i]+new_infections-new_recoveries
      R[i+1] = R[i]+new_recoveries
            
   # Return lists of values at each time step
   return S, I, R