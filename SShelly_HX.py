'''
Date of last edit: June 14th, 2020
Author(s): Ryan McGuire*    Lane Carasik^
*Virginia CommonWealth University
*FAST Research Group
Equations for Shell-Side Reactor Heat Exchangers: 
Python Script for calculation of Shell-Side heat transfer coefficient using correctional factors 
for Shell-and-Tube bundles based on models found in literature.
'''

'''
Revision Points
1) A_o_cr
2) A_o_tb
3) N_b
'''

##Basic Imports and Function Representatives

import math
pi = math.pi
e = math.e
acos = math.acos
sin = math.sin

##Measured Information used in calculations

#Tube Material admiralty (70% Cu 30% Ni)
D_s = 0.336            #Shell-side inside diameter (m)
d_o = 0.019            #Tube-side outer diameter (m)
d_i = 0.0166           #Tube-side inside diameter (m)
d_l = 0.019749         #Tube whole diameter (m)
p_t = 0.020            #Tube pitch (m)
Tbl = 45               #Tube bundle Layout (degrees)
L_bc = 0.279           #Central Baffle Spacing (m)
L_bi = 0.318           #Inlet Baffle Spacing (m)
L_bo = 0.318           #Outlet baffle spacing (m)
l_c = 0.0867           #Baffle Cut (m) or 25.8%
N_ss = 1               #Number of sealing strip pairs
D_baffle = 0.333054    #Baffle distance (m)
N_t = 102              #Total number of tubes
L = 4.3                #Tube lenght (m)
w_p = 0.019            #Width of bypass lane (m)
n_p = 2                #Number of tube passes
N_rcc = 9              #Number of effective rube rows crossed durring flow
N_p = 2                #Number of pass partitions
D_otl = 0.321          #Diameter of the outer tube limit
d_tb = 0.000794        #Tube-to-baffle hole diametral clearance (m)
d_sb = 0.002946        #Shell-to-baffle diametral clearance (m)
k_w = 111              #Thermal conductivity of tube wall (W/m * K)
m_s = 36.3             #Oil Flow rate (kg/s)
T_si = 65.6            #Oil inlet temperature (degrees C)
R_of = 0.000176        #Oil side fouling factor (m^2*W/K)
X_t = 0.0354           #Tranverse tube pitch (m)
R_if = 0.000088        #Water side fouling factor (m^2*W/K)
m_t = 18.1             #Water flow rate (kg/s)
T_ti = 32.2            #Water inlet temperature (degrees C)
k = 0.140              #Thermal cunductivity for fluid
m_s = 36.3             #Fluid mass flow rate (kg/s)
u_s = 0.00646          #Fluid dynamic viscosity
Laminar_flow = 0       #Boolean to see if the flow is laminar, 1 is Laminar, 0 is Non-Laminar
Nu_s = 125             #Nusselt number
Bundle_layout = 30     #Tube bundle Layout degree
##Correctional Factor Calculations

#Baffle configuration correctional factor
D_ctl = D_otl-d_o
theta_ctl = 2*acos((D_s-2*l_c)/D_ctl)
F_c = 1-(theta_ctl/pi)+(sin(theta_ctl)/pi)
J_c = round((0.55+0.72*F_c),5)
#print(J_c)

#Bundle Leakage effects correctional factor
theta_b = 2*acos(1-(2*l_c/D_s))
F_w = (theta_ctl/(2*pi))-(sin(theta_ctl)/(2*pi))
d_sb = D_s-D_baffle
d_tb = d_l-d_o
A_o_sb = pi*D_s*(d_sb/2)*(1-(theta_b/(2*pi)))
A_o_tb = (pi*d_o*d_tb*N_t*(1-F_w))/2
A_o_cr = 0.03275
'''
if Bundle_layout in (30,90):
    A_o_cr = (D_s-D_otl+(D_ctl/X_t)*(X_t-d_o))*L_bc
elif Bundle_layout in (45,60):
    print("no1")
elif Bundle_layout in (30,90):
    print("no2")
else Bundle_layout in (30,90):
    print("no3")
'''
r_s = (A_o_sb)/(A_o_sb+A_o_tb)
r_lm = (A_o_sb+A_o_tb)/(A_o_cr)
J_l = round((0.44*(1-r_s)+(1-0.44*(1-r_s))*math.e**(-2.2*r_lm)),5)
#print(J_l)

#Bundle and pass partition bypass correctional factor
G_s = (m_s/A_o_cr)
A_o_bp = L_bc*(D_s-D_otl+0.5*N_p*w_p)
Re_s = (G_s*d_o)/u_s
N_ss_plus = N_ss/N_rcc
r_b = A_o_bp/A_o_cr
if Re_s<=100:
    C = 1.35
else:
    C = 1.25
if N_ss_plus >= 1/2:
    J_b = 1
else:
    J_b = round((e**(-C*r_b*(1-(2*N_ss_plus)**(1/3)))),5)
#print(J_b)

#Larger baffle spacing correctional factor
if Laminar_flow == 1:
    n = 1/3
else:
    n = 3/5
L_i_plus = L_bi/L_bc
L_o_plus = L_bo/L_bc
N_b = ((L-L_bi-L_bo)/L_bc)+1
J_s = round(((N_b-1+L_i_plus**(1-n)+L_o_plus**(1-n))/(N_b-1+L_i_plus+L_o_plus)),5)
#print(J_s)

#Adverse Temperature Gradient Buildup in Laminar Flow correctional factor
N_rcw = (0.8/X_t)*(l_c-(1/2)*(D_s-D_ctl))
N_rc = N_rcc+N_rcw
if Re_s >= 100:
    J_r = 1
elif Re_s <= 20:
    J_r = (10/N_rc)**0.18
else:
    J_r = 1+(((10/N_rc)**0.18)-1)/(20-100)*(Re_s-100)
#print(J_r)

##Heat Transfer Calculations - Core

#Perfect Heat Transfer Coefficent 
h_id = round((((Nu_s*k)/d_o)*(1)**-0.14),3)
#print(h_id)

#Calculated Heating Coefficent
h_s = round((h_id*J_c*J_l*J_b*J_s*J_r),3)
print("The Calculated Shell Side Heat Transfer Coefficient is " +str(h_s)+ " W/m^2 * K")
