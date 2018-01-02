## ============= thesis.py ===================================
# This file includes the functions for calculations reltated to 
# drive system effectiveness, mission energy, and battery utilization
#
# Function List:
#       
#
## ===========================================================

def missionrange(eta, mass, omega, torque):
    """This function determines the total mission range capability given 
a certain motor efficiency and mass. The other aspects are as follows:
    - Battery Parameters (except mass)
    - VMTOW
    - Structural/Payload Mass
    - IGBT/MOSFET Parameters
but these aspects are considered system constants and are not changed 
for simplicity. 

INPUTS: 
    eta     NP.MATRIX   Motor efficiency within range 0 -> 1
                        n x m matrix

    mass    FLOAT       Motor mass in kg

    omega   FLOAT       Range of motor rotational speeds
                        n x 1 vector

    torque  FLOAT       Range of motor torques
                        m x 1 vector
   
OUTPUTS:
    range   FLOAT   Mission range in nm"""
    
    # Imports this repository for helper functions
    import matplotlib.pyplot as plt
    import numpy as np
    import thesis as th
    import pandas as pd
    from scipy import integrate

    # System Constants
    NM_TO_KM = 1.85         # [km/NM]
    CRUISE_SPEED = 110      # [NM/h]
    CLIMB_SPEED = 80        # [NM/h]
    DESC_SPEED = 100        # [NM/h]

    VMTOW = 1633            # [kg]
    STRUCT_FRAC = 1/3       # [-]    
    PAYLOAD_FRAC = 1/3      # [-]
    EPPS_FRAC = 1/3         # [-]

    P_TO = 335              # [kW]
    P_CLIMB = 150           # [kW]
    P_CRUISE = 69           # [kW]
    P_DEC = 20              # [kW]
    P_SLOW = 80             # [kW]
    P_STOP = 300            # [kW]
    P_LAND = 300            # [kW]

    T_TO = 37               # [s]
    T_CLIMB = 300           # [s]
    T_DEC = 300             # [s]
    T_SLOW = 60             # [s]
    T_TRAN = 15             # [s] Uses P_CLIMB
    T_STOP = 7              # [s]
    T_LAND = 30             # [s] Uses P_TO
    
    LGMJ1 = pd.read_csv('/Users/Asmodeus/Documents/Danmarks_Tekniske_Universitet/S161568MSC/LGMJ1.csv', header=None,delimiter = ";") # Pandas dataframe
    CELLS_S = 195           # [-]

    # Derived constants
    cruiseS = CRUISE_SPEED * NM_TO_KM   # [km/h]
    climbS = CLIMB_SPEED * NM_TO_KM     # [km/h]
    descS = DESC_SPEED * NM_TO_KM       # [km/h]

    EPPSmass = VMTOW * EPPS_FRAC    # [kg]

    P = [P_TO,P_CLIMB,P_CRUISE,P_DEC,P_SLOW,P_CLIMB,P_STOP,P_LAND]
    T = [T_TO,T_CLIMB, 0,      T_DEC,T_SLOW,T_TRAN, T_STOP,T_LAND]
    S = [0,   climbS, cruiseS, descS,0,     0,      0,      0]
    
    # Find value of non-cruise energy and distance
    time = range(0,sum(T)) 
    power = [0] * (max(time)+1) # one second time step
    speed = [0] * (max(time)+1)

    for x in range(len(P)):
        if x == 0:
            tlast = T[x]
            for y in range(T[x]):
                power[y] = P[x]
                speed[y] = S[x]
        else:
            for y in range(tlast,tlast+T[x]):
                power[y] = P[x]
                speed[y] = S[x]
            tlast += T[x]
     
    energy0 = np.trapz(power, time)/3600    # [kWh] Total energy of non-cruise 
    dist0 = np.trapz(speed, time)/3600      # [km]  Total dist non-cruise
   
    print("The energy usage in the non-cruise portions is {0:0.2f} [kWh]".format(energy0))
    print("                 and the distance traveled is {0:.2f} [km].".format(dist0))

    # Find battery mass
    battMass = EPPSmass - mass  # [kg] total EPPS mass minus motor system mass
    
    

    # Flight Power Plot
    T_CRUISE = 2600         # [s]   Assumed from initial trace, edit 
                            # depending on final plot
    T[2] = T_CRUISE

    time = range(0,60+T_TO+T_CLIMB+T_CRUISE+T_DEC+T_SLOW+T_TRAN+T_STOP+T_LAND)

    print("The total mission duration is", max(time), "[s].")

    power = [0] * (max(time)+1) # one second time step
    speed = [0] * (max(time)+1)

    for x in range(len(P)):
        if x == 0:
            tlast = 30 + T[x]
            for y in range(30,30+T[x]):
                power[y] = P[x]
                speed[y] = S[x]
        else:
            for y in range(tlast,tlast+T[x]):
                power[y] = P[x]
                speed[y] = S[x]
            tlast += T[x]
    
    energy = integrate.cumtrapz(power, time)
    energy = np.append(energy,energy[-1])/3600
    
    distPlot = integrate.cumtrapz(speed, time)
    distPlot = np.append(distPlot, distPlot[-1])/3600
    
    dist = np.trapz(speed,time)/3600

    print("The total distance traveled during the mission is {0:.2f} [km].".format(dist))

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(time,power,'b')
    ax2.plot(time,energy,'r')
    ax1.set_xlabel('Time [s]')
    ax1.set_ylabel('Power [kW]',color='b')
    ax2.set_ylabel('Energy [kWh]',color='r')
    th.nicefig('Mission0.png',title="Mission 0 Baseline")


    # Battery Mass Determinant
    CELLS_P = 24            # [-] CHANGES DEPENDING ON MASS
                            # Important when considering Pmax,disc due to resistance


# End MissionRange() ================================================================


def nicefig(name, w=18, h=10,**kwargs):
    """## ======== Nice figure function ==============================================
# Creates a nice looking figure with 
#       - gridlines 
#       - LaTeX interpreter fonts
#       - Nice labels with equal digits
#       - Adjusts axes limits
#       - Resizes figure to A4 standard
#       - Removes external whitespace when printing figure
# 
# Inputs: 
#       fig             INT     Figure object to modify
#       name            STRING  Name of output file, set type by suffix fx .pdf or .eps
#                                   Supports eps, pdf, png, svg, and a few others.
#                                   See matplotlib.pyplot.savefig() documentation.
#       * w             INT     Width in cm. Default 18
#       * h             INT     Height in cm. Default 10
#       ** title        STRING  Title of plot in latex format
#       ** xlabel       STRING  Xlabel in latex format
#       ** ylabel       STRING  Ylabel in latex format
#       ** filepath     STRING  Filepath for saving output figure
# 
## ============================================================================ 
   """
    import matplotlib.pyplot as plt

    # ========== Optional args handling ============
    # Sets up kwarg default values, or input values 
    # if kwargs exist
    prop = {} # Empty dict 
    for x in ['title', 'xlabel','ylabel','filepath','fignum']:
        if x in kwargs:
            prop[x] = kwargs[x]
        elif x == 'filepath' and x not in kwargs:
            prop['filepath'] = '/Users/Asmodeus/Desktop/'
        elif x != 'filepath' and x not in kwargs:
            prop[x] = 'none'
    # End optional args
    
    print('============== Figure Properties ===============')
    print('Print to: {}'.format(prop['filepath'] + name))
    print('Width: {} cm'.format(w))
    print('Height: {} cm'.format(h))
    print('Title: {}'.format(prop['title']))
    print('X label: {}'.format(prop['xlabel']))
    print('Y label: {}'.format(prop['ylabel']))

    # ========== setup =============
    # change input cm to in for plt use
    IN_TO_CM = 1/2.54
    w *= IN_TO_CM
    h *= IN_TO_CM

    # change fig size
    if prop['fignum'] != 'none':
        figout = plt.figure(prop['fignum'], figsize=(w,h))
    else:
        figout = plt.gcf()
        figout.set_size_inches(w,h)
    
    # Use LaTeX and add figure labels
    # Check if all labels exist then set interpreter
        # if not any(x is 'none' for x in prop.itervalues()):
    plt.rc('text', usetex=True)    
    plt.rc('font', family='serif')

    # Set titles
    if prop['title'] != 'none':
        plt.suptitle(prop['title'])
    if prop['xlabel'] != 'none':
        plt.xlabel(prop['xlabel'])
    if prop['ylabel'] != 'none':
        plt.ylabel(prop['ylabel'])
    
    # add grid
    plt.grid(True)
    
    # show figure
    plt.draw()
    plt.show()
    
    # save figure
    nameandpath = prop['filepath'] + name
    figout.savefig(nameandpath, bbox_inches='tight')
# End function














