"""# ============= thesis.py ===================================
This file includes the functions for calculations reltated to 
drive system effectiveness, mission energy, and battery utilization

#Class List:
    mission()       Contains mission constants and some derived values
    cell(battMass)  Contains cell consants and some derived pack values

#Function List:
    missionrange(eta,mass,omega,torque,name='Mission0.png',title='Baseline Power Profile') 
                                Plots a flight profile for one motor design
    heatrange(eta,mass,name='heatrange.png',title='Total Range vs. Motor Mass, \eta')
                                Plots a heat plot of maximum range for many motor designs
    nicefig(name, w,h,**kwargs)         Cleans up plot 
# ==========================================================="""

def missionrange(M, name='Mission0.png', title='Baseline Power Profile'):
    """## =========== Mission Range function ======================================
This function determines the total mission range capability given 
a certain motor efficiency and mass. The other aspects are as follows:
    - Battery Parameters (except mass)
    - VMTOW
    - Structural/Payload Mass
    - IGBT/MOSFET Parameters
but these aspects are considered system constants and are not changed 
for simplicity. 

INPUTS: 
    M       MISSION Class with info from battrange(9
OUTPUTS:
    range   FLOAT   Mission range in nm
## =============================================================================="""
    
    # Imports this repository for helper functions
    import matplotlib.pyplot as plt
    import numpy as np
    import thesis as th
    import pandas as pd
    from scipy import integrate

    # Find value of non-cruise energy and distance
    energy0 = np.trapz(M.power, M.time)/3600    # [kWh] Total energy of non-cruise 
    dist0 = np.trapz(M.speed, M.time)/3600      # [km]  Total dist non-cruise
   
    print("The energy usage in the non-cruise portions is {0:0.2f} [kWh]".format(energy0))
    print("                 and the distance traveled is {0:.2f} [km].".format(dist0))

    time = range(0,60+sum(M.T))

    print("The total mission duration is", max(time), "[s].")

    power = [0] * (max(time)+1) # one second time step
    speed = [0] * (max(time)+1)

    for x in range(len(M.P)):
        if x == 0:
            tlast = 30 + M.T[x]
            for y in range(30,30+M.T[x]):
                power[y] = M.P[x]
                speed[y] = M.S[x]
        else:
            for y in range(tlast,tlast+M.T[x]):
                power[y] = M.P[x]
                speed[y] = M.S[x]
            tlast += M.T[x]
    
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
    th.nicefig(name,title)

    # Battery Mass Determinant
    CELLS_P = 24            # [-] CHANGES DEPENDING ON MASS
                            # Important when considering Pmax,disc due to resistance


# End MissionRange() ================================================================

def  heatrange(eta,mass,name='heatrange.png',title='Total Range vs. Motor Mass, \eta'):
    """## ======== Heat range function ==============================================
Creates a 3D plot to display variations in range given many possible motor designs, with axes:
        x - motor mass
        y - motor efficiency
        z - aircraft range
 
#Inputs:
        eta         LIST    List of efficiencies. MUST CORRESPOND WITH mass
        mass        LIST    List of motor masses.
        name        STRING  Name to save output figure
        title       STRING  Title of plot
## ==========================================================================="""


# End heatrange() =================================================================

def battrange(eta,mass,omega,torque,**kwargs):
    """## ======== Battery Range function ========================================
Function that takes known battery params and simulates the battery response given a
specific motor efficiency. Assumes a full battery and ideal mission conditions.

Uses th.mission and th.cell classes. Reference th.mission.__doc__ or th.cell.__doc__
for more information.

#Inputs: 
    eta     NP.MATRIX   n by m 2D matrix of efficiency data for the motor. MUST CORRESPOND TO omega AND torque.
                                    In range of (0,1)
    mass    FLOAT       Motor system mass, used to calculate battery mass
    omega   LIST        n long list of motor output rotational speed [rad/s]
                        Solely used for plotting, not relevant for power calculations
    torque  LIST        m long list of motor output torques [Nm]
                        Solely used for plotting, not relevant for power calculations
    **filename STRING   Optional output parameter to change filename, set type by 
                                    suffix. Supports eps, pdf, png, svg, and a few 
                                    others. See matplotlib.pyplot.savefig() 
                                    documentation.
#Outputs:
    range   FLOAT       Maximum theoretical range of the system

## ========================================================================"""

    import numpy as np
    import pandas as pd
    from scipy import integrate
    import scipy.interpolate
    import thesis as th
    import matplotlib.pyplot as plt

    ## System contants, see th.mission.__doc__ for detailed info
    M = th.mission()

    # Battery mass from motor system mass and total EPPS mass
    battMass = M.EPPSmass - mass
    
    # Battery info, see th.cell.__doc__ for detailed info
    C = th.cell(battMass)
    
    # interpolate the eta matrix for use in sim
    etaInterp = scipy.interpolate.interp2d(omega,torque,eta)

    ## Create simulation
    # initialize values to hold sim outputs
    S = th.simData(C)

    # Set cruise length flag
    landingPower = True

    while (S.t < 30 + M.T[0]): # do TO (hover)
        S = sim(M,C,S,etaInterp,0)         
        Pmax = (C.Voc(S.SOC[S.t])-C.Vmin)/C.Res(S.SOC[S.t])*C.Vmin/1000 # [kW]
    
    while (S.t < (30 + sum(M.T[0:2]))): # do phase 2 (climb)
        S = sim(M,C,S,etaInterp,1)
        Pmax = (C.Voc(S.SOC[S.t])-C.Vmin)/C.Res(S.SOC[S.t])*C.Vmin/1000 # [kW]

    # while cruise THIS ONE IS VARIABLE LENGTH
    while landingPower is True:
        S = sim(M,C,S,etaInterp,2) # returns tCruise
        tCruise = S.t - 30 - sum(M.T[0:2])
        # Calculate current max possible power draw and set flag if under requirements
        Pmax = (C.Voc(S.SOC[S.t])-C.Vmin)/C.Res(S.SOC[S.t])*C.Vmin/1000 # [kW]
        if Pmax < M.power[2]:
            landingPower = False
            S.removeLast()  # Removes elements that would place landing reqs over max batt output     
    print("Cruise duration is",tCruise,"seconds")
    print("Current time is", S.t)

    # while desc,slow,stop,tran,land
    for x in range(3,8):
        while S.t < (30 + sum(M.T[0:x+1]) + tCruise):
            S = sim(M,C,S,etaInterp,x)
    
    S.finalValues()
    
    print("Simulation complete!")
    
    if "filename" in kwargs:
        S.plot(C,kwargs["filename"])
    else:
        S.plot(C,"SimulationResults.png")
    return S
# End battrange() ====================================================================

def sim(M,cell,S,eta,phase):
    import numpy as np
    
    battPower = S.battPower
    SOC = S.SOC
    current = S.current
    t = S.t

    S.battPower.append(M.P[phase] / eta(M.omega[phase],M.torque[phase])) # may need to swap index order
    
    S.current.append((cell.Voc(SOC[t]) - np.sqrt(max(0.01, cell.Voc(SOC[t])**2-4*cell.Res(SOC[t])*battPower[t]*1000)))/(2*cell.Res(SOC[t])))
    
    S.Vterm.append(cell.Voc(SOC[t]) - S.current[t] * cell.Res(SOC[t]))
    S.SOC.append(max(0, min(1,SOC[t] - np.trapz([float(current[t-1]),float(current[t])],np.divide([t-1,t],3600))/cell.packCap)))
          
    S.t += 1
    S.time.append(t)
    
    return S

# End sim() ==========================================================================
    

class simData():
    """##============ simData ==========================================================
Class that contains simulation results. Used as single argument to sim function.

#FIELDS:
    battpower   LIST    contains the power drawn from the battery
    current     LIST    contains the current drawn from the battery
    SOC         LIST    contains the SOC of the battery
    t           INT     contains the present time of the simulation
    time        LIST    contains the overall time vector
    Vterm       LIST    contains the terminal voltage of the battery

#METHODS:
    removeLast()        removes the last element of all LIST fields
                        and updates t to the last timestep
    plot(cell)          plots simulation data
    finalValues()       appends final values to complete plot       
"""
    def __init__(self,cell):
        self.battPower = [0] * 60
        self.current = [0] * 60 
        self.SOC = [1] * 60
        self.t = 59
        self.time = list(range(60))
        self.Vterm = [cell.Voc(1)] * 60

    def removeLast(self):
        """##========= removeLast Method ===============================================
Method that resets to the previous time step."""
        self.battPower = self.battPower[:-1]
        self.current = self.current[:-1]
        self.SOC = self.SOC[:-1]
        self.time = self.time[:-1]
        self.t -= 1
        self.Vterm = self.Vterm[:-1]

    def plot(self,cell,filename):
        """##=========== simPlot Method =================================================
Method that plots the contents of the simData class. 

#INPUTS:
    cell    OBJECT  class object that contains the cell parameters
"""
        import matplotlib.pyplot as plt
        import thesis as th
        from matplotlib.offsetbox import TextArea, VPacker, AnchoredOffsetbox, HPacker, VPacker

        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        l1 = ax1.plot(self.time,self.battPower,'b',label="P_{batt}")
        l2 = ax1.plot(self.time,cell.Voc(self.SOC),'g:',label="V_{oc}")
        l3 = ax1.plot(self.time,self.Vterm,'g',label="V_{term}")
        l4 = ax2.plot(self.time,self.SOC,'r',label="SOC")
        ax2.set_ylim([0,1])
        ax1.set_xlim([0,max(self.time)])
        ax1.set_xlabel('Time [s]')
        
        # Draw multicolored axis label
        ybox1 = TextArea("Battery Power [kW]", textprops=dict(color='b',size=10, rotation=90,ha='left',va='bottom'))
        ybox2 = TextArea("Voltage [V]", textprops=dict(color='g',size=10, rotation=90,ha='left',va='bottom'))
        ybox = VPacker(children=[ybox1,ybox2],align="right",pad=5,sep=10)
        anchored_ybox = AnchoredOffsetbox(loc=8,child=ybox,pad=0.,frameon=False,bbox_to_anchor=(-0.08,0.2),bbox_transform=ax1.transAxes, borderpad=0.)
        ax1.add_artist(anchored_ybox)   

        ax2.set_ylabel('State of Charge [-]',color='r')

        lines = l1 + l2 + l3 + l4
        labels = [l.get_label() for l in lines]
        plt.legend(lines, labels, loc=8)

        th.nicefig(filename,title="Mission Simulation Results") 
        
    def finalValues(self):
        for x in range(60):
            self.battPower.append(0)
            self.current.append(0)
            self.SOC.append(self.SOC[-1])
            self.time.append(max(self.time)+1)
            self.Vterm.append(self.Vterm[-1])

## End sim class ========================================================================
class cell():
    """#=========== cell(battMass = 250) class ===============================================
Class that contains the constant values for the battery cell, and some derived values.
Default battery mass of 250 [kg] corresponds to 2 packs of 195 series by 13 parallel cells.

#FILEDS:
    #Constants:
        LGMJ1       Pandas dataframe containing cell data. Print dataframe for more info
        CELLS_S     number of cells in series                       [-]
        cellMass    mass of a single cell                           [kg]
        SOCinit     LIST of state of chagres for interpolation      [-]
        cellRes     LIST of cell resistances at corresponding SOC   [Ohm]
        cellVoc     LIST of cell open-circuit voltages at SOC       [V]
        cellCap     cell capacity                                   [Ah]
        Vmin        Minimum pulse voltage                           [V]

    #Derived Values:
        noCells     number of cells in pack                         [-]
        cellsP      number of parallel cells in pack                [-]
        packCap     pack capacity                                   [Ah]
        packNomEnergy   nominal pack energy (at Vnom)               [kWh]
        Res         Scipy interpolation. Pack discharge resistance. [Ohm]
                    Call this as cell.Res(SOC)
        Voc         Scipy interpolation. Pack open circuit voltage. [V]
                    Call this as cell.Voc(SOC)
"""
    def __init__(self,battMass = 250): # default of 250kg. 195s13p pack x2 
        import pandas as pd
        import numpy as np
        import scipy.interpolate
        # import data
        self.LGMJ1 = pd.read_csv('/Users/Asmodeus/Documents/Danmarks_Tekniske_Universitet/S161568MSC/LGMJ1.csv', header=None,delimiter = ";") # Pandas dataframe
        self.CELLS_S = 195              # [-]
        self.cellMass = float(self.LGMJ1[1][7].replace(',','.')) # [kg] Gets mass from dataframe and changes to decimal radix

        # Approximate cells in parallel
        self.noCells = battMass / self.cellMass   # Wont be INT
        self.cellsP = self.noCells / self.CELLS_S

        # Calculate capacity
        self.cellCap = float(self.LGMJ1[1][11].replace(',','.')) # [Ah] Gets cell capacity and changes radix 
        self.packCap = self.cellsP * self.cellCap # [Ah]
        self.packNomEnergy = self.CELLS_S*float(self.LGMJ1[1][12].replace(',','.'))*self.packCap/1000 # [kWh] NOMINAL
    
        ## Interpolate Resistance and Voc 
        self.SOCinit = [1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0]    # initial data is every 10% 

        # Preallocate
        self.cellRes = [0]*len(self.SOCinit)
        self.cellVoc = [0]*len(self.SOCinit)

        for x in range(len(self.SOCinit)):
            self.cellVoc[x] = float(self.LGMJ1[1][28-x].replace(',','.')) # This flips the data, since in DFit is 0-100
            self.cellRes[x] = float(self.LGMJ1[1][39-x].replace(',','.'))

        # Pack values
        self.Voc = np.multiply(self.cellVoc, self.CELLS_S)
        self.Res = np.divide(np.multiply(self.cellRes, self.CELLS_S),self.cellsP)
        #  interpolate
        self.Res = scipy.interpolate.interp1d(self.SOCinit,self.Res)   # call this as Res(SOC)
        self.Voc = scipy.interpolate.interp1d(self.SOCinit,self.Voc)

        # pack minimum
        self.Vmin = 2.6 * self.CELLS_S  # [V]

class mission():
    """#============ mission() class ====================================================
Class that contains the constant values of the mission, and some derived values.

#FIELDS:
    #Constants:
    NM_TO_CRUISE    conversion factor from NMi to km            [km/NMi]
    CRUISE_SPEED    groundspeed of the aircraft at cruise       [NMi/h]
    CLIMB_SPEED     groundspeed of the aircraft in climb        [NMi/h]
    DESC_SPEED      groundspeed of the aircraft in desc         [NMi/h]
    VMTOW           vehicle maximum takeoff weight              [kg]
    STRUCT_FRAC     fraction of VMTOW for structural components [-]
    PAYLOAD_FRAC    fraction of VMTOW for payload               [-]
    EPPS_FRAC       fraction of VMTOW for EPPS                  [-]
    P_TO            power required for takeoff (hover)          [kW]
    P_CLIMB         power required for climb                    [kW]
    P_CRUISE        power required for cruise                   [kW]
    P_DESC          power required for descent                  [kW]
    P_SLOW          power required to slow after descent        [kW]
    P_STOP          power required to stop after slow           [kW]
    P_LAND          power required to land (hover)              [kW]
    T_TO            duration of takeoff                         [s]
    T_CLIMB         duration of climb                           [s]
    T_DESC          duration of descent                         [s]
    T_SLOW          duration of slow                            [s]
    T_TRAN          duration of transition. Uses P_CLIMB        [s]
    T_STOP          duration of stop                            [s]
    T_LAND          duration of landing.                        [s] 
    RPM_TO          rotational speed of takeoff                 [RPM]
    RPM_CLIMB       rotational speed of climb                   [RPM]
    RPM_CRUISE      rotational speed of cruise                  [RPM]
    RPM_DESC        rotational speed of descent                 [RPM]
    RPM_SLOW        rotational speed of slow                    [RPM]
    RPM_TRAN        rotational speed of transition              [RPM]
    RPM_STOP        rotational speed of stop                    [RPM]
    RPM_LAND        rotational speed of landing                 [RPM]

    #Derived constants:
    cruiseS         cruise groundspeed converted to             [km/h]
    climbS          climb groundspeed converted to              [km/h]
    descS           descent groundspeed converted to            [km/h]
    EPPSmass        absolute mass of EPPS                       [kg]
    P = [P_TO,P_CLIMB,P_CRUISE,P_DEC,P_SLOW,P_CLIMB,P_STOP,P_LAND]
    T = [T_TO,T_CLIMB,0,       T_DEC,T_SLOW,T_TRAN, T_STOP,T_LAND]
                    as T_CRUISE is what changes as battery mass and 
                    motor efficiency changes, it is not set here.
    S = [0,   climbS, cruiseS, descS,0,     0,      0,     0]
                    assume negligible groundspeed during hover portions
    RPM = [RPM_TO,RPM_CLIMB,RPM_CRUISE,RPM_DESC,RPM_SLOW,RPM_STOP,RPM_LAND]
    omega           LIST is equivalent to RPM, converted to     [rad/s]
    time            LIST of times 1s time step, 60s longer than sum of T
                    30s header and 30s tail
    power           LIST of powers at 1s time step
    speed           LIST of speeds at 1s time step
    rotspeed        LIST of rotor rotational speeds at 1s time step
    """
    def __init__(self):    # System Constants
        import numpy as np
        import scipy.interpolate
        
        self.NM_TO_KM = 1.85         # [km/NMi]
        self.CRUISE_SPEED = 110      # [NMi/h]
        self.CLIMB_SPEED = 80        # [NMi/h]
        self.DESC_SPEED = 100        # [NMi/h]

        self.VMTOW = 1633            # [kg]
        self.STRUCT_FRAC = 1/3       # [-]    
        self.PAYLOAD_FRAC = 1/3      # [-]
        self.EPPS_FRAC = 1/3         # [-]

        self.P_TO = 335              # [kW]
        self.P_CLIMB = 150           # [kW]
        self.P_CRUISE = 69           # [kW]
        self.P_DESC = 20             # [kW]
        self.P_SLOW = 80             # [kW]
        self.P_STOP = 300            # [kW]
        self.P_LAND = 300            # [kW]

        self.T_TO = 37               # [s]
        self.T_CLIMB = 300           # [s]
        self.T_DESC = 300            # [s]
        self.T_SLOW = 60             # [s]
        self.T_TRAN = 15             # [s] Uses P_CLIMB
        self.T_STOP = 7              # [s]
        self.T_LAND = 30             # [s] Uses P_TO

        self.RPM_TO = 1200           # [RPM]
        self.RPM_CLIMB = 1050        # [RPM]
        self.RPM_CRUISE = 1000       # [RPM]
        self.RPM_DESC = 850          # [RPM]
        self.RPM_SLOW = 1100         # [RPM]
        self.RPM_TRAN = 1000         # [RPM]
        self.RPM_STOP = 1200         # [RPM]
        self.RPM_LAND = 1200         # [RPM]
    
        # Derived constants
        self.cruiseS = self.CRUISE_SPEED * self.NM_TO_KM   # [km/h]
        self.climbS = self.CLIMB_SPEED * self.NM_TO_KM     # [km/h]
        self.descS = self.DESC_SPEED * self.NM_TO_KM       # [km/h]

        self.EPPSmass = self.VMTOW * self.EPPS_FRAC    # [kg]

        self.P = [self.P_TO,self.P_CLIMB,self.P_CRUISE,self.P_DESC,self.P_SLOW,self.P_CLIMB,self.P_STOP,self.P_LAND]
        self.T = [self.T_TO,self.T_CLIMB, 0,self.T_DESC,self.T_SLOW,self.T_TRAN,self.T_STOP,self.T_LAND]
        self.S = [0,   self.climbS, self.cruiseS, self.descS,0,     0,      0,      0]
        self.RPM=[self.RPM_TO,self.RPM_CLIMB,self.RPM_CRUISE,self.RPM_DESC,self.RPM_SLOW,self.RPM_TRAN,self.RPM_STOP,self.RPM_LAND]
        self.omega = np.multiply(np.divide(self.RPM,60),2*np.pi)
        self.Torque = np.divide(np.multiply(self.P,1000),self.omega)

        time = range(0,sum(self.T)) 
        power = [0] * (max(time)+1) # one second time step
        speed = [0] * (max(time)+1)
        rotspeed = [0] * (max(time)+1)
        torque = [0] * (max(time)+1)

        for x in range(len(self.P)):
            if x == 0:
                tlast = self.T[x]
                for y in range(self.T[x]):
                    power[y] = self.P[x]
                    speed[y] = self.S[x]
                    rotspeed[y] = self.omega[x]
                    torque[y] = self.Torque[x]
            else:
                for y in range(tlast,tlast+self.T[x]):
                    power[y] = self.P[x]
                    speed[y] = self.S[x]
                    rotspeed[y] = self.omega[x]
                    torque[y] = self.Torque[x]
                tlast += self.T[x]
        
        self.torque = torque
        self.time = time
        self.power = power
        self.speed = speed
        self.rotspeed = rotspeed
        self.PProp = scipy.interpolate.interp1d(self.time, self.power)
    
# End class mission() =========================================================================


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
    for x in ['title', 'xlabel','ylabel','filepath','fignum','legend']:
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

    if prop["legend"] is True: 
        plt.legend(loc=0)
        print("Legend: On")
    # show figure
    plt.draw()
    plt.show()
    
    # save figure
    nameandpath = prop['filepath'] + name
    figout.savefig(nameandpath, bbox_inches='tight')
# End function














