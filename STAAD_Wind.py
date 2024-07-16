# -------------------------------------------------------------
# Based on: NS-EN 1991-1-4 Actions on Structures
# -------------------------------------------------------------

# glob and os don't need a virtual environment
import glob
import os
import matplotlib.pyplot as plt
import numpy as np

# Functions: --------------------------------------------------
def ln(x):
    return 99999*((x**(1/99999))-1)
# -------------------------------------------------------------

# Find all files with the extension '.inp':
files = glob.glob("*.inp")

print('\n\n')

# Read and execute each input file:
for file in files:
    name = os.path.splitext(os.path.basename(file))[0]
    with open(file, 'r') as input:
        code = input.read()
        exec(code)
        
    print("From '{}.inp'".format(name))
    print()

    # Several setups cases are allowed in same input file
    for i in range(0, len(case)):
        # ------ Get un-scaled wind profile -------------------
        # Wind statistics:
        k = 0.2
        n = 0.5
        p = 0.2 if stage[i] == 2 else 0.9999
        Cprob = ((1-k*ln(-ln(1-p)))/(1-k*ln(-ln(0.98))))**n
        # Wind factors:
        Cdir = 1
        Calt = 1
        Cseason = 1
        # Bse wind speed:
        Vb = Vref[i]*Cprob*Cdir*Calt*Cseason
        # Terrain rougness
        Kr = 0.215
        # Turbulence
        Kl = 1
        # Rougness length (Table 4.1)
        Z0 = 0.3
        # Min. height (Table 4.1)
        Zmin = 5
        # Terrain form factor
        C0 = 1
        # Air density
        Density = 0.00125
        
        Sigmay = Vb*Kr*Kl
        
        # ------ Get scale factor -----------------------------
        Ffact = 1.2
        Efact = 1.0
        Zref = Hw[i]-Y0[i]+H1[i]
        ce_zref = 1+7*Kl/(C0*ln(Zref/Z0))
        cr_z_ref = Kr*ln(Zref/Z0)
        vm_zref = cr_z_ref*C0*Vb
        qpref = ce_zref*Density*vm_zref**2/2
        S = round(Ffact*Efact*qpref*Mfact[i]*(W2[i]*H2[i]+\
                  H2[i]*(W1[i]-W2[i])/2)/Fwind[i], 2)
        
        _qb_ = []
        _Z_ = []
        XX = np.array([])
        YY = np.array([])
        j = 0
        while j-3 < H2[i]:
            j = j+3
            cr = Kr*ln(j/Z0)
            vm = cr*C0*Vb
            lv = Sigmay/vm
            ce = 1+7*lv
            qb = 0.5*Density*vm**2
            qblv = ce*qb
            Re = vm*1/(1.5*10**-5)/10**6
            Z = Y0[i]-H1[i]+j
            _qb_.append(round(qblv, 2))
            _Z_.append(round(Y0[i]-H1[i]+j, 2))
            XX = np.append(XX, qblv*S)
            YY = np.append(YY, Z-Y0[i])
            
        # Add profiles to plot
        plt. plot(XX, YY, label=case[i])
        
        # Print into terminal:
        print('   ', 32*'-', case[i], 37*'-')
        print('    TYPE', stage[i], 'Storage' if stage[i] == 2\
                                     else 'Load-in')
        print('    INT', ' '.join(map(str,_qb_)),'HEIG', ' '.join(map(str,_Z_)))
        #print('    HEIG', ' '.join(map(str,_Z_)))
        print()
        print('    LOAD # LOADTYPE Wind TITLE WIND', case[i])
        print('    WIND LOAD', head[i], S, 'TYPE', stage[i],\
              'XR', xrng[i][0], xrng[i][1],\
              'YR', yrng[i][0], yrng[i][1],\
              'ZR', zrng[i][0], zrng[i][1], 'OPEN')
        print()
    # end for i in range(0, len(case))
    print()
    
    # Display plot
    plt.title(structure)
    plt.xlabel("Wind pressure [kN/m^2]")
    plt.ylabel("Elevation above ground level [m]")
    plt.legend()
    plt.show()
    
# end for file in files
print()
