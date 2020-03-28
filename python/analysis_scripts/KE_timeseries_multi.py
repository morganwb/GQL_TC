def idkrn(t, KE):
    import numpy as np
    import matplotlib.pyplot as plt
    #find growth region
    first_found = False
    error = 10**-20
    first = 0
    last = 0
    delta_ts = []
    
    #Find max delta_t
    for i in range(5,len(t)):
        delta_ts.append(t[i]-t[i-1])
        max_delta_t = max(delta_ts)

    #fit curve
    for i in range(5,len(t)):
        delta_t = t[i]-t[i-1]
        if not first_found and delta_t>=max_delta_t:
            first = i
            first_found = True
        if first_found and delta_t<max_delta_t:
            last = i
            break

    if last == 0 or last-first==1:
        last = first+4
        first = first
    #last = first + 10
    print(first,last)

    polfit = np.polyfit(t[first:last],np.log(KE)[first:last],1)
    fit_equation = lambda t: t*polfit[0]+polfit[1]
    plt.plot(t,fit_equation(t),label='Exponential growth fit \n gamma ='+str(round(polfit[0],3)))
    pass

def plot_energies(folder_prefix):

    import matplotlib.pyplot as plt
    import h5py
    import glob
    import numpy as np
    from pathlib import Path
    import sys

    # get list of folders to pull data from, and sort by Lambda value
    folders = sorted(glob.glob(folder_prefix), key=lambda folder: int(folder.split("_")[-1]))
    folders.insert(0, folders[-1])
    del folders[-1]
    plt.figure()
    
    for folder in folders:
        filename = folder + "/scalar/scalar_s1.h5"
        df = h5py.File(filename,"r")
        t = df["scales/sim_time"][:]
        KE = df["tasks/KE"][:].ravel()
        
        
        if "GQL" not in folder:
            idkrn(t,KE)
            label = "DNS"
        else:
            label = "GQL, $\Lambda = " + folder.split("_")[-1] + "$"
        plt.plot(t, np.log(KE), label=label)
    plt.legend()
    plt.xlabel('Time ')
    plt.ylabel('ln(Energy)') 
    plt.savefig("KE_multifig.png")

if __name__ == '__main__':
    from pathlib import Path
    import sys
    scalar_file = Path(sys.argv[-1])
    plot_energies("TC_3d_re_100*")

