import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import matplotlib.ticker as ticker

# USER MUST ADD FOLDER PATH TO FILES
folder_path='FOLDER_PATH//FOLDER_PATH//FOLDER_PATH//'
folder_path='/Users/pedrosantanagiacon/Documents/Documents/DU4/Diss/Data/2_Evaluated//'

# calculates average of a list
def average(lst): 
    return sum(lst) / len(lst) 

# function for single exponential decay
def exp_decay(t, A, LT):
    intensity = A * np.exp(-t/LT)
    return intensity

# function for bi-exponential decay
def bi_exp_decay(t, A_PF, LT_PF, A_DF, LT_DF):
    intensity = A_PF * np.exp(-t/LT_PF) + A_DF * np.exp(-t/LT_DF)
    return intensity

# calculates the chi square value between an experimental dataset and the fitted expression
def chi_squared_(experiment_y, fitted_y):
    try:
        chi = sum(((experiment_y - fitted_y)**2)/ experiment_y)
    except:
        # if function fails chi square is assumed to be big
        chi = 10**10
    return chi



directory = os.getcwd()
sample = []
# gets files in directory
for file in os.listdir(str(folder_path)):
    if file.endswith(".txt"):
        sample.append(os.path.join(file))
    elif file.endswith(".TXT"):
        sample.append(os.path.join(file))
for z in range(0,len(sample)):
    print(sample[z][:])
    array2D = []
    file = str(folder_path+str(sample[z]))
    with open(file, "r") as data:
        for line in data.readlines():
            array2D.append(line.split())
    # obtain times and normalized intensity for the decay
    times = np.array(array2D[0][1:])
    normalized = np.array(array2D[-3][1:])
    times = times.astype(float)
    normalized = normalized.astype(float)
    array2D = np.array(array2D)

    # get rid of measurements dominated by noise by comparing their average to their standard deviation
    clean_columns = []
    for i in range(0, len(array2D[0])-4):
        column = [float(row[i]) for row in array2D[0:-6]]
        if average(column) < np.std(column):
            clean_columns.append(i)
    first_column = 1
    last_column = clean_columns[-1]

    # plots decay
    plt.scatter(times[first_column:last_column], normalized[first_column:last_column], color = 'k', marker = 's', s = 10)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Delay Time (s)')
    plt.ylabel('Integrated Area (a.u.)')
    plt.show()

    # intial guesses for fittings
    initialA = 0.00002
    initialC = 2.39*10**-3
    initial = np.array([initialA, initialC])

    # obtains the two fittings of the bi-exponential decay
    # different fittings were obtained
    chi = 10**10
    for i in range(5, last_column-first_column-5):
        parameters1, covariance1 = curve_fit(exp_decay, times[first_column:i], normalized[first_column:i], initial, bounds=([0, 0],[1,1]))
        parameters2, covariance2 = curve_fit(exp_decay, times[i:last_column], normalized[i:last_column], initial, bounds=([0, 0],[1,1]))

        y_values1 = exp_decay(times[first_column:
                                    i], parameters1[0], parameters1[1])
        y_values2 = exp_decay(times[i:last_column], parameters2[0], parameters2[1])
        y_values = bi_exp_decay(times, parameters1[0], parameters1[1], parameters2[0], parameters2[1])

        fit1 = chi_squared_(normalized[first_column:i],y_values1)
        fit2 = chi_squared_(normalized[i:last_column], y_values2)
        fit = chi_squared_(normalized, y_values)

        # checks for best fitting out of all fittings attempted
        if fit < chi:
            chi = fit
            # extracts parameters from fitting (Amplitudes and Lifetimes)
            lifetime_PF = parameters1[1]
            lifetime_DF = parameters2[1]
            A_PF = parameters1[0]
            A_DF = parameters2[0]
            A_PF_uncertainty = np.sqrt(np.diag(covariance1))[0]
            A_DF_uncertainty = np.sqrt(np.diag(covariance2))[0]
            lifetime_PF_uncertainty = np.sqrt(np.diag(covariance1))[1]
            lifetime_DF_uncertainty = np.sqrt(np.diag(covariance2))[1]

    # plots experimental data and obtained fitting
    plt.scatter(times[first_column:last_column], normalized[first_column:last_column], color = 'k', marker = 's', s = 10)
    plt.plot(times[first_column:last_column], y_values[first_column:last_column], color = 'red')
    plt.xlabel('Delay Time (s)')
    plt.ylabel('Integrated Area (a.u.)')
    plt.xscale('log')
    plt.yscale('log')
    plt.show()

    # prints parameters from fitting
    print('Lifetimes and amplitudes')
    print(lifetime_PF)
    print(lifetime_PF_uncertainty)
    print(lifetime_DF)
    print(lifetime_DF_uncertainty)
    print(A_PF)
    print(A_PF_uncertainty)
    print(A_DF)
    print(A_DF_uncertainty)

    # User is asked to input the measured aerated PLQY and aerated-degassed ratio
    input_PF_aer = float(input('Aerated PLQY: '))
    input_RATIO = float(input('Ratio: '))

    # PLQY value calculated from ratio and aerated PLQY, as well as its uncertainty
    PLQY = (input_RATIO + 1) * input_PF_aer
    PLQY_UNC = PLQY * np.sqrt((1/1728)**2+(1/329)**2)

    # calculates DF/PF ratio for bi-exponential model (DIFFERENT THAN AERATED-DEGASSED RATIO!)
    ratio = (lifetime_DF*A_DF)/(lifetime_PF*A_PF)
    ratio_uncertainty = ratio * np.sqrt((lifetime_PF_uncertainty/lifetime_PF)**2+(lifetime_DF_uncertainty/lifetime_DF)**2+(A_PF_uncertainty/A_PF)**2+(A_DF_uncertainty/A_DF)**2)

    # calculates PF and DF individually, and their uncertainties
    efficiency_PF = PLQY/(1+ratio)
    efficiency_DF = PLQY - efficiency_PF
    efficiency_PF_uncertainty = efficiency_PF * np.sqrt((PLQY_UNC/PLQY)**2+(ratio_uncertainty/ratio)**2)
    efficiency_DF_uncertainty = np.sqrt(PLQY_UNC**2 + efficiency_PF_uncertainty**2)

    # calculates maximum and minimum values for ISC*RISC expression
    ISC_RISC_MAX = (ratio + ratio_uncertainty)/(1+ ratio - ratio_uncertainty)
    ISC_RISC_MIN = (ratio - ratio_uncertainty)/(1+ ratio + ratio_uncertainty)

    # outputs PLQY, PF, DF and RATIO values
    print('PLQY = ' + str(PLQY) + ' + ' + str(PLQY_UNC))
    print('PF = ' + str(efficiency_PF) + ' + ' + str(efficiency_PF_uncertainty))
    print('DF = ' + str(efficiency_DF) + ' + ' + str(efficiency_DF_uncertainty))
    print('RATIO = ' + str(ratio) + ' + ' + str(ratio_uncertainty))

    # obtains all possible ISC and RISC values within the range
    values = []
    for ISC in np.arange(0,1,0.001):
        for RISC in np.arange(0,1,0.001):
            if ISC*RISC > ISC_RISC_MIN and ISC*RISC < ISC_RISC_MAX:
                if ISC <= 1 - efficiency_PF:
                    values.append([ISC, RISC, ISC*RISC])

    # creates arrays for storing all decay rate constants calculated for each ISC RISC pair
    k_r_S_values = []
    k_r_S_err_values = []
    k_nr_S_values = []
    k_nr_S_err_values = []
    k_isc_values = []
    k_isc_err_values = []
    k_ic_S_values = []
    k_ic_S_err_values = []
    k_risc_values = []
    k_risc_err_values = []
    k_ict_values = []
    ISCs = []
    ICs = []
    RISCs = []
    ICTs = []

    # calculates all rate constants and efficiencies, as well as their respective uncertainties
    for value in values:
        ISC = value[0]
        RISC = value[1]
        IC = 1 - ISC - efficiency_PF
        ICT = 1 - RISC

        k_r_S = efficiency_PF / lifetime_PF
        k_r_S_err = np.sqrt((efficiency_PF_uncertainty/efficiency_PF)**2+(efficiency_DF_uncertainty/efficiency_DF)**2) * k_r_S

        k_isc = ISC/lifetime_PF
        k_isc_err = np.sqrt((0.01/ISC)**2+(lifetime_PF_uncertainty/lifetime_PF)**2) * k_isc

        k_ic_S = (1/lifetime_PF) - k_r_S - k_isc
        k_ic_S_err = np.sqrt(k_r_S_err**2 + k_isc_err**2 + (lifetime_PF_uncertainty/(lifetime_PF**2))**2)
        
        k_risc = (RISC/lifetime_DF)*((efficiency_PF + efficiency_DF)/efficiency_PF)
        k_risc_err = np.sqrt((0.01/RISC)**2+(lifetime_DF_uncertainty/lifetime_DF)**2 + (efficiency_PF_uncertainty/efficiency_PF)**2 + (np.sqrt(efficiency_PF_uncertainty**2 + efficiency_DF_uncertainty**2)/(efficiency_PF + efficiency_DF))**2) * k_risc

        k_ict = (ICT/lifetime_DF)*((efficiency_PF + efficiency_DF)/efficiency_PF)

        k_r_S_values.append(k_r_S)
        k_r_S_err_values.append(k_r_S_err)
        k_isc_values.append(k_isc)
        k_isc_err_values.append(k_isc_err)
        k_ic_S_values.append(k_ic_S)
        k_ic_S_err_values.append(k_ic_S_err)
        k_risc_values.append(k_risc)
        k_risc_err_values.append(k_risc_err)
        k_ict_values.append(k_ict)

        ISCs.append(ISC)
        ICs.append(IC)
        RISCs.append(RISC)
        ICTs.append(ICT)

    # calculates rate constants and efficiencies from average, and uncertainty with various methods (highest uncertainty will be quoted)
        
    k_r_S = average(k_r_S_values)
    k_r_S_err3 = np.std(k_r_S_values) / np.sqrt(len(k_r_S_values))
    k_r_S_err1 = max(k_r_S_values) - k_r_S
    k_r_S_err2 = min(k_r_S_values) - k_r_S

    k_isc = average(k_isc_values)
    k_isc_err3 = np.std(k_isc_values) / np.sqrt(len(k_isc_values))
    k_isc_err1 = max(k_isc_values) - k_isc
    k_isc_err2 = min(k_isc_values) - k_isc

    k_ic_S = average(k_ic_S_values)
    k_ic_S_err3 = np.std(k_ic_S_values) / np.sqrt(len(k_ic_S_values))
    k_ic_S_err1 = max(k_ic_S_values) - k_ic_S
    k_ic_S_err2 = min(k_ic_S_values) - k_ic_S

    k_risc = average(k_risc_values)
    k_risc_err3 = np.std(k_risc_values) / np.sqrt(len(k_risc_values))
    k_risc_err1 = max(k_risc_values) - k_risc
    k_risc_err2 = min(k_risc_values) - k_risc


    k_ict = average(k_ict_values)
    k_ict_err1 = max(k_ict_values) - k_ict
    k_ict_err2 = min(k_ict_values) - k_ict
    k_ict_err3 = np.std(k_ict_values) / np.sqrt(len(k_ict_values))

    EFF_ISC = average(ISCs)
    EFF_err1 = max(ISCs) - EFF_ISC
    EFF_err2 = min(ISCs) - EFF_ISC
    EFF_err3 = np.std(ISCs) / np.sqrt(len(ISCs))

    EFF_IC = average(ICs)
    IC_err1 = max(ICs) - EFF_IC
    IC_err2 = min(ICs) - EFF_IC
    IC_err3 = np.std(ICs) / np.sqrt(len(ICs))

    EFF_RISC = average(RISCs)
    RISC_err1 = max(RISCs) - EFF_RISC
    RISC_err2 = min(RISCs) - EFF_RISC
    RISC_err3 = np.std(RISCs) / np.sqrt(len(RISCs))

    EFF_ICT = average(ICTs)
    ICT_err1 = max(ICTs) - EFF_ICT
    ICT_err2 = min(ICTs) - EFF_ICT
    ICT_err3 = np.std(ICTs) / np.sqrt(len(ICTs))




    # returns all calculated values
    print('\n\n\n\n\n\n')
    print('k_r_S = ' + str(k_r_S) + ' +- ' + str(k_r_S_err) + ' +- ' + str(k_r_S_err2) + '+-' + str(k_r_S_err3))
    print('k_ic_S = ' + str(k_ic_S) + ' +- ' + str(k_ic_S_err) + ' +- ' + str(k_ic_S_err2) + ' +- ' + str(k_ic_S_err3))
    print('k_isc = ' + str(k_isc) + ' +- ' + str(k_isc_err) + ' +- ' + str(k_isc_err2) + ' +- ' + str(k_isc_err3))
    print('k_risc = ' + str(k_risc) + ' +- ' + str(k_risc_err) + ' +- ' + str(k_risc_err2) + ' +- ' + str(k_risc_err3))
    print('k_nr = ' + str(k_ict) + ' +- ' + str(k_ict_err1) + ' +- ' + str(k_ict_err2) + ' +- ' + str(k_ict_err3))
    print('EFF_r_S = ' + str(efficiency_PF) + ' +- ' + str(efficiency_PF_uncertainty))
    print('EFF_ic_S = ' + str(EFF_IC) + ' +- ' + str(IC_err1) + ' +- ' + str(IC_err2)+ ' +- ' + str(IC_err3))
    print('EFF_isc_S = ' + str(EFF_ISC) + ' +- ' + str(EFF_err1) + ' +- ' + str(EFF_err2)+ ' +- ' + str(EFF_err3))
    print('EFF_risc_S = ' + str(EFF_RISC) + ' +- ' + str(RISC_err1) + ' +- ' + str(RISC_err2)+ ' +- ' + str(RISC_err3))
    print('EFF_ic_T = ' + str(EFF_ICT) + ' +- ' + str(ICT_err1) + ' +- ' + str(ICT_err2) + ' +- ' + str(ICT_err3))