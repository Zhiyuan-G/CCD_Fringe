import scipy.interpolate
import astropy.table
import astropy.units as u
import csv
import pandas as pd
import tmm
import numpy as np

def load_refraction_data(Temp,Epoxy_ind = 1.5,kind='linear'):
    
    interpolators = {}
    # Vacuum.
    interpolators['Vacuum'] = lambda wlen: np.ones_like(wlen)
    
    
    ### Silicon (Temperature independent)
    abs_table = astropy.table.Table.read('data/Si-absorption.csv', format='ascii.csv', names=('wlen', 'k', 'a'))
    table = astropy.table.Table.read('data/Si-index.csv', format='ascii.csv', names=('wlen', 'n'))
    wlen = abs_table['wlen']
    k = abs_table['k']
    
    # Interpolate n onto the finer grid used for k.
    n = np.interp(wlen, table['wlen'], table['n'])
    # Build a linear interpolator of the complex index of refraction vs wavelength [nm].
    n = n + 1.j * k
    interpolators['Si'] = scipy.interpolate.interp1d(wlen, n, copy=True, kind=kind)
    
    
    ### Silicon (Temperature dependent)
    Green_table = astropy.table.Table.read('data/Si_index_Green.csv',format = 'ascii.csv')
    wlen_Green = Green_table['wlen']*1e3
    c_n = Green_table['C_n']
    c_k = Green_table['C_k']
    
    def Temp_Si_index (index, tem_coeff, T = 173):
        T0 = 300.0
        b = tem_coeff * 1e-4 * T0
        index_temp = index * (T / T0)**b
        return(index_temp)
    
    n_temp = Temp_Si_index(Green_table['n'],c_n,T = Temp)
    k_temp = Temp_Si_index(np.imag(interpolators['Si'](wlen_Green)),c_k,T = Temp)
    # Interpolate n onto the finer grid used for k.
    
    n_temp = np.interp(wlen_Green, Green_table['wlen']*1e3, n_temp)

    # Build a linear interpolator of the complex index of refraction vs wavelength [nm].
    n_temp = n_temp + 1.j * k_temp
    interpolators['Si_Temp'] = scipy.interpolate.interp1d(wlen_Green, n_temp, copy=True, kind= kind)
    ###################################
    
    abs_table = astropy.table.Table.read('data/Si_index_Green.csv',format = 'ascii.csv')
    wlen = abs_table['wlen']*1e3
    k = abs_table['k']
    c_n = abs_table['C_n']
    c_k = abs_table['C_k']
    def Temp_Si_index (index, tem_coeff, T = Temp):
        
        T0 = 300.0
        b = tem_coeff * 1e-4 * T0
        index_temp = index * (T / T0)**b
    
        return(index_temp)
    n_173 = Temp_Si_index(abs_table['n'],c_n,T = Temp)
    k_173 = Temp_Si_index(abs_table['k'],c_k,T = Temp)
    # Interpolate n onto the finer grid used for k.
    n_173 = np.interp(wlen, abs_table['wlen']*1e3, n_173)
    # Build a linear interpolator of the complex index of refraction vs wavelength [nm].
    n_173 = n_173 + 1.j * k_173
    interpolators['Si_Green'] = scipy.interpolate.interp1d(wlen, n_173, copy=True, kind= kind)
    ###################################
    
    # Read tabulated Si3N4 data from
    # http://refractiveindex.info/?shelf=main&book=Si3N4&page=Philipp
    table = astropy.table.Table.read('data/Si3N4-index.csv', format='ascii.csv', names=('wlen', 'n'))
    # Convert from um to nm.
    wlen = 1e3 * table['wlen']
    n = table['n']    
    interpolators['Si3N4'] = scipy.interpolate.interp1d(wlen, n, copy=True, kind=kind)

    # Read SiO2 tabulated data from
    # http://refractiveindex.info/?shelf=main&book=SiO2&page=Malitson
    table = astropy.table.Table.read('data/SiO2-index.csv', format='ascii.csv', names=('wlen', 'n'))
    # Convert from um to nm.
    wlen = 1e3 * table['wlen']
    n = table['n']    
    interpolators['SiO2'] = scipy.interpolate.interp1d(wlen, n, copy=True, kind=kind)

    # Read MgF2 tabulated data from
    # http://refractiveindex.info/?shelf=main&book=MgF2&page=Li-o
    table = astropy.table.Table.read('data/MgF2-index.csv', format='ascii.csv', names=('wlen', 'n'))
    # Convert from um to nm.
    wlen = 1e3 * table['wlen']
    n = table['n']    
    interpolators['MgF2'] = scipy.interpolate.interp1d(wlen, n, copy=True, kind=kind)
    
    #Epoxy
    n = np.full_like(n,1)*Epoxy_ind
    interpolators['Epoxy'] = scipy.interpolate.interp1d(wlen, n, copy=True, kind=kind)
    
    
    # Read Ta2O5 tabulated data from
    # https://refractiveindex.info/?shelf=main&book=Ta2O5&page=Rodriguez-de_Marcos
    table = astropy.table.Table.read('data/Ta2O5-index.csv', format='ascii.csv', names=('wlen', 'n','k'))
    # Convert from um to nm.
    wlen = 1e3 * table['wlen']
    n = table['n']    
    k = table['k']
    n = n + 1.j*k
    interpolators['Ta2O5'] = scipy.interpolate.interp1d(wlen, n, copy=True, kind=kind)
    

    
    return interpolators