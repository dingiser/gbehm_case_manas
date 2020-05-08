import numpy as np
import pandas as pd
from datetime import datetime
from netCDF4 import Dataset,num2date

def etorh(qe, airtk, pres):
    '''
    This code comes from forcing_class.py by Li Hongyi,it is used
    to convert the specific humidity to relative humidity

    :param qe:
    :param airtk:
    :param pres:
    :return:
    '''
    airtc = airtk - 273.15
    es = 611. * np.exp(17.27 * airtc / (237.3 + airtc))
    es1 = 611. * np.exp(21.87 * airtc / (265.5 + airtc))
    es[airtc < 0.] = es1[airtc < 0.]

    rh = qe * pres / (0.622 * es)
    rh = np.where(rh > 1., 0.95, rh)
    rh = np.where(rh < 0., 0.05, rh)

    qe = np.where(rh > 1., 0.622 * es * 0.95 / pres, qe)
    qe = np.where(rh < 0., 0.622 * es * 0.05 / pres, qe)


    return rh, qe

def rhtoe(rh,airtk,pres):
    '''

    :param rh: 0-1
    :param airtk: C
    :param pres: Pa
    :return: rh 0-1
            qe kg/kg
    '''
    airtc = airtk - 273.15    #units C
    es = 611. * np.exp(17.27 * airtc / (237.3 + airtc)) #units Pa
    es1 = 611. * np.exp(21.87 * airtc / (265.5 + airtc))
    es[airtc < 0.] = es1[airtc < 0.]

    rh = np.where(rh > 1., 0.95, rh)
    rh = np.where(rh < 0., 0.05, rh)

    qe=0.622*es*rh/pres #units

    return rh,qe

def Tstamp2idate(ts,deltim):
    '''

    :param ts:
    :param deltim:
    :return: year, julian day, seconds of the starting time
    '''
    idc=pd.to_datetime(ts).dayofyear
    ihc=(idc-1)*24
    idate = np.array([ts.year, idc, ihc * deltim])

    return idc,ihc,idate

def _driverXJ(times,MeteoPath):
    '''

    :param times: list
        [iyear,imonth]
    :param MeteoPath:
    :return:
    '''
    iyear, imonth = times[0], times[1]
    meteoNf = Dataset(MeteoPath[0] + 'Manas_drivers_XJ_' + str(iyear) + '%02d' % (imonth) + '.nc')
    glw_nc = meteoNf['lwr'][:]
    psfc_nc = meteoNf['PSFC'][:]
    swd_nc = meteoNf['swr'][:]
    airt_nc = meteoNf['T2'][:]
    U10m_nc = meteoNf['U10'][:]
    V10m_nc = meteoNf['V10'][:]
    rh_nc = meteoNf['2mrh2'][:]
    rh_nc, Q2_nc = rhtoe(rh_nc, airt_nc, psfc_nc)
    prec_nc = meteoNf['RAINC'][:]
    #ncTime = num2date(meteoNf['time'][:], units=meteoNf['time'].units)
    ncTime = pd.date_range(datetime(iyear, imonth, 1),
                           periods=glw_nc.shape[0], freq='H')
    ncTime.to_pydatetime()
    meteoNf.close()
    return glw_nc, psfc_nc, swd_nc, airt_nc, U10m_nc, V10m_nc, rh_nc, Q2_nc, prec_nc, ncTime

def _driverWRF(times,MeteoPath):
    iyear, imonth = times[0], times[1]
    ncfile = Dataset(MeteoPath[0] + 'Manas_CGBM_' + str(iyear) + '-%02d' % (imonth) + "_v4" + '.nc',
                     format='NETCDF4')
    glw_nc = ncfile['forc_frl'][:]
    psfc_nc = ncfile['forc_psrf'][:]
    swd_nc = ncfile['forc_solarin'][:]
    airt_nc = ncfile['forc_t'][:]
    U10m_nc = ncfile['forc_us'][:]
    V10m_nc = ncfile['forc_vs'][:]
    Q2_nc = ncfile['forc_q'][:]
    rh_nc = ncfile.variables['forc_rh'][:]
    prec_nc = np.load(MeteoPath[1]+"Manas_WRF_Prec_"+str(iyear)+"-%02d"%(imonth)+".npy")/3600.
    ncTime=pd.date_range(datetime(iyear,imonth,1),
                         periods=glw_nc.shape[0],freq='H')
    ncTime=ncTime.to_pydatetime()
    return glw_nc, psfc_nc, swd_nc, airt_nc, U10m_nc, \
           V10m_nc, rh_nc, Q2_nc, prec_nc, ncTime

def _driverERA(times,MeteoPath):
    iyear, imonth = times[0], times[1]
    ncfile = Dataset(MeteoPath[0] + 'Manas_driver_' + str(iyear) + '%-02d' % (imonth) + '.nc')
    glw_nc=ncfile['lwr'][:]
    psfc_nc=ncfile['pres'][:]
    swd_nc=ncfile['swr'][:]
    airt_nc=ncfile['temp2m'][:]
    U10m_nc=ncfile['U10m'][:]
    V10m_nc=ncfile['V10m'][:]
    Q2_nc=ncfile['q2m'][:]
    rh_nc=ncfile['rh2m'][:]
    prec_nc=ncfile['prec'][:]/3600.
    ncTime = pd.date_range(datetime(iyear, imonth, 1),
                           periods=glw_nc.shape[0], freq='H')
    ncTime = ncTime.to_pydatetime()
    return glw_nc, psfc_nc, swd_nc, airt_nc, U10m_nc, \
           V10m_nc, rh_nc, Q2_nc, prec_nc, ncTime


def ReadMeteoDeriver(times,MeteoPath,name):
    '''

    :param itime: list
        [iyear,imonth]
    :param MeteoPath: lsit
        [str,...]path to store meteological data.
    :return:
    '''
    if name=='XJ':
        return _driverXJ(times,MeteoPath)
    elif name == 'WRF':
        return _driverWRF(times,MeteoPath)
    elif name == 'ERA_Interim':
        return _driverERA(times,MeteoPath)
    else:
        raise KeyError('No meteological driver be found.')
