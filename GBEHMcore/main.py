import numpy as np
import pandas as pd
from datetime import datetime

from netCDF4 import Dataset,num2date,date2num
import pygbhm,pyice
from AuxTools import *


path1="./gisdata/"
path2="./parameter/"
out_path="./result/"

meteoDriverPath="../data/meteo/WRF_result/"


#-------constant define------------------------------------------------
newcols = 111              # columns of the input raster
newrows = 107              # rows of the input raster
nl = 7
nl_soil  = 7               # number of soil layers,default:7
nl_lake  = 10              # number of lake layers
maxsnl   = -5              # max number of snow layers
deltim   = 3600.         # seconds in a time-step
custom_soillayer = 0
record = 1
run_pbsm = 0               #wind blow sonw
run_gbhm=1
builtmod=1
nc=17                      #total number of the subbasin
nx=300
startTime="2004-01-01"
endTime="2009-12-31"
endPHTime="2004-12-31"
PHTimes=4

endletter='_v2'



#-----------Set preheating date----------------------------------------
tPH=pd.date_range(startTime,endPHTime)
ta=pd.date_range(startTime,endTime)

startyear = ta[0].year
startmonth= ta[0].month
startday  = ta[0].day
endyear   = ta[-1].year
endmonth  = ta[-1].month
endday    = ta[-1].day

start_sub = 1
end_sub   = nc

dates=PHTimes*tPH.tolist()+ta.tolist()

#------------Read GIS data--------------------------------------------------
nf=Dataset(path1+'GisPara_lupdate.nc')
longxy=nf['longitude'][:]
latixy=nf['latitude'][:]
elev=nf['elevation'][:].data
mask=nf['mask'][:]
slope=nf['slope'][:].data
length=nf['slope_length'][:].data
ivt=nf['landuse'][:]
ivt[ivt==15]=7
area=mask[:]
#------------------read soil data--------------------------------------------
ncf=Dataset(path1+"SOIL_Para.nc")
wsat=ncf['THSCH'][0]
wrsd=ncf['THR'][0]
watern=ncf['N'][0]
alpha=ncf['ALPHA'][0]

BD = ncf.variables['BD'][:]
CL = ncf.variables['CL'][:]
SA = ncf.variables['SA'][:]
SI = ncf.variables['SI'][:]
OC = ncf.variables['OC'][:]
SOM = ncf.variables['SOM'][:]
GRAV = ncf.variables['GRAV'][:]
thsgm = ncf.variables['THSGM'][:]  # 0-1
ksvg = ncf.variables['K_SVG'][:]  # cm.d-1
psis = ncf.variables['PSI_S'][:]  # cm
wfld = ncf.variables['TH33'][:]
lamd = ncf.variables['LAMBDA'][:]

ncf.close()

#------------------extent the soil data to nl layers--------------
gravel_grid = np.zeros((nl,newrows,newcols))
sand_grid   = np.zeros((nl,newrows,newcols))
clay_grid   = np.zeros((nl,newrows,newcols))
soc_grid    = np.zeros((nl,newrows,newcols))
bd_grid     = np.zeros((nl,newrows,newcols))
ks_grid     = np.zeros((nl,newrows,newcols))
ths_grid    = np.zeros((nl,newrows,newcols))
psi_grid    = np.zeros((nl,newrows,newcols))
wfld_grid   = np.zeros((nl,newrows,newcols))
lamd_grid   = np.zeros((nl,newrows,newcols))

if custom_soillayer == 0:
    if nl_soil <= 7:
        gravel_grid[0:nl_soil, :, :] = GRAV[0:nl_soil, :, :]
        sand_grid[0:nl_soil, :, :] = SA[0:nl_soil, :, :]
        clay_grid[0:nl_soil, :, :] = CL[0:nl_soil, :, :]
        soc_grid[0:nl_soil, :, :] = OC[0:nl_soil, :, :]
        bd_grid[0:nl_soil, :, :] = BD[0:nl_soil, :, :]
        ks_grid[0:nl_soil, :, :] = ksvg[0:nl_soil, :, :]
        ths_grid[0:nl_soil, :, :] = thsgm[0:nl_soil, :, :]
        psi_grid[0:nl_soil, :, :] = psis[0:nl_soil, :, :]
        wfld_grid[0:nl_soil, :, :] = wfld[0:nl_soil, :, :]
        lamd_grid[0:nl_soil, :, :] = lamd[0:nl_soil, :, :]
    if nl_soil > 7:
        gravel_grid[0:7, :, :] = GRAV[:, :, :]
        sand_grid[0:7, :, :] = SA[:, :, :]
        clay_grid[0:7, :, :] = CL[:, :, :]
        soc_grid[0:7, :, :] = OC[:, :, :]
        bd_grid[0:7, :, :] = BD[:, :, :]
        ks_grid[0:7, :, :] = ksvg[:, :, :]
        ths_grid[0:7, :, :] = thsgm[:, :, :]
        psi_grid[0:7, :, :] = psis[:, :, :]
        wfld_grid[0:7, :, :] = wfld[:, :, :]
        lamd_grid[0:7, :, :] = lamd[:, :, :]
        for itmp in range(7, nl_soil):
            gravel_grid[itmp, :, :] = gravel_grid[6, :, :]
            sand_grid[itmp, :, :] = sand_grid[6, :, :]
            clay_grid[itmp, :, :] = clay_grid[6, :, :]
            soc_grid[itmp, :, :] = soc_grid[6, :, :]
            bd_grid[itmp, :, :] = bd_grid[6, :, :]
            ks_grid[itmp, :, :] = ks_grid[6, :, :]
            ths_grid[itmp, :, :] = ths_grid[6, :, :]
            psi_grid[itmp, :, :] = psi_grid[6, :, :]
            wfld_grid[itmp, :, :] = wfld_grid[6, :, :]
            lamd_grid[itmp, :, :] = lamd_grid[6, :, :]
elif custom_soillayer == 1:
    # [0-45,45-91,91-166,166-289,289-493,493-829,829-1383] (/mm)
    for itmp in range(nl):
        if itmp == 0: jtmp = itmp + 0
        if itmp > 0 and itmp <= 3: jtmp = itmp - 1
        if itmp >= 4 and itmp < 6: jtmp = 3
        if itmp >= 6 and itmp < 8: jtmp = 4
        if itmp >= 8 and itmp < 11: jtmp = 5
        if itmp >= 11: jtmp = 6

        gravel_grid[itmp, :, :] = GRAV[jtmp, :, :]
        sand_grid[itmp, :, :] = SA[jtmp, :, :]
        clay_grid[itmp, :, :] = CL[jtmp, :, :]
        soc_grid[itmp, :, :] = OC[jtmp, :, :]
        bd_grid[itmp, :, :] = BD[jtmp, :, :]
        ks_grid[itmp, :, :] = ksvg[jtmp, :, :]
        ths_grid[itmp, :, :] = thsgm[jtmp, :, :]
        psi_grid[itmp, :, :] = psis[jtmp, :, :]
        wfld_grid[itmp, :, :] = wfld[jtmp, :, :]
        lamd_grid[itmp, :, :] = lamd[jtmp, :, :]


#-----------------initialize GBHM --------------------
subbasin,nsub,nflow,dx,dr_gbhm,s0,wriver,mnroughness,ngrid,grid_row,grid_col,\
                                psubbasin,nbasinup,pbasinup,Dr_grid,basinmap\
                                = pygbhm.read_basin(newrows,newcols,nc)
drw_gbhm = 0.5*dr_gbhm
Dr_grid  = np.flipud(Dr_grid)
s0[s0==0.]    =0.00001


#---------------------???-----------------------------
q1    = np.zeros((nc,nx))
qr1   = np.zeros((nc,nx))
qr2   = np.zeros((nc,nx))
qlin1 = np.zeros((nc,nx))

qd    = np.zeros((nc,366) )
qm    = np.zeros((nc,12)  )
sst   = np.zeros((newrows,newcols) )
testeva = np.zeros((newrows,newcols) )
testrnf = np.zeros((newrows,newcols) )
testp   = np.zeros((newrows,newcols) )
tmpzero2d  = np.zeros((newrows,newcols))
tmpzero3d  = np.zeros((nl_soil,newrows,newcols))





#--------------pyice initialize--------------------
snow_d_grid = np.zeros((newrows, newcols))
lakedepth = np.zeros((newrows, newcols))
soil_t_grid = np.zeros((nl, newrows, newcols)) + 280.
soil_w_grid = np.zeros((nl, newrows, newcols)) + 0.3
_,_,idate0=Tstamp2idate(dates[0],deltim)

lons, dlon, dlat, itypwat, dz_lake, \
soil_s_v_alb, soil_d_v_alb, soil_s_n_alb, soil_d_n_alb, \
porsl, psi0, bsw, hksati, csol, dksatu, dkdry, rootfr, \
z0m, displa, sqrtdi, effcon, vmax25, slti, hlti, shti, hhti, \
trda, trdm, trop, gradm, binter, extkn, chil, ref, tran, \
z_soisno, dz_soisno, t_soisno, wliq_soisno, wice_soisno, smf_soisno, smf_ice, t_grnd, \
tlsun, tlsha, ldew, sag, scv, snowdp, fveg, fsno, sigf, green, lai, sai, \
coszen, albg, albv, alb, ssun, ssha, \
thermk, extkb, extkd, zwt, wa, t_lake, lake_icefrac, qref, rst, trad, tref, \
zlnd, zsno, csoilc, dewmx, wtfact, capr, cnfac, ssi, \
wimp, pondmx, smpmax, smpmin, trsmx0, tcrit, \
emis, z0ma, zol, rib, ustar, qstar, tstar, fm, fh, fq \
    = pyice.initial.initialize(idate0, newrows, newcols, nl, nl_soil, maxsnl, nl_lake, ivt,
                                latixy, longxy, snow_d_grid, lakedepth, gravel_grid, sand_grid, clay_grid,
                                soc_grid, bd_grid, soil_t_grid, soil_w_grid)

ds = z_soisno[nl_soil - maxsnl - 1, :, :] + 0.5 * dz_soisno[nl_soil - maxsnl - 1, :, :]  # depth of topsoil(m)
dg = 5 * ds + 17  # depth of unconfined acquifer (m)
wa = dg * 0.12 * 0.5 * 1000  # water storage in aquifer [mm]

# Suppose water table is half depth of dg,
# and the porosity is 0.1.

zwt = ds + dg - (wa * 0.001) / 0.12  # the depth to water table [m]``
# zwt    = copy.deepcopy(ds)
sst = wa * 0.
porsl = ths_grid + 0.
hksati = ks_grid * 10. / 24. / 3600.  # cm/day -> mm/s
# hksati = ks_grid * 10. / 24. / 3600.
# hksati[0,:,:]=20.*hksati[0,:,:]


psi0 = psi_grid * 10.  # cm -mm
bsw = 1. / lamd_grid
# print np.shape(lamd_grid)

tcrit = 2.  # changed by HY. 2016-7-16
zsno = 0.024  # changed by HY. 2016-8-7
# zlnd  = 0.01



#-------------set parameter-------------------------------------------
kground=np.zeros([newrows,newcols])+0.8*0.01/24./3600.0
ssf = np.zeros((newrows,newcols)) + 0.06
adjfac = 0.1


temp_year=0
temp_month=0
start=1

RWData=True

for num,i_date in enumerate(dates):
    print (i_date)
    iyear=i_date.year
    imonth=i_date.month
    iday=i_date.day
    idays=ta[(ta.year==iyear)&(ta.month==imonth)][0].day
    idaye=ta[(ta.year==iyear)&(ta.month==imonth)][-1].day

    idc, ihc, idate=Tstamp2idate(i_date,deltim)

    if (iyear % 4==0):
        dayinmonth=np.array([31,29,31,30,31,30,31,31,30,31,30,31])
    else:
        dayinmonth=np.array([31,28,31,30,31,30,31,31,30,31,30,31])


    if  RWData:
        if record == 1:

            glw_nc, psfc_nc, swd_nc, airt_nc, U10m_nc, \
            V10m_nc, rh_nc, Q2_nc, prec_nc, ncTime=\
                ReadMeteoDeriver([iyear,imonth],[meteoDriverPath],name='ERA_Interim')

            newnc = Dataset(out_path + str(iyear) + '%02d' % (imonth) + '_result' + endletter + '.nc', \
                            'w', format='NETCDF4_CLASSIC')
            newnc.description = 'Manas River basin'
            newnc.creattime = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

            # dimensions
            newnc.createDimension('time', None)
            newnc.createDimension('layer', nl_soil + 1)
            newnc.createDimension('LAT', newrows)
            newnc.createDimension('LON', newcols)

            # variables
            latgrid = newnc.createVariable('latitude', 'f4', ('LAT', 'LON'))
            longrid = newnc.createVariable('longitude', 'f4', ('LAT', 'LON'))
            timegrid = newnc.createVariable('time', 'f4', ('time'))

            snowdep = newnc.createVariable('snowd', 'f4', ('time', 'LAT', 'LON',))
            swe = newnc.createVariable('swe', 'f4', ('time', 'LAT', 'LON',))
            runoff_sur = newnc.createVariable('rsfc', 'f4', ('time', 'LAT', 'LON',))
            runoff_tot = newnc.createVariable('rtot', 'f4', ('time', 'LAT', 'LON',))
            runoff_grd = newnc.createVariable('rgrd', 'f4', ('time', 'LAT', 'LON',))
            runoff_lat = newnc.createVariable('qlat', 'f4', ('time', 'LAT', 'LON',))
            qinflrate = newnc.createVariable('qinfl', 'f4', ('time', 'LAT', 'LON',))
            runoff_sm = newnc.createVariable('qsmelt', 'f4', ('time', 'LAT', 'LON',))
            sublrate = newnc.createVariable('subl', 'f4', ('time', 'LAT', 'LON',))
            frosrate = newnc.createVariable('fros', 'f4', ('time', 'LAT', 'LON',))
            evaporate = newnc.createVariable('evpa', 'f4', ('time', 'LAT', 'LON',))
            etrate = newnc.createVariable('et', 'f4', ('time', 'LAT', 'LON',))
            evplrate = newnc.createVariable('evpl', 'f4', ('time', 'LAT', 'LON',))
            evpgrate = newnc.createVariable('evpg', 'f4', ('time', 'LAT', 'LON',))
            sengrate = newnc.createVariable('seng', 'f4', ('time', 'LAT', 'LON',))
            senlrate = newnc.createVariable('senl', 'f4', ('time', 'LAT', 'LON',))
            senarate = newnc.createVariable('sena', 'f4', ('time', 'LAT', 'LON',))
            lenarate = newnc.createVariable('lena', 'f4', ('time', 'LAT', 'LON',))
            fluxgrnd = newnc.createVariable('fgrnd', 'f4', ('time', 'LAT', 'LON',))
            tempgrnd = newnc.createVariable('tgrnd', 'f4', ('time', 'LAT', 'LON',))
            ref_srad = newnc.createVariable('refsrad', 'f4', ('time', 'LAT', 'LON',))
            forcrain = newnc.createVariable('rain', 'f4', ('time', 'LAT', 'LON',))
            gstorage = newnc.createVariable('gstroage', 'f4', ('time', 'LAT', 'LON',))
            prec = newnc.createVariable('prec', 'f4', ('time', 'LAT', 'LON',))
            forcsnow = newnc.createVariable('snow', 'f4', ('time', 'LAT', 'LON',))
            snowpixel = newnc.createVariable('snowpixel', 'f4', ('time', 'LAT', 'LON',))
            snowevp = newnc.createVariable('snowevp', 'f4', ('time', 'LAT', 'LON',))
            qchargegrid = newnc.createVariable('qchargegrid', 'f4', ('time', 'LAT', 'LON',))
            smsoil = newnc.createVariable('smsoil', 'f4', ('time', 'LAT', 'LON',))
            smice = newnc.createVariable('smice', 'f4', ('time', 'LAT', 'LON',))
            soilw = newnc.createVariable('soilw', 'f4', ('time', 'layer', 'LAT', 'LON',))
            soilt = newnc.createVariable('soilt', 'f4', ('time', 'layer', 'LAT', 'LON',))
            snowsmf = newnc.createVariable('snowsmf', 'f4', ('time', 'layer', 'LAT', 'LON',))
            icesmf = newnc.createVariable('icesmf', 'f4', ('time', 'layer', 'LAT', 'LON',))
            if (run_pbsm == 1):
                blow_snow = newnc.createVariable('blow', 'f4', ('time', 'LAT', 'LON',))
                blow_sub = newnc.createVariable('blowsub', 'f4', ('time', 'LAT', 'LON',))

            days=ta[ta.month==imonth]
            timegrid[:] =date2num(days.to_pydatetime(),units='days since 2010-01-01 00:00:00')
            timegrid.units="days since 2010-01-01 00:00:00"
            timegrid.calendar='standard'
            longrid[:]=longxy
            latgrid[:]=latixy


        RWData=False


    if record == 1:
        dailysnowdp = tmpzero2d * 0.
        dailyscv = tmpzero2d * 0.
        dailyrsur = tmpzero2d * 0.
        dailyrnof = tmpzero2d * 0.
        dailyqlat = tmpzero2d * 0.
        dailyqsmelt = tmpzero2d * 0.
        dailysubl = tmpzero2d * 0.
        dailyfros = tmpzero2d * 0.
        dailyqground = tmpzero2d * 0.
        dailyfevpa = tmpzero2d * 0.
        dailyfevpl = tmpzero2d * 0.
        dailyetr = tmpzero2d * 0.
        dailyfevpg = tmpzero2d * 0.
        dailyseng = tmpzero2d * 0.
        dailysena = tmpzero2d * 0.
        dailysenl = tmpzero2d * 0.
        dailysmsoil = tmpzero2d * 0.
        dailysmice = tmpzero2d * 0.
        dailyforc_rain = tmpzero2d * 0.
        dailywa = tmpzero2d * 0.
        dailyforc_prec = tmpzero2d * 0.
        dailyforc_snow = tmpzero2d * 0.
        dailysoilw = tmpzero3d * 0.
        dailysoilt = tmpzero3d * 0.
        dailyfgrnd = tmpzero2d * 0.
        dailytgrnd = tmpzero2d * 0.
        dailylena = tmpzero2d * 0.
        dailysr = tmpzero2d * 0.
        dailysnowevp = tmpzero2d * 0.
        dailyqchargegrid = tmpzero2d * 0.
        dailyqinfl = tmpzero2d * 0.


    for ihour in range(0, 24):
        ihc = ihc + 1
        idate = np.array([iyear, idc, (ihc - 0.5) * deltim])

        nchour=int((datetime(iyear,imonth,iday,ihour)-ncTime[0]).total_seconds()/deltim)
        oro = np.ones((newrows, newcols))
        rssca = np.ones([newrows, newcols])

        forc_t = airt_nc[nchour]
        forc_q = Q2_nc[nchour]
        forc_psrf = psfc_nc[nchour]
        forc_solarin = swd_nc[nchour]
        forc_frl = glw_nc[nchour]
        forc_prec = prec_nc[nchour]
        forc_us = U10m_nc[nchour]
        forc_vs = V10m_nc[nchour]
        forc_rh = rh_nc[nchour]

        oro, z_soisno, dz_soisno, t_soisno, wliq_soisno, wice_soisno, sst, smf_soisno, smf_ice, t_grnd, tlsun, tlsha, \
        ldew, sag, scv, snowdp, fveg, fsno, sigf, green, lai, sai, coszen, albg, albv, alb, ssun, ssha, \
        thermk, extkb, extkd, zwt, wa, t_lake, lake_icefrac, laisun, laisha, region_qsub, net_blow, \
        blowing_add_region, rstfac, h2osoi, wat, qlat, forc_rain, forc_snow, taux, tauy, fsena, \
        fevpa, lfevpa, fsenl, fevpl, etr, fseng, fevpg, olrg, fgrnd, trad, tref, qref, rsur, rnof, qintr, \
        qinfl, qdrip, rst, assim, respc, sabvsun, sabvsha, sabg, sr, solvd, solvi, solnd, solni, srvd, \
        srvi, srnd, srni, solvdln, solviln, solndln, solniln, srvdln, srviln, srndln, srniln, qcharge, \
        xerr, zerr, drw_gbhm, drw, qground, qsmelt, qsubl, qfros, snoweva, start, qr1, qlin1, qd, qm, q1, qr2, q2, qh \
        = pyice.core(
            mask, run_pbsm, nl, nl_soil, maxsnl, nl_lake, newrows, newcols, deltim, oro,\
            lons, dlon, dlat, ivt, itypwat, lakedepth, dz_lake, soil_s_v_alb, soil_d_v_alb, soil_s_n_alb,\
            soil_d_n_alb, wsat, wrsd, watern, alpha, slope, length, elev, ds, Dr_grid, dg, kground,\
            porsl, psi0, bsw,\
            hksati, wfld_grid, csol, dksatu, dkdry, rootfr, z0m, displa, sqrtdi, effcon, vmax25, slti, hlti,\
            shti, hhti, trda, trdm, trop, gradm, binter, extkn, chil, ref, tran, forc_rh, forc_t, forc_q,\
            forc_psrf, forc_solarin, forc_frl, forc_prec, forc_us, forc_vs, idate, z_soisno, dz_soisno,\
            t_soisno, wliq_soisno, wice_soisno, sst, smf_soisno, smf_ice, t_grnd, tlsun, tlsha, ldew, sag, scv, snowdp,\
            fveg,\
            fsno, sigf, green, lai, sai, coszen, albg, albv, alb, ssun, ssha, thermk, extkb, extkd, zwt,\
            wa, t_lake, lake_icefrac, zlnd, zsno, csoilc, dewmx, wtfact, capr, cnfac, ssi, wimp, pondmx,\
            smpmax, smpmin, trsmx0, tcrit, adjfac, drw_gbhm, run_gbhm, subbasin, area, dx, dr_gbhm,\
            wriver, s0, mnroughness, nflow, psubbasin, nbasinup, nsub, ngrid, grid_row, grid_col,\
            start, iyear, imonth, iday, ihour + 1, nc, startyear, endyear, dayinmonth, startmonth, endmonth,\
            startday, endday, idc, ihc, start_sub, end_sub, qr1, qlin1, qd, qm, q1, qr2, ssf, rssca)

        if record == 1:
            dailysnowdp  = snowdp/24. + dailysnowdp
            dailyscv     = scv/24.    + dailyscv
            dailyrsur    = rsur/24.   + dailyrsur
            dailyrnof    = rnof/24.   + dailyrnof
            dailyqlat    = qlat/24.   + dailyqlat
            dailyqsmelt  = qsmelt/24. + dailyqsmelt
            dailysubl    = qsubl/24.  + dailysubl
            dailyfros    = qfros/24.  + dailyfros
            dailyqground = qground/24.+ dailyqground
            dailyfevpa   = fevpa/24.  + dailyfevpa
            dailyfevpl   = fevpl/24.  + dailyfevpl
            dailyetr     = etr/24.    + dailyetr
            dailyfevpg   = fevpg/24.  + dailyfevpg
            dailyseng    = fseng/24.  + dailyseng
            dailysena    = fsena/24.  + dailysena
            dailysenl    = fsenl/24.  + dailysenl
            dailyfgrnd   = fgrnd/24.  + dailyfgrnd
            dailytgrnd   = t_grnd/24. + dailytgrnd
            dailylena    = lfevpa/24. + dailylena
            dailysr      = sr/24. + dailysr

            dailyforc_rain = forc_rain / 24. + dailyforc_rain
            dailywa = wa / 24. + dailywa
            dailyforc_prec = forc_prec / 24. + dailyforc_prec
            dailyforc_snow = forc_snow / 24. + dailyforc_snow

            tmpsmf = smf_soisno[6:nl_soil + 7, :, :]
            meltsoil = np.sum(tmpsmf[0:nl_soil, :, :] * wliq_soisno[5:nl_soil + 5, :, :], axis=0) \
                       + tmpsmf[nl_soil, :, :] * wa
            meltice = np.sum(smf_ice[0:nl_soil, :, :] * wice_soisno[5:nl_soil + 5, :, :], axis=0)
            dailysmsoil = meltsoil / 24. + dailysmsoil
            dailysmice =  meltice / 24. + dailysmice

            dailysnowevp = snoweva / 24. + dailysnowevp
            dailyqchargegrid = qcharge / 24. + dailyqchargegrid
            if ihour == 10:
                dailysca = snowdp * 0.
                dailysca = np.where(snowdp > 0.005, 1., 0.)
            dailysoilt = dailysoilt + t_soisno[5:nl_soil + 5, :, :] / 24.
            dailysoilw = dailysoilw + wliq_soisno[5:nl_soil + 5, :, :] * 0.001 / dz_soisno[5:nl_soil + 5, :, :] / 24.

        testp = testp + forc_prec
        testeva = testeva + fevpa
        testrnf = testrnf + rnof
        # if ihour==22: print wliq_soisno[5,78,110],snowdp[78,110]
        # if ihour==23: print fevpa[78,110] *3600
        maskeva = testeva[mask > 0]
        maskp = testp[mask > 0]
        maskrnf = testrnf[mask > 0]
        maskqchargegrid = qcharge[mask > 0]
        maskqground = qground[mask > 0]
        maskqlat = qlat[mask > 0]
        maskrsur = rsur[mask > 0]
        maskrnof = rnof[mask > 0]

    if record == 1:
        snowdep   [iday-idays,:,:] = dailysnowdp
        swe       [iday-idays,:,:] = dailyscv
        runoff_sur[iday-idays,:,:] = dailyrsur
        runoff_tot[iday-idays,:,:] = dailyrnof
        runoff_lat[iday-idays,:,:] = dailyqlat
        qinflrate [iday-idays,:,:] = dailyqinfl
        sublrate  [iday-idays,:,:] = dailysubl
        frosrate  [iday-idays,:,:] = dailyfros
        runoff_sm [iday-idays,:,:] = dailyqsmelt
        runoff_grd[iday-idays,:,:] = dailyqground
        evaporate [iday-idays,:,:] = dailyfevpa
        evplrate  [iday-idays,:,:] = dailyfevpl
        etrate    [iday-idays,:,:] = dailyetr



        evpgrate  [iday-idays,:,:] = dailyfevpg
        sengrate  [iday-idays,:,:] = dailyseng
        senlrate  [iday-idays,:,:] = dailysenl
        senarate  [iday-idays,:,:] = dailysena
        snowpixel [iday-idays,:,:] = dailysca
        lenarate  [iday-idays,:,:] = dailylena
        ref_srad  [iday-idays,:,:] = dailysr
        fluxgrnd  [iday-idays,:,:] = dailyfgrnd
        tempgrnd  [iday-idays,:,:] = dailytgrnd
        forcrain  [iday-idays,:,:] = dailyforc_rain
        gstorage  [iday-idays,:,:] = dailywa
        prec      [iday-idays,:,:] = dailyforc_prec
        forcsnow  [iday-idays,:,:] = dailyforc_snow
        qchargegrid   [iday-idays,:,:] = dailyqchargegrid
        snowevp   [iday-idays,:,:] = dailysnowevp
        smsoil    [iday-idays,:,:] = dailysmsoil
        smice     [iday-idays,:,:] = dailysmice
        soilt [iday-idays,0:nl_soil,:,:] = dailysoilt
        soilw [iday-idays,0:nl_soil,:,:] = dailysoilw
        icesmf[iday-idays,0:nl_soil,:,:] = smf_ice[0:nl_soil,:,:]
        snowsmf[iday-idays,:,:,:] = tmpsmf
        if (run_pbsm ==idays):
            blow_snow [iday-idays,:,:] = dailynet_blow
            blow_sub  [iday-idays,:,:] = dailyregion_qsub


        if iday == idaye:
            newnc.close()
            RWData=True
