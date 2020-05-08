module blowing_varpar
!----------------------------------------------------------------------------------------
! Variables and parameters of snow pack.
!
! Author: Hongyi Li (lihongyi@lzb.ac.cn)
! Created Date: 2014-10-15;
! Updated Date(version 1.0): 2014-10-30
! Updated Date(version 1.1): 2014-11-22
!----------------------------------------------------------------------------------------
use precision
    implicit none
    save

! !PUBLIC DATA MEMBERS:
    public
    real(r8), parameter :: grid_size = 1000. !m
    !_____________________
    !physical constants
    ! real(8), parameter :: Lf = 333600.       !Latent heat of fusion
    ! real(8), parameter :: Lv = 2510400.      !Latent heat of evaporation
    real(r8), parameter :: Liv = 2844000.     !Latent heat of sublimation

    !_____________________
    !observation instruments parameters
    ! real(8), parameter :: zt = 2.0   !reference height for air temperature observation, /m
    ! real(8), parameter :: zq = 2.0   !reference height for humidity observation, /m
    ! real(8), parameter :: zw = 10.0  !reference height for wind speed observation, /m
    !_____________________
    !adjustable parameters

    ! real(8), parameter :: time_step = 3600.!*2*24 !half-hourly, /sec
    ! real(8), parameter :: ci = 2102. !heat capacity of ice, J/kg/K
    ! real(8), parameter :: cl = 4188. !heat capacity of liquid, J/kg/K
    ! real(8), parameter :: cv = 1870. !heat capacity of vapor, J/kg/K

    ! real(8), parameter :: liquid_density = 1000. !kg.m-3
    ! real(8), parameter :: ice_density = 917.

    ! real(8), parameter :: emissivity = 0.98
    ! real(8), parameter :: sbc = 5.6704E-8            !Stefan-Boltzman constant

    ! real(8), parameter :: sr = 0.04
    ! real(8), parameter :: ul = 1.787/1000.

    ! real(8), parameter :: min_value = 10.E-10
    ! real(8), parameter :: Tf = 273.15
    !calibrated parameters
!    real(8), parameter :: zo = 0.00002   !roughness length for snow surface, /m [0.00001,0.004]
    ! real(8)  :: zo=0.008  !roughness length for snow surface, /m [0.00001,0.004]
    ! real(8)  :: Krho=1.0    !1.   !Adjustable proportionality factor for new snow density, [0.8,1.2]
!    real(8), parameter :: Krho = 1.   !Adjustable proportionality factor for new snow density, [0.8,1.2]
    ! real(8)  :: adj_crp = 2.           !Critical temperature for defining snow/rain
    ! real(8), parameter :: KKS = 1.   !Adjustable proportionality factor for KS, [0.5,1.5]
    ! real(8), parameter :: VIR0 = 0.98 !reflective of new snow at VIR band, [0.85,1.0],0.98
    ! real(8), parameter :: NIR0 = 0.70 !reflective of new snow at NIR band, [0.5,0.8],0.70
    ! real(8), parameter :: local_lon = 120.  !longtitude of local time zone, such as 120 for beijing time
    ! integer, parameter :: max_layer = 5 ! set 5 layers as the maximum layer number
end module blowing_varpar

