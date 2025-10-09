These initialization files add options for setting up ideal cases in WRF

WRFv3.9x:
  module_initialize_quarter_ss.F (obviously only this one ideal setup)

WRFv4:
  module_initialize_ideal.F 
  

Namelist options (for WRFv4, notes indicate which specific ideal cases are applicable):

   NAMELIST /init_ideal/ nxc,nyc,delt,use_eta_levels, &
        ugrid,vgrid,setccn,use_snd_plevels,use_snd_zlevels,numbub,nxb,nyb, &
        xradbub,yradbub,zradbub,moisten,mtn_ht,ido_terrain,centlat, centlon, lm, &
        lu_input, vegfra, canwat, isltyp, sst_input, &
        real_soil,dz_bot, dz_top, nbndlyr, config_stretch_z

nxc,nyc : For single bubble, sets x,y location (in gridpoints). Default location is
      nxc = (ide-ids)/3 + 1
      nyc = (jde-jds)/2 + 1
      Use nxc = -2 , nyc = -2 to set to center of domain

ugrid, vgrid : Set grid motion (default values of 0)

use_eta_levels : If true, Sets vertical levels according to the eta_levels namelist item (default false; only for quarter_ss and squall2d)

use_snd_zlevels : If true, uses the z levels of the input sounding (number of levels must match; only for quarter_ss and squall2d) 
                  NOTE: The sounding must have a z=0 level. If not, the read subroutine will add one.

use_snd_plevels : If true, uses the pressure levels of the input sounding (number of levels must match) (only for quarter_ss and squall2d)

config_stretch_z = 0 ! (WRF4 only) 0 for 'equal spacing' (wrf default), 1 for geometric stretch using the following parameters:
    ! interpolates input sounding to the internally-generated z levels (close to them, anyway) 
    dz_bot = 200. ! grid spacing at z=0 
    dz_top = 500., ! = max stretched dz
    nbndlyr = 0 ! number of layers with dz_bot before stretching starts
    ! note that ztop = config_flags%ztop

numbub : number of thermal bubbles (default 1, up to maxbub=10) (only for quarter_ss)

nxb,nyb : gridpoint locations of bubbles (up to numbub) (only for quarter_ss)
    example:
    &init_ideal
      numbub = 3
      nxb = 60, 60, 60,
      nyb =  5, 15, 25,
    /


xradbub,yradbub,zradbub : spheroid radii (meters) (default values 10000, 10000, 1500.) (only for quarter_ss)

moisten : if true, add water vapor to bubble to preserve RH (default false) (only for quarter_ss)

mtn_hgt : Height of ideal mountain (default 500 m)

ido_terrain : default ido_terrain = 0, otherwise 1 = mountain, 2 = EW Ridge, 3 = NS Ridge

sst_input : set sst (default 300 K, for convrad only, used instead of theta_surf from input_sounding)

centlat, centlon : set latitude, longitude (affects radiation, coriolis, if turned on)

lm : sets half-width of land section for convrad (default = 0) and seabreeze2d_x (default = 25)

setccn (logical, defaualt .false.) : Set true to use loop to set the qnn field to a given 
           value or vertical profile (user-supplied code)

