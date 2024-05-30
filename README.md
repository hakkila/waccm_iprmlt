# IPRMLT modification for WACCM
MLT-dependency coding for the Whole Atmosphere Community Climate Model (WACCM),
version 6.

Read magnetic local time (MLT) dependent forcing into WACCM-D.

Created by Tuomas Häkkilä, Finnish Meteorological Institute.
Based on previous modifications by Christine Smith-Johnsen (replacing aurora), and
Pekka Verronen and Daniel Marsh (MLT-dependency for WACCM 4).

## Versions
The IPRMLT code comes in two versions: (1) the default version, which can be used
to add MLT-dependent ionization forcing input to WACCM simulations, and (2) the
aurora version, which turns off the Kp-driven parameterization of auroral forcing
within WACCM, so that the MLT input can be used to replace it.

## Contents
The following modified files are included:
- chemistry.F90
    - Modified to include MLT dependent forcing. Largely based on the default
      chemistry.F90 file, as well as modifications by Christine Smith-Johnsen.
    - Calls iprmlt_readnl, iprmlt_init, iprmlt_adv from mo_iprmlt.F90.
- mo_setext.F90
    - Modified to include MLT dependent forcing. Largely based on the default
      chemistry.F90 file, as well as modifications by Christine Smith-Johnsen
      and by Daniel Marsh and Pekka Verronen.
    - The aurora-version excludes the default Kp-index parameterized aurora
      (to be replaced by the MLT dependent forcing) 
    - Calls iprmlt_ionization_noxhox from mo_iprmlt.F90
- mo_iprmlt.F90 (code that reads the MLT dependent forcing)
    - Code that reads the MLT dependent forcing. Updated, based on work by
      Daniel Marsh and Pekka Verronen.
    - Includes the following subroutines: iprmlt_readnl, iprmlt_init,
      iprmlt_adv, iprmlt_ionization_noxhox.
- mo_iprmlt.F90 (code that reads the MLT dependent forcing)
    - Code that reads the MLT dependent forcing. Updated, based on work by
      Daniel Marsh and Pekka Verronen.
    - Includes the following subroutines: iprmlt_readnl, iprmlt_init,
      iprmlt_adv, iprmlt_ionization_noxhox.
- mo_aurora.F90 (makes sure there is no Kp aurora)
    - aurora-version only, used to make sure there is no Kp-driven auroral
      forcing. From Christine Smith-Johnsen.
- namelist_definition.xml
        Modifications to the namelist file to include IPRMLT fields:
        iprmlt_ionization_datapath, iprmlt_ionization_filename,
        iprmlt_ionization_filelist, iprmlt_ionization_fldname.
These files should be placed in SourceMods/src.cam/

## Usage
### Namelist definitions
namelist_definition.xml

A new group called 'iprmlt_ionization_nl' is created for user_nl_cam. The
following fields within this group can be defined:
- 'iprmlt_ionization_datapath': Directory of the forcing data files.
- 'iprmlt_ionization_filename': Name of the first forcing file. The filename
  is relative to iprmlt_ionization_datapath.
- 'iprmlt_ionization_filelist': Name of the list of forcing files in
  chronological order (day-by-day). The filename is relative to
  iprmlt_ionization_datapath.
- 'iprmlt_ionization_fldname': Name of the field within the forcing datafiles
  that contains the forcing. This is set to 'ipr' by default, and need
  not be included in user_nl_cam.

Note that other modifications (e.g. cyclical year) are not supported by this
modification.

### Forcing data
mo_iprmlt.F90
    
The forcing files must be such that each file includes the forcing for only
one day. The forcing can have dimensions (lev, lshell, mlt, time) or
(lev, glat, mlt, time) (this is the order in which Fortran (and eg. Matlab)
reads the NetCDF dimensions. The 'default' NetCDF dimension order is
reversed). If each file only includes one timestep, the forcing does not
need to have a time dimension.

- MLT: required fields:
    - mlt: Magnetic local time
        - units: seconds
    - mlt_width: MLT bin widths
        - units: seconds
        - size: (length(mlt))
- LSHELL: required fields:
    - lshell: McIlwain L shell number
        - units: Earth radii
    - lshell_width: L-shell bin widths
        - units: Earth radii
        - size: (length(lshell))
- GLAT: required fields:
    - glat: Geomagnetic latitude
        - units: degrees north
    - glat_bnds: Geomagnetic latitude bin bounds
        - units: degrees north
        - size: (length(glat),2) (again, this is the order of dimensions in Fortran!)
- LEV: required fields:
    - lev: Pressure level
        - units: hPa
        - note: The code DOES NOT interpolate vertically, so this should
          match the WACCM grid.
- TIME: required fields:
    - datesec: Time of day
        - units: seconds (since beginning of day)
    - date: Date of file as int, format yyyymmdd
        - note: Can have same length as datesec, only first value will be read.

- IPR: the MLT dependent forcing data
    - field name: 'ipr', or defined in user_nl_cam with
      'iprmlt_ionization_fldname'
        - dimensions: (lev,lshell,mlt,time) OR (lev,glat,mlt,time) (OR, if only
          one timestep is included per day, (lev,lshell,mlt) OR
          (lev,glat,mlt))
        - units: ions/cm^3/s OR ions/g/s (conversion to ions/cm^3/s is included
          in the code) (other spellings of these units also work: /cm^3/s,
          /g/s, cm^-3s^-1, g^-1s^-1)

### Output

The modification enables the following output variables for WACCM:
- MLT: magnetic local time (h)
- MLTI: MLT index
- IPRMLT: MLT dependent ionization (ions/cm3/s)
- LSHELL: L-shell (Earth radii), if the input is on L-shell grid
- LSHELLI: L-shell index, if the input is on L-shell grid
- IPRgs: MLT dependent ionization (ions/g/s), if the input is in units ions/g/s
