!  GREGKW -- keyword-layer control codes for Greg Johnson hook family.
!  INCLUDE in any routine that needs to read or set hook-driver codes.
!  Initialised to -1 (unset) by BLOCK DATA GREGKWBD in greghghg.f90;
!  set by keyword handlers in initre.f90 (DGDRIVER, HGDRIVER, CROWNDRIVER).
!
!  IDGDRV   -- DG site driver code, -1 unset; 0 none,1 elev,2 cspi,3 bgi,4 esi,5 emt
!  IMORTDRV -- mortality coefficient TIER code, -1 unset; 0 crown,1 size,2 size+BGI
!  IHGDRV   -- HG arm code, -1 unset; 1 Bayes-calibrated regional,2 Greg ORGANON CONUS,
!              3 CONUS species-dependent,4 CONUS species-independent
!  ICRWDRV  -- crown recession driver code, -1 unset; 1 Greg CONUS recession model
INTEGER IDGDRV, IMORTDRV, IHGDRV, ICRWDRV
COMMON /GREGKW/ IDGDRV, IMORTDRV, IHGDRV, ICRWDRV
