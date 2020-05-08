MODULE topoConst

!=======================================================================
! Topographical constants 
!=======================================================================

  use precision
  IMPLICIT NONE

  public
  real(r8), parameter :: &
  aniks(25)=(/   1.0,&        ! 0  Ocean
             1.0,&        ! 1  Urban and Built-Up Land
             1.0,&        ! 2  Dryland Cropland and Pasture
             1.0,&        ! 3  Irrigated Cropland and Pasture
             1.0,&        ! 4  Mixed Dryland/Irrigated Cropland and Pasture
             5.0,&        ! 5  Cropland/Grassland Mosaic
             5.0,&        ! 6  Cropland/Woodland Mosaic
             5.0,&        ! 7  Grassland
             3.0,&        ! 8  Shrubland
             5.0,&        ! 9  Mixed Shrubland/Grassland
             5.0,&        !10  Savanna
             5.0,&        !11  Deciduous Broadleaf Forest 
             5.0,&        !12  Deciduous Needleleaf Forest 
             5.0,&        !13  Evergreen Broadleaf Forest
             5.0,&        !14  Evergreen Needleleaf Forest
             5.0,&        !15  Mixed Forest
             1.0,&        !16  Inland Water
             1.0,&        !17  Herbaceous Wetland
             1.0,&        !18  Wooded Wetland
             3.0,&        !19  Barren or Sparsely Vegetated
             3.0,&        !20  Herbaceous Tundra
             3.0,&        !21  Wooded Tundra
             3.0,&        !22  Mixed Tundra
             3.0,&        !23  Bare Ground Tundra
             3.0  /)      !24  Snow or Ice
             
  real(r8), parameter :: &           
  sstmaxs(25)=(/  0.0,&        ! 0  Ocean
             10.0,&        ! 1  Urban and Built-Up Land
             15.0,&        ! 2  Dryland Cropland and Pasture
             15.0,&        ! 3  Irrigated Cropland and Pasture
             15.0,&        ! 4  Mixed Dryland/Irrigated Cropland and Pasture
             25.0,&        ! 5  Cropland/Grassland Mosaic
             25.0,&        ! 6  Cropland/Woodland Mosaic
             20.0,&        ! 7  Grassland
             20.0,&        ! 8  Shrubland
             20.0,&        ! 9  Mixed Shrubland/Grassland
             15.0,&        !10  Savanna
             25.0,&        !11  Deciduous Broadleaf Forest 
             25.0,&        !12  Deciduous Needleleaf Forest 
             25.0,&        !13  Evergreen Broadleaf Forest
             25.0,&        !14  Evergreen Needleleaf Forest
             25.0,&        !15  Mixed Forest
              0.0,&        !16  Inland Water
             10.0,&        !17  Herbaceous Wetland
             10.0,&        !18  Wooded Wetland
             15.0,&        !19  Barren or Sparsely Vegetated
             15.0,&        !20  Herbaceous Tundra
             15.0,&        !21  Wooded Tundra
             15.0,&        !22  Mixed Tundra
             15.0,&        !23  Bare Ground Tundra
              0.0  /)      !24  Snow or Ice  
              
  real(r8), parameter :: &            
  surfns(25) =(/ 0.02,&        ! 0  Ocean
              1.0,&        ! 1  Urban and Built-Up Land
              0.1,&        ! 2  Dryland Cropland and Pasture
              0.1,&        ! 3  Irrigated Cropland and Pasture
              0.1,&        ! 4  Mixed Dryland/Irrigated Cropland and Pasture
              0.1,&        ! 5  Cropland/Grassland Mosaic
              0.1,&        ! 6  Cropland/Woodland Mosaic
              0.1,&        ! 7  Grassland
              0.2,&        ! 8  Shrubland
             0.15,&        ! 9  Mixed Shrubland/Grassland
              0.1,&        !10  Savanna
              0.2,&        !11  Deciduous Broadleaf Forest 
              0.2,&        !12  Deciduous Needleleaf Forest 
              0.2,&        !13  Evergreen Broadleaf Forest
              0.2,&        !14  Evergreen Needleleaf Forest
              0.2,&        !15  Mixed Forest
             0.02,&        !16  Inland Water
              0.1,&        !17  Herbaceous Wetland
              0.1,&        !18  Wooded Wetland
              0.2,&        !19  Barren or Sparsely Vegetated
              0.1,&        !20  Herbaceous Tundra
              0.1,&        !21  Wooded Tundra
              0.1,&        !22  Mixed Tundra
              0.1,&        !23  Bare Ground Tundra
              0.1  /)      !24  Snow or Ice   
END MODULE topoConst
