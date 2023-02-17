      subroutine rdaer1

c  reads aerosol extinctions
c  on dynamics model lat grid and radiation pressure grid
c  (month,lat,press)

c  FOR HEATING RATE COMPUTATION

      logical aeread
      real lataer

      common/aerheat/nzaer,nltaer,paer(45),lataer(36),backext(12,36,45),
     * volcext(96,36,45),aeread

      data aeread/.false./

c  read background and volcanic aerosol
      open(unit=100,file='aero_heat2d_dbc_2.dat',status='old',
     * convert="big_endian", form='unformatted')
      read(100) fnz,fnlt
      nzaer = int(fnz)
      nltaer = int(fnlt)
      read(100) paer
      read(100) lataer
      read(100) volcext
      close(100)
      
cdbc use volcanic 1990 for background aerosol distribution

      do 10 i=1,12
      do 10 j=1,36
      do 10 k=1,45
         backext(i,j,k)=volcext(i,j,k)
 10   continue

      write(6,*) 'read aerosol 1 micron extinctions for heating'

      aeread = .true.

      return
      end

      subroutine rdaer2

c  reads aerosol extinctions
c  on chemistry model latitude and pressure grid, with extra level at surface

c  FOR MULTIPLE SCATTERING COMPUTATION

      real lataer
      common/aerphot/nzaer,nltaer,paer(47),lataer(18),backext(12,18,47),
     * volcext(96,18,47)
      open(unit=105, file='aero_phot2d_dbc.dat', status='old',
     * convert="big_endian", form='unformatted')
      read(105) fnz,fnlt
      nzaer=int(fnz)
      nltaer=int(fnlt)
      read(105) paer
      read(105) lataer
      read(105) volcext
      close(105)
      
cdbc use volcanic 1990 for background aerosol distribution

      do 10 i=1,12
      do 10 j=1,18
      do 10 k=1,47
         backext(i,j,k)=volcext(i,j,k)
 10   continue

      write(6,*) 'read aerosol 1 micron extinctions for photolysis'
      
      return
      end
