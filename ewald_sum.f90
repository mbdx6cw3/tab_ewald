      PROGRAM elec_pot
      USE header_file
!     ***********************************************************************
!     ***********************************************************************
!     Description:
!       Program to calculate the electrostatic potential at a number of grid
!       points near a surface. Currently setup to calculate potential along z-
!       axis (x=0, y=0) through a porous system with 3D periodic boundaries. 
!       User must to fiddle with the write_tab subroutine to change this.
!
!     Notes:
!       1. The atomic coordinates of the adsorbent from a DL_POLY 'CONFIG' output 
!          file.
!       2. Ewald sum parameters and atomic partial charges are provided in the
!          input.dat file. 
!       3. Output from simulation is written to 'monitor.dat'.
!
!     Contains the subroutines:
!       1. kspace_ewald     - tabulates the electrostatic potential at grid points
!       2. write_tab        - calculates the reciprocal-space part of Ewald sum (surface-probe)
!       3. rspace_ewald     - calculates the real-space part of Ewald sum (surface-probe)
!       4. apply_mic        - does periodic boundary conditions
!       5. cell_dim         - returns the dimensional properties of a simulation cell
!
!     Contains the function:
!       1. erfc          - returns the complimentary error function.
!
!     Christopher D. Williams (2014)
!     ***********************************************************************
!     ***********************************************************************
      ! read in input parameters and conditios
      open(unit = 1, file = 'input.dat', status = 'old')
      read(unit = 1, nml = input_deck)
      close(unit = 1)

      open(unit = 10, file = 'monitor.dat', status = 'unknown')
      open(unit = 20, file = 'stats.dat', status = 'unknown')
      write(10,*) 'Input parameters read -----------------------------'
      write(10,*) '---------------------------------------------------'
      write(10,*) 'Adsorbent parameters'
      write(10,*)
      write(10,*) 'Oxygen atoms: '
      write(10,*) 'Charge = ',qO
      write(10,*)
      write(10,*) 'Silicon atoms: '
      write(10,*) 'Charge = ',qSi
      write(10,*)
      write(10,*) 'Silanol oxygen atoms: '
      write(10,*) 'Charge = ',qO_s
      write(10,*) 
      write(10,*) 'Nitrogen atoms: '
      write(10,*) 'Charge = ',qN
      write(10,*)
      write(10,*) 'EDA C atoms: '
      write(10,*) 'Charge = ',qC
      write(10,*)
      write(10,*) 'Tether C atoms: '
      write(10,*) 'Charge = ',qC_1
      write(10,*)
      write(10,*) 'Nitrogen hydrogen atoms: '
      write(10,*) 'Charge = ',qH_N
      write(10,*)
      write(10,*) 'Carbon hydrogen atoms: '
      write(10,*) 'Charge = ',qH_C
      write(10,*) '---------------------------------------------------'

      ! read in starting configuration from DL_POLY CONFIG
      open(unit = 3, file = "CONFIG", status = "old")
      read(3,*) name
      read(3,*) CK, BK, n, u

      allocate (cell(10),atype(n),aindex(n),rx(n),ry(n),rz(n),q(n))
      allocate (rxtab(n+1),rytab(n+1),rztab(n+1),qtab(n+1))     
 
      cell(:) = 0
      aindex(:) = 0
      rx(:) = 0.0d0
      ry(:) = 0.0d0
      rz(:) = 0.0d0
      q(:) = 0.0d0
      rxtab(:) = 0.0d0
      rytab(:) = 0.0d0
      rztab(:) = 0.0d0
      qtab(:) = 0.0d0 

      read(3,*) ax, ay, az
      read(3,*) bx, by, bz
      read(3,*) cx, cy, cz

      lx = ax  
      ly = by
      lz = cz

      lx2 = lx/2
      ly2 = ly/2
      lz2 = lz/2     
      
      ! read in atomic configuration and assign charges
      do i = 1, n
        read(3,*) atype(i), aindex(i)
        read(3,*) rx(i), ry(i), rz(i)
        !read(3,*) vxi, vyi, vzi
        !read(3,*) fxi, fyi, fzi
        if (atype(i) == 'O') then
          q(i) = qO
        else if (atype(i) == 'Si') then
          q(i) = qSi
        else if (atype(i) == 'sO') then
          q(i) = qO_s
        else if (atype(i) == 'H') then
          q(i) = qH
        end if
      end do
      close(unit=3)

      ! do something to adjust N and H_N charges

      ! calculate dimensional properties of simulation cell
      call cell_dim

      ! calculate volume
      vol = cell(10) 

      ! calculate total adsorbent charge
      qtot = sum(q)

      ! calculate real-space cutoff
      rcsq = rc * rc

      ! reset origin
      do i = 1, n
        rx(i) = rx(i) + lx2
        ry(i) = ry(i) + ly2
        rz(i) = rz(i) + lz2
      end do

      write(10,*) '---------------------------------------------------'
      write(10,*) 'CONFIG file read ----------------------------------'
      write(10,*) 'Total number of adsorbent atoms = ', n
      write(10,*) 'Volume of box = ', vol, 'ang**3'
      write(10,*) 'Total adsorbent charge = ', qtot
      write(10,*) '---------------------------------------------------'
      write(10,*)

      ! setup Ewald sum
      write(10,*) '--------------------------------------------------'
      write(10,*) 'Ewald Parameters ---------------------------------'

      ! determine k-vector cutoff
      kcut = min(dble(kmax1) * (1/cell(7)), dble(kmax2) * (1/cell(8)), dble(kmax3) * (1/cell(9)))
      kcut = 1.05D0 * twopi * kcut
      ksqcut = kcut * kcut
      write(10,*) 'Reciprocal Space Cutoff =', ksqcut
      write(10,*) 'Real Space Cutoff       =', rc
      write(10,*) 'kvectors in x,y,z       =', kmax1, kmax2, kmax3
      write(10,*) 'Gaussian width          =', kappa
      write(10,*)
      write(10,*) '-------------------------------------------------'
      write(10,*) 'Calculating Ewald Sum ---------------------------'
      write(10,*)

      ! reciprocal-space contribution to potential
      call kspace_ewald(q,rx,ry,rz)

      ! calculate self-interactions
      do i = 1, n
        ewald_self = ewald_self + (q(i) * q(i))
      end do
      ewald_self = -kappa  * ewald_self / root_pi  ! Angstrom ** -1

      ! calculate neutralising background charge correction
      if (qtot.ne.0.0d0) then
        write(10,*) 'System not neutral - use neutralising correction'
        ewald_neut = -1.0d0 * qtot * qtot * pi / 2.0d0 / vol / kappa / kappa   ! Angstrom ** -1
      end if

      ! calculate real-space part
      call rspace_ewald(q,rx,ry,rz)

      ! calculate total ewald sum 
      ewald_tot = real_sum + recip_sum + ewald_self + ewald_neut   ! Angstrom ** -1

      write(10,*) 'Total Ewald Sum               =', ewald_tot * 1.0d10 * qesq * r4pie0,  'J'
      write(10,*) 'Real Space Contribution       =', real_sum * 1.0d10 * qesq * r4pie0,   'J'
      write(10,*) 'Reciprocal Space Contribution =', recip_sum * 1.0d10 * qesq * r4pie0,  'J'      
      write(10,*) 'Self-Interaction Correction   =', ewald_self * 1.0d10 * qesq * r4pie0, 'J'
      write(10,*) 'Fuchs Correction              =', ewald_neut * 1.0d10 * qesq * r4pie0, 'J'
      write(10,*)
      write(10,*) '-------------------------------------------------'

      ! write new array for adsorbent plus a charge
      do i = 1, n
        qtab(i)  = q(i)
        rxtab(i) = rx(i)
        rytab(i) = ry(i)
        rztab(i) = rz(i)
      end do
 
      ! new total charge   
      qtot = qtot + 1.0d0

      write(10,*) 'Tabulating Ewald Sum'
      write(10,*)

      ! put this in a loop eventually
      !something to determine new positions

      qtab(n+1) = 1.0D0
      rxtab(n+1) = 0.0D0
      rytab(n+1) = 0.0D0
      rztab(n+1) = 0.0D0

      ! reciprocal-space contribution to potential
      call kspace_ewald(qtab,rxtab,rytab,rztab)

      ! calculate self-interactions
      ewald_self = 0.0D0
      do i = 1, n
        ewald_self = ewald_self + (qtab(i) * qtab(i))
      end do
      ewald_self = -kappa  * ewald_self / root_pi  ! Angstrom ** -1

      ! calculate neutralising background charge correction
      if (qtot.ne.0.0D0) then
        ewald_neut = -1.0d0 * qtot * qtot * pi / 2.0d0 / vol / kappa / kappa   ! Angstrom ** -1
      end if      

      ! calculate real-space part
      call rspace_ewald(qtab,rxtab,rytab,rztab)

      ! calculate total ewald sum
      ewald_tot_plus_tab = real_sum + recip_sum + ewald_self + ewald_neut   ! Angstrom ** -1

      write(10,*) 'Total Ewald Sum               =', ewald_tot_plus_tab * 1.0d10 * qesq * r4pie0,  'J'
      write(10,*) 'Real Space Contribution       =', real_sum * 1.0d10 * qesq * r4pie0,   'J'
      write(10,*) 'Reciprocal Space Contribution =', recip_sum * 1.0d10 * qesq * r4pie0,  'J'
      write(10,*) 'Self-Interaction Correction   =', ewald_self * 1.0d10 * qesq * r4pie0, 'J'
      write(10,*) 'Fuchs Correction              =', ewald_neut * 1.0d10 * qesq * r4pie0, 'J'

      ! calculate difference between this and orginal calculation
      ewald_tab = ewald_tot_plus_tab - ewald_tot

      write(10,*) 'Ewald Difference              =', ewald_tab * 1.0d10 * qesq * r4pie0, 'J'

      ! calculate potential from energy - may need to adjust to get correct units
      ewald_pot = ewald_tab 

      write(10,*) 'Ewald Potential              =', ewald_pot, 'V'

      ! write to monitor.dat
      write(10,*) rztab, ewald_pot

!     ***********************************************************************
      END PROGRAM      
!     -----------------------------------------------------------------------      



!     1.---------------------------------------------------------------------
      SUBROUTINE kspace_ewald(q,rx,ry,rz)
      USE header_file, ONLY: kappa,kmax1,kmax2,kmax3,maxk,ksqcut,&
     & vol,lx,ly,lz,n,dp,twopi
!     ***********************************************************************
!     Subroutine to calculate the reciprocal space part of the Ewald sum.
!   
!     Comments:
!     Wavevectors must fit into a box of unit length.
!     
!     k is a vector in reciprocal space
!     calculate and store exp(-rksq/b)/rksq terms for each k-vector
!     compute k = 2 * pi / L (l,m,n)
!     kx sum starts at zero to exploit symmetry
!
!     ewald_recip is the reciprocal space energy (in Joules) due to surface
!
!     Christopher D. Williams (Apr. 2013)
!     ***********************************************************************
      implicit none
      integer :: i, j, kx, ky, kz, totk
      real(dp) :: rkx, rky, rkz, r4kappa2, kappa2, d 
      real(dp) :: rksq, arg, rrksq, ksq
      real(dp) :: eikr, kvec_sum, kvec, factor
      real(dp) :: sum_real, sum_imag, eikr_real, eikr_imag
      real(dp) :: recip_sum     

      kappa2 = kappa * kappa
      r4kappa2 = 1.0D0 / (4.0d0 * kappa2)

      ! zero accumulators/counters      
      recip_sum = 0.0D0
      totk = 0
      
      ! loop over all k-vectors
      do kx = 0, kmax1
        
        rkx = twopi * dble(kx)/lx

        do ky = -kmax2, kmax2
        
          rky = twopi * dble(ky)/ly

          do kz = -kmax3, kmax3
              
            rkz = twopi * dble(kz)/lz
            
            rksq = rkx * rkx + rky * rky + rkz * rkz
            ksq = kx * kx + ky * ky + kz * kz            
              
            ! restrict sum to within cutoff and do not consider central cell
            if ((rksq.le.ksqcut).and.(ksq.ne.0.0D0)) then
              totk = totk + 1
              kvec_sum = 0.0D0
              sum_real = 0.0D0
              sum_imag = 0.0D0

              if (totk.gt.maxk) then
                write(10,*) 'kvec is too small - error!'
                stop

              else
                  
                ! factor means only half of k-space need be considered
                if(kx.eq.0) then
                  factor = 1.0d0
                else
                  factor = 2.0d0
                end if

                rrksq = 1.0D0 / rksq
                kvec = exp(-rksq * r4kappa2) * rrksq

                ! for this k-vector loop over all ions
                do i = 1, n
 
                  arg = rkx * rx(i) + rky * ry(i) + rkz * rz(i)
                  eikr_real = q(i) * cos(arg)         ! sum real part
                  eikr_imag = q(i) * sin(arg)
                  sum_real = sum_real + eikr_real
                  sum_imag = sum_imag + eikr_imag

                end do

              end if

              kvec_sum = factor * kvec * ((sum_real * sum_real) + (sum_imag * sum_imag))
              recip_sum = recip_sum + kvec_sum

            end if
            
          end do

        end do

      end do

      recip_sum = (twopi / vol) * recip_sum ! Angstrom ** -1
!     ************************************************************************
      END SUBROUTINE
!-----------------------------------------------------------------------------


!-----2.---------------------------------------------------------------------
      SUBROUTINE write_tab
      USE header_file, ONLY: lx,ly,lz,lx2,ly2,lz2,ugrid,zint,&
     & grid_res,recip_sum,dp,r4pie0,qe,ewald_self,ewald_neut,real_sum
!     ***********************************************************************
!     Subroutine that writes a 3D mesh for the simulation cell and at each
!     mesh point calculates the electrostatic potential due to the adsorbent
!     using the Ewald summation.
!
!     Christopher D. Williams (Feb. 2013)
!     ***********************************************************************
!     variable specification
      implicit none
      integer :: i, j, k
      real(kind=dp) :: rxg, ryg, rzg, ewald_tot
      integer :: ngridx, ngridy, ngridz

      write(10,*)
      write(10,*) 'Tabulating Ewald -----------------------------------'
      write(10,*) 'z-coordinate (ang.)         Electrostatic Potential'

      ngridx = 0
      ngridy = 0
      ngridz = ceiling(lz/grid_res)

      allocate(ugrid(0:ngridz))
      ugrid = 0.0d0

      zint = lz/ngridz

      do i = 0, ngridx

        rxg = 0.0d0

        do j = 0, ngridy

          ryg = 0.0d0

          ! calculate electric potential at this grid point
          do k = 0, ngridz

            rzg = -lz2 + k * zint

            ! real-space contribution to potential
            real_sum = 0.0d0
            !call rspace_ewald(rxg,ryg,rzg)

            ! total energy of +1 point charge at this site
            ewald_tot = real_sum + recip_sum + ewald_self + ewald_neut   ! Angstrom ** -1
            write(10,*)
            write(10,*) real_sum, recip_sum, ewald_self, ewald_neut
            ! work back to find potential at that point
            ugrid(i) = ewald_tot * 1.0d10  * qe * r4pie0      ! V
            write(10,*) rzg, ugrid(i)

          end do

        end do

      end do

      write(10,*) 'Ewald Tabulation Complete ---------------------------'
      write(10,*) '-----------------------------------------------------'
!     ***********************************************************************
      END SUBROUTINE
!----------------------------------------------------------------------------


!-----3.----------------------------------------------------------------------
      SUBROUTINE rspace_ewald(q,rx,ry,rz)
      USE header_file, ONLY: lx,ly,lz,kappa,real_sum,rcsq,n
!     ************************************************************************
!     Subroutine to calculate the real-space part of the Ewald sum. Requires 
!     the funtion 'ERFC' to calculate the complimentary error function.
!
!     Calculates the probe-adsorbent real space interaction only. Adsorbent -  
!     adsorbent real-space sum was calculated earlier.
!
!     Christopher D. Williams (Apr. 2013)
!     ************************************************************************
      implicit none
      integer :: i, j
      real*8 :: rxij, ryij, rzij, rijsq, rij, vij, krij

      ! zero accumulator
      real_sum = 0.0d0

      ! loop over all ions 
      do i = 1, n - 1
        
        do j = i + 1, n    

          rxij = rx(i) - rx(j)
          ryij = ry(i) - ry(j)
          rzij = rz(i) - rz(j)

          ! minimum image
          rxij = rxij - lx * nint(rxij/lx)
          ryij = ryij - ly * nint(ryij/ly)
          rzij = rzij - lz * nint(rzij/lz)

          rijsq = rxij * rxij + ryij * ryij + rzij * rzij

          ! accumulate real-space sum
          if (rijsq.lt.rcsq) then
            rij = sqrt(rijsq)
            krij = kappa * rij
            vij = q(i) * q(j) * erfc(krij) / rij
            real_sum = real_sum + vij  ! Angstrom ** -1
          end if

        end do

      end do

!     ************************************************************************
      END SUBROUTINE
!-----------------------------------------------------------------------------


!----4-.---------------------------------------------------------------------
     SUBROUTINE cell_dim
     USE header_file, ONLY: ax,ay,az,bx,by,bz,cx,cy,cz,cell
!    ************************************************************************
!    Calculates the dimensional properties of a simulation cell specified by
!    the input parameters ax - cz, returned in the cell matrix.
!
!    cell(1 to 3) - lengths of cell vectors
!    cell(4 to 6) - cosines of cell angles
!    cell(7 to 9) - perpendicular cell widths
!    cell(10)     - cell volume
!
!    Christopher D. Williams May. 2013
!
!    Adapted from DL_POLY Classic subroutine dcell
!    ************************************************************************
     implicit none
     real*8 :: axb1,axb2,axb3,bxc1,bxc2,bxc3,cxa1,cxa2,cxa3

     ! calculate lengths of cell vectors 
     cell(1)=sqrt(ax*ax+ay*ay+az*az)
     cell(2)=sqrt(bx*bx+by*by+bz*bz)
     cell(3)=sqrt(cx*cx+cy*cy+cz*cz)
      
     ! calculate cosines of cell angles
     cell(4)=(ax*bx+ay*by+az*bz)/(cell(1)*cell(2))
     cell(5)=(ax*cx+ay*cy+az*cz)/(cell(1)*cell(3))
     cell(6)=(bx*cx+by*cy+bz*cz)/(cell(2)*cell(3))
      
     ! calculate vector products of cell vectors
     axb1=ay*bz-az*by
     axb2=az*bx-ax*bz
     axb3=ax*by-ay*bx
     bxc1=by*cz-bz*cy
     bxc2=bz*cx-bx*cz
     bxc3=bx*cy-by*cx
     cxa1=cy*az-ay*cz
     cxa2=ax*cz-az*cx
     cxa3=ay*cx-ax*cy
      
     ! calculate volume of cell
     cell(10)=abs(ax*bxc1+ay*bxc2+az*bxc3)
      
     ! calculate cell perpendicular widths
     cell(7)=cell(10)/sqrt(bxc1*bxc1+bxc2*bxc2+bxc3*bxc3)
     cell(8)=cell(10)/sqrt(cxa1*cxa1+cxa2*cxa2+cxa3*cxa3)
     cell(9)=cell(10)/sqrt(axb1*axb1+axb2*axb2+axb3*axb3)
!    ************************************************************************
     END SUBROUTINE
!----------------------------------------------------------------------------



!-----FUNCTIONS--------------------------------------------------------------
!----------------------------------------------------------------------------
!-----1.---------------------------------------------------------------------
      REAL FUNCTION erfc(x)
!     ***********************************************************************
!     Approximation to the complimentary error function.
!
!     Reference: Abramowitz and Stegun, Handbook of Mathematical Functions,
!                National Bureau of Standards, Formula 7.1.26
!     ***********************************************************************
      implicit none
      real*8, parameter :: a1 = 0.254829592, a2 = -0.284496736 
      real*8, parameter :: a3 = 1.421413741, a4 = -1.453152027
      real*8, parameter :: a5 = 1.061405429, p = 0.3275911
      real*8 :: t, xsq, tp
      real*8, intent(in) :: x 

      t = 1.0 / (1.0 + p * x)
      xsq = x * x
      tp = t * (a1 + t * (a2 + t * (a3 + t * (a4 + t * a5))))
      erfc = tp * exp(-xsq)
!     ***********************************************************************
      END FUNCTION
!----------------------------------------------------------------------------
       


