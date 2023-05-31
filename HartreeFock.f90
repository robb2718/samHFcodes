!  HartreeFock.f90 
!
!  FUNCTIONS:
!  HartreeFock - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: HartreeFock
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

program HartreeFock

implicit none

!####################################################################################################################################################
!####################################################################################################################################################

!read in inputs 
       
    double precision, parameter :: rho_inf = 0.1d0
    
    integer, parameter :: nmax = 18
    integer, parameter :: partnum = 66
    
    integer, parameter :: strmin = 25
    integer, parameter :: strmax = 25
    integer, parameter :: strstep = 5
       
    integer, parameter :: periodmin = 2
    integer, parameter :: periodmax = 2
    integer, parameter :: periodstep = 1
    
    integer, parameter :: accuracy = 1000
    
    double precision, parameter :: errortol = 1.D-6
    integer, parameter :: maxevals = 20
    
!####################################################################################################################################################
!####################################################################################################################################################    
    
!    !skyrme parameters from SLy4 model
!    double precision, parameter :: t0 = -2488.91d0 !MeV fm^3
!    double precision, parameter :: t1 = 486.82d0 !MeV fm^5
!    double precision, parameter :: t2 = -546.39d0 !MeV fm^5
!    double precision, parameter :: t3 = 13777.0d0 !MeV fm^(3+3*sigma)
!    double precision, parameter :: x0 = 0.834d0
!    double precision, parameter :: x1 = -0.344d0
!    double precision, parameter :: x2 = -1.0d0
!    double precision, parameter :: x3 = 1.354d0
!    double precision, parameter :: sigma = 1.0d0/6.0d0
!    double precision, parameter :: W0 = 123.0d0 !MeV fm^5

    !skyrme parameters for KDE0v1 model
    double precision, parameter :: t0 = -2553.1d0 !MeV fm^3
    double precision, parameter :: t1 =  411.7d0 !MeV fm^5
    double precision, parameter :: t2 = -419.9d0 !MeV fm^5
    double precision, parameter :: t3 =  14603.6d0 !MeV fm^(3+3*sigma)
    double precision, parameter :: x0 =  0.65d0
    double precision, parameter :: x1 = -0.35d0
    double precision, parameter :: x2 = -0.93d0
    double precision, parameter :: x3 =  0.95d0
    double precision, parameter :: sigma = 1.0d0/6.0d0

   !skyrme parameters for NRAPR
   ! double precision, parameter :: t0 = -2719.7d0 !MeV fm^3
   ! double precision, parameter :: t1 =  417.6d0 !MeV fm^5
   ! double precision, parameter :: t2 = -66.7d0 !MeV fm^5
   ! double precision, parameter :: t3 =  15042.0d0 !MeV fm^(3+3*sigma)
   ! double precision, parameter :: x0 =  0.16d0
   ! double precision, parameter :: x1 = -0.05d0
   ! double precision, parameter :: x2 =  0.03d0
   ! double precision, parameter :: x3 =  0.14d0
   ! double precision, parameter :: sigma = 0.14d0

!    !skyrme parameters for SkM*
!    double precision, parameter :: t0 = -2645.d0 !MeV fm^3
!    double precision, parameter :: t1 = 410.d0 !MeV fm^5
!    double precision, parameter :: t2 = -135.d0 !MeV fm^5
!    double precision, parameter :: t3 = 15595.d0 !MeV fm^(3+3*sigma)
!    double precision, parameter :: x0 = 0.09d0
!    double precision, parameter :: x1 = 0.d0
!    double precision, parameter :: x2 = 0.d0
!    double precision, parameter :: x3 = 0.d0
!    double precision, parameter :: sigma = 1.d0/6.d0
!    double precision, parameter :: W0 = 130.d0 !MeV fm^5

    !!skyrme parameters from SKRA model
    !double precision, parameter :: t0 = -2895.4 !MeV fm^3
    !double precision, parameter :: t1 = 405.5 !MeV fm^5
    !double precision, parameter :: t2 = -89.1 !MeV fm^5
    !double precision, parameter :: t3 = 16660.0 !MeV fm^(3+3*sigma)
    !double precision, parameter :: x0 = 0.08
    !double precision, parameter :: x1 = 0.00
    !double precision, parameter :: x2 = 0.20
    !double precision, parameter :: x3 = 0.0
    !double precision, parameter :: sigma = 0.14
    
    !parameters for mathieu function solutions
!    double precision, parameter :: t0 = 0.0 !MeV fm^3
!    double precision, parameter :: t1 = 0.0 !MeV fm^5
!    double precision, parameter :: t2 = 0.0 !MeV fm^5
!    double precision, parameter :: t3 = 0.0 !MeV fm^(3+3*sigma)
!    double precision, parameter :: x0 = 0.0
!    double precision, parameter :: x1 = 0.0
!    double precision, parameter :: x2 = 0.0
!    double precision, parameter :: x3 = 0.0
!    double precision, parameter :: sigma = 0.0
    
    
!##############################################################################################################################################################################################    
!##############################################################################################################################################################################################  
!##############################################################################################################################################################################################  
!##############################################################################################################################################################################################  
!##############################################################################################################################################################################################  

!constants

    !definition of Pi
    double precision, PARAMETER :: pi = 4*atan(1.0_16)
    
    !hbar
    double precision, parameter :: hbar = 6.58212D-22 !MeV*s
    !m (neutron mass)
    double precision, parameter :: m = 939.564133d0 !MeV/c^2
    !c (speed of light)
    double precision, parameter :: c = 2.99792D23 !fm/s
    !m0 (m/c^2)
    double precision, parameter :: m0 = m/(c**2)
    
    !hbar^2/2m
    double precision, parameter :: dispersion = 20.72125d0 !MeV fm^2 
    
    !kinetic density
    double precision :: tau_inf 
    
    !thermodynamic limit E/A (energy per particle)
    double precision :: ea
      
!##############################################################################################################################################################################################  
!##############################################################################################################################################################################################  
    
!accuracy parameter initializations

    !accuracy of the solver (also the size of the finite difference matrix)
    integer, parameter :: stencil = 5

!####################################################################################################################################################
!####################################################################################################################################################

!physical system variables 
    
    !box length and stepsize
    double precision :: length_0
    double precision :: length
    double precision :: stepsize

!####################################################################################################################################################
!####################################################################################################################################################
    
!external (cosinusoidal) potential parameters
    
    !fermi energy
    double precision :: E_f 
    
    !strength is ~0.5xFermi Energy
    double precision :: strengthfactor
    double precision :: strength 
    integer :: periodicity
    
!####################################################################################################################################################
!####################################################################################################################################################

!derived type definitions

    type :: states
    !stores eigenvalues and eigenvectors
        double precision :: eval
        double precision, dimension(accuracy) :: evec
    !stores k_x and k_y values (for use in kinetic energy density)
        integer :: nx
        integer :: ny
    !order at which the eigenvalue appears for given nx/ny pairs
        integer :: evalindex
    !IMPORTANT! this is not the norm of the wavefunction, this is the squared norm of the function
        double precision :: norm
    !multiplicity stores the number of times a wavefunction will be repeated (1 time for nx=ny=0, 4 times for nx=ny or nx/=0 and ny=0 or vice versa, and 8 times for nx/=ny/=0)
        integer :: multiplicity
    !check if the states form a closed shell
        logical :: closed
    end type states
 
!####################################################################################################################################################
!####################################################################################################################################################

!variable initilizations
    
    !finite difference matrix used at each step
    double precision, dimension(accuracy, accuracy) :: findiffmtx
    
    !array storing eigenpairs
    type(states), dimension(accuracy) :: eigenpairs
    type(states) :: Matt1, Matt2 !QOL variables
    !array storing all states to be used in calculation (has size numparticles/2 because of double multiplicity of particles)
    integer, parameter :: searchspacesize = 3 !determines the size of the search space to be used (over 7000 particles requires searchspacesize >=2, for 4224 particles searchspacesize = 1
    integer, parameter :: sortsize = searchspacesize*accuracy
    type(states), dimension(accuracy*searchspacesize) :: sortedstates, MattCopy
    integer :: appendloops !counting parameter for append loops to sortedstates array
    integer :: appendlooppos !position in append to start
    
    !nx and ny parameters as well as looping parameters
    integer :: nx, ny
    integer :: iMatt, jMatt, MattCount, MattUnique ! Looping and counting parameters for Matt's additions
    !number of times the code will be run 
    integer :: evals
    
    !ranges from strmin to strmax
    integer :: strengthindex
    
    !boolean variable checking if the calculation is the first one of the loop (used for the eigenvalue appending)
    logical :: appendcheck
    
    !densities
    double precision, dimension(accuracy) :: rho
    double precision, dimension(accuracy) :: rhoder
    double precision, dimension(accuracy) :: rho2der
    double precision, dimension(accuracy) :: tau
    
    !array reading parameter
    integer :: a_j
    
    !writing parameters to find where closed shells exist
    integer :: dummywrite
    integer :: dummypartnum
    
    !derivative of integrand function, change of variable array for change of variable
    double precision, dimension(accuracy) :: coeffIarray
    double precision, dimension(accuracy) :: coeffIprime
    double precision, dimension(accuracy) :: covarray
    integer :: changeindex
    
    !total energy
    double precision :: totalenergy
    double precision, dimension(accuracy) :: H_array
    
    !convergence testing
    double precision :: totalenergy_0
    double precision :: reldiff, reldiff_0
    double precision :: minima
    integer :: finaleval
    integer :: minimacheck  !checks two iterations after minima is found to make sure minima is true
    integer :: reldiffcheck  !checks two iterations after reasonable error range is found to make sure minima is true
    
!####################################################################################################################################################
!####################################################################################################################################################

    !file writing formats
    character(200) :: filename
    character(200) :: energyformat
    character(200) :: rhoformat
    character(200) :: tauformat
    
!####################################################################################################################################################
!####################################################################################################################################################

    !nmax only needs to be 5 for rho = 0.02,0.04,0.08 at 0.5EF strength for 600 particles at 1 period
    !IMPORTANT: nmax must be greater or equal the maximum n value on the grid (note k = 2pi*n/length is the wavenumber)
    !nmax = 5
    !tune this parameter to specify the maximum number of iterations used per loop as a cutoff
    !only need maxevals = 1 for mathieu solutions
    
    
!##############################################################################################################################################################################################    
!##############################################################################################################################################################################################  
!##############################################################################################################################################################################################  
!##############################################################################################################################################################################################  
!##############################################################################################################################################################################################  
    
    
!main code body

open(17,file="output.dat") ! Matt Addition

    do strengthindex = strmin, strmax, strstep
        do periodicity = periodmin, periodmax, periodstep
            
            !checks at specific particle numbers
            if(partnum == 8250)then   
                !checking specific periodicites (matches for 4224 particles)
                if((periodicity .ne. 1) .and. (periodicity .ne. 5) .and. (periodicity .ne. 10) .and. (periodicity .ne. 15) &
.and. (periodicity .ne. 20) .and. (periodicity .ne. 30) .and. (periodicity .ne. 40) .and. (periodicity .ne. 50))then
                    cycle
                end if  
            else if(partnum == 66)then             
                !checking specific periodicites (matches for 66 particles)
                if((periodicity .ne. 1) .and. (periodicity .ne. 2) .and. (periodicity .ne. 3) .and. (periodicity .ne. 4) &
 .and. (periodicity .ne. 6) .and. (periodicity .ne. 8) .and. (periodicity .ne. 10))then
                    cycle
                end if                  
            end if

            !skips all repeats at 0 strength
            if((strengthindex == 0) .and. (periodicity .ne. 1))then
                cycle
            end if
            
            !parameter initializations            
            tau_inf = (3.0d0/5.0d0)*(3.0d0*pi**2)**(2.0d0/3.0d0)*rho_inf**(5.0d0/3.0d0)
            E_f = dispersion*((3.0d0*(pi**2)*rho_inf)**(2.0d0/3.0d0))
        
            strengthfactor = strengthindex*0.01d0
            strength = strengthindex*0.01d0*E_f
            
            length_0 = (dble(partnum)/rho_inf)**(dble(1.0d0/3.0d0))
            length = length_0*(accuracy - 1)/(real(accuracy))
            stepsize = length/(real(accuracy-1))
        
            !array initializations
            rho = rho_inf
            rhoder = 0.0d0
            rho2der = 0.0d0
            tau = tau_inf
            findiffmtx = 0.0d0
        
            eigenpairs(:)%eval = 0.0d0
            do a_j = 1, accuracy
                eigenpairs(a_j)%evec(:) = 0.0d0
            end do
            eigenpairs(:)%nx = 0
            eigenpairs(:)%ny = 0
            eigenpairs(:)%evalindex = 0
            eigenpairs(:)%norm = 0.0d0
            eigenpairs(:)%multiplicity = 0.0d0
        
            sortedstates(:)%eval = 0.0d0
            do a_j = 1, accuracy
                sortedstates(a_j)%evec(:) = 0.0d0
            end do
            sortedstates(:)%nx = 0
            sortedstates(:)%ny = 0
            sortedstates(:)%evalindex = 0
            sortedstates(:)%norm = 0.0d0
            sortedstates(:)%multiplicity = 0.0d0
        
            covarray = 0.0d0
            changeindex = 0
    
            appendloops = 0
            appendlooppos = 0
        
            coeffIarray = 0.0d0
            coeffIprime = 0.0d0
        
            ea = (3.0d0/5.0d0)*dispersion*(3.0d0*(pi**2))**(2.0d0/3.0d0)*rho_inf**(2.0d0/3.0d0) &
            + (1.0d0/4.0d0)*t0*(1.0d0-x0)*rho_inf &
            +(3.0d0/40.0d0)*((1.0d0-x1)*t1+3.0d0*(1.0d0+x2)*t2)*(3.0d0*pi**2)**(2.0d0/3.0d0)*rho_inf**(5.0d0/3.0d0) &
            +(1.0d0/24.0d0)*t3*(1.0d0-x3)*rho_inf**(sigma+1.0d0)
                
            write(*,*) "particle number", partnum
            write(*,*) "density", rho_inf, ",", "strength", strengthfactor, ",", "periodicity", periodicity
            !write(*,*) ea
                
            evalloop: do evals = 1, maxevals

                rhoformat = "(a7, I5.5, a4, F4.2, a3, F4.2, a7, I2.2, a10, I3.2, a4)"
                write(filename, rhoformat) "partnum", partnum, "rho", rho_inf, "_Ef", strengthfactor, "_period",  periodicity, &
                "_iteration", evals, "_rho"
                open(unit=16, file=trim(adjustl(filename))//'.dat', action="write", status="replace", RECL=200)
                DO a_j = 1, accuracy
                  WRITE(16,*) (a_j-1)*stepsize ,",", rho(a_j)
                END DO
                CLOSE(16)

                appendloops = 0
                call coeffIarraygenerate(coeffIarray)
                call firstder(coeffIarray, coeffIprime, stencil)
                !finds eigenpairs in one eighth of the nx-ny plane and assigns multiplicities to each pair    
                do nx = 0, nmax
                    do ny = 0, nx
                        call finitediffmatrix(findiffmtx, stencil)
                        call eigenvalvec(findiffmtx, eigenpairs)
                        !multiplicity checks
                        if(nx==0)then
                            eigenpairs(:)%multiplicity = 1*2
                        else
                            if((ny==0) .or. (nx==ny))then
                                eigenpairs(:)%multiplicity = 4*2
                            else
                                eigenpairs(:)%multiplicity = 8*2
                            end if   
                        end if
                        eigenpairs(:)%nx = nx
                        eigenpairs(:)%ny = ny
                        if(appendloops < searchspacesize)then
                            call quicksort(eigenpairs, 1, accuracy)
                            appendlooppos = accuracy*appendloops
                            sortedstates((appendlooppos+1):(appendlooppos+accuracy)) = eigenpairs
                            appendloops = appendloops + 1
                        else
                            call eigenvalappend(sortedstates, eigenpairs)
                        end if
                    end do
                end do

               !MattCopy=sortedstates 

                call density(sortedstates, rho, tau)
                call firstder(rho, rhoder, stencil)
                call secondder(rho, rho2der, stencil)
            
                !calculating the total energy
                H_array = 0.0d0
                do dummywrite = 1, accuracy
                    H_array(dummywrite) = H(dummywrite)
                end do
                
                !integerates energy density functional to find total energy
                totalenergy = (length_0**2.0d0)*gaulegint(H_array)
                    
                !checks for minima, then convergence and exits once minima or minimum converged value is found  
                if(evals == 1)then
                    totalenergy_0 = totalenergy
                    minima = totalenergy
                    finaleval = 1
                    minimacheck = 0
                    reldiffcheck = 0
                else if(evals == maxevals)then
                    if(totalenergy < minima)then
                        minima = totalenergy
                        finaleval = maxevals
                    end if
                    write(*,*) "loop number", evals, ",", "energy/particle", minima/real(partnum)
                    write(*,*) "completed" 
                    exit evalloop
                else
                    if(minimacheck <= 1)then
                        !checks for minima
                        if(totalenergy > minima)then
                            totalenergy_0 = totalenergy
                            minimacheck = minimacheck + 1
                            reldiffcheck = 0
                        else if(totalenergy <= minima)then  
                            !checks for congergence
                            reldiff = abs(totalenergy - totalenergy_0)/(0.5d0*(abs(totalenergy_0) + abs(totalenergy)))
                            if(reldiff <= errortol)then
                                if(reldiffcheck <= 1)then
                                    totalenergy_0 = totalenergy
                                    minima = totalenergy
                                    finaleval = evals
                                    minimacheck = 0
                                    reldiffcheck = reldiffcheck + 1
                                else if(reldiffcheck == 2)then
                                    minima = totalenergy
                                    finaleval = maxevals
                                    write(*,*) "loop number", evals, ",", "energy/particle", minima/real(partnum)
                                    write(*,*) "completed" 
                                    exit evalloop
                                end if    
                            else if(reldiff > errortol)then
                                totalenergy_0 = totalenergy
                                minima = totalenergy
                                finaleval = evals
                                minimacheck = 0
                                reldiffcheck = 0
                            end if
                        end if
                    else if(minimacheck == 2)then !first 2 searches for minima completed, exits loop
                        write(*,*)  "loop number", evals, ",", "energy/particle", minima/real(partnum)
                        write(*,*) "completed" 
                        exit evalloop     
                    end if
                end if
                
                write(*,*) "loop number", evals, ",", "energy/particle", totalenergy/real(partnum)

            end do evalloop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Counts the number of distinct states
!            MattCount=0
!            do iMatt=1,sortsize
!              write(*,*) iMatt, MattCopy(iMatt)%multiplicity
!              MattCount=MattCount+MattCopy(iMatt)%multiplicity
!              if (MattCount>=partnum) then
!                MattUnique=iMatt
!                exit
!              end if
!            end do
  
!            write(*,*) "Number of states =", MattUnique 

! Normalizes the z-orbitals
 !           do iMatt = 1, MattUnique
 !              MattCopy(iMatt)%norm = gaulegint(MattCopy(iMatt)%evec**2)
 !              MattCopy(iMatt)%evec = MattCopy(iMatt)%evec/sqrt(MattCopy(iMatt)%norm)
 !              write(*,*) iMatt, "total eigenvalue:", MattCopy(iMatt)%eval,"nx=",MattCopy(iMatt)%nx,"ny=",MattCopy(iMatt)%ny,&
!"e_i-(hbar^2/2m)(kx^2+ky^2)=",MattCopy(iMatt)%eval-(MattCopy(iMatt)%nx**2+MattCopy(iMatt)%ny**2)*dispersion*((2*pi)/length_0)**2
!            end do

! Writes the normalized z-orbitals to file
          !  open(19,file="phi.dat")

           ! do iMatt=1,10
           !   write(19,*) iMatt, MattCopy(iMatt)%nx,MattCopy(iMatt)%ny
           !   do jMatt=1,accuracy
           !     write(19,*) (jMatt-1)*stepsize, MattCopy(iMatt)%evec(jMatt)
           !   end do
           ! end do
           ! close(19)

! Computes the inner products
!            do iMatt = 1, MattUnique
!              Matt1=MattCopy(iMatt)
!              do jMatt = 1, iMatt
!                Matt2=MattCopy(jMatt)

!                if(Matt1%nx==Matt2%nx .and. Matt1%ny==Matt2%ny) then
            !      write(*,*) 1-(iMatt-jMatt)/max((iMatt-jMatt),1), gaulegint(Matt1%evec*Matt2%evec)
            !      write(*,*) 1-(iMatt-jMatt)/max((iMatt-jMatt),1), trapint(Matt1%evec*Matt2%evec)
!                  write(*,*) iMatt,jMatt, trapint(Matt1%evec*Matt2%evec)
!                end if
!              end do
!            end do
           
           ! write(*,*) covarray 

   ! write quantity to a file
!            open(19,file="density.dat")
!            do iMatt=1,accuracy
!              write(19,*) (iMatt-1)*stepsize, rho(iMatt)
!            end do
!            close(19)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !file writing
        

!!        !writing energies
        energyformat = "(a7, I5.5, a4, F4.2, a3, F4.2, a7, I2.2, a7)"
        write(filename, energyformat) "partnum", partnum, "_rho", rho_inf, "_Ef", strengthfactor, "_period",&
  periodicity,  "_energy"
        open(unit=16, file=trim(adjustl(filename))//'.dat', action="write", status="replace", RECL=200)
            write(16,*) minima, ",", minima/real(partnum)
        close(16)

        write(17,*) strengthfactor,",",periodicity,",", minima/real(partnum) 

       !writing densities
        rhoformat = "(a7, I5.5, a4, F4.2, a3, F4.2, a7, I2.2, a4)"
        write(filename, rhoformat) "partnum", partnum, "rho", rho_inf, "_Ef", strengthfactor, "_period",  periodicity,  "_rho"
        open(unit=16, file=trim(adjustl(filename))//'.dat', action="write", status="replace", RECL=200)
            DO a_j = 1, accuracy
                WRITE(16,*) (a_j-1)*stepsize ,",", rho(a_j) 
            END DO
        CLOSE(16)
        
        !writing kinetic densities
        tauformat = "(a7, I5.5, a4, F4.2, a3, F4.2, a7, I2.2, a4)"
        write(filename, tauformat) "partnum", partnum, "rho", rho_inf, "_Ef", strengthfactor, "_period",  periodicity,  "_tau"
        open(unit=16, file=trim(adjustl(filename))//'.dat', action="write", status="replace", RECL=200)
            DO a_j = 1, accuracy
                WRITE(16,*) (a_j-1)*stepsize ,",", tau(a_j) 
            END DO
        CLOSE(16)
            
        end do
    end do
    close(17) 
!####################################################################################################################################################
!####################################################################################################################################################

!subroutines and functions

    contains
 
    subroutine closurecheck(closurestates)
    !this subroutine checks whether a shell closure occurs at a given particle number, and updates the number to keep the shell closed if not
        implicit none
        
        !inputs
        type(states), dimension(accuracy), intent(inout) :: closurestates
        
        !looping variables
        integer closureloopindex
        integer closurelooppartnum
        
        !array storing the shell closure data
        integer, dimension(0:nmax,0:nmax) :: shellclosurearray
        integer :: shellclosurenx, shellclosureny
        
        shellclosurearray = 0
        closurelooppartnum = 0
        
        do closureloopindex = 1, accuracy
            shellclosurearray(closurestates(closureloopindex)%nx, closurestates(closureloopindex)%ny) &
 = shellclosurearray(closurestates(closureloopindex)%nx, closurestates(closureloopindex)%ny) + 1  
            
            checkloop: do shellclosurenx = 0, nmax
                do shellclosureny = 0, nmax
                    if((mod(shellclosurearray(shellclosurenx, shellclosureny),2) == 0) .and. &
(shellclosurearray(shellclosurenx, shellclosureny) /= 0))then
                        closurestates(closureloopindex)%closed = .false.
                        exit checkloop
                    else if((shellclosurenx==nmax) .and. (shellclosureny==nmax))then
                        closurestates(closureloopindex)%closed = .true.
                    end if
                end do
            end do checkloop
        
        end do
        
    end subroutine closurecheck

    subroutine varchange(varchangearray)
    !this subroutine calculates the chane of variable array
        implicit none

        !outputs
        double precision, dimension(accuracy), intent(out) :: varchangearray

        !looping parameters
        integer :: varchangex
        double precision :: varchangeint

        varchangearray = 0.0d0
        do varchangex = 2, accuracy
            varchangearray(varchangex) = varchangearray(varchangex-1) + stepsize*(coeffI(varchangex)+coeffI(varchangex-1))/2.0d0
        end do

        varchangearray = exp(-0.5d0*varchangearray)

    end subroutine varchange


!####################################################################################################################################################


    subroutine eigenvalvec(A, pairs)
    !this subroutine computes the eigenvalues ER and eigenvectors VSR for a given input NxN Matrix A
        implicit none
    
        !inputs

        !matrices to find eigenvalues of (matrix B is the identity)
        double precision, dimension(accuracy,accuracy), intent(inout) :: A
        !vector storing eigenpairs
        type(states), dimension(accuracy), intent(out) :: pairs
        
        !process variables
        
        !array storing the right eigenvalues
        double precision, dimension(accuracy) :: ER
        !arrays storing right eigenvectors 
        double precision, dimension(accuracy, accuracy) :: VSR
        !column index of vsr
        integer :: vsr_col

        !order of matrix A (and B)
        integer, parameter :: N = accuracy
        !dimensions of matrix A (and B)
        integer, parameter :: lda = accuracy, ldb = accuracy
        !dimensions of the matrices VL and VR (LDVL is not used)
        integer, parameter :: ldvsl = accuracy, ldvsr = accuracy
        !workspace size
        integer, parameter :: lwork = 8*N + 16
        !arrays used in the eigenvalue finding process (alphai should be = 0 since energies should be positive)
        double precision, dimension(accuracy) :: alphar, alphai, beta
        !arrays storing left eigenvectors (not used)
        double precision, dimension(accuracy, accuracy) :: vsl
        !identity matrix
        double precision, dimension(accuracy, accuracy) :: B
        !workspace 
        double precision, dimension(lwork) :: work
        !eigenvalue sorting
        integer :: sdim
        double precision, dimension(N) :: bwork
        !process information
        integer :: info

        !counter to generate identity matrix
        integer :: B_i, B_j

        !creating identity matrix
        do B_i = 1, ldb
    	    do B_j = 1, ldb
                if(B_i == B_j)then
                    B(B_i, B_j) = 1.0d0
                else
                    B(B_i, B_j) = 0.0d0
                end if
            end do
        end do

        !LAPACK subroutine finding right eigenvalues and eigenvectors
        call dgges('N', 'V', 'N', SELCTG, N, A, LDA, B, LDB, SDIM, ALPHAR, ALPHAI, BETA, VSL,&
 LDVSL, VSR, LDVSR, WORK, LWORK, BWORK, INFO)

        ER = ALPHAR/BETA
        
        !bundles the eigenvalues and eigenvectors into a paired data type
        pairs(:)%eval = ER
        do vsr_col = 1, accuracy
            pairs(vsr_col)%evec(:) = VSR(:,vsr_col)
        end do
        
    end subroutine eigenvalvec
    
!####################################################################################################################################################

    logical function SELCTG(alphar, alphai, beta)
    !dummy sorting function for the matrix eigenvalue solver
        double precision, intent(in) :: alphar, alphai, beta
        selctg = .true.
        if(abs(alphai) >= 10**8)then
            selctg = .false.
        end if
    end function SELCTG

!####################################################################################################################################################

    subroutine finitediffmatrix(A, stenciltype)
    !generates the finite difference matrix using 3 point stencil
        implicit none
        
        !difference matrix
        double precision, dimension(accuracy, accuracy), intent(out) :: A
        
        !defines the stencil type to be used (3 or 5 point)
        integer, intent(in) :: stenciltype
        
        !rows and column indices of A
        integer :: a_row, a_col
        
        !clearing the matrix
        A = 0.0d0
        
        if(stenciltype == 3)then 
            do a_row = 1, accuracy
                do a_col = 1, accuracy
                    if(a_row == a_col)then
                        !diagonal elements
                        A(a_row, a_col) = coeffC(a_row) - 2.0d0*coeffA(a_row)/(stepsize**2)
                    else if((a_col == a_row + 1) .or. ((a_row == accuracy) .and. (a_col == 1)))then
                        !upper off diagonal elements and lower left hand corner element (offset 1 right from corner for periodicity)
                        A(a_row, a_col) = coeffA(a_row)/(stepsize**2) + coeffB(a_row)/(2.0d0*stepsize)
                    else if((a_col == a_row - 1) .or. ((a_row == 1) .and. (a_col == accuracy)))then
                        !lower off diagonal elements and upper right hand corner element (offset 1 left from corner for periodicity)
                        A(a_row, a_col) = coeffA(a_row)/(stepsize**2) - coeffB(a_row)/(2.0d0*stepsize)
                    else
                        !all other elements are zero
                        A(a_row, a_col) = 0.0d0
                    end if
                end do
            end do
        else if(stenciltype == 5)then
            do a_row = 1, accuracy
                do a_col = 1, accuracy
                    if(a_row == a_col)then
                        !diagonal elements
                        A(a_row, a_col) = coeffCold(a_row) - 5.0d0*coeffA(a_row)/(2.0d0*stepsize**2)
                    else if((a_col == a_row + 1) .or. ((a_row == accuracy) .and. (a_col == 1)))then
                        !first upper off diagonal elements and lower left hand corner element (offest 1 right from corner for periodicity)
                        A(a_row, a_col) = 4.0d0*coeffA(a_row)/(3.0d0*stepsize**2) + 2.0d0*coeffBold(a_row)/(3.0d0*stepsize)
                    else if((a_col == a_row + 2) .or. ((a_row == accuracy - 1) .and. (a_col == 1)) &
.or. ((a_row == accuracy) .and. (a_col == 2)))then
                        !second upper off diagonal elements and lower left hand corner elements 
                        A(a_row, a_col) = - coeffA(a_row)/(12.0d0*stepsize**2) - coeffBold(a_row)/(12.0d0*stepsize)
                    else if((a_col == a_row - 1) .or. ((a_row == 1) .and. (a_col == accuracy)))then
                        !first lower off diagonal elements and upper right hand corner (offset 1 left from corner for periodicity)
                        A(a_row, a_col) = 4.0d0*coeffA(a_row)/(3.0d0*stepsize**2) - 2.0d0*coeffBold(a_row)/(3.0d0*stepsize)
                    else if((a_col == a_row - 2) .or. ((a_row == 1) .and. (a_col  == accuracy - 1)) &
.or. ((a_row == 2) .and. (a_col == accuracy)))then
                        !second lower off diagonal elements and upper right hand corner elements 
                        A(a_row, a_col) = - coeffA(a_row)/(12.0d0*stepsize**2) + coeffBold(a_row)/(12.0d0*stepsize)
                    else 
                        A(a_row, a_col) = 0.0d0
                    end if
                end do
            end do
        end if    
    end subroutine finitediffmatrix

!####################################################################################################################################################
    
    !gauss-legendre integration subroutines and functions
    
    SUBROUTINE gauleg(x_initial, x_final, x_array,w_array,n)
    !calculates the gauss-legendre weights for a given accuracy and range
        INTEGER n
        DOUBLE PRECISION x_initial,x_final,x_array(n),w_array(n)
        DOUBLE PRECISION EPS
        PARAMETER (EPS=1.d-14)
        INTEGER i,j,m
        DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
        m=(n+1)/2
        xm=0.5d0*(x_final+x_initial)
        xl=0.5d0*(x_final-x_initial)
        do i=1,m
            z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1           continue
                p1=1.d0
                p2=0.d0
                do j=1,n
                    p3=p2
                    p2=p1
                    p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11              end do
                pp=n*(z*p1-p2)/(z*z-1.d0)
                z1=z
                z=z1-p1/pp
            if(abs(z-z1).gt.EPS)goto 1
            x_array(i)=xm-xl*z
            x_array(n+1-i)=xm+xl*z
            w_array(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
            w_array(n+1-i)=w_array(i)
12      end do
        return
    END SUBROUTINE gauleg
    
    double precision function gaulegint(A)
        !calculates the integral of a function using the gauss-legendre quadrature
        double precision, dimension(accuracy) :: A
        !index of the array
        integer :: a_index
        !range of the integration
        double precision :: x_init, x_fin
        !arrays containing the abscissa and the gauss-legendre weights
        double precision, dimension(accuracy) :: x_array, w_array
        
        !linear interpolation points
        integer :: interlow, interhigh
        double precision :: interweight
        
        x_init = 0.0d0
        x_fin = length_0
        
        call gauleg(x_init, x_fin, x_array, w_array, accuracy)
        
        gaulegint = 0.0d0
        do a_index = 1, accuracy
            !linear interpolation between points is used for better accuracy
            interlow = floor(x_array(a_index)/stepsize)+1
            interhigh = interlow + 1
            interweight = (x_array(a_index)/stepsize) - floor(x_array(a_index)/stepsize)
            
            if(interlow == accuracy)then
                gaulegint = gaulegint +  w_array(a_index)*(A(accuracy)*(1.0d0 - interweight) + A(1)*interweight)
            else
                gaulegint = gaulegint +  w_array(a_index)*(A(interlow)*(1.0d0 - interweight) + A(interhigh)*interweight)
            end if
        end do
    end function gaulegint
    
    
    
!####################################################################################################################################################
    double precision function trapint(A) 
        !calculates the integral of a function using trapezoid rule
        double precision, dimension(accuracy) :: A
        !index of the array
        integer :: a_index
        
        trapint = 0.0d0
        do a_index = 1, accuracy - 1
            trapint = trapint + stepsize*(A(a_index)+A(a_index+1))/2.0d0
        end do  
        trapint = trapint + stepsize*(A(accuracy) + A(1))/2.0d0
    end function trapint

!####################################################################################################################################################

    !quicksort subroutines/functions
    
    recursive subroutine quicksort(arraysort, low, high)
        implicit none
        type(states), dimension(:), intent(inout) :: arraysort
        integer, intent(in) :: low, high
        integer :: pivot_location
        
        if(low < high)then
            call partition(arraysort, low, high, pivot_location)
            call quicksort(arraysort, low, pivot_location)
            call quicksort(arraysort, pivot_location + 1, high)
        end if 
    end subroutine quicksort
 
    subroutine partition(arraysort, low, high, pivot_location)
        implicit none
        type(states), dimension(:), intent(inout) :: arraysort
        type(states) :: swap
        integer, intent(in) :: low, high
        integer, intent(out) :: pivot_location
        double precision :: pivot
        
        integer :: counter
        
        pivot = arraysort(low)%eval
        pivot_location = low
        
        do counter = low + 1, high
            if(arraysort(counter)%eval < pivot)then
                swap = arraysort(counter)
                arraysort(counter) = arraysort(pivot_location + 1)
                arraysort(pivot_location + 1) = swap
                pivot_location = pivot_location + 1
            end if
        end do
        
        swap = arraysort(low)
        arraysort(low) = arraysort(pivot_location)
        arraysort(pivot_location) = swap   
    end subroutine partition

!####################################################################################################################################################
    
    subroutine eigenvalappend(mainarray, appendarray)
        !appends and sorts in new eigenstates into states used for calculations (discards unused states)
        implicit none
        !main array holding used states
        type(states), dimension(sortsize), intent(inout) :: mainarray
        !states that have not yet been appended
        type(states), dimension(accuracy), intent(in) :: appendarray
        !temporary array used for sorting
        type(states), dimension(sortsize+accuracy) :: sortingarray
        
        !adds already known eigenstates to sorting array
        sortingarray(1:sortsize) = mainarray
        !adds new eigenstates to sorting array, and sorts them into the known states
        sortingarray((1+sortsize):(sortsize+accuracy)) = appendarray
        call quicksort(sortingarray, 1, (sortsize+accuracy))
        !takes the first (accuracy) eigenstates from the new sorted states and modifies the main array to reflect this
        mainarray = sortingarray(1:sortsize)
    end subroutine eigenvalappend
    
!####################################################################################################################################################
    
    subroutine density(eigenstates, rho_out, tau_out)
    !calculates particle number density and kinetic energy density
        implicit none
        !input list of sorted and selected eigenstates
        type(states), dimension(sortsize), intent(inout) :: eigenstates
        !output density
        double precision, dimension(accuracy), intent(out) :: rho_out
        !output kinetic density
        double precision, dimension(accuracy), intent(out) :: tau_out
        
        !temporary wavefunction storing the derivative of the wavefunctions
        double precision, dimension(accuracy) :: wavefunder_temp
        !counter to run through all the states
        integer :: eigenstates_i
        
        !counters for the multiplicites of each state and the total number of particles run through
        integer :: multiplicitycount
        integer :: particlecount
        
        !zero the output array
        rho_out = 0.0d0
        tau_out = 0.0d0
        
        particlecount = 0
        outerloop1: do eigenstates_i = 1, sortsize
            do multiplicitycount = 1, eigenstates(eigenstates_i)%multiplicity
                eigenstates(eigenstates_i)%norm = gaulegint(eigenstates(eigenstates_i)%evec**2)
                rho_out = rho_out + (eigenstates(eigenstates_i)%evec**2)/eigenstates(eigenstates_i)%norm
                particlecount = particlecount + 1
                if(particlecount == partnum)then
                    exit outerloop1
                end if
            end do
        end do outerloop1
        if(particlecount < partnum)then
            write(*,*) "the search space size is too small, searchspacesize must be increased"
            read(*,*)
        end if
        !normalization in x/y gives factor of 1/(length_0**2)
        rho_out = rho_out/(length_0**2)
        
        particlecount = 0
        outerloop2: do eigenstates_i = 1, sortsize
            do multiplicitycount = 1, eigenstates(eigenstates_i)%multiplicity
                wavefunder_temp = 0.0d0
                call firstder(eigenstates(eigenstates_i)%evec, wavefunder_temp, stencil)
                tau_out = tau_out + (((2.0d0*pi/length_0)**2)*(eigenstates(eigenstates_i)%nx**2 +&
 eigenstates(eigenstates_i)%ny**2)*(eigenstates(eigenstates_i)%evec**2) + (wavefunder_temp**2))/eigenstates(eigenstates_i)%norm
                particlecount = particlecount + 1
                if(particlecount == partnum)then
                    exit outerloop2
                end if
            end do
        end do outerloop2
        !normalization in x/y gives factor of 1/(length_0**2)
        tau_out = tau_out/(length_0**2)
        
        
    end subroutine density
    

!####################################################################################################################################################
    
    !derivative subroutines
    
    subroutine firstder(A,Ader,stenciltype)
        !computes first derivative using 3 or 5 point stencil        
        implicit none
        !incoming array (change bounds to allow for using modulo at boundaries)
        double precision, dimension(0:accuracy-1), intent(in) :: A
        !derivative matrix
        double precision, dimension(accuracy), intent(out) :: Ader
        !stencil type (3 or 5 point)
        integer, intent(in) :: stenciltype
        integer :: a_index
        
        if(stenciltype == 3)then 
            do a_index = 0, accuracy - 1
                Ader(a_index + 1) = (A(modulo(a_index + 1, accuracy)) - A(modulo(a_index - 1, accuracy)))/(2.0d0*stepsize)
            end do
        else if(stenciltype == 5)then
            do a_index = 0, accuracy - 1           
                Ader(a_index+1) = (-A(modulo(a_index + 2, accuracy)) + 8.0d0*A(modulo(a_index + 1, accuracy)) &
- 8.0d0*A(modulo(a_index - 1, accuracy)) + A(modulo(a_index - 2, accuracy)))/(12.0d0*stepsize)
            end do
        end if  
    end subroutine firstder
    
    subroutine secondder(A,A2der,stenciltype)
        !computes second derivative using 3 or 5 point stencil        
        implicit none
        !incoming array
        double precision, dimension(0:accuracy-1), intent(in) :: A
        !derivative matrix
        double precision, dimension(accuracy), intent(out) :: A2der
        !stencil type (3 or 5 point)
        integer, intent(in) :: stenciltype
        integer :: a_index
        
        if(stenciltype == 3)then
            do a_index = 0, accuracy - 1
                A2der(a_index + 1) = (A(modulo(a_index + 1, accuracy)) - 2.0d0*A(a_index) + &
A(modulo(a_index - 1, accuracy)))/(stepsize**2)
            end do
        else if(stenciltype == 5)then
            do a_index = 0, accuracy - 1
                A2der(a_index+1) = (-A(modulo(a_index + 2, accuracy)) + 16.0d0*A(modulo(a_index + 1, accuracy)) &
- 30.0d0*A(a_index) + 16.0d0*A(modulo(a_index - 1, accuracy)) - A(modulo(a_index - 2, accuracy)))/(12.0d0*(stepsize**2))
            end do
        end if
    end subroutine secondder

!####################################################################################################################################################
   
    !these are the coefficients used for the transformed problem with the first derivative term removed

    !Writing a 2nd order ODE eigenvalue problem as A(x)y''(x) + B(x)y'(x) + C(x)y(x) = Ey(x) these functions represent each of the A,B,C coefficients in this equation
    double precision function coeffA(N)
        integer, intent(in) :: N
        coeffA = coeffAold(N)
    end function coeffA
    
    double precision function coeffB(N)
        integer, intent(in) :: N
        coeffB = 0.d0
    end function coeffB
    
    double precision function coeffC(N)
        integer, intent(in) :: N
        coeffC = coeffCold(N) - coeffAold(N)*(0.5d0*coeffIprime(N)+0.25d0*((coeffI(N))**2))
        !converting from nx -> kx gives factor of (2.0*pi*/length_0)**2
    end function coeffC   

    !integrand used for the change of variables
    double precision function coeffI(N)
        integer, intent(in) :: N
        coeffI = coeffBold(N)/coeffAold(N)
        !converting from nx -> kx gives factor of (2.0*pi*/length_0)**2
    end function coeffI 
    
    subroutine coeffIarraygenerate(Iarray)
    !subroutine calculates the coefficient I array
        implicit none
        
        double precision, dimension(accuracy), intent(inout) :: Iarray
        integer :: Iarrayindex
        
        do Iarrayindex = 1, accuracy
            Iarray(Iarrayindex) = coeffBold(Iarrayindex)/coeffAold(Iarrayindex)
        end do
        
    end subroutine coeffiarraygenerate

    !these are the untransformed coefficents for the neutron matter problem with the first derivative term intact

    !Writing a 2nd order ODE eigenvalue problem as A(x)y''(x) + B(x)y'(x) + C(x)y(x) = Ey(x) these functions represent each of the A,B,C coefficients in this equation
    double precision function coeffAold(N)
        integer, intent(in) :: N
        coeffAold = -RM(N)
    end function coeffAold
    
    double precision function coeffBold(N)
        integer, intent(in) :: N
        coeffBold = -RMder(N)
    end function coeffBold
    
    double precision function coeffCold(N)
        integer, intent(in) :: N
        coeffCold = U(N) + V(N) + RM(N)*(nx**2+ny**2)*((2.0d0*pi/length_0)**2)
        !converting from nx -> kx gives factor of (2.0*pi*/length_0)**2
    end function coeffCold


    
!####################################################################################################################################################
    
    !functions for skyrme interactions    

    !reduced mass term
    double precision function RM(N)
        integer, intent(in) :: N  
        RM = dispersion + (1.0d0/8.0d0)*((1.0d0-x1)*t1+3.0d0*(1.0d0+x2)*t2)*rho(N)
    end function RM
    
    !reduced mass derivative term
    double precision function RMder(N)
        integer, intent(in) :: N  
        RMder = (1.0d0/8.0d0)*((1.0d0-x1)*t1+3.0d0*(1.0d0+x2)*t2)*rhoder(N)
    end function RMder
    
    !skyrme interaction energy
    double precision function U(N)
        integer, intent(in) :: N
        U = (1.0d0/2.0d0)*t0*(1.0d0-x0)*rho(N) &
        + (1.0d0/24.0d0)*t3*(1.0d0-x3)*(2.0d0+sigma)*(rho(N)**(sigma+1.0d0)) &
        + (1.0d0/8.0d0)*((1.0d0-x1)*t1+3.0d0*(1.0d0+x2)*t2)*tau(N) &
        + (3.0d0/16.0d0)*((x1-1.0d0)*t1+(1.0d0+x2)*t2)*rho2der(N)
    end function U
    
    !external potential
    double precision function V(N)
        integer, intent(in) :: N
        !start at N-1 since arrays are indexed to begin at 1 -> V(N=1) = strength*cos(0) = strength
        V = strength*cos(2.0d0*Pi*periodicity*(N-1)*stepsize/length_0)
    end function V
    
    !hamiltonian density functional
    double precision function H(N)
        integer, intent(in) :: N
        H = dispersion*tau(N) &
            + (1.0d0/4.0d0)*t0*(1.0d0-x0)*rho(N)**2 &
            + (1.0d0/24.0d0)*t3*(1.0d0-x3)*rho(N)**(2.0d0+sigma) &
            + (1.0d0/8.0d0)*(t1*(1.0d0-x1)+3.0d0*t2*(1.0d0+x2))*tau(N)*rho(N) &
            + (3.0d0/32.0d0)*(t1*(1.0d0-x1)-t2*(1.0d0+x2))*(rhoder(N)**2) &
            + V(N)*rho(N)
    end function H

!####################################################################################################################################################

 end program HartreeFock

