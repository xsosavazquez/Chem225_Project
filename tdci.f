      Program tdci

      Real::NAtoms,NStates,NMOs,NSim,EMax,Maxstep,Delta,w,chirp
      Real::ncyc,phase,pi=3.14159,
      Character(len=10)::envelope

      NAtom=2           !Number of atoms in the system
      NStates=3         !Number of excited states in the system
      NMOs=2            !Number of molecular orbitals in the CAS active space
      NSim=3            !Number of states in the simulation
      EMax=0.025        !Field Strength in AU
      Maxstep=1000      !Number of iteration steps
      Delta=0.05        !Stepsize in AU
      w=0.147           !Frequency of the field in AU
      envelope="cosine" !An envelope that could be created (trapezoidal, gaussian or cosine)
      chirp=0.0         !chirp in AU
      ncyc=7            !Number of cycles in the envelope
      phase=0*pi/2

!PUT INPUT DATA HERE!!!!!!!!!!!!!!!!!!!!!!!

!checks to see if the file exists, if it does not exist, then exit the
!program
      open(unit=6,File=tdm_tdcis.txt,IOStat=IError)
      if(IError.ne.0) then
        write(*,'(A)')' Error Opening tdm_tdcis.txt!'
        fail = .True.
      endif      


!Begin TDCI Code
      period=2*pi/(w*delta)
      real,dimension(:),allocate::y=0
      allocate(y(maxstep))
      do 20 t = 0, maxstep
        if(envelope.eq.'cosine') then
          if(t.lt.ncyc*period) then
            y(t)=0.5-cos(2*pi*t/(ncyc*period))/2
          else
            y(t)=0
          endif
        endif
   20 continue
   
      t=0 !Define time variable
	real, dimension(:), allocatable :: direction,norm_direc
	direction= (/CoordsX(1)-CoordsX(0),CoordsY(1)-CoordsY(0),CoordsZ(1)-CoordsZ)(0)/)
	norm_direc = norm(direction) !may have to create "real function norm(x)" function if I don't find the intrinsic function of this
	
	Write(*,*) direction
	
	real, dimension(:) :: dip = 0, norm = 0, maxstep = 0
	real, dimension(:,:) :: v = 0
	
	allocate (dip(maxstep), norm(maxstep), v(maxstep,NSim), efield(maxstep))
	
	norm(0)=1
	norm(1)=1
	v(0,0) =1
	v(1,0) =1
	
	do 10 t=1, maxstep
	
	 efield(t) = emax * (sin(w*delta*t + phase))	
	 
10 continue	
	
	integer :: i = 0,j = 0 , k=0
!	real, dimension(:) :: 
	real, dimension(:,:) :: h0 = 0 , d0 = 0
	
	allocate (h0(NSim,NSim), d0(NSim,NSim))
	
	do 20 k=0 , NSim
	
	h0(i,i) = Energies(i)
	i = i + 1 
	
20 Continue	

	do 30 i=0, NSim
		do 40 j=0, NSim
			d0(i,j)= dot_product(norm_direc,TransDipoles(i,j))
	40 Continue
30 Continue
   
      End Program tdci
