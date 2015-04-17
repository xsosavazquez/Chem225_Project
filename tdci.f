      Program tdci

      Real::NAtoms,NStates,NMOs,NSim,EMax,Maxstep,Delta,w,chirp
      Real::ncyc,phase,pi=3.14159,
      Character(len=10)::envelope

      NAtom=2           !Number of atoms in the system
      NStates=3         !Number of excited states in the system
      NMOs=2            !Number of molecular orbitals in the CAS active space
      NSim=2            !Number of states in the simulation
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
      float,dimension(:),allocate::y=0
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
   
   testing
   
      End Program tdci
