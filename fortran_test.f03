Program Fortran_Test 

	Implicit None
	
	type :: vector_trans
		real, dimension(3) :: coord
	end type vector_trans

	integer :: maxstep=10000, NSim=3, NStates=3
	real:: emax=0.0025, w=0.147, delta=0.05, phase=0, mag_direc=0, sum_direc=0, pi=3.14159265359
	
!anything before this is already defined before

	real, dimension(4) :: Energies=0
	real, dimension(2) :: CoordsX, CoordsY, CoordsZ

	integer :: i=1, j=1, k=1
	real, dimension(:,:),Allocatable :: h0, d0
	real, dimension(:,:),Allocatable :: h
	complex, dimension(:),Allocatable :: ctmp,ctmp1,ctmp2,comp_imag

	integer :: t = 1 
	real, dimension(3):: direction=0, norm_direc=0
	
	complex :: imag=(0.0,1.0)
	real, dimension(:), Allocatable :: dip, norm, efield
	real, dimension (:,:), Allocatable :: v
	
	type(vector_trans),dimension(:,:), Allocatable :: TransDipoles
	allocate(TransDipoles(NStates+2,NStates+2))

	phase=0*pi/2
	
!	write(*,*) 'Trans Dipoles empty', TransDipoles

	TransDipoles(1,2)%coord=(/0.0,0.0,1.1568909999999999/)
	TransDipoles(2,1)%coord=(/0.0,0.0,1.1568909999999999/)
	TransDipoles(2,3)%coord=(/0.0,0.0,1.448285/)
	TransDipoles(3,2)%coord=(/0.0,0.0,1.448285/)
	Energies=(/0.0,0.97666132,1.63684807,0.0/)

!	write(*,*) 'Trans Dipoles', TransDipoles
!	write(*,*) 'Energies', Energies
	
	CoordsX=(/0.0,0.0/)
	CoordsY=(/0.0,0.0/)
	CoordsZ=(/-0.36656,0.36656/)
	

	direction=(/CoordsX(2)-CoordsX(1),CoordsY(2)-CoordsY(1),CoordsZ(2)-CoordsZ(1)/)
	do i=1,3	
	sum_direc=sum_direc+(direction(i)*direction(i))
	end do
	mag_direc=sqrt(sum_direc)
	norm_direc=direction/mag_direc
	 
!	write(*,*) 'Direction', direction

	allocate(dip(maxstep),norm(maxstep),v(maxstep,NSim),efield(maxstep))
	
	dip=0.0
	norm=0.0
	v=0.0
	efield=0.0
	
	norm(1)=1.0
	norm(2)=1.0
	v(1,1)=1.0
	v(2,1)=1.0

	do t=1,maxstep
	   efield(t+1)=emax * sin(w*delta*t+phase)
	end do
	
!	write(*,*) 'efield', efield

	allocate(h0(NSim,NSim), d0(NSim, NSim))
	h0=0.0
	d0=0.0

!Erase when done
!	d0(1,:)=(/0.0,1.156891,0.0/)
!	d0(2,:)=(/1.156891,0.0,1.448285/)
!	d0(3,:)=(/0.0,1.448285,0.0/)
!	data d0 / 0.0,1.156891,0.0,1.156891,0.0,1.448285,0.0,1.448285,0.0 /
!	write(*,*) d0
!Erase when done

	do k=1,NSim
	h0(k,k) = energies(k)
	end do
!	write(*,*) h0

	do i=1,NSim
	   do j=1,NSim
		  d0(i,j)=dot_product(norm_direc,TransDipoles(i,j)%coord)
	   end do
	end do
!	write(*,*) 'd0', d0

	allocate(ctmp(NSim),ctmp1(NSim),ctmp2(NSim),comp_imag(NSim),h(NSim,NSim)) !Initialize as zero

	ctmp=(0.0,0.0)
	ctmp1=(0.0,0.0)
	ctmp2=(0.0,0.0)
	comp_imag=(0.0,0.0)
	h=0.0

1349 FORMAT (F15.10,F15.10,'i')
535 FORMAT (ES14.5,ES14.5,'i')
!	write(*,1349) ctmp
!	write(*,1349) ctmp2

	
	h = h0+efield(2)*d0
10 FORMAT ('h',F10.8)
!	write(*,10) h	
	
	ctmp1(1)=(1.0,0.0)
	comp_imag=MatMul(imag*delta*h,ctmp1)
!	write(*,1349) comp_imag
	ctmp = ctmp1+comp_imag
!	write(*,535) ctmp
	ctmp2=ctmp/sqrt(abs(dot_product(conjg(ctmp),ctmp)))
!	write(*,*) 'ctmp2'
!	write(*,535) ctmp2
	
	do i =3,maxstep
		h=h0 + efield(i) * d0
		comp_imag = MatMul(2*imag*delta*h,ctmp2)
		ctmp=ctmp1+comp_imag
!	write(*,*) 'ctmp1'
!	write(*,535) ctmp1
!	write(*,*) 'comp_imag'
!	write(*,535) comp_imag
!	write(*,*) 'ctmp'
!	write(*,535) ctmp
		ctmp1=ctmp2
		norm(i)=sqrt(abs(dot_product(conjg(ctmp),ctmp)))
		ctmp2=ctmp/norm(i)
!	write(*,*) 'ctmp2'
!	write(*,535) ctmp2
		v(i,:)=ctmp2
		dip(i)=real(dot_product(matmul(d0,conjg(ctmp2)),ctmp2))
!	write(*,*) 'dip', dip(i)
	end do

End Program Fortran_Test
