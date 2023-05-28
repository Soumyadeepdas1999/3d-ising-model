program ising
  implicit none
  
  integer :: i,j,L,p,a,b,c,d,niter,time,mm,nn,N
  real :: r,q,E,M,mag,Ei,Ef,dE,h,u
  
  real :: T=2.0,J_ising=1.0 ! assigning value to relevent parameters: k_B T=1
  
  integer, dimension(:,:),allocatable :: spin
  integer :: seed 
  character(len=30) :: charac_a,charac_b ! charac_b stores the name dump_pos
  
  seed=44859
  charac_b = 'store_config'
  
  print *,'Enter the number of lattice points in one dimension'
  read *,L
  print *,'Enter the number of iterations'
  read *,niter
  
  allocate(spin(L,L))
  E=0.0  ! Instantaneous Energy of the lattice
  M=0.0  ! Instantaneous magnetization of the lattice.
  N=L*L  ! Total number of spins in lattice.
  
  call random_seed 
 
!Iinitialize your lattice
  open(71,file='initial_isisng.dat')
  p=0
  do i=1,L
     do j=1,L
        call RANDOM_NUMBER(r)
        spin(j,i) = -1
     
!       if(r<0.5)then
!          spin(j,i)=-1
!       else
!          spin(j,i)=+1
!       end if
        ! WRITING DOWN INITIAL CONFIGURATION.
        write(71,fmt='(4g10.8)') float(i),float(j),float(p),float(spin(j,i))
      end do
   end do
close(71) 

!Calculate initial magnetization and energy
do i=1,L
   do j=1,L 
      a=i+1;b=i-1;c=j+1;d=j-1  ! identifying the 4 neighbours of spin(j,i)
      
      if(i==L)a=1  ! PBC: periodic boundary conditions.
      if(i==1)b=L
      if(j==1)d=L
      if(j==L)c=1
     
      M=M+spin(i,j)
      E=E-J_ising*float((spin(j,i))*(spin(a,j)+spin(b,j)+spin(i,c)+spin(i,d)))
      
   end do
end do

mag=M/(float(N)) ! N=L*L: magnetization (instantaneous) per spin.
E=E*0.5d0
print *,'initial energy E, E per spin =',E, E/float(N) 
print *,'initial magnetization M, M per spin =',M, mag

!INITIALIZATION COMPLETE.
!______________________________________________
!Evolve it to reach equilibrium

open(10,file='ising_T2_L20_init_random.dat')
do time=1,niter ! loop over no of MCS

   do mm=1,L
      do nn=1,L
      
         call RANDOM_NUMBER(r);   i=int(r*float(L))+1 ! CHOOSING a LATTICE SITE.
         call RANDOM_NUMBER(r);   j=int(r*float(L))+1
         
         a=i+1;b=i-1;c=j+1;d=j-1 ! identifying neighbours of spin(i,j)
         
         if(i==L)a=1; if(i==1)b=L;     if(j==1)d=L; if(j==L)c=1 ! PBC
         Ei=-J_ising*float((spin(i,j))*(spin(a,j)+spin(b,j)+spin(i,c)+spin(i,d)))! BEFORE TRIAL FLIP
         
         spin(i,j)=-spin(i,j) ! TRIAL FLIP
         
         Ef=-J_ising*float((spin(i,j))*(spin(a,j)+spin(b,j)+spin(i,c)+spin(i,d)))! AFTER TRIAL FLIP
         
         dE=Ef-Ei ! DIFFERENCE in ENERGIES.
         
         iF(dE<=0.0)then
            E=E+dE    ! UPDATING ENERGY AND MAGNETIZATION OF LATTICE.
            M=M+(2.0*float(spin(i,j)))
         else
            u=exp(-dE/(T))
            call RANDOM_NUMBER(h)
            if(h<u)then
               E=E+dE
               M=M+(2.0*float(spin(i,j))) ! INSTANEOUS MAG. of ENTIRE LATTICE.
            else
               spin(i,j)=-spin(i,j) ! TRIAL FLIP NOT ACCEPTED; E and M NOT UPDATED.
            end if
         end if
    
      end do    
   end do 




!Calculate initial magnetization and energy
do i=1,L
   do j=1,L
   
      a=i+1;b=i-1;c=j+1;d=j-1 ! identifying the 4 neighbours of spin(j,i)
      
      if(i==L)a=1 ! PBC: periodic boundary conditions.
      if(i==1)b=L
      if(j==1)d=L
      if(j==L)c=1
      
      M=M+spin(i,j)
      E=E-J_ising*dfloat((spin(j,i))*(spin(a,j)+spin(b,j)+spin(i,c)+spin(i,d)))
      
   end do
end do

mag=M/(dfloat(N)) ! N=L*L: magnetization (instantaneous) per spin.
E=E*0.5d0
print *,'initial energy E, E per spin =',E, E/dfloat(N)
print *,'initial magnetization M, M per spin =',M, mag

! INITIALIZATION COMPLETE.
 
 open(10,file='ising_L40_heating_T30_01_dT1.dat')
 
 do T_temp=99,300,3 ! TEMPERATURE LOOP
 
      
       T =dfloat(T_temp)/100.0d0 ! FIX T
       
       av_m=0.0d0; av_E=0.0d0
       av_m_N = 0.0d0	;	av_e_N = 0.0d0 ! AV. E, M of ENTIRE LATTICE
       av_m2  = 0.0d0	;	av_e2  = 0.0d0 ! <E2>, <M2> of ENTIRE LATTICE  
       
       
       
do time=1,niter ! loop over no of MCS

   do mm=1,L ! 1 MCS
      do nn=1,L
      
         call RANDOM_NUMBER(r);    i=int(r*float(L))+1 ! CHOOSING A LATTICE SITE.
         call RANDOM_NUMBER(r);    j=int(r*float(L))+1
         
         a=i+1;b=i-1;c=j+1;d=j-1 ! identifying neighbours of spin(i,j)         
         if(i==L)a=1; if(i==1)b=L;     if(j==1)d=L; if(j==L)c=1 ! PBC
         
         Ei=-J_ising*float((spin(i,j))*(spin(a,j)+spin(b,j)+spin(i,c)+spin(i,d)))! BEFORE TRIAL FLIP
         spin(i,j)=-spin(i,j) ! TRIAL FLIP         
         Ef=-J_ising*float((spin(i,j))*(spin(a,j)+spin(b,j)+spin(i,c)+spin(i,d)))! AFTER TRIAL FLIP
         
         dE=Ef-Ei ! DIFFERENCE in ENERGIES.
         
! METROPOLIS ALGORITHM.
 	 iF(dE<=0.0)then
            E=E+dE    ! UPDATING ENERGY AND MAGNETIZATION OF LATTICE.
            M=M+(2.0*float(spin(i,j)))
         else
            u=exp(-dE/(T))
            call RANDOM_NUMBER(h)
            if(h<u)then
               E=E+dE  
               M=M+(2.0*float(spin(i,j))) ! INSTANEOUS MAG. of ENTIRE LATTICE.
            else
               spin(i,j)=-spin(i,j) ! TRIAL FLIP NOT ACCEPTED; E and M NOT UPDATED.
            end if
         end if
    
      end do    
   end do 
      
!__________________________________________________________
!AFTER reaching equilibrium, COLLECT STATISTICAL DATA

   if(time.gt.n_equil) then
!     if(mod(time,n_stat).eq.0) then
   
      mag=abs(M)/(dfloat(N)) ! N=L*L: magnetization (instantaneous) per spin.
      av_m = av_m + mag ; av_e = av_e +E/dfloat(N) ! PER SPIN.
      
      av_m_N = av_m_N + abs(M)	  ;	av_e_N = av_e_N + E ! AV. E, M of ENTIRE LATTICE
      av_m2 = av_m2 + (M*M)  ;	   av_e2 = av_e2 + (E*E) ! <E2>, <M2> of ENTIRE LATTICE
      
!     endif            
   endif
   
enddo ! do time=1,niter ! loop over no of MCS

   av_m = av_m/dfloat(niter - n_equil); av_e = av_e/dfloat(niter - n_equil)
   
   av_e2= av_e2/dfloat(niter - n_equil); av_e_N=av_e_N/dfloat(niter - n_equil)
   cv = (av_e2 - av_e_N*av_e_N)/(T*T)
   av_m2= av_m2/dfloat(niter - n_equil); av_m_N=av_m_N/dfloat(niter - n_equil)
   chi = (av_m2 - av_m_N*av_m_N)/(T)

   write(10,*)T,av_M,av_e,cv,chi ! writing down E and M with no. of iterations.
   
enddo ! do T_temp=300,150

close(10)

deallocate(spin)

end program

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
