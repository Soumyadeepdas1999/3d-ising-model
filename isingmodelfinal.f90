PROGRAM isingmodelfinal
IMPLICIT NONE
REAL:: E, M, r, avg_en, avg_sq_en, avg_sq_mag, mag_sus, sums, dE, rnd
REAL, DIMENSION(4):: cum
REAl:: heat_capacity
INTEGER, PARAMETER:: J=-1
INTEGER, PARAMETER:: L=40
INTEGER, PARAMETER:: totstep=20000
INTEGER, PARAMETER:: equil=2000
INTEGER:: x,y,z, up, T, right, accept, N,i,k, left, down,w, front, back

INTEGER,DIMENSION(L,L,L):: spin
accept=0
N=L*L*L

call RANDOM_SEED()

E = 0.0
M=0.0
cum(1)=0.0
cum(2)=0.0
cum(3)=0.0
cum(4)=0.0

do x=1,L
  do y=1,L
    do z=1,L
    call RANDOM_NUMBER(r)
    if(r .lt. 0.5) then
      spin(x,y,z) = -1
    else
      spin(x,y,z) = +1
    end if
    M=M+spin(x,y,z)
    end do
  end do
end do

do x=1,L
  do y=1,L
    do z=1,L
    !impose cyclic boundary conditions and give values to variables
    up = y + 1
    right = x + 1
    front = z + 1
    !cyclic conditions
    if(x.eq.L) right = 1
    if(y .eq. L) up = 1
    if(z .eq. L) front = 1
    sums = Spin(x,up,z) + Spin(right,y,z) + Spin(x,y,front)
    E = E + J*Spin(x,y,z)*sums
    end do
  end do
end do

PRINT *, 'E', 'M'


do T=6,3,-1
   avg_en = 0.0
   avg_sq_en = 0.0
   avg_sq_mag = 0.0
   
  do i = 1,equil !bring system to equilibrium
    do x = 1,L
      do y = 1,L
        do z = 1,L
           if (x == 1) then
              left = spin(L,y,z)
              right = spin(2,y,z)
           else if (x == L) then
              left = spin(L-1,y,z)
              right = spin(1,y,z)
           else
              left = spin(x-1,y,z)
              right = spin(x+1,y,z)
           end if
           
           if (y == 1) then
              up = spin(x,2,z)
              down = spin(x,L,z)
           else if (y == L) then
              up = spin(x,1,z)
              down = spin(x,L-1,z)
           else
              up = spin(x,y+1,z)
              down = spin(x,y-1,z)
           end if
           
           if (z == 1) then
              back = spin(x,y,L)
              front = spin(x,y,2)
           else if (z == L) then
              back = spin(x,y,L-1)
              front = spin(x,y,1)
           else
              back = spin(x,y,z-1)
              front = spin(x,y,z+1)
           end if
           
           dE = 2*spin(x,y,z)*(left + right + up + down + front + back)
           call RANDOM_NUMBER(rnd)
           if (rnd .lt. exp(-dE/T)) then
              spin(x,y,z) = -spin(x,y,z)
              accept = accept + 1
              M = M + 2*spin(x,y,z)
              E = E + dE
           end if
        end do !for z
      end do !for y
    end do !for x
  end do !for equil
  do k = 1,totstep
    do x = 1,L
      do y = 1,L
        do z = 1,L
           if (x == 1) then
              left = spin(L,y,z)
              right = spin(2,y,z)
           else if (x == L) then
              left = spin(L-1,y,z)
              right = spin(1,y,z)
           else
              left = spin(x-1,y,z)
              right = spin(x+1,y,z)
           end if
           
           if (y == 1) then
              up = spin(x,2,z)
              down = spin(x,L,z)
           else if (y == L) then
              up = spin(x,1,z)
              down = spin(x,L-1,z)
           else
              up = spin(x,y+1,z)
              down = spin(x,y-1,z)
           end if
           
           if (z == 1) then
              back = spin(x,y,L)
              front = spin(x,y,2)
           else if (z == L) then
              back = spin(x,y,L-1)
              front = spin(x,y,1)
           else
              back = spin(x,y,z-1)
              front = spin(x,y,z+1)
           end if
           
           dE = 2*spin(x,y,z)*(left + right + up + down + front + back)
           call RANDOM_NUMBER(rnd)
           if (rnd .lt. exp(-dE/T)) then
              spin(x,y,z) = -spin(x,y,z)
              accept = accept + 1
              M = M + 2*spin(x,y,z)
              E = E + dE
           end if
        end do !for z
      end do !for y
    end do !for x
    cum(1) = cum(1) + E
    cum(2) = cum(2) + E*E
    cum(3) = cum(3) + M
    cum(4) = cum(4) + M*M
  end do
  avg_en = cum(1)/totstep
  avg_sq_en = cum(2)/totstep
  avg_sq_mag = cum(4)/totstep
  heat_capacity = (avg_sq_en - avg_en*avg_en)/(T*T)
  
  PRINT *, 'T',  'avg_en/N', 'avg_sq_en/N', 'avg_mag/N', 'avg_sq_mag/N'
  open(unit=1,file='latticesize40.dat')
  write(1,*)T,avg_en/N,heat_capacity/N
  
  do w=1,4
    cum(w)=0.0
  end do
  
end do !for T

close(unit=1)

end program isingmodelfinal
  
