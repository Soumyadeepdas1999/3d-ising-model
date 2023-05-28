program externalfield
implicit none
integer(8), parameter :: dp=kind(0.d0)
real(dp):: T 
real:: E, M, r, avg_en, avg_sq_en, avg_mag, avg_sq_mag, mag_sus, sums, dE, rnd
real,dimension(4)::cum
real:: heat_capacity
integer, parameter :: J = -1
integer, parameter :: L = 15
integer, parameter :: totstep = 20000
integer, parameter :: equil = 2000
integer, parameter :: B = 0.1
integer :: x,y,z, up, right, N,i,k, left, down,w, front, back
integer,dimension(L,L,L):: spin
N = L*L*L
call random_seed()
E = 0.0
M = 0.0
cum(1) =0.0
cum(2) =0.0
cum(3) =0.0
cum(4) =0.0

do x=1,L
  do y=1,L
    do z=1,L
    call RANDOM_NUMBER(r)
    if(r .lt. 0.5) then
      Spin(x,y,z) = -1
    else
      Spin(x,y,z) = +1
    end if
    M = M + Spin(x,y,z)
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
    E = E + J*Spin(x,y,z)*(sums + B)
    end do
  end do
end do

print *, 'E','M'


do T=6,35E(-1),-5E(-2)
  avg_en = 0.0
  avg_sq_en = 0.0
  avg_mag = 0.0
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
        
        dE = 2*spin(x,y,z)*(left + right + up + down + front + back + B)
        call RANDOM_NUMBER(rnd)
        if (rnd .lt. exp(-dE/T)) then
        spin(x,y,z) = -spin(x,y,z)
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
       dE = 2*spin(x,y,z)*(left + right + up + down + front + back + B)
       call RANDOM_NUMBER(rnd)
       if (rnd .lt. exp(-dE/T)) then
       spin(x,y,z) = -spin(x,y,z)
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
 avg_mag = cum(3)/totstep
 avg_sq_mag = cum(4)/totstep
 heat_capacity = (avg_sq_en - avg_en*avg_en)/(T*T)
 mag_sus = (avg_sq_mag - avg_mag*avg_mag)/T
 print *, 'T', 'avg_en/N', 'avg_sq_en/N', 'avg_mag/N', 'avg_sq_mag/n'
 open(unit=1,file='externalfield01.dat')
 write(1,*) T,avg_en/n,avg_mag/n,heat_capacity/n,mag_sus/n
 do w=1,4
   cum(w)=0.0
 end do
end do !for T

close(unit=1)

end program externalfield
