Parameter(N=100)

real L,dE,TM,TE,E,M,X,E1,C1,C2,Esq,TEsq,C,AvgE,AvgM,TMsq,Msq

Integer S(0:N+1,0:N+1,0:N+1)

do while (.true.)

   read(*,*,end=10)T

   do i=0,N+1

      do j=0,N+1

         do k=0,N+1

            S(i,j,k)=1

         enddo

      enddo

    enddo

    AJ=1.0

    TM=0.0

    TE=0.0

    TEsq=0.0

    TMsq=0.0

    R=rand(23123)

    do kmc=1,1100

       M=0.0

       E=0.0

       Esq=0.0

       Msq=0.0

      do i=1, N

         do j=1, N

            do k=1, N

               If (i==1) S(N+1,j,k) = S(i,j,k)

               If (i==N) S(0,j,k) = S(i,j,k)

               If (j==1) S(i,N+1,k) = S(i,j,k)

               If (j==N) S(i,0,k) = S(i,j,k)

               If (k==1) S(i,j,N+1) = S(i,j,k)

               If (k==N) S(i,j,0) = S(i,j,k)

               dE=2*AJ*S(i,j,k)*(S(i+1,j,k)+S(i-1,j,k)+S(i,j+1,k)+S(i,j-1,k)+S(i,j,k+1)+S(i,j,k-1))

               if (dE<=0) then

                  S(i,j,k)=(-1)*S(i,j,k)

                  else

                     L=exp(-dE/T)

                     R=rand()

                     if (R < L) then

                       S(i,j,k)=(-1)*S(i,j,k)

                     endif

                endif

             enddo

          enddo

      enddo

      if (kmc>1000) then

      do i=1, N

         do j=1, N

            do k=1, N

               M=M+S(i,j,k)

               E1=-AJ*S(i,j,k)*(S(i+1,j,k)+S(i-1,j,k)+S(i,j+1,k)+S(i,j-1,k)+S(i,j,k+1)+S(i,j,k-1))

               E=E+E1

            enddo

          enddo

       enddo

       E = E/2.0

       Esq= E*E

       Msq= M*M

       TM=TM+M

       TE=TE+E

       TESq=TEsq+Esq

       TMsq=TMsq+Msq

       else

    endif

enddo

Vinv = 1./(1.*N*N*N)

C2=TEsq*Vinv*Vinv/100.

AvgE=TE*Vinv/100.

C1=AvgE*AvgE

C=(1./T*T)*((1.*C2)-(1.*C1))

AvgM=Vinv*TM/100.

M2=AvgM*AvgM

M1=TMsq*Vinv*Vinv/100.

X=(1./T)*(M1-M2)

write(7,*)T,AvgE

write(8,*)T,AvgM

write(9,*)T,C

write(10,*)T,X

enddo

10 continue

stop

end program
