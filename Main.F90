!----------------------------------------------------------
! This Program Calculates 2D Laminar Incompressible Flow and Heat Transfer in Channel
! using SIMPLE and SIMPLER Algorithms and Staggered Uniform Grid
! Programmer  : Reza Foroozani
! Last change : 7/9/2024 - 00:54 PM
!----------------------------------------------------------
PROGRAM  laminar_incompressible_staggered_grid_and_heat_transfer
implicit none  
INTEGER ,PARAMETER :: M=500 , N=100

integer imethod,i_unsteady,inerLoopCounter,timeStepCounter,Nmax,i,j,steady,frequency
integer iconv,ienergy,i_simple,iter_show_and_write,max_iter,anim
real(8) AEu(m-1,n),AWu(m-1,n),ANu(m-1,n),ASu(m-1,n),APu(m-1,n),bu(m-1,n)
real(8) AEv(m,n-1),AWv(m,n-1),ANv(m,n-1),ASv(m,n-1),APv(m,n-1),bv(m,n-1)
real(8) AET(m,n),AWT(m,n),ANT(m,n),AST(m,n),APT(m,n),bT(m,n)     
real(8) x(0:m+1),y(0:n+1),xv(0:m+1),yv(0:n),xu(0:m),yu(0:n+1)
real(8) Wp,Wu,Gama,Re,Utop,flux,WT,GamaT,Pr,flux_T
real(8)u(0:m,0:n+1),v(0:m+1,0:n),up(0:m+1,0:n+1),vp(0:m+1,0:n+1),Ubar(0:m,0:n+1),Vbar(0:m+1,0:n)
real(8) P(0:m+1,0:n+1),Pprim(0:m+1,0:n+1)
real(8) T(0:m+1,0:n+1),T0(0:m+1,0:n+1)
real(8) error_v,error_u,errorp,errorbp,error,norm_u, norm_v,error_T        
real(8) dx,dy,lx,ly
real(8) dt,t_final,time,u_ref
real(8) P_error,u_error,v_error,UnSteadyerror_v,T_error ,norm_T ,UnSteadyerror_T        
real(8) P0(0:m+1,0:n+1),u0(0:m,0:n+1),v0(0:m+1,0:n) 
real(8) p_init,pp_init,u_init,v_init,T_init,T_top,T_down
real(8)    mn,un_err,iter_err
real(8) i_u(4),i_v(4),i_p(4),i_T(4)
real(8) T_left,T_right,top_u_velocity,down_u_velocity ,right_u_velocity,left_u_velocity ,top_v_velocity ,down_v_velocity ,right_v_velocity,left_v_velocity
real(8) start , finish
CHARACTER (LEN=10) :: TIME1, DATE1, DATE2, ZONE1,TIME0
INTEGER :: TIME2(8) 
!
common / unsteady / i_unsteady
common / length /dx,dy,lx,ly
common / index_of_boundary / i_u,i_v,i_p,i_T
common / time/dt
common / MeshVariable / x,y,xv,yv,xu,yu
common / SolveMomentom/Wu,error_u,error_v
common / VelocityVariable / AEu,AWu,ANu,ASu,bu,AEv,AWv,ANv,ASv,bv
common / VelocityComponent / u,v,u0,v0
common / DummyVelocity / Ubar,Vbar
common / PressureComponent / P,Pprim
common / coefficient / Gama,GamaT
common / pressureVariable/APu,APv,Wp,errorbp,errorp
common / ResultVelocity/up,vp
common / SolveEnergy/ WT,error_T
common / Energyvariable/ AET,AWT,ANT,AST,APT,bT 
common  /energycomponent /T,T0
common / convective / iconv
common / unsteady / anim
!
call cpu_time(start)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read input file
open (101,file="input.inp")
read (101,*)
read (101,*)
read (101,*) imethod
read(101,*) i_unsteady

if (i_unsteady==0) then
    print*,"Please Enter Maximum Iteration for Steady Flow Calculations :  "
    read(*,*) max_iter
    
    print*,"Please Enter Show and Write Iterations for This Caclulations : "
    read(*,*) iter_show_and_write
end if

read(101,*) ienergy
read(101,*) iconv
read(101,*) RE
read(101,*) pr 
read(101,*) dt
read(101,*) lx 
read(101,*) ly

read(101,*)  Wp
read(101,*)  Wu
read(101,*)  WT     

do i=1,4
read(101,*) i_u(i)    ! 1=top 2-down 3-left 4-right
end do

do i=1,4
read(101,*) i_v(i)
end do

do i=1,4        
read(101,*) i_p(i)
end do

do i=1,4
read(101,*) i_T(i)
end do


read(101,*)  u_ref

read(101,*)  top_u_velocity    
read(101,*)  down_u_velocity 
read(101,*)  right_u_velocity    
read(101,*)  left_u_velocity 
read(101,*)  top_v_velocity    
read(101,*)  down_v_velocity 
read(101,*)  right_v_velocity    
read(101,*)  left_v_velocity


read(101,*) p_init
read(101,*) pp_init
read(101,*) u_init
read(101,*) v_init
read(101,*) T_init

read(101,*) T_top
read(101,*) T_down
read(101,*) T_right
read(101,*) T_left

read(101,*) un_err
read(101,*) iter_err

read(101,*) anim

close (101)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
time=0
error=1.
inerLoopCounter=0
timeStepCounter=0


!generate mesh
call mesh

CALL DATE_AND_TIME (DATE1,TIME0,ZONE1,TIME2)


Gama=u_ref/Re

GamaT=u_ref/(Re*pr)
		
P(:,:)=p_init
Pprim(:,:)=pp_init
u(:,:)=u_init
v(:,:)=v_init
!INITIALIZING Temperature
if (ienergy==1) then
T(:,:)=T_init
end if 

! bcs velocities
u(:,n+1)=top_u_velocity    
u(:,0)=down_u_velocity 
u(m,:)=right_u_velocity    
u(0,:)=left_u_velocity 
v(:,n)=top_v_velocity    
v(:,0)=down_v_velocity 
v(0,:)=right_v_velocity    
v(m+1,:)=left_v_velocity

!bcs Temperature :  
if (ienergy==1) then
T(:,n+1)=T_top
T(:,0)=T_down 
T(0,:)=T_left
T(m+1,:)=T_right
end if





!loop for transit solution (external loop)
!**********************Unsteady loop****************************
!                      _____________		

timeStepCounter=0
UnSteadyerror_v=1.0
UnSteadyerror_T=1.0
steady=0

! Starting external transient loop
if (i_unsteady==1) then
frequency=1/dt
write(*,*)frequency
write(*,'(5x,a,/)')"timeStepCounter  inerLoopCounter      error_u          error_v          error_p"
write(*,'(5x,a)')"(number of repeat (number of repeat  error in calculated error in calculated"
write(*,'(5x,a,/)')" external loop)    simple loop)     U velocity filed   V velocity filed   P pressure filed"
    do while ( UnSteadyerror_v.gt.un_err .or. Unsteadyerror_T.gt.un_err )   !External transient Loop
	time=time+dt
    timeStepCounter=timeStepCounter+1 
    inerLoopCounter=0
	error=1.
    u0(:,:)=u(:,:)
	v0(:,:)=v(:,:)
    T0(:,:)=T(:,:)

	do while(error>iter_err)              !Begining of Internal Loop  

		inerLoopCounter=inerLoopCounter+1
!**********************SIMPLE Algorithm****************************
!                      _________________		
	
		!calculate the coefficient of momentum Eq.
		 CALL coefficient_momentom

        !calculate the coefficient of pressure Eq. & solve it
		if (imethod==1) CALL pressure 		 

		!solve the momentum Eq.
		 CALL solve_momentum

		!solve the pressure_correction Eq. 
		 CALL pressure_correction 

        !correct the velocity field
		 CALL CorrectVelocity (APu,APv)


        !calculate the coefficient of energy Eq.
        if (ienergy==1)        CALL coefficient_energy
        
        !solve the energy Eq.(To find temperature distribution)
        if (ienergy==1)        CALL solve_energy        

        
        
		error=max ( error_u,error_v,error_T)          
		
		write(*,1) timeStepCounter,inerLoopCounter,error_u,error_v,error_T
	 1	format (5x,I10,8x,I10,5x,3(6x,D12.6))


		
end do  !End of Inner Loop or simple loop      




        P_error=0.0
        u_error=0.0
        v_error=0.0
        T_error=0.0
		flux=0.0
		UnSteadyerror_v=0.0
        UnSteadyerror_T=0.0


do i=0,m
	   do J=0,n
	      u_error=u_error+(u(i,j)-u0(i,j))**2.0
		  v_error=v_error+(v(i,j)-v0(i,j))**2.0
          T_error=T_error+(T(i,j)-T0(i,j))**2.0
		  norm_u=sqrt(u_error)
		  norm_v=sqrt(v_error)
          norm_T=sqrt(T_error)
	    end do
	end do

    UnSteadyerror_v=max(norm_u, norm_v)
    UnSteadyerror_T=norm_T

    u0(:,:)=u(:,:)
	v0(:,:)=v(:,:)
	P0(:,:)=P(:,:)
    T0(:,:)=T(:,:)
	
  

  write(*,'(5x,a,f10.5)')'time of soulution is:',time
  write(*,'(5x,a,f10.7,//)')"UnSteady error in this tim step is",UnSteadyerror_v
  write(*,'(5x,a,f10.7,//)')"UnSteady error_T in this tim step is",UnSteadyerror_T

  mn=1
  steady=1
  

  
end do  !unsteady loop
call result(mn,steady)
end if    !unsteady loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!steady loop 
!
if (i_unsteady == 0 ) then 
    inerLoopCounter=0
	error=1.
    u0(:,:)=u(:,:)
	v0(:,:)=v(:,:)
    T0(:,:)=T(:,:)

    
open (103,file="residual_steady.plt")    
	do while(error>iter_err)              !Begining of Internal Loop  

		inerLoopCounter=inerLoopCounter+1
        
!        if (inerLoopCounter>10000 ) exit
!**********************SIMPLE Algorithm****************************
!                      _________________		
	
		!calculate the coefficient of momentum Eq. 
		 CALL coefficient_momentom
        
         !calculate the coefficient of pressure Eq. & solve it
		if (imethod==1) CALL pressure 		 
		 

		!solve the momentum Eq.
		 CALL solve_momentum

		!solve the pressure_correction Eq. 
		CALL pressure_correction 

        !correct the velocity field
		 CALL CorrectVelocity (APu,APv)

        
 
        !calculate the coefficient of energy Eq.
        if (ienergy==1)        CALL coefficient_energy
        
        !solve the energy Eq.(To find temperature distribution)
        if (ienergy==1)        CALL solve_energy       

        
        
		error=max ( error_u,error_v,error_T)          
		

        if (mod(inerloopcounter,iter_show_and_write)==0) then
              mn=1
             steady=1
                    CALL DATE_AND_TIME(DATE2,TIME1,ZONE1,TIME2)
 
                    PRINT*,'_______________________________'
                    PRINT*,'Iter No.:',inerloopcounter,"steady"
                    PRINT*,'P Err:   ',errorbp
                    PRINT*,'U Err:   ',error_u
                    PRINT*,'V Err:   ',error_v
                    PRINT*,'T Err:   ',error_T                    
                    PRINT*,'    Start Date: ',DATE1,' , Last Date: ',DATE2
                    PRINT*,'    Start Time: ',TIME0,' , Last Time:',TIME1
                    PRINT*
                    print*,"start writing ..."
            call result(mn,steady)
            print*,"end writing ..."

            end if 
        
		if (inerLoopCounter==max_iter ) then
            exit 
        end if
        
        write(103,*) inerloopcounter,log(max(error_u,error_v,error_T))

		
    end do  !End of Inner Loop or simple loop  
close(103)    
end if ! steady loop

!
  mn=1
  steady=1
  
call result(mn,steady)

call cpu_time(finish)
print*,"cpu time is " , finish-start
print*,"press enter to finish"
read(*,*)


end PROGRAM  laminar_incompressible_staggered_grid_and_heat_transfer



!******************** subroutine mesh**********************
!                     _______________
subroutine mesh
implicit none 
INTEGER ,PARAMETER :: M=500 , N=100
integer i,j
REAL(8) x(0:m+1),y(0:n+1),xv(0:m+1),yv(0:n),xu(0:m),yu(0:n+1)
real(8) dx,dy,lx,ly

common / length /dx,dy,lx,ly
common / MeshVariable / x,y,xv,yv,xu,yu

!GEOMETRY OF VARIABLE (Pressure and Velocity) 
dx=lx/(m+1)
dy=ly/(n+1)
! Constructing Scalar cells
!
!position of scaler variable cells
!__________________________________
!X
!***************

! X positions of U cells
!____________________
!X
xu(0)=0
do i=1,m
    xu(i)=xu(i-1)+dx

end do 

! Y positions of V cells
!____________________
!Y
yv(0)=0
do j=1,n
       yv(j)=yv(j-1)+dy

end do
!
! X positions of Scalar cells(forward staggered grids)
!
x(0)=0
do i=1,m
   x(i)=(xu(i-1)+xu(i))/2.
end do 
x(m+1)=x(m)+dx/2
!
! Y positions of Scalar cells(forward staggered grids)
!
y(0)=0
do i=1,n
   y(i)=(yv(i-1)+yv(i))/2.
end do 
y(n+1)=y(n)+(y(n-1)-y(n-2))

!
!Y position of U cells
!
do i=0,n+1
  yu(i)=y(i)
end do
!
!
!X position of V cells
!
do i=0,m+1
   xv(i)=x(i)
end do 
!*****************
return 
end
   
   


!******************** coefficient_momentum**********************
!                     ____________________

!calculate the coefficient of momentum Eq.

subroutine  coefficient_momentom 
implicit none 
INTEGER ,PARAMETER :: M=500 , N=100

integer i,j,i_unsteady,iconv
real(8),dimension(0:m+1,0:n+1) ::  De,Dw,Dn,Ds,Fe,Fw,Fn,Fs
real(8) AEu(m-1,n),AWu(m-1,n),ANu(m-1,n),ASu(m-1,n),APu(m-1,n),bu(m-1,n),A
real(8) AEv(m,n-1),AWv(m,n-1),ANv(m,n-1),ASv(m,n-1),APv(m,n-1),bv(m,n-1)
real(8) x(0:m+1),y(0:n+1),xv(0:m+1),yv(0:n),xu(0:m),yu(0:n+1)
real(8) u(0:m,0:n+1),v(0:m+1,0:n),u0(0:m,0:n+1),v0(0:m+1,0:n)
real(8) Ubar(0:m,0:n+1),Vbar(0:m+1,0:n)
real(8) T(0:m+1,0:n+1)
real(8) AP0u(m-1,n),su_u(m-1,n),sp_u(m-1,n),AP0v(m,n-1),su_v(m,n-1),sp_v(m,n-1)
real(8) dx,dy,dt,Gama,Wp,errorbp, ALPHA_e, ALPHA_w, ALPHA_n, ALPHA_s
real(8) i_u(4),i_v(4),i_p(4),i_T(4),pee,pew,pen,pes,sorc1(0:m+1,0:n+1),sorc2(0:m+1,0:n+1)
real(8) alphae,alphaw,alphan,alphas

common / convective / iconv
common / unsteady / i_unsteady
common / length /dx,dy
common / index_of_boundary / i_u,i_v,i_p,i_T
common / time/dt
common / VelocityVariable/ AEu,AWu,ANu,ASu,bu,AEv,AWv,ANv,ASv,bv
common / MeshVariable / x,y,xv,yv,xu,yu
common / VelocityComponent / u,v,u0,v0
common / DummyVelocity / Ubar,Vbar
common / coefficient / Gama
common / pressureVariable/APu,APv,Wp,errorbp


!COMPUTING COEFFICIENT OF x-MOMENTUM
do j=1,n
	do i=1,m-1
		
		Fe(i,j)=(yv(j)-yv(j-1))*(u(i,j)+u(i+1,j))/2.       
		Fw(i,j)=(yv(j)-yv(j-1))*(u(i,j)+u(i-1,j))/2.      
		Fn(i,j)=(xu(i)-xu(i-1))*(v(i,j)+v(i+1,j))/2.       
		Fs(i,j)=(xu(i)-xu(i-1))*(v(i,j-1)+v(i+1,j-1))/2.  
        
		De(i,j)=Gama*(yu(j)-yu(j-1))/(x(i)-x(i-1))
        Dw(i,j)=Gama*(yu(j)-yu(j-1))/(x(i)-x(i-1))
		Dn(i,j)=Gama*(xv(i)-xv(i-1))/(y(j)-y(j-1))
		Ds(i,j)=Gama*(xv(i)-xv(i-1))/(y(j)-y(j-1))
        
        if (iconv==1) then !upwind
            
        AEu(i,j)=de(i,j)+max(-fe(i,j),0.0)
		AWu(i,j)=dw(i,j)+max(fw(i,j),0.0)
		ANu(i,j)=dn(i,j)+max(-fn(i,j),0.0)
		ASu(i,j)=ds(i,j)+max(fs(i,j),0.0)
		AP0u(i,j)=dx*dy/dt*i_unsteady        
        sorc1(i,j)=fe(i,j)-fw(i,j)+fn(i,j)-fs(i,j)
        
        else if (iconv==2) then !hybrid
            
        AEu(i,j)=max(-fe(i,j),(de(i,j)-(fe(i,j)/2)),0.0)
		AWu(i,j)=max(fw(i,j),(dw(i,j)+(fw(i,j)/2)),0.0)
		ANu(i,j)=max(-fn(i,j),(dn(i,j)-(fn(i,j)/2)),0.0)
		ASu(i,j)=max(fs(i,j),(ds(i,j)+(fs(i,j)/2)),0.0)
		AP0u(i,j)=dx*dy/dt*i_unsteady        
        sorc1(i,j)=fe(i,j)-fw(i,j)+fn(i,j)-fs(i,j)
        
        else if (iconv==3) then !powerlaw    
            pee=fe(i,j)/de(i,j)
            pew=fw(i,j)/dw(i,j)
            pes=fs(i,j)/ds(i,j)
            pen=fn(i,j)/dn(i,j)
            
        AEu(i,j)=   de(i,j)*max(0.0,(1.0-0.1*pee)**5) + max(-fe(i,j),0.0)
		AWu(i,j)=  dw(i,j)*max(0.0,(1.0-0.1*pew)**5) + max(fw(i,j),0.0)
		ANu(i,j)=  dn(i,j)*max(0.0,(1.0-0.1*pen)**5) +max(-fn(i,j),0.0)
		ASu(i,j)=   ds(i,j)*max(0.0,(1.0-0.1*pes)**5) + max(fs(i,j),0.0) 
		AP0u(i,j)=dx*dy/dt*i_unsteady        
        sorc1(i,j)=fe(i,j)-fw(i,j)+fn(i,j)-fs(i,j)

        else if (iconv==4) then ! Quick Hayase
            

        IF (Fw(i,j) > 0) THEN
            ALPHAw = 1
          ELSE
            ALPHAw = 0
        END IF
        IF (Fe(i,j) > 0) THEN
            ALPHAe = 1
          ELSE
            ALPHAe = 0
        END IF
        IF (Fs(i,j) > 0) THEN
            ALPHAs = 1
          ELSE
            ALPHAs = 0
        END IF
        IF (Fn(i,j) > 0) THEN
            ALPHAn = 1
          ELSE
            ALPHAn = 0
        END IF
    
		AEu(i,j)=De(i,j)-(1-ALPHAe)*Fe(i,j)
		AWu(i,j)=Dw(i,j)+ALPHAw*Fw(i,j)
		ANu(i,j)=Dn(i,j)-(1-ALPHAn)*Fn(i,j)
		ASu(i,j)=Ds(i,j)+ALPHAs*Fs(i,j)
		AP0u(i,j)=dx*dy/dt*i_unsteady        
        sorc1(i,j)=fe(i,j)-fw(i,j)+fn(i,j)-fs(i,j)
        end if             
            
               
	end do
end do



!COMPUTING COEFFICIENT OF y-MOMENTUM
do i=1,m
	do j=1,n-1
	
		Fe(i,j)=(yv(j)-yv(j-1))*(u(i,j)+u(i,j+1))/2.       
		Fw(i,j)=(yv(j)-yv(j-1))*(u(i-1,j)+u(i-1,j+1))/2.      
		Fn(i,j)=(xu(i)-xu(i-1))*(v(i,j)+v(i,j+1))/2.       
		Fs(i,j)=(xu(i)-xu(i-1))*(v(i,j)+v(i,j-1))/2.  
        
		De(i,j)=Gama*(yu(j)-yu(j-1))/(x(i)-x(i-1))
        Dw(i,j)=Gama*(yu(j)-yu(j-1))/(x(i)-x(i-1))
		Dn(i,j)=Gama*(xv(i)-xv(i-1))/(y(j)-y(j-1))
		Ds(i,j)=Gama*(xv(i)-xv(i-1))/(y(j)-y(j-1))

        if (iconv==1) then !upwind
            
        AEv(i,j)=de(i,j)+max(-fe(i,j),0.0)
		AWv(i,j)=dw(i,j)+max(fw(i,j),0.0)
		ANv(i,j)=dn(i,j)+max(-fn(i,j),0.0)
		ASv(i,j)=ds(i,j)+max(fs(i,j),0.0)
		AP0v(i,j)=dx*dy/dt*i_unsteady        
        sorc2(i,j)=fe(i,j)-fw(i,j)+fn(i,j)-fs(i,j)
        
        else if (iconv==2) then !hybrid
            
        AEv(i,j)=max(-fe(i,j),(de(i,j)-(fe(i,j)/2)),0.0)
		AWv(i,j)=max(fw(i,j),(dw(i,j)+(fw(i,j)/2)),0.0)
		ANv(i,j)=max(-fn(i,j),(dn(i,j)-(fn(i,j)/2)),0.0)
		ASv(i,j)=max(fs(i,j),(ds(i,j)+(fs(i,j)/2)),0.0)
		AP0v(i,j)=dx*dy/dt*i_unsteady        
        sorc2(i,j)=fe(i,j)-fw(i,j)+fn(i,j)-fs(i,j)
        
        else if (iconv==3) then !powerlaw    
            pee=fe(i,j)/de(i,j)
            pew=fw(i,j)/dw(i,j)
            pes=fs(i,j)/ds(i,j)
            pen=fn(i,j)/dn(i,j)
            
        AEv(i,j)=   de(i,j)*max(0.0,(1.0-0.1*pee)**5) + max(-fe(i,j),0.0)
		AWv(i,j)=  dw(i,j)*max(0.0,(1.0-0.1*pew)**5) + max(fw(i,j),0.0)
		ANv(i,j)=  dn(i,j)*max(0.0,(1.0-0.1*pen)**5) +max(-fn(i,j),0.0)
		ASv(i,j)=   ds(i,j)*max(0.0,(1.0-0.1*pes)**5) + max(fs(i,j),0.0) 
		AP0v(i,j)=dx*dy/dt*i_unsteady        
        sorc2(i,j)=fe(i,j)-fw(i,j)+fn(i,j)-fs(i,j)

        else if (iconv==4) then ! Quick Hayase
            

        IF (Fw(i,j) > 0) THEN
            ALPHAw = 1
          ELSE
            ALPHAw = 0
        END IF
        IF (Fe(i,j) > 0) THEN
            ALPHAe = 1
          ELSE
            ALPHAe = 0
        END IF
        IF (Fs(i,j) > 0) THEN
            ALPHAs = 1
          ELSE
            ALPHAs = 0
        END IF
        IF (Fn(i,j) > 0) THEN
            ALPHAn = 1
          ELSE
            ALPHAn = 0
        END IF
    
		AEv(i,j)=De(i,j)-(1-ALPHAe)*Fe(i,j)
		AWv(i,j)=Dw(i,j)+ALPHAw*Fw(i,j)
		ANv(i,j)=Dn(i,j)-(1-ALPHAn)*Fn(i,j)
		ASv(i,j)=Ds(i,j)+ALPHAs*Fs(i,j)
		AP0v(i,j)=dx*dy/dt*i_unsteady        
        sorc2(i,j)=fe(i,j)-fw(i,j)+fn(i,j)-fs(i,j)
        end if     
    end do
end do


!SOURCE TERM OF y-COMPONENT
su_v(:,:)=0. ;sp_v(:,:)=0.


!SOURCE TERM OF x-COMPONENT
su_u(:,:)=0. ;sp_u(:,:)=0.



do j=1,n
    do i=1,m-1

        if (iconv==4) then
            
            IF (Fw(i,j) > 0) THEN
                ALPHAw = 1
            ELSE
                ALPHAw = 0
            END IF
            IF (Fe(i,j) > 0) THEN
                ALPHAe = 1
            ELSE
                ALPHAe = 0
            END IF
            IF (Fs(i,j) > 0) THEN
                ALPHAs = 1
            ELSE
                ALPHAs = 0
            END IF
            IF (Fn(i,j) > 0) THEN
                ALPHAn = 1
            ELSE
                ALPHAn = 0
            END IF

            IF(i==1 .and. j==1) THEN
                bu(i,j)=sorc1(i,j)+AP0u(i,j)*u0(i,j)+su_u(i,j)+ &
                        1/8*Fw(i,j)*ALPHAw*(3*u(i,j)-2*u(i-1,j)-(-u(i,j)))+ &
                        1/8*Fe(i,j)*ALPHAe*(u(i-1,j)+2*u(i,j)-3*u(i+1,j))+ &
                        1/8*Fw(i,j)*(1-ALPHAw)*(3*u(i-1,j)-2*u(i,j)-u(i+1,j))+ &
                        1/8*Fe(i,j)*(1-ALPHAe)*(2*u(i+1,j)-u(i+2,j)-3*u(i,j))+ &
                        1/8*Fs(i,j)*ALPHAs*(3*u(i,j)-2*u(i,j-1)-u(i,j-2))+ &
                        1/8*Fn(i,j)*ALPHAn*(u(i,j-1)+2*u(i,j)-3*u(i,j+1))
            ELSE
                bu(i,j)=sorc1(i,j)+AP0u(i,j)*u0(i,j)+su_u(i,j)+ &
                        1/8*Fw(i,j)*ALPHAw*(3*u(i,j)-2*u(i-1,j)-u(i-2,j))+ &
                        1/8*Fe(i,j)*ALPHAe*(u(i-1,j)+2*u(i,j)-3*u(i+1,j))+ &
                        1/8*Fw(i,j)*(1-ALPHAw)*(3*u(i-1,j)-2*u(i,j)-u(i+1,j))+ &
                        1/8*Fe(i,j)*(1-ALPHAe)*(2*T(i+1,j)+(2*u(i+1,j)-u(i,j))-3*T(i,j))+ &
                        1/8*Fs(i,j)*ALPHAs*(3*u(i,j)-2*u(i,j-1)-u(i,j-2))+ &
                        1/8*Fn(i,j)*ALPHAn*(u(i,j-1)+2*u(i,j)-3*u(i,j+1))
            END IF
        end if
        
        APu(i,j)=AEu(i,j)+AWu(i,j)+ANu(i,j)+ASu(i,j)+AP0u(i,j)-sp_u(i,j) +sorc1(i,j)

        if (iconv/=4) then
            bu(i,j)=sorc1(i,j)+AP0u(i,j)*u0(i,j)
        end if
    end do
end do

do i=1,m
    do j=1,n-1
        if (iconv==4) then
            
            IF (Fw(i,j) > 0) THEN
                ALPHAw = 1
            ELSE
                ALPHAw = 0
            END IF
            IF (Fe(i,j) > 0) THEN
                ALPHAe = 1
            ELSE
                ALPHAe = 0
            END IF
            IF (Fs(i,j) > 0) THEN
                ALPHAs = 1
            ELSE
                ALPHAs = 0
            END IF
            IF (Fn(i,j) > 0) THEN
                ALPHAn = 1
            ELSE
                ALPHAn = 0
            END IF

            IF(i==1 .and. j==1) THEN
                bv(i,j)=sorc2(i,j)+AP0v(i,j)*v0(i,j)+su_v(i,j)+ &
                        1/8*Fw(i,j)*ALPHAw*(3*v(i,j)-2*v(i-1,j)-(-v(i,j)))+ &
                        1/8*Fe(i,j)*ALPHAe*(v(i-1,j)+2*v(i,j)-3*v(i+1,j))+ &
                        1/8*Fw(i,j)*(1-ALPHAw)*(3*v(i-1,j)-2*v(i,j)-v(i+1,j))+ &
                        1/8*Fe(i,j)*(1-ALPHAe)*(2*v(i+1,j)-v(i+2,j)-3*v(i,j))+ &
                        1/8*Fs(i,j)*ALPHAs*(3*v(i,j)-2*v(i,j-1)-v(i,j-2))+ &
                        1/8*Fn(i,j)*ALPHAn*(v(i,j-1)+2*v(i,j)-3*v(i,j+1))
            ELSE
                bv(i,j)=sorc2(i,j)+AP0v(i,j)*v0(i,j)+su_v(i,j)+ &
                        1/8*Fw(i,j)*ALPHAw*(3*v(i,j)-2*v(i-1,j)-v(i-2,j))+ &
                        1/8*Fe(i,j)*ALPHAe*(v(i-1,j)+2*v(i,j)-3*v(i+1,j))+ &
                        1/8*Fw(i,j)*(1-ALPHAw)*(3*v(i-1,j)-2*v(i,j)-v(i+1,j))+ &
                        1/8*Fe(i,j)*(1-ALPHAe)*(2*T(i+1,j)+(2*v(i+1,j)-v(i,j))-3*T(i,j))+ &
                        1/8*Fs(i,j)*ALPHAs*(3*v(i,j)-2*v(i,j-1)-v(i,j-2))+ &
                        1/8*Fn(i,j)*ALPHAn*(v(i,j-1)+2*v(i,j)-3*v(i,j+1))
            END IF
        end if
        
        if (iconv/=4) then
            bv(i,j)=sorc2(i,j)+AP0v(i,j)*v0(i,j)
        end if
        
        APv(i,j)=AEv(i,j)+AWv(i,j)+ANv(i,j)+ASv(i,j)+AP0v(i,j)-sp_v(i,j)+sorc2(i,j)
    end do
end do
     
!calculate dummy velocity
Ubar(:,:)=u(:,:)		
Vbar(:,:)=v(:,:)		
do j=1,n
	do i=1,m-1
		Ubar(i,j)=(AEu(i,j)*u(i+1,j)+AWu(i,j)*u(i-1,j)+ANu(i,j)*u(i,j+1)+ASu(i,j)*u(i,j-1)+bu(i,j))/APu(i,j)
	end do
end do
do i=1,m
	do j=1,n-1
		Vbar(i,j)=(AEv(i,j)*v(i+1,j)+AWv(i,j)*v(i-1,j)+ANv(i,j)*v(i,j+1)+ASv(i,j)*v(i,j-1)+bv(i,j))/APv(i,j)
	end do
end do

RETURN
end
!********************pressure**********************
!                   ___________
!calculate the coefficient of pressure Eq. & solve it

subroutine pressure 
implicit none 
INTEGER ,PARAMETER :: M=500 , N=100

real(8) AE(m,n),AW(m,n),AN(m,n),AS(m,n),AP(m,n),b(m,n),su(m,n),sp(m,n)
real(8) AX(m),BX(m-1),CX(2:m),QX(m),AY(n),BY(n-1),CY(2:n),QY(n)
integer i,j,iconv

real(8) Ubar(0:m,0:n+1),Vbar(0:m+1,0:n),P(0:m+1,0:n+1),Pprim(0:m+1,0:n+1)
real(8) errorp,errorbp,Wp
real(8) APu(m-1,n),APv(m,n-1)
real(8) dx,dy

common / convective / iconv
common / length /dx,dy
common / PressureComponent / P,Pprim
common / DummyVelocity / Ubar,Vbar
common / pressureVariable/APu,APv,Wp,errorbp,errorp

!COMPUTING COEFFICIENT OF PRESSURE EQ.
do j=1,n
	do i=1,m-1
		AE(i,j)=dy**2.D0/APu(i,j)
	end do
end do
do j=1,n
	do i=2,m
		AW(i,j)=dy**2.D0/APu(i-1,j)
	end do
end do
do i=1,m
	do j=1,n-1
		AN(i,j)=dx**2.D0/APv(i,j)
	end do
end do
do i=1,m
	do j=2,n
		AS(i,j)=dx**2.D0/APv(i,j-1)
	end do
end do

! Allocating boundary side neighboring coefficients to zerro
AS(:,1)=0.
AN(:,n)=0.
AW(1,:)=0.
AE(m,:)=0.

!computing source term, bP(i,j) and central cofficient, AP(i,j)
su(:,:)=0.
sp(:,:)=0.

do j=1,n
	do i=1,m
		b(i,j)=(Ubar(i-1,j)-Ubar(i,j))*dy+(Vbar(i,j-1)-Vbar(i,j))*dx+su(i,j)
		AP(i,j)=AE(i,j)+AW(i,j)+AN(i,j)+AS(i,j)-sp(i,j)
	end do
end do

!SOLVE THE PRESSURE EQ. WITH AADI METHOD
AX(:)=0. ;BX(:)=0. ;CX(:)=0. ;QX(:)=0.
AY(:)=0. ;BY(:)=0. ;CY(:)=0. ;QY(:)=0.
!SWEEP Y
do j=1,n
	do i=1,m
		AX(i)=AP(i,j)/Wp  
		QX(i)=AN(i,j)*P(i,j+1)+AS(i,j)*P(i,j-1)+b(i,j)+(1.D0-Wp)/Wp*AP(i,j)*P(i,j)
	end do
	QX(1)=QX(1)+AW(1,j)*P(0,j)   !AW(1,:) already was set to zero so it is not necessary to change QX(1)!!
	QX(m)=QX(m)+AE(m,j)*P(m+1,j) !AE(m,j) already was set to zero so it is not necessary to change QX(m)!!
	do i=1,m-1
		BX(i)=-AE(i,j)
	end do
	do i=2,m
		CX(i)=-AW(i,j)
	end do
	CALL Crute_method(m,AX,BX,CX,QX)
	do i=1,m
		P(i,j)=QX(i)
	end do
	AX(:)=0. ;BX(:)=0. ;CX(:)=0. ;QX(:)=0.
end do

!SWEEP X
do i=1,m
	do j=1,n
		AY(j)=AP(i,j)/Wp
		QY(j)=AE(i,j)*P(i+1,j)+AW(i,j)*P(i-1,j)+b(i,j)+(1.D0-Wp)/Wp*AP(i,j)*P(i,j)
	end do
	QY(1)=QY(1)+AS(i,1)*P(i,0)    !AS(:,1) already was set to zero so it is not necessary to change QY(1)!!
	QY(n)=QY(n)+AN(i,n)*P(i,n+1)  !AN(:,n) already was set to zero so it is not necessary to change QY(1)!!
	do j=1,n-1
		BY(j)=-AN(i,j)
	end do
	do j=2,n
		CY(j)=-AS(i,j)
	end do
	CALL Crute_method(n,AY,BY,CY,QY)
	do j=1,n
		P(i,j)=QY(j)
	end do
	AY(:)=0. ;BY(:)=0. ;CY(:)=0. ;QY(:)=0.
end do
!ONE STEP 
!

RETURN
end

!********************solve_momentum*********
!                    _______________

!solve the momentum Eqs.

subroutine  solve_momentum
implicit none 
INTEGER ,PARAMETER :: M=500 , N=100
real(8) AX(m-1),BX(m-2),CX(2:m-1),QX(m-1),AY(n),BY(n-1),CY(2:n),QY(n)
real(8) AXY(m),BXY(m-1),CXY(2:m),QXY(m),AYY(n-1),BYY(n-2),CYY(2:n-1),QYY(n-1)
integer i,j,ARVE1(2),ARVE2(2)
real(8) AEu(m-1,n),AWu(m-1,n),ANu(m-1,n),ASu(m-1,n),APu(m-1,n),bu(m-1,n)
real(8) AEv(m,n-1),AWv(m,n-1),ANv(m,n-1),ASv(m,n-1),APv(m,n-1),bv(m,n-1)
real(8) x(0:m+1),y(0:n+1)
real(8) u(0:m,0:n+1),v(0:m+1,0:n),u0(0:m,0:n+1),v0(0:m+1,0:n)
!real(8) u_2(0:m,0:n+1),v_2(0:m+1,0:n)
real(8) P(0:m+1,0:n+1)
real(8) error_u,error_v,Wu
real(8) Flux_u,Flux_v,Residu_u,Residu_v
real(8) dx,dy,Wp,errorbp
real(8) i_u(4),i_v(4),i_p(4),i_T(4)

common / SolveMomentom/Wu,error_u,error_v
common / length /dx,dy
common / index_of_boundary / i_u,i_v,i_p,i_T
common / VelocityVariable/ AEu,AWu,ANu,ASu,bu,AEv,AWv,ANv,ASv,bv
common / MeshVariable / x,y
common / VelocityComponent / u,v,u0,v0
common / PressureComponent / P
common / pressureVariable/APu,APv,Wp,errorbp

AX(:)=0. ;BX(:)=0. ;CX(:)=0. ;QX(:)=0.
AY(:)=0. ;BY(:)=0. ;CY(:)=0. ;QY(:)=0.

!SOLVE THE x-MUMENTUM EQ. WITH AADI METHOD
!SWEEP Y
do j=1,n
	do i=1,m-1
        if (i_u(3)==0) then
        u(0,:)=u(1,:)
        end if 
        if (i_u(4)==0) then
        u(m,:)=u(m-1,:)
        end if
        if (i_u(1)==0) then
        u(:,n+1)=u(:,n)
        end if 
        if (i_u(2)==0) then
        u(:,0)=u(:,1)
        end if
        
		AX(i)=APu(i,j)/Wu
		QX(i)=ANu(i,j)*u(i,j+1)+ASu(i,j)*u(i,j-1)+bu(i,j)+(1.D0-Wu)/Wu*APu(i,j)*u(i,j)+(P(i,j)-P(i+1,j))*dy
	end do
	QX(1)=QX(1)+AWu(1,j)*u(0,j)          !u(0,j) is the U velocity of the left vertical wall
	QX(m-1)=QX(m-1)+AEu(m-1,j)*u(m,j)    !u(m,j) is the U velocity of the right vertical wall
	do i=1,m-2
		BX(i)=-AEu(i,j)
	end do
	do i=2,m-1
		CX(i)=-AWu(i,j)
	end do
	CALL Crute_method(m-1,AX,BX,CX,QX)
	do i=1,m-1
		u(i,j)=QX(i)
	end do
	AX(:)=0. ;BX(:)=0. ;CX(:)=0. ;QX(:)=0.
end do

!SWEEP X
do i=1,m-1
	do j=1,n
        if (i_u(3)==0) then
        u(0,:)=u(1,:)
        end if 
        if (i_u(4)==0) then
        u(m,:)=u(m-1,:)
        end if
        if (i_u(1)==0) then
        u(:,n+1)=u(:,n)
        end if 
        if (i_u(2)==0) then
        u(:,0)=u(:,1)
        end if
        AY(j)=APu(i,j)/Wu
		QY(j)=AEu(i,j)*u(i+1,j)+AWu(i,j)*u(i-1,j)+bu(i,j)+(1.D0-Wu)/Wu*APu(i,j)*u(i,j)+(P(i,j)-P(i+1,j))*dy
	end do
	QY(1)=QY(1)+ASu(i,1)*u(i,0)      !u(i,0) is the U velocity of the bottom horizontal wall
	QY(n)=QY(n)+ANu(i,n)*u(i,n+1)    !u(i,n+1) is the U velocity of the top horizontal wall
	do j=1,n-1
		BY(j)=-ANu(i,j)
	end do
	do j=2,n
		CY(j)=-ASu(i,j)
	end do
	CALL Crute_method(n,AY,BY,CY,QY)
	do j=1,n
		u(i,j)=QY(j)
	end do
	AY(:)=0. ;BY(:)=0. ;CY(:)=0. ;QY(:)=0.
end do
!ONE STEP

!SOLVE THE y-MUMENTUM EQ. WITH AADI METHOD
!SWEEP Y
do j=1,n-1
	do i=1,m
        if (i_v(3)==0) then
        v(0,:)=v(1,:)
        end if 
        if (i_v(4)==0) then
        v(m+1,:)=v(m,:)
        end if
        if (i_v(1)==0) then
        v(:,n)=v(:,n-1)
        end if 
        if (i_v(2)==0) then
        v(:,0)=v(:,1)
        end if
        
		AXY(i)=APv(i,j)/Wu
		QXY(i)=ANv(i,j)*v(i,j+1)+ASv(i,j)*v(i,j-1)+bv(i,j)+(1.D0-Wu)/Wu*APv(i,j)*v(i,j)+(P(i,j)-P(i,j+1))*dx
	end do
	QXY(1)=QXY(1)+AWv(1,j)*v(0,j)      !v(0,j) is the V velocity of the left vertical wall
	QXY(m)=QXY(m)+AEv(m,j)*v(m+1,j)    !v(m+1,j) is the V velocity of the right vertical wall
	do i=1,m-1
		BXY(i)=-AEv(i,j)
	end do
	do i=2,m
		CXY(i)=-AWv(i,j)
	end do
	CALL Crute_method(m,AXY,BXY,CXY,QXY)
	do i=1,m
		v(i,j)=QXY(i)
	end do
	AXY(:)=0. ;BXY(:)=0. ;CXY(:)=0. ;QXY(:)=0.
end do

!SWEEP X
do i=1,m
	do j=1,n-1
        if (i_v(3)==0) then
        v(0,:)=v(1,:)
        end if 
        if (i_v(4)==0) then
        v(m+1,:)=v(m,:)
        end if
        if (i_v(1)==0) then
        v(:,n)=v(:,n-1)
        end if 
        if (i_v(2)==0) then
        v(:,0)=v(:,1)
        end if
        AYY(j)=APv(i,j)/Wu
		QYY(j)=AEv(i,j)*v(i+1,j)+AWv(i,j)*v(i-1,j)+bv(i,j)+(1.D0-Wu)/Wu*APv(i,j)*v(i,j)+(P(i,j)-P(i,j+1))*dx
	end do
	QYY(1)=QYY(1)+ASv(i,1)*v(i,0)        !v(i,0) is the V velocity of the bottom horizontal wall
	QYY(n-1)=QYY(n-1)+ANv(i,n-1)*v(i,n)  !u(i,n) is the V velocity of the top horizontal wall
	do j=1,n-2
		BYY(j)=-ANv(i,j)
	end do
	do j=2,n-1
		CYY(j)=-ASv(i,j)
	end do
	CALL Crute_method(n-1,AYY,BYY,CYY,QYY)
	do j=1,n-1
		v(i,j)=QYY(j)
	end do
	AYY(:)=0. ;BYY(:)=0. ;CYY(:)=0. ;QYY(:)=0.
end do
! ONE STEP

Flux_u=0.D0 ;Flux_v=0.D0 ;Residu_u=0.D0 ;Residu_v=0.D0
do i=1,m-1
	do j=1,n
		Residu_u=Residu_u+DABS(ANu(i,j)*u(i,j+1)+ASu(i,j)*u(i,j-1)+AEu(i,j)*u(i+1,j)+AWu(i,j)*u(i-1,j)+bu(i,j)+(P(i,j)-P(i+1,j))*dy-APu(i,j)*u(i,j))
		Flux_u=Flux_u+DABS(APu(i,j)*u(i,j))
	end do
end do
do i=1,m
	do j=1,n-1
		Residu_v=Residu_v+DABS(ANv(i,j)*v(i,j+1)+ASv(i,j)*v(i,j-1)+AEv(i,j)*v(i+1,j)+AWv(i,j)*v(i-1,j)+bv(i,j)+(P(i,j)-P(i,j+1))*dx-APv(i,j)*v(i,j))
		Flux_v=Flux_v+DABS(APv(i,j)*v(i,j))
	end do
end do
error_u=Residu_u/Flux_u
error_v=Residu_v/Flux_v


RETURN
end







!********************solve_pressure _correction equation*********
!                    _______________
!solve the pressure_correction Eq.

subroutine PRESSURE_CORRECTION 
implicit none 
INTEGER ,PARAMETER :: M=500 , N=100
real(8) AE(m,n),AW(m,n),AN(m,n),AS(m,n),AP(m,n),bprim(m,n),su(m,n),sp(m,n)
real(8) AX(m),BX(m-1),CX(2:m),QX(m),AY(n),BY(n-1),CY(2:n),QY(n)

integer i,j,iconv
real(8) u(0:m,0:n+1),v(0:m+1,0:n),P(0:m+1,0:n+1),Pold(0:m+1,0:n+1),Pprim(0:m+1,0:n+1),u0(0:m,0:n+1),v0(0:m+1,0:n)
real(8) errorbp
real(8) APu(m-1,n),APv(m,n-1)
real(8) Flux,Sum
real(8) dx,dy,wp,wpp
real(8) i_u(4) , i_v(4) , i_T(4) , i_p(4)
REAL(8) x(0:m+1),y(0:n+1),xv(0:m+1),yv(0:n),xu(0:m),yu(0:n+1)

common / convective / iconv
common / MeshVariable / x,y,xv,yv,xu,yu
common / length /dx,dy
common / index_of_boundary / i_u,i_v,i_p,i_T
common / VelocityComponent / u,v,u0,v0
common / PressureComponent / P,Pprim
common / pressureVariable/APu,APv,wp,errorbp


wpp=0.2
pold=p

!COMPUTING COEFFICIENT OF PRESSURE-CORRECTION EQ.
do j=1,n
	do i=1,m-1
		AE(i,j)=((Y(J)-Y(J-1))**2)/APu(i,j)
	end do
end do
do j=1,n
	do i=2,m
		AW(i,j)=((Y(J)-Y(J-1))**2)/APu(i-1,j)
	end do
end do
do i=1,m
	do j=1,n-1
		AN(i,j)=((X(I)-X(I-1))**2)/APv(i,j)
	end do
end do
do i=1,m
	do j=2,n
		AS(i,j)=((X(I)-X(I-1))**2)/APv(i,j-1)
	end do
end do


su(:,:)=0.
sp(:,:)=0.

!SOLVE THE PRESSURE-CORRECTION EQ.
do j=1,n
	do i=1,m
		bprim(i,j)=u(i-1,j)*(Y(J)-Y(J-1))-u(i,j)*(Y(J)-Y(J-1))+v(i,j-1)*(X(I)-X(I-1))-v(i,j)*(X(I)-X(I-1))+su(i,j)
		AP(i,j)=AE(i,j)+AW(i,j)+AN(i,j)+AS(i,j)-sp(i,j)
	end do
end do
!SOLVE THE PRESSURE-CORRECTION EQ. WITH AADI METHOD

Pprim(:,:)=1e-16
AX(:)=0. ;BX(:)=0. ;CX(:)=0. ;QX(:)=0.
AY(:)=0. ;BY(:)=0. ;CY(:)=0. ;QY(:)=0.
!SWEEP X
do j=1,n
	do i=1,m
	    if (i_p(3)==0) then
        pprim(0,:)=pprim(1,:)
        end if 
        if (i_p(4)==0) then
        pprim(m+1,:)=pprim(m,:)
        end if
        if (i_p(1)==0) then
        pprim(:,n+1)=pprim(:,n)
        end if 
        if (i_p(2)==0) then
        pprim(:,0)=pprim(:,1)
        end if

        AX(i)=AP(i,j)/Wpp
		QX(i)=AN(i,j)*Pprim(i,j+1)+AS(i,j)*Pprim(i,j-1)+bprim(i,j)+(1.D0-Wpp)/Wpp*AP(i,j)*Pprim(i,j)
	end do
	QX(1)=QX(1)+AW(1,j)*Pprim(0,j)     !AW(1,:) already was set to zero so it is not necessary to change QX(1)!!
	QX(m)=QX(m)+AE(m,j)*Pprim(m+1,j)   !AE(m,j) already was set to zero so it is not necessary to change QX(m)!!
	do i=1,m-1
		BX(i)=-AE(i,j)
	end do
	do i=2,m
		CX(i)=-AW(i,j)
	end do
	CALL Crute_method(m,AX,BX,CX,QX)
	do i=1,m
		Pprim(i,j)=QX(i)
	end do
	AX(:)=0. ;BX(:)=0. ;CX(:)=0. ;QX(:)=0.
end do

!SWEEP y
do i=1,m
	do j=1,n
	    if (i_p(3)==0) then
        pprim(0,:)=pprim(1,:)
        end if 
        if (i_p(4)==0) then
        pprim(m+1,:)=pprim(m,:)
        end if
        if (i_p(1)==0) then
        pprim(:,n+1)=pprim(:,n)
        end if 
        if (i_p(2)==0) then
        pprim(:,0)=pprim(:,1)
        end if


        AY(j)=AP(i,j)/Wpp
		QY(j)=AE(i,j)*Pprim(i+1,j)+AW(i,j)*Pprim(i-1,j)+bprim(i,j)+(1.D0-Wpp)/Wpp*AP(i,j)*Pprim(i,j)
	end do
	QY(1)=QY(1)+AS(i,1)*Pprim(i,0)
	QY(n)=QY(n)+AN(i,n)*Pprim(i,n+1)
	do j=1,n-1
		BY(j)=-AN(i,j)
	end do
	do j=2,n
		CY(j)=-AS(i,j)
	end do
	CALL Crute_method(n,AY,BY,CY,QY)
	do j=1,n
		Pprim(i,j)=QY(j)
	end do
	AY(:)=0. ;BY(:)=0. ;CY(:)=0. ;QY(:)=0.
end do

errorbp=MAXVAL(DABS(bprim)) 

!for simple algorithm
p=(1-wp)*pold+wp*(p+pprim)

RETURN
end


!********************CorrectVelocity*********
!                    _______________
!!correct the velocity field
subroutine CorrectVelocity (APu,APv)
implicit none 
INTEGER ,PARAMETER :: M=500 , N=100

integer i,j
real(8) APu(m-1,n),APv(m,n-1)
real(8) u(0:m,0:n+1),v(0:m+1,0:n),u0(0:m,0:n+1),v0(0:m+1,0:n)
real(8) P(0:m+1,0:n+1),Pprim(0:m+1,0:n+1) 
real(8) dx,dy
real(8) i_u(4),i_v(4),i_p(4),i_T(4)
 REAL(8) x(0:m+1),y(0:n+1),xv(0:m+1),yv(0:n),xu(0:m),yu(0:n+1)

common / MeshVariable / x,y,xv,yv,xu,yu
common / length / dx,dy
common / index_of_boundary / i_u,i_v,i_p,i_T
common / VelocityComponent / u,v,u0,v0
common / PressureComponent / P,Pprim

!CORRECT THE U Velocity Fields(u) 
do j=1,n
	do i=1,m-1
		u(i,j)=u(i,j)+(YU(J)-YU(J-1))/APu(i,j)*(Pprim(i,j)-Pprim(i+1,j))
	end do
end do
!CORRECT THE V Velocity Fields(v)
do i=1,m
	do j=1,n-1
		v(i,j)=v(i,j)+(XV(I)-XV(I-1))/APv(i,j)*(Pprim(i,j)-Pprim(i,j+1))
	end do
end do

RETURN
end



!********************Crute_method *********
!                    _______________
!solve  the 3diagonal matrix by Compact_LU Crute method
subroutine Crute_method (n,a,b,c,d)
implicit none 
integer n,i,j,k
real(8) a(n),b(n-1),c(2:n),d(n),x(n)


	b(1)=b(1)/a(1)
	do i=2,n-1
		a(i)=a(i)-b(i-1)*c(i)
		b(i)=b(i)/a(i)	
	end do
	a(n)=a(n)-c(n)*b(n-1)

	x(1)=d(1)/a(1)
	do j=2,n
		x(j)=(d(j)-c(j)*x(j-1))/a(j)
	end do
	
	do k=n-1,1,-1
		x(k)=x(k)-b(k)*x(k+1)
	end do
	d(:)=x(:)

RETURN
end

!******************** subroutine coefficient energy**********
!                                __________________
SUBROUTINE coefficient_energy
IMPLICIT NONE
INTEGER ,PARAMETER :: M=500 , N=100
INTEGER i,j,iconv,i_unsteady
real(8),dimension(0:m+1,0:n+1) ::  DeT,DwT,DnT,DsT,FeT,FwT,FnT,FsT 
real(8) AET(m,n),AWT(m,n),ANT(m,n),AST(m,n),APT(m,n),bT(m,n)  
real(8) x(0:m+1),y(0:n+1),xv(0:m+1),yv(0:n),xu(0:m),yu(0:n+1)
real(8) u(0:m,0:n+1),v(0:m+1,0:n)
real(8) T(0:m+1,0:n+1),T0(0:m+1,0:n+1)             
real(8) AP0T(m,n),su_T(m,n),sp_T(m,n)
real(8) dx,dy,dt
real(8) GamaT,WT  
real(8) ALPHAeT, ALPHAwT, ALPHAnT, ALPHAsT
real(8) i_u(4),i_v(4),i_p(4),i_T(4)
real(8) pee,pew,pen,pes
real(8) sorc(0:m+1,0:n+1)


common /convective / iconv
common / unsteady / i_unsteady
common / length /dx,dy
common / index_of_boundary / i_u,i_v,i_p,i_T
common / time/dt
common / MeshVariable / x,y,xv,yv,xu,yu
common / VelocityComponent / u,v
common / coefficient / GamaT 
common / EnergyComponent / T,T0 
common / Energyvariable/ AET,AWT,ANT,AST,APT,bT,AP0T,su_T,sp_T,WT  

!COMPUTING COEFFICIENT OF Energy Eq.
do j=1,n
	do i=1,m
		
		FeT(i,j)=dy*u(i,j)
		FwT(i,j)=dy*u(i-1,j)
		FnT(i,j)=dx*v(i,j)
		FsT(i,j)=dx*v(i,j-1)
        
		DeT(i,j)=GamaT*(y(j)-y(j-1))/(x(i)-x(i-1))
        DwT(i,j)=GamaT*(y(j)-y(j-1))/(x(i)-x(i-1))
		DnT(i,j)=GamaT*(x(i)-x(i-1))/(y(j)-y(j-1))
		DsT(i,j)=GamaT*(x(i)-x(i-1))/(y(j)-y(j-1))
		
        if (iconv==1) then ! upwind

        AET(i,j)=det(i,j)+max(-fet(i,j),0.0)
		AWT(i,j)=dwt(i,j)+max(fwt(i,j),0.0)
		ANT(i,j)=dnt(i,j)+max(-fnt(i,j),0.0)
		AST(i,j)= dst(i,j)+max(fst(i,j),0.0)
		AP0T(i,j)=dx*dy/dt*i_unsteady        
        sorc(i,j)=fet(i,j)-fwt(i,j)+fnt(i,j)-fst(i,j)
            
        end if 
         if (iconv==2) then ! hybrid

        AET(i,j)=max(-feT(i,j),(deT(i,j)-(feT(i,j)/2)),0.0)
		AWT(i,j)=max(fwT(i,j),(dwT(i,j)+(fwT(i,j)/2)),0.0)
		ANT(i,j)=max(-fnT(i,j),(dnT(i,j)-(fnT(i,j)/2)),0.0)
		AST(i,j)= max(fsT(i,j),(dsT(i,j)+(fsT(i,j)/2)),0.0)
		AP0T(i,j)=dx*dy/dt*i_unsteady        
        sorc(i,j)=fet(i,j)-fwt(i,j)+fnt(i,j)-fst(i,j)    
         end if 
         
         if (iconv==3) then ! power_law

         pee=feT(i,j)/deT(i,j)
         pew=fwT(i,j)/dwT(i,j)
         pen=fnT(i,j)/dnT(i,j)
         pes=fsT(i,j)/dsT(i,j)
             
        AET(i,j)=det(i,j)*max(0.0,(1.0-0.1*pee)**5) + max(-feT(i,j),0.0)
		AWT(i,j)=dwt(i,j)*max(0.0,(1.0-0.1*pew)**5) + max(fwT(i,j),0.0)
		ANT(i,j)=dnt(i,j)*max(0.0,(1.0-0.1*pen)**5) + max(-fnT(i,j),0.0)
		AST(i,j)=dst(i,j)*max(0.0,(1.0-0.1*pes)**5) + max(fsT(i,j),0.0) 
		AP0T(i,j)=dx*dy/dt*i_unsteady        
        sorc(i,j)=fet(i,j)-fwt(i,j)+fnt(i,j)-fst(i,j)            
         end if 
        
        if (iconv==4) then     !quick
        
		 !QUICK_Hayase method
		 
        IF (FwT(i,j) > 0) THEN
            ALPHAwT = 1
          ELSE
            ALPHAwT = 0
        END IF
        IF (FeT(i,j) > 0) THEN
            ALPHAeT = 1
          ELSE
            ALPHAeT = 0
        END IF
        IF (FsT(i,j) > 0) THEN
            ALPHAsT = 1
          ELSE
            ALPHAsT = 0
        END IF
        IF (FnT(i,j) > 0) THEN
            ALPHAnT = 1
          ELSE
            ALPHAnT = 0
        END IF
    
		AET(i,j)=DeT(i,j)-(1-ALPHAeT)*FeT(i,j)
		AWT(i,j)=DwT(i,j)+ALPHAwT*FwT(i,j)
		ANT(i,j)=DnT(i,j)-(1-ALPHAnT)*FnT(i,j)
		AST(i,j)=DsT(i,j)+ALPHAsT*FsT(i,j)
		AP0T(i,j)=dx*dy/dt*i_unsteady        
        sorc(i,j)=fet(i,j)-fwt(i,j)+fnt(i,j)-fst(i,j)
        end if 
        
        
	end do
end do

!SOURCE TERM OF Energy-COMPONENT
su_T(:,:)=0. ;sp_T(:,:)=0.

do j=1,n
    do i=1,m
        if (iconv==4) then
            
            IF (FwT(i,j) > 0) THEN
                ALPHAwT = 1
            ELSE
                ALPHAwT = 0
            END IF
            IF (FeT(i,j) > 0) THEN
                ALPHAeT = 1
            ELSE
                ALPHAeT = 0
            END IF
            IF (FsT(i,j) > 0) THEN
                ALPHAsT = 1
            ELSE
                ALPHAsT = 0
            END IF
            IF (FnT(i,j) > 0) THEN
                ALPHAnT = 1
            ELSE
                ALPHAnT = 0
            END IF

            IF(i==1 .and. j==1) THEN
                bT(i,j)=sorc(i,j)+AP0T(i,j)*T0(i,j)+su_T(i,j)+ &
                        1/8*FwT(i,j)*ALPHAwT*(3*T(i,j)-2*T(i-1,j)-(-T(i,j)))+ &
                        1/8*FeT(i,j)*ALPHAeT*(T(i-1,j)+2*T(i,j)-3*T(i+1,j))+ &
                        1/8*FwT(i,j)*(1-ALPHAwT)*(3*T(i-1,j)-2*T(i,j)-T(i+1,j))+ &
                        1/8*FeT(i,j)*(1-ALPHAeT)*(2*T(i+1,j)+T(i+2,j)-3*T(i,j))+ &
                        1/8*FsT(i,j)*ALPHAsT*(3*T(i,j)-2*T(i,j-1)-T(i,j-2))+ &
                        1/8*FnT(i,j)*ALPHAnT*(T(i,j-1)+2*T(i,j)-3*T(i,j+1))
            ELSE
                bT(i,j)=sorc(i,j)+AP0T(i,j)*T0(i,j)+su_T(i,j)+ &
                        1/8*FwT(i,j)*ALPHAwT*(3*T(i,j)-2*T(i-1,j)-T(i-2,j))+ &
                        1/8*FeT(i,j)*ALPHAeT*(T(i-1,j)+2*T(i,j)-3*T(i+1,j))+ &
                        1/8*FwT(i,j)*(1-ALPHAwT)*(3*T(i-1,j)-2*T(i,j)-T(i+1,j))+ &
                        1/8*FeT(i,j)*(1-ALPHAeT)*(2*T(i+1,j)+(2*T(i+1,j)-T(i,j))-3*T(i,j))+ &
                        1/8*FsT(i,j)*ALPHAsT*(3*T(i,j)-2*T(i,j-1)-T(i,j-2))+ &
                        1/8*FnT(i,j)*ALPHAnT*(T(i,j-1)+2*T(i,j)-3*T(i,j+1))
            END IF
        end if
        
        if (iconv/=4) then
            bT(i,j)=sorc(i,j)+AP0T(i,j)*T0(i,j)
        end if

        APT(i,j)=AET(i,j)+AWT(i,j)+ANT(i,j)+AST(i,j)+AP0T(i,j)-sp_T(i,j)
    end do
end do

		
RETURN
END

!******************** subroutine solve_energy****************
!                     _______________________

subroutine solve_energy
implicit none 
INTEGER ,PARAMETER :: M=500 , N=100
real(8) AX(m),BX(m-1),CX(2:m),QX(m),AY(n),BY(n-1),CY(2:n),QY(n) 
integer i,j
real(8) x(0:m+1),y(0:n+1)
real(8) dx,dy,dt
real(8) u(0:m,0:n+1),v(0:m+1,0:n)
real(8) AET(m,n),AWT(m,n),ANT(m,n),AST(m,n),APT(m,n),bT(m,n)  
real(8) T(0:m+1,0:n+1),T0(0:m+1,0:n+1)             
real(8) WT,Flux_T,Residu_T ,error_T                          
real(8) DeT,DwT,DnT,DsT,FeT,FwT,FnT,FsT,GamaT         
real(8) AP0T(m,n),su_T(m,n),sp_T(m,n)  
real(8) i_u(4),i_v(4),i_p(4),i_t(4)

common / length /dx,dy
common / index_of_boundary / i_u,i_v,i_p,i_T
common / Time/ dt
common / MeshVariable / x,y
common / VelocityComponent / u,v   
common / coefficient / GamaT
common / Energyvariable/ AET,AWT,ANT,AST,APT,bT,AP0T,su_T,sp_T 
common / EnergyComponent / T,T0 
common / SolveEnergy/ WT,error_T    


 AX(:)=0. ;BX(:)=0. ;CX(:)=0. ;QX(:)=0.
 AY(:)=0. ;BY(:)=0. ;CY(:)=0. ;QY(:)=0.
 
  !SOLVE THE Energy Eq. WITH AADI METHOD
  !SWEEP Y
 do j=1,n
	do i=1,m
        if (i_T(3)==0) then
        T(0,:)=T(1,:)
        end if 
        if (i_T(4)==0) then
        T(m+1,:)=T(m,:)
        end if
        if (i_T(1)==0) then
        T(:,n+1)=T(:,n)
        end if 
        if (i_T(2)==0) then
        T(:,0)=T(:,1)
        end if

		AX(i)=APT(i,j)/WT
		QX(i)=ANT(i,j)*T(i,j+1)+AST(i,j)*T(i,j-1)+bT(i,j)+(1.D0-WT)/WT*APT(i,j)*T(i,j)   !R.H.S. of the Eq.
	end do
	QX(1)=QX(1)+AWT(1,j)*T(0,j)          
	QX(m)=QX(m)+AET(m,j)*T(m+1,j)  
	do i=1,m-1
		BX(i)=-AET(i,j)
	end do
	do i=2,m
		CX(i)=-AWT(i,j)
	end do
	CALL Crute_method(m,AX,BX,CX,QX)
	do i=1,m
		T(i,j)=QX(i)
	end do
	AX(:)=0. ;BX(:)=0. ;CX(:)=0. ;QX(:)=0.
 end do

 !SWEEP X
 do i=1,m
        if (i_T(3)==0) then
        T(0,:)=T(1,:)
        end if 
        if (i_T(4)==0) then
        T(m+1,:)=T(m,:)
        end if
        if (i_T(1)==0) then
        T(:,n+1)=T(:,n)
        end if 
        if (i_T(2)==0) then
        T(:,0)=T(:,1)
        end if
        
        do j=1,n
		AY(j)=APT(i,j)/WT
		QY(j)=AET(i,j)*T(i+1,j)+AWT(i,j)*T(i-1,j)+bT(i,j)+(1.D0-WT)/WT*APT(i,j)*T(i,j)
	end do
	QY(1)=QY(1)+AST(i,1)*T(i,0)      
	QY(n)=QY(n)+ANT(i,n)*T(i,n+1)
	do j=1,n-1
		BY(j)=-ANT(i,j)
	end do
	do j=2,n
		CY(j)=-AST(i,j)
	end do
	CALL Crute_method(n,AY,BY,CY,QY)
	do j=1,n
		T(i,j)=QY(j)
	end do
	AY(:)=0. ;BY(:)=0. ;CY(:)=0. ;QY(:)=0.
 end do

!One step.
Flux_T=0.D0 ;Residu_T=0.D0
do i=1,m
   do j=1,n
      Residu_T=Residu_T+DABS(ANT(i,j)*T(i,j+1)+AST(i,j)*T(i,j-1)+AET(i,j)*T(i+1,j)+AWT(i,j)*T(i-1,j)+bT(i,j)-APT(i,j)*T(i,j))
      Flux_T=Flux_T+DABS(APT(i,j)*T(i,j))
   end do
end do
error_T=Residu_T/Flux_T


RETURN
    END

!******************** subroutine result**********************
!                     _________________
subroutine      result(time,steady)
implicit none 
INTEGER ,PARAMETER :: M=500 , N=100

integer i,j,steady,anim
real(8) x(0:m+1),y(0:n+1),xv(0:m+1),yv(0:n),xu(0:m),yu(0:n+1)
real(8) up(0:m+1,0:n+1),vp(0:m+1,0:n+1),p(0:m+1,0:n+1)
real(8) u(0:m,0:n+1),v(0:m+1,0:n),u0(0:m,0:n+1),v0(0:m+1,0:n)
real(8) time
real(8) T(0:m+1,0:n+1)

common / unsteady / anim
common / MeshVariable / x,y,xv,yv,xu,yu
common / ResultVelocity/up,vp
common / VelocityComponent / u,v,u0,v0
common / PressureComponent / P
common /EnergyComponent/T

if (anim==1) then  !unsteady animation 

   open(1000,File='animation.plt',position='APPEND')                      
write(1000,*) 'VARIABLES=',"X",',',"Y",',',"u",',',"v",',',"p",',',"T"
write(1000,*) 'ZONE I=',m+2,' J=',n+2
  write(1000,*)"SOLUTIONTIME=",time
   up=0.0
   vp=0.0

   

       do i=1,m+1
          do j=1,n+1
		     up(i,j)=(u(i,j)+u(i-1,j))/2.0  
			 vp(i,j)=(v(i,j)+v(i,j-1))/2.0    
	      end do
       end do


       do i=0,m+1
          do j=0,n+1
		    write(1000,*)x(i),y(j),up(i,j),vp(i,j),p(i,j),T(i,j)
	

	      end do
       end do  
     close(1000)
     
     end if !unsteady animation 
     
        open(1001,File='final_contour.plt')

write(1001,*) 'VARIABLES=',"X",',',"Y",',',"u",',',"v",',',"p",',',"T"
write(1001,*) 'ZONE I=',m+2,' J=',n+2
   up=0.0
   vp=0.0

   

       do i=1,m+1
          do j=1,n+1
		     up(i,j)=(u(i,j)+u(i-1,j))/2.0  
			 vp(i,j)=(v(i,j)+v(i,j-1))/2.0    
	      end do
       end do


          do j=0,n+1
       do i=0,m+1
		    write(1001,*)x(i),y(j),up(i,j),vp(i,j),p(i,j),T(i,j)
	
	      end do
       end do  
     close(1001)
     open(1002,File='U @ end.plt')
     write(1002,*) 'VARIABLES=',"U",',',"Y"
     do j=0,n+1
     write(1002,*) up(m,j),y(j)
     end do
    close(1002)
    
     open(1002,File='T @ end.plt')
     write(1002,*) 'VARIABLES=',"T",',',"Y"
     do j=0,n+1
     write(1002,*) T(m,j),y(j)
     end do
    close(1002)
    
     return
    end