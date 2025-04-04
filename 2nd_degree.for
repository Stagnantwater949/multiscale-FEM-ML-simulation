      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, 
     &                STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, 
     &                CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS, COORDS, 
     &                DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, 
     &                KSPT, KSTEP, KINC)

      INCLUDE 'ABA_PARAM.INC'

	!step-0
      ! Declare variables
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS,NTENS)
      DIMENSION STRAN(NTENS), DSTRAN(NTENS), TIME(2), PREDEF(1)
      DIMENSION PROPS(NPROPS), COORDS(3), DROT(3,3), DFGRD0(3,3)
      DIMENSION DDSDDT(NTENS), DRPLDE(NTENS), SSE(1), SPD(1), SCD(1)
      DIMENSION DRPLDT(1), DPRED(1), DFGRD1(3,3),JSTEP(4)

      ! Material properties
      REAL*8 E, NU, SIGY, G, K, STR, DEQPL, PLASTIC_MULTIPLIER, YIELD_FUNC
      INTEGER U, J,P,O,I
      

      !implicit real(a-h o-z)
      
      dimension eelas(ntens), eplas (ntens), flow(ntens),olds (ntens), oldpl(ntens)

      parameter(toler=1.d-6,newton=20)

      ddsdde=0.d0
      E = props(1)
      xnu = props(2)
      Sy = props(3)
	  Xn=props(4)
	  SR=props(5)
	  GrS=props(6)
      !Young's modulus
      !Poisson's ratio
      !Yield stress
      !xn=props(4) strain hardening exponent
	  !SR=props(5) strain rate
	  
	  !Recover elastic and plastic strain and rotate forward
	  !recover equivalent plastic strain from statevariable
      call rotsig(statev(1), drot,eelas,2,ndi,nshr)
      call rotsig(statev(ntens+1),drot,eplas,2,ndi,nshr)
      eqplas=statev(1+2*ntens)
      
	  !assign previous stress and plastic strain tensor to old stress and old plastic strain
      olds=stress
      oldpl=eplas
      
    !Step-1: Elastic Predictor
		!Build elastic stiffness matrix
      eg=E/(1.d0+xnu)/2.d0
      elam=(E/(1.d0-2.d0*xnu)-2.d0*eg)/3.d0
      do i=1,3
          do j=1,3
              ddsdde(j,i)=elam
          end do
          ddsdde(i,i)=2.d0*eg+elam
      end do
      do i=4,ntens
          ddsdde(i,i)=eg
      end do
	  
		!calculate the predictor stress and elastic strain
      stress=stress+matmul (ddsdde,dstran)
      eelas=eelas+dstran
      
	  !Calculate equivalent vonMises stress
      Smises=(stress(1)-stress (2))**2+(stress (2)-stress (3))**2+(stress(3)-stress(1))**2
      do i=4,ntens
          Smises=Smises+6.d0*stress(i)**2
      end do
      Smises=sqrt(Smises/2.d0)
      
	  !Get yield stress from the specified hardening curve  
      !Sf=178.7707+(0.0000)+(847.7242*eqplas)+(568800.6245*SR)-(1.1593*GrS)-(1695.4791*eqplas**2)+(382357.1240*eqplas*SR)+(0.1552*(eqplas*GrS))-(506963633.4680*SR**2)-(109.2859*SR*GrS)+(0.0165*GrS**2)-(842.0720*eqplas**3)-(106326.2052*eqplas**2*SR)-(3.5993*eqplas**2*GrS)-(31480885.0535*eqplas*SR**2)+(2.4799*eqplas*SR*GrS)+(0.0002*eqplas*GrS**2)+(45189533554.1358*SR**3)+(7198.6076*SR**2*GrS)+(0.2116*SR*GrS**2)-(0.0001*GrS**3)
	  Sf = 194.7439 + (0.0000*1) + (95440.2556*eqplas) + (975510.7493*SR) - (0.6449*GrS) - (2227.2416*eqplas**2) + (545657.1438*eqplas*SR) + (0.2969*eqplas*GrS) - (8201273.1963*SR**2) - (12.3624*SR*GrS) + (0.0041*GrS**2)



		
	  !Determine if active yielding is there or not
      if (Smises.gt.(1.d0+toler)*Sf) then
	  !Calculate the flow direction
      Sh=(stress(1)+stress(2)+stress(3))/3.d0
      flow(1:3)=(stress(1:3)-Sh)/Smises
      flow(4:ntens)=stress(4:ntens)/Smises
	!STEP:2 COMPUTE INCREMENTAL EQUIVALENT PLASTIC STRAIN
		!Solve for Smises and deqpl using Newton's method
      deqpl=0.d0
      ET = 950.2556 - (2 * 2227.2416 * eqplas) + (54657.1438 * SR) + (0.2969 * GrS)

	        


      do kewton=1,newton 
          rhs=Smises-(3.d0*eg)*deqpl-Sf
          deqpl=deqpl+rhs/((3.d0*eg)+Et)
		  Sf = 194.7439 + (0.0000*1) + (95440.2556*eqplas) + (975510.7493*SR) - (0.6449*GrS) - (2227.2416*(eqplas+deqpl)**2) + (545657.1438*(eqplas+deqpl)*SR) + (0.2969*(eqplas+deqpl)*GrS) - (8201273.1963*SR**2) - (12.3624*SR*GrS) + (0.0041*GrS**2)

		  
		  !Sf=178.7707+(0.0000)+(847.7242*(eqplas+deqpl))+(568800.6245*SR)-(1.1593*GrS)-(1695.4791*(eqplas+deqpl)**2)+(382357.1240*(eqplas+deqpl)*SR)+(0.1552*(eqplas+deqpl)*GrS)-(506963633.4680*SR**2)-(109.2859*SR*GrS)+(0.0165*GrS**2)-(842.0720*(eqplas+deqpl)**3)-(106326.2052*(eqplas+deqpl)**2*SR)-(3.5993*(eqplas+deqpl)**2*GrS)-(31480885.0535*(eqplas+deqpl)*SR**2)+(2.4799*(eqplas+deqpl)*SR*GrS)+(0.0002*(eqplas+deqpl)*GrS**2)+(45189533554.1358*SR**3)+(7198.6076*SR**2*GrS)+(0.2116*SR*GrS**2)-(0.0001*GrS**3)
		  !Sf=10.2356+(320.3661*(eqplas+deqpl))+(164834.6640*SR)+(-274.1557*(eqplas+deqpl)**2)+(19843.3736*(eqplas+deqpl)*SR)+(-159123458.1398*(SR**2))
		  !Et=320.3661 + (-548.3114*(eqplas+deqpl)) + (19843.3736*SR)
		  ET = 95440.2556 - (2 * 2227.2416 * (eqplas+deqpl)) + (545657.1438 * SR) + (0.2969 * GrS)
		  !Et=847.7242+(-3390.9582*(eqplas+deqpl))+(382357.1240*SR)+(0.1552*GrS)+(-2526.2160*(eqplas+deqpl)**2)+(-212652.4104*(eqplas+deqpl)*SR)+(-7.1986*(eqplas+deqpl)*GrS)-31480885.0535*SR**2+(2.4799*SR*GrS)+(0.0002*GrS**2)


          if(abs(rhs).lt.toler*Sy) exit
      end do
      
      if (kewton.eq.newton) write(7,*)'WARNING: plasticity loop failed'
    !
	!STEP-3: UPDATE STRESS, ELASTIC STRAIN AND PLASTIC STRAIN TENSORS 
      ! update stresses and strains
      stress(1:3)=flow(1:3)*Sf+Sh
      eplas(1:3)=eplas(1:3)+3.d0/2.d0*flow(1:3)*deqpl
      eelas(1:3)=eelas(1:3)-3.d0/2.d0*flow(1:3)*deqpl
      stress(4:ntens)=flow(4:ntens)*Sf
      eplas(4:ntens)=eplas(4:ntens)+3.d0*flow(4:ntens)*deqpl
      eelas(4:ntens)=eelas(4:ntens)-3.d0*flow(4:ntens)*deqpl
      eqplas=eqplas+deqpl
      ! Calculate the plastic strain energy density
      do i=1,ntens
          spd=spd+(stress(i)+olds(i))*(eplas(i)-oldpl(i))/2.d0
      end do
    !STEP-4: COMPUTATION OF CONSISTENT MATERIAL JACOBIAN  
      ! Formulate the jacobian (material tangent)
      effg=eg*Sf/Smises
      efflam=(E/(1.d0-2.d0*xnu)-2.d0* effg)/3.d0
      effhrd=3.d0*eg*Et/(3.d0*eg+Et)-3.d0*effg
      do i=1,3
          do j=1,3
              ddsdde(j, i)=efflam
          end do
          ddsdde(i,i)=2.d0*effg+efflam
      end do
      do i=4,ntens
          ddsdde(i,i)=effg
      end do
      do i=1,ntens
          do j=1,ntens
            ddsdde(j,i)=ddsdde(j, i)+effhrd*flow(j)*flow(i)
          end do
      end do
      endif
      
	  !Store the strains in state variable array
      statev(1:ntens)=eelas
      statev((ntens+1):2*ntens)=eplas
      statev(1+2*ntens)=eqplas
      return
      end