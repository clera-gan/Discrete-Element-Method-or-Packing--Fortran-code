C1   ******************************************************************
C    ******************************************************************
C    *                 PROGRAM MAIN				      *
C    ******************************************************************
       INCLUDE  'common.for'
       DIMENSION
     & LFIRST(NN),LLAST(NN),LLIST(MAXIL),
     & LFIRST_O(NN),LLAST_O(NN),LLIST_O(MAXIL),
     & NPZON(NZSX,NZSY,NZSZ),IPZON(NZSX,NZSY,NZSZ,NZSP),
     & IZONX(NN),IZONY(NN),IZONZ(NN),icolor(nn),
     & fcontact(nn),ncontact(nn),NCN(nn),ivancontact(nn),
     & fcontactx(nn),fcontacty(nn),fcontactz(nn),fvanforcex(nn),
     & fvanforcey(nn),fvanforcez(nn),incontact(nn),
     & Fvanforce(nn),NXangpij(200),
     & NXang1(1000),NXang2(1000),NXang3(1000),NXangpp(200)
       DOUBLE PRECISION
     & QX(NN),QY(NN),QZ(NN),FRICP(MAXIL),FRICP_O(MAXIL),
     & ANGPX(NN),ANGVX(NN),RMOMX(NN),DISPTX(MAXIL),DISPTX_O(MAXIL),
     & ANGPY(NN),ANGVY(NN),RMOMY(NN),DISPTY(MAXIL),DISPTY_O(MAXIL),
     & ANGPZ(NN),ANGVZ(NN),RMOMZ(NN),DISPTZ(MAXIL),DISPTZ_O(MAXIL),
     & ANGVXB(NN),ANGVYB(NN),ANGVZB(NN),
     & DISPTWX(NN),DISPTWY(NN),DISPTWZ(NN),FRICW(NN),
     & RMASS(NN),XXI(NN),YYI(NN),ZZI(NN),
     & RRING(1000),XP0(1000),YP0(1000),yring(1000),
     & rpx(5000000),rpy(5000000),rpz(5000000)
       DOUBLE PRECISION
     & EMOD,VPOIS,DiamA,DiamB,DiamC,DENP,EMODW,VPOISW,
     & DiamA5,DiamB5,DiamC5,CUA,CUB,CUA2,CUB2,CUAW,ESTARP,ESTARW,
     & VEL0,UFRICP,UFRICW,DAMPNP,DAMPTP,DAMPNW,DAMPTW,CDAMP,
     & stiffness,DampCoeffN,DampCoeffT,UPP,UPW,
     & TwoTH,gg,DAMPNP2,DAMPTP2,DAMPNW2,DAMPTW2,
     & RNO,ANN2,DT2,DZONE,DRINGS,ZP0,DXY,ANG,
     & ANG1,ANG2,ANG3,Ea01,Ea02,Ea03,Eb01,Eb02,Eb03,Ec01,Ec02,Ec03,
     & COORDA,COORDB,RXI,RYI,RZI,X,Y,Z,RR,RxCont1,RyCont1,RzCont1,
     & RxCont2,RyCont2,RzCont2,RxCont,RyCont,RzCont,DSNORM,
     & A1,B1,C1,F1,G1,H1,P1,Q1,R1,D1,A2,B2,C2,F2,G2,H2,P2,Q2,R2,D2,
c------Normal and tangential vectors--
     & VctNXi,VctNYi,VctNZi,VCtNXj,VctNYj,VctNZj,VctNXYZi,VctNXYZj,
     & VctTXi,VctTYi,VctTZi,VCtTXj,VctTYj,VctTZj,
c-----Curvature-----------------------
     & curveradi,UPi,DNi,SINangLi,UPj,DNj,SINangLj,EccI2,WWWI,WWWJ,
     & RcurveIX,RcurveIY,RcurveI,RcurveJX,RcurveJY,RcurveJ, Rstar,EccJ2,
c----- contact---------------------- 
     & Fvander,cosgama,theta,fac1,fac2,fac3,fac4,Hamaker,
     & Closedis,Cutoffdis,
c-----Forces--------------------------
     & Rmav,HERTZC,FNORMT,UFNORMT,DSMAXFP,DSMAX,DTLIM,DISPM,
     & RCI,XNUI,YNUI,ZNUI,RXCONI,RYCONI,RZCONI,DXROT,DYROT,DZROT,
     & RCJ,XNUJ,YNUJ,ZNUJ,RXCONJ,RYCONJ,RZCONJ,DXCONI,DYCONI,DZCONI,
     & DXCONJ,DYCONJ,DZCONJ,DXCONIJ,DYCONIJ,DZCONIJ,DXTIJ,DYTIJ,DZTIJ,
     & DD,FFRICP,FFRICW,DISPT,FFRICX,FFRICY,FFRICZ,VREL,PKN,CN,FDAMPN,
     & FDAMPX,FDAMPY,FDAMPZ,PKT,CT,FDTX,FDTY,FDTZ,FDX,FDY,FDZ,FDN,
     & FWALL,DSMAXFW,UFWALL,XNU,YNU,ZNU,RYYP,RXXP,
c-----Torques
     & FDNX,FDNY,FDNZ,distNi,XpositNi,YpositNi,ZpositNi,XNUNi,YNUNi,
     & ZNUNi,distNj,XpositNj,YpositNj,ZpositNj,XNUNj,YNUNj,ZNUNj,
     & ffricxyz,distTi,XpositTi,YpositTi,ZpositTi,XNUTi,YNUTi,
     & ZNUTi,distTj,XpositTj,YpositTj,ZpositTj,XNUTj,YNUTj,ZNUTj,
     & DMOMiXN,DMOMiYN,DMOMiZN,DMOMjXN,DMOMjYN,DMOMjZN,
     & DMOMiXT,DMOMiYT,DMOMiZT,DMOMjXT,DMOMjYT,DMOMjZT,
c-----Update
     & ALPHA,ADAMP,BDAMP, RMOV2,ANGVXB0,ANGVYB0,ANGVZB0,
     & RMOMXB,RMOMYB,RMOMZB,YZI,coeffM,coeffP,angvxBt,angvzBt,
     & cbot,cmd1,cmd2,DTQ00,DTQ01,DTQ02,DTQ03,QQQ,
     & wxold,wyold,wzold,wxnew,wynew,wznew,Cridis,SAI0,SBI0,SCI0,
     & dffx,dfpressx,dfvisx,dffz,dfpressz,dfvisz,
     & ANGIJX1,ANGIJY1,ANGIJZ1,angijx,angijy,angijz,TAANG,TAANGW,rijx,
     & rijy,rijz,UPPO,UPWO,ROLLPX,ROLLPY,ROLLPZ
      character*20::str1,str2,str3
C    ******************************************************************
C    *	         read initial values to parameters	              *
C    ******************************************************************
       OPEN(UNIT=1,FILE='hopp3d.inp',STATUS='OLD')
       open(unit=6,file='view.dat')
C-----Define the Youngs modulus, major semi-axis, and density
       READ(1,*)
       READ(1,*)IforceMode,Iforcemodel
       READ(1,*)
       READ(1,*)EMOD,VPOIS,DENP
       READ(1,*)
       READ(1,*)apratio,DiamA,DiamB,DiamC      
      
       DIAM=DiamA
       if(apratio.lt.1.0)then
c       diamA*apratio=diamB*apratio=diamC
c       diamA=diamC/apratio
c       diamB=diamC/apratio
       DiamC=(diam*diam*diam*apratio*apratio)**(1.D0/3.D0)
       DiamA=DiamC/apratio
       DiamB=DiamC/apratio
       endif
       if(apratio.eq.1.0)then
       DiamC=diam*0.999999999D0
       DiamA=diam
       DiamB=diam
       endif
       if(apratio.gt.1.0)then
c       diamA=diamB*apratio=diamC*apratio
c       diamB=diamA/apratio
c       diamC=diamA/apratio
       DiamA=(diam*diam*diam*apratio*apratio)**(1.D0/3.D0)
       DiamB=DiamA/apratio
       DiamC=DiamA/apratio
       endif
c-------------------------------------
        EMODW=EMOD
        VPOISW=VPOIS
        DIAM=DiamA
        DiamA=DiamA/Diam
        DiamB=DiamB/Diam
        DiamC=DiamC/Diam
        DiamA5=diamA*0.5D0
        DiamB5=diamB*0.5D0
        DiamC5=diamC*0.5D0
c-------------------------------------
        EMOD=6.D0*EMOD/(3.14159265D0*DENP*9.81D0*DIAM)
        EMODW=6.D0*EMODW/(3.14159265D0*DENP*9.81D0*DIAM)
        CUA=1.0D0
        CUAW=1.0D0
        ESTARP=1.33333333D0*EMOD/(2.D0*(1.D0-VPOIS*VPOIS))
        ESTARW=1.33333333D0/((1.D0-VPOISW*VPOISW)/EMOD+
     *         (1.D0-VPOISW*VPOISW)/EMODW)
c-------------------------------------
       READ(1,*)
       READ(1,*)NTOT
       READ(1,*)
       READ(1,*)DT,TSTOP
       READ(1,*)
       READ(1,*)DMT1,DMT2,thick
       READ(1,*)
       READ(1,*)ZHT1,ZHT2
       READ(1,*)
       READ(1,*)VEL0,NRINGS
       READ(1,*)
       READ(1,*)UFRICP,UFRICW
       READ(1,*)
       READ(1,*)DAMPNP,DAMPTP,DAMPNW,DAMPTW,CDAMP
       READ(1,*)
       READ(1,*)stiffness,DampCoeffN,DampCoeffT,Hamaker
       READ(1,*)
       READ(1,*)UPP,UPW
       CLOSE(UNIT=1)
c---------------------
       REAL_HEIGHT=ZHT2
       REAL_D2=DMT2
       REAL_DT=DT*sqrt(DIAM/9.81D0)
       REAL_THICK=thick
c--------------------------------------------------------------
c       call gasflow to calculate the flow filed
c------------------------------------------------
C    ******************************************************************
C    *        Set Some Constants that will be used in the code        *
C    ******************************************************************
       If(diamA.eq.diamB.and.DiamB.eq.diamC)Ishape=0
       If(diamA.gt.diamB.and.DiamB.eq.diamC)Ishape=1
       If(diamA.eq.diamB.and.DiamB.gt.diamC)Ishape=2
       If(diamA.gt.diamB.and.DiamB.gt.diamC)Ishape=3
c---------------------------
       PI=ATAN(1.D0)*4.D0
       TWOTH=2.D0/3.D0
       realtime=1200.0
       itime=1
       gg=9.81D0
       xsmall=1.0D-30
       fac=(pi*denp*gg*diam*diam*diam)/6.  
       ALPHA=0.D0
       ADAMP=(2.D0-ALPHA*DT)/(2.D0+ALPHA*DT)
       BDAMP=2.D0/(2.D0+ALPHA*DT)
       DAMPNP2=2.D0*DAMPNP
       DAMPTP2=2.D0*DAMPTP
       DAMPNW2=2.D0*DAMPNW
       DAMPTW2=2.D0*DAMPTW
c      effective diameter and sphericity
       tops=(DiamA5*DiamB5*DiamC5)**0.666667
       palfa=1.6075
       bots=(DiamA5*DiamB5)**palfa+(DiamA5*DiamC5)**palfa
     &       +(DiamB5*DiamC5)**palfa
       sphety=tops/((bots/3.)**(1./palfa))
       Ediam=diam*(DiamA*DiamB*DiamC)**0.333333
       fctxz=13.0*sphety**4-33.31*sphety**3+34.22*sphety**2-
     &       19.14*sphety+6.235
C    ******************************************************************
C    *       Set the maxmium particle diameter 1, and the reduced     *
C    *       units were used                                          *
C    ******************************************************************
       DO I=1,NTOT
        RNO=random()
c---------------------
        DiaA(I)=diamA
        DiaB(I)=diamB
        DiaC(I)=diamC
        RadA(I)=diamA5
        RadB(I)=diamB5
        RadC(I)=diamC5
        Q00(I)=0.D0
        Q01(I)=0.D0
        Q02(I)=0.D0
        Q03(I)=0.D0
        RMASS(I)=DiaA(I)*DiaB(I)*DiaC(I)
        XXI(i)=0.2D0*RMass(i)*(RadB(i)*RadB(i)+RadC(i)*RadC(i))
        YYI(i)=0.2D0*RMass(i)*(RadA(i)*RadA(i)+RadC(i)*RadC(i))
        ZZI(i)=0.2D0*RMass(i)*(RadA(i)*RadA(i)+RadB(i)*RadB(i))
       ENDDO
c---------------------
        THICK=Real_thick/DIAM
        THICK5=THICK*0.5D0
        DMT1=DMT1/DIAM
        DMT2=DMT2/DIAM
        ZHT1=ZHT1/DIAM
        ZHT2=ZHT2/DIAM
c        UPP=UPP/DIAM
c        UPW=UPW/DIAM
        DMT15=DMT1*0.5D0
        DMT25=DMT2*0.5D0
        HGT1=ZHT1
        HGT2=HGT1+ZHT2
        HOPZT=HGT2
        Cutoffdis=1.0D-9/DIAM
       write(6,*)'Number of particles',NTOT
       write(6,*)'Three diameters',DiamA,DiamB,DiamC
       write(6,*)'Time step', DT,DMT1,DMT2
       write(6,*)'Glabal damping',CDAMP,HGT1,HGT2
C    ******************************************************************
C    *    definitions of cua2, cub, cub2, dt2, and the zones 	      *
C    *     in the whole calculation domain: NZONX,NZONY,NZONZ         *
C    ******************************************************************
       CUB=CUA*1.55D0
       CUA2=CUA*CUA
       CUB2=CUB*CUB
       ANN2=(0.5D0*(CUB-CUA))**2
       DT2=DT*DT
       DZONE=CUB
       NZONX=INT(dmt2/DZONE)+1
       NZONY=int(thick/dzone)+1
       NZONZ=INT(HOPZT/DZONE)+1
       IF(NZONY.LE.0.OR.NZONZ.LE.0) THEN
        WRITE(6,*) 'NZONY OR NZONZ <0'
        STOP ' NZONY OR NZONZ < 0 '
       ENDIF
C    ******************************************************************
C    *    set initial positions and velocity for each particle        *
C    ******************************************************************
       i=0
       ip=0
       n=0
       Nthickness=int(dmt2)-1
       Yring(1)=0.0D0
       DO II=1,Nthickness
       yring(2*II)=yring(2*II-2)+1.0D0
       yring(2*II+1)=yring(2*II-1)-1.0D0
       ENDDO
 
       nrings=int(dmt25)-1
       IF(dmt25.GT.1.99D0)THEN
       DRINGS=(dmt25-0.5D0)/DBLE(NRINGS)
       DO 8 IR=1,NRINGS
       RRING(IR)=(DBLE(IR)-0.5D0)*DRINGS
       rring(ir+nrings)=-rring(ir)
   8   CONTINUE

        DO 10 ij=1,NRINGS*2
        DO 10 IR=1,NRINGS*2
          IP=IP+1
          XP0(IP)=RRING(IR)
          YP0(IP)=RRING(ij)
  10    CONTINUE
          NPLAYER=IP
       ELSE
        NPLAYER=1
       ENDIF

       NLAYERS=INT(NTOT/NPLAYER)+1
       ZP0=HOPZT-5.0D0
       DO 11 ILAY=1,NLAYERS
       DO 11 IP=1,NPLAYER
         I=I+1
         IF(I.GT.NTOT)GOTO 11
           RX(I)=XP0(IP)
           RY(I)=YP0(IP)
           RZ(I)=ZP0
           ANGPX(I)=0.D0
           ANGPY(I)=0.D0
           ANGPZ(I)=0.D0
           RNO=RANDOM()
           ANG=PI*(RNO*.33333D0+.33333D0)
           DZ(I)=-VEL0*SIN(ANG)*DT
           DXY=-VEL0*COS(ANG)*DT
           RNO=RANDOM()
           ANG=2.D0*PI*RNO
           DX(I)=DXY*COS(ANG)
           DY(I)=DXY*SIN(ANG)
           ANGVX(I)=0.D0
           ANGVY(I)=0.D0
           ANGVZ(I)=0.D0
           ANGVXB(I)=0.D0
           ANGVYB(I)=0.D0
           ANGVZB(I)=0.D0
           RNO=1.0D0*Random()
           ANG1=pi*RNO/2.D0*1.0D0
           ANG2=pi*RNO/2.D0*1.0D0
           ANG3=pi*RNO/2.D0*1.0D0
           Q00(i)=cos(0.5D0*ANG2)*cos(0.5D0*ANG1+0.5D0*ANG3)
           Q01(i)=sin(0.5D0*ANG2)*cos(0.5D0*ANG1-0.5D0*ANG3)
           Q02(i)=sin(0.5D0*ANG2)*sin(0.5D0*ANG1-0.5D0*ANG3)
           Q03(i)=cos(0.5D0*ANG2)*sin(0.5D0*ANG1+0.5D0*ANG3)
  11   CONTINUE
C    ******************************************************************
C    *    set neighbour list initial value and other variable 	      *
C    *	  initial value                                               *
C    ******************************************************************
       DO 2003 I=1,NTOT
        lfirst(i)=0
        llast(i)=0
        disptwx(i)=0.D0
        disptwy(i)=0.D0
        disptwz(i)=0.D0
        fricw(i)=0.D0
        icolor(i)=i
2003   CONTINUE
       do k=1,maxil
        llist(k)=0
        disptx(k)=0.D0
        dispty(k)=0.D0
        disptz(k)=0.D0
        fricp(k)=0.D0
        enddo

       TIME=0.0
       NALL=0
       NLIST=0
       IALL=1
       it=0
       itt=0
       olap_max1=0.0
       olap_max2=0.0
       olap_min1=100.
       olap_min2=100.
       inext=0
       ippp=0
       Nparticledat=0
       XmaxRstar=0.
       XminRstar=1.
       height0=0.0D0
C**********************************************************************
C    *    IF the value of NRESTART is ZERO, then read data from the   *
C    *     file of preflow.dat                                        *
C**********************************************************************
       open(9,file='restart.dat',status='old')
       read(9,*)nrestart,nparticledat
       if(nrestart.eq.0)then
          open(2,file = 'preflow.dat', status = 'old')
          read(2,*)it,n,inext,ippp,itime
          if (n.ne.ntot) then
          write(6,*) 'n is unequal to ntot !'
c          stop  
          endif
          do i=1,n
          read (2,*)rx(i),ry(i),rz(i),dx(i),dy(i),dz(i),
     &         angpx(i),angpy(i),angpz(i),angvx(i),angvy(i),angvz(i),
     &         angvxB(i),angvyB(i),angvzB(i),velo_x(i),velo_z(i),
     &         orenxa(i),orenya(i),orenza(i)
	     hmaxrz=max(hmaxrz,rz(i)*1.0)
          enddo
          DO 3003 I=1,N
          read (2,*) lfirst(i),llast(i),disptwx(i),disptwy(i),
     &         disptwz(i),fricw(i),Q00(i),Q01(i),Q02(i),Q03(i),icolor(i)
3003      CONTINUE
		goto 110

         do k=1,maxil
            read(2,*)llist(k),disptx(k),dispty(k),disptz(k),fricp(k)
          enddo
          close(2)
          endif
       close(9)
       close(6)

c    ******************************************************************
C    ******************************************************************
C    *             Main Program Starts Here			      *
C    ******************************************************************
6     continue
       open(unit=6,file='view.dat',position='append')
       it=it+1
       TIME=DT*iT
       if(n.lt.ntot)then
       TNEWPAR=0.5D0
         if(time.gt.float(n)*tnewpar/float(nplayer))then
            n=n+nplayer
            if(n.gt.ntot)n=ntot
            iall=1
         endif
       endif
       DO 2008 I=1,N
        fx(i)=0.D0
        fy(i)=0.D0
        fz(i)=0.D0
        rmomx(i)=0.D0
        rmomy(i)=0.D0
        rmomz(i)=0.D0
        fcontact(i)=0.0
        fcontactx(i)=0.0
        fcontacty(i)=0.0
        fcontactz(i)=0.0
        Fvanforce(i)=0.0
        Fvanforcex(i)=0.0
        Fvanforcey(i)=0.0
        Fvanforcez(i)=0.0
        ncontact(i)=0
        ncn(i)=0
2008   CONTINUE
C    ******************************************************************
C    *             Renew Neighbour List  			      *
C    ******************************************************************
       IF (IALL.EQ.1) THEN
        NALL=NALL+1
        COORDA=0.D0
        COORDB=0.D0
        maxpzon=0
        DO 452 IX=1,NZONX
        DO 452 IY=1,NZONY
        DO 452 IZ=1,NZONZ
         NPZON(IX,IY,IZ)=0
  452   CONTINUE
        DO 453 I=1,N
         IZONX(I)=INT((RX(I)+dmt25) /DZONE)+1
         IZONY(I)=INT((RY(I)+thick5)/DZONE)+1
         IZONZ(I)=INT( RZ(I)        /DZONE)+1
         if(rz(i).le.0.0) izonz(i)=0
          IX=IZONX(I)
          IY=IZONY(I)
          IZ=IZONZ(I)
          IF(IZ.GT.NZONZ) THEN
           WRITE(6,*)I,IZ,NZONZ,rz(i),' IZ > NZONZ '
           STOP ' IZ > NZONZ '
          ENDIF
          IF(IZ.NE.0)THEN
           NPZON(IX,IY,IZ)=NPZON(IX,IY,IZ)+1
           IF(NPZON(IX,IY,IZ).GT.NZSP) THEN
           WRITE(6,*) 'ERROR npzon > nzsp'
           STOP 'ERROR npzon > nzsp'
           ENDIF
           IPZON(IX,IY,IZ,NPZON(IX,IY,IZ))=I
           maxpzon=max(maxpzon,npzon(ix,iy,iz))
          ENDIF
  453   CONTINUE
        DO 7532 I=1,N
         QX(I)=0.D0
         QY(I)=0.D0
         QZ(I)=0.D0
         LFIRST_O(I)=LFIRST(I)
         LLAST_O(I)=LLAST(I)
7532    CONTINUE
        DO 7533 L=1,MAXIL
         LLIST_O(L)=LLIST(L)
         DISPTX_O(L)=DISPTX(L)
         DISPTY_O(L)=DISPTY(L)
         DISPTZ_O(L)=DISPTZ(L)
         FRICP_O(L)=FRICP(L)
7533    CONTINUE
        L=1      
        DO 456 I=1,N-1
         LFIRST(I)=0
         LLAST(I)=0
         INX=IZONX(I)
         INY=IZONY(I)
         INZ=IZONZ(I)
         IF(INZ.NE.0)THEN
          LFIRST(I)=L
          RXI=RX(I)
          RYI=RY(I)
          RZI=RZ(I)
          INXA=INX-1
          IF(INXA.EQ.0)INXA=1
          INXB=INX+1
          IF(INXB.GT.NZONX)INXB=NZONX 
          INYA=INY-1
          IF(INYA.EQ.0)INYA=1
          INYB=INY+1
          IF(INYB.GT.NZONY)INYB=NZONY 
          INZA=INZ-1
          IF(INZA.EQ.0)INZA=1
          INZB=INZ+1
          IF(INZB.GT.NZONZ)INZB=NZONZ 
          DO 455 IX=1,NZONX
          DO 455 IY=1,NZONY
          DO 455 IZ=INZA,INZB
            NP=NPZON(IX,IY,IZ)
            IF(NP.EQ.0)GOTO 455
            DO 454 II=1,NP
              J=IPZON(IX,IY,IZ,II)
              IF(J.LE.I)GOTO 454
              X=RXI-RX(J)
              Y=RYI-RY(J)
              Z=RZI-RZ(J)
              IF(abs(y).gt.thick*0.5)then
                if(y.gt.0.)then
                y=y-thick
                goto 19
                endif
                if(y.lt.0.)then
                y=y+thick
                goto 19
                endif
              ENDIF
19            IF(abs(x).gt.dmt25)then
                if(x.gt.0.)then
                x=x-dmt2
                goto 18
                endif
                if(x.lt.0.)then
                x=x+dmt2
                goto 18
                endif
              ENDIF
18            RR=X*X+Y*Y+Z*Z
              IF (RR.LT.CUB2) THEN
                COORDB=COORDB+1.D0
                LLIST(L)=J
                DISPTX(L)=0.
                DISPTY(L)=0.
                DISPTZ(L)=0.
                FRICP(L)=0.
                IF(LFIRST_O(I).NE.0)THEN
                  DO 66 LO=LFIRST_O(I),LLAST_O(I)
                    JO=LLIST_O(LO)
                    IF(JO.EQ.J)THEN
                      DISPTX(L)=DISPTX_O(LO)
                      DISPTY(L)=DISPTY_O(LO)
                      DISPTZ(L)=DISPTZ_O(LO)
                      FRICP(L)=FRICP_O(LO)
                    ENDIF
  66              CONTINUE
                ENDIF
                L=L+1
                ENDIF
 454        CONTINUE
 455      CONTINUE
         ENDIF
         LLAST(I)=L-1
         IF(LLAST(I).LT.LFIRST(I))LFIRST(I)=0
456     CONTINUE
        NINL=L-1
        IF (NINL.GT.MAXIL) THEN
          WRITE(6,7644)NINL,MAXIL
7644      FORMAT('!!!! NUMBER IN LIST = ',I8,' MAX. IN LIST = ',I8)
          STOP ' ERROR NINL > MAXIL'
        ENDIF
        IALL=0
       ELSE
        NLIST=NLIST+1
       ENDIF
C    ******************************************************************
C    *             All Kinds of Forces Calculation  		      *
C    ******************************************************************
C    Calculating forces of Particle-Particle 
       MF=1
       DO 4756 I=1,N-1
        RXI=RX(I)
        RYI=RY(I)
        RZI=RZ(I)
        ML=LLAST(I)
C    There are no neighbours if MF>ML
        IF (MF.LE.ML) THEN
          DO 4757 K=MF,ML
            J=LLIST(K)
            X=RXI-RX(J)
            Y=RYI-RY(J)
            Z=RZI-RZ(J)
            RYYP=RY(J)
            RXXP=RX(J)
            IFVDW=0
            IF(abs(y).gt.thick*0.5)then
              if(y.gt.0.)then
              RYYP=RY(J)+thick
              goto 17
              endif
              if(y.lt.0.)then
              RYYP=RY(J)-thick
              goto 17
              endif
            ENDIF
17          IF(abs(x).gt.dmt25)then
              if(x.gt.0.)then
              RXXP=RX(J)+dmt2
              goto 16
              endif
              if(x.lt.0.)then
              RXXP=RX(J)-dmt2
              goto 16
              endif
            ENDIF
16         Y=RYI-RYYP
           X=RXI-RXXP
           RR=X*X+Y*Y+Z*Z
           IF(RR.ge.1.5)goto 4757
c------------------------------------------------------------
c     This is for spherical particles if Ishpae=0
c     This is for non-spherical particles if ishape=\0
c---------------------------------------
       IF(Ishape.eq.0)THEN
         DSNORM=RadA(i)+RadA(j)-sqrt(RR)
         RxCont=0.5D0*(Rx(i)+RXXP)
         RyCont=0.5D0*(Ry(i)+RYYP)
         RzCont=0.5D0*(Rz(i)+Rz(j))
         RxCont1=Rxcont
         RyCont1=Rycont
         RzCont1=Rzcont
         RxCont2=Rxcont
         RyCont2=Rycont
         RzCont2=Rzcont
         Icontact1=1
         Icontact2=1
         if(DSNORM.ge.0)then
         Icollision=1
         IFvdw=-1
         Closedis=cutoffdis
         goto 476
         else
         Icollision=0
         IFvdw=1
          IF(Closedis.lt.Cutoffdis)then
          Closedis=cutoffdis
          ENDIF
         goto 476 
         endif
         ELSE
c----------------------------
         CALL DETECT(RadA(i),RadB(i),RadC(i),RadA(j),RadB(j),RadC(j),
     &               RX(I),RY(I),RZ(I),Rxxp,Ryyp,RZ(J),
     &               Q00(I),Q01(I),Q02(I),Q03(I),
     &               Q00(J),Q01(J),Q02(J),Q03(J),
     &               RxCont1,RyCont1,RzCont1,Icontact1)
         IF(Icontact1.eq.0)goto 475
         CALL DETECT(RadA(j),RadB(j),RadC(j),RadA(i),RadB(i),RadC(i),
     &               Rxxp,Ryyp,RZ(J),RX(I),RY(I),RZ(I),
     &               Q00(J),Q01(J),Q02(J),Q03(J),
     &               Q00(I),Q01(I),Q02(I),Q03(I),
     &               RxCont2,RyCont2,RzCont2,Icontact2)
         IF(Icontact2.eq.0)goto 4757
         RxCont=0.5D0*(RxCont1+RxCont2)
         RyCont=0.5D0*(RyCont1+RyCont2)
         RzCont=0.5D0*(RzCont1+RzCont2)
         DSNORM=(RxCont1-RxCont2)**2+(RyCont1-RyCont2)**2+
     &           (RzCont1-RzCont2)**2+1.D-30
         DSNORM=SQRT(DSNORM) 
         IF(Icontact1.eq.1.and.Icontact2.eq.1)THEN
         Icollision=1
         IFvdw=-1
         Closedis=cutoffdis
         goto 476
         ENDIF
c----------------------------------------
c----------------------------------------
475     Icollision=0
        CriDis=0.5D0*Ediam/diam
        SAI0=RadA(i)+CriDis
        SBI0=RadB(i)+CriDis
        SCI0=RadC(i)+CriDis
            CALL DETECT(SAI0,SBI0,SCI0,RadA(j),RadB(j),RadC(j),
     &               RX(I),RY(I),RZ(I),Rxxp,Ryyp,RZ(J),
     &               Q00(I),Q01(I),Q02(I),Q03(I),
     &               Q00(J),Q01(J),Q02(J),Q03(J),
     &               RxCont1,RyCont1,RzCont1,Icontact1)
           IF(Icontact1.eq.0)goto 4757
            CALL DETECT(RadA(j),RadB(j),RadC(j),SAI0,SBI0,SCI0,
     &               Rxxp,Ryyp,RZ(J),RX(I),RY(I),RZ(I),
     &               Q00(J),Q01(J),Q02(J),Q03(J),
     &               Q00(I),Q01(I),Q02(I),Q03(I),
     &               RxCont2,RyCont2,RzCont2,Icontact2)
          IF(Icontact2.eq.0)goto 4757
          Icontact1=1
          Icontact2=1
         SlopXYZ=(RxCont1-RxCont2)**2+(RyCont1-RyCont2)**2+
     &           (RzCont1-RzCont2)**2+1.D-30
          SlopXYZ= sqrt(SlopXYZ)
c-------- unit vector, point 2 to 1' ,that is from point 1 to 2
          SlopX=(RxCont1-RxCont2)/SlopXYZ
          SlopY=(RyCont1-RyCont2)/SlopXYZ
          SlopZ=(RzCont1-RzCont2)/SlopXYZ
          Closedis=abs(Cridis-SlopXYZ)
          RxCont1=RxCont1-Cridis*SlopX
          RyCont1=RyCont1-Cridis*SlopY
          RzCont1=RzCont1-Cridis*SlopZ
          IFvdw=1
          IF(Closedis.lt.Cutoffdis)then
          Closedis=cutoffdis
          endif
          goto 476 
       ENDIF
c------------------------------------------------------------
c------------------------------------------------------------
c    If two particels are contacting
c-------------------------------------
476       IF(Icontact1.eq.1.and.Icontact2.eq.1)THEN
c----------------------------------------------------------------
c    Calculation of Normal and Tangential vectors and Curvature
c----------------------------------------------------------------
       IF(Ishape.eq.0)THEN
         VctNXi=(rx(i)-rxxp)/sqrt(RR)
         VctNYi=(ry(i)-ryyp)/sqrt(RR)
         VctNZi=(rz(i)-rz(j))/sqrt(RR)
         VctNXj=-VctNXi
         VctNYj=-VctNYi
         VctNZj=-VctNZi
         RcurveIX=0.5
         RcurveIY=0.5
         RcurveJX=0.5
         RcurveJY=0.5
         Rstar=0.25D0
       ELSE
       CALL PartPara(RadA(i),RadB(i),RadC(i),Q00(i),Q01(i),Q02(i),
     &               Q03(i),Rx(i),Ry(i),Rz(i),A1,B1,C1,F1,G1,H1,
     &               P1,Q1,R1,D1)
       CALL PartPara(RadA(j),RadB(j),RadC(j),Q00(j),Q01(j),Q02(j),
     &               Q03(j),Rxxp,Ryyp,Rz(j),A2,B2,C2,F2,G2,H2,
     &               P2,Q2,R2,D2)
C**********************************************************
       VctNXi=-2.D0*(A1*Rxcont1+G1*Rzcont1+H1*Rycont1+P1)
       VctNYi=-2.D0*(B1*Rycont1+F1*Rzcont1+H1*Rxcont1+Q1)
       VctNZi=-2.D0*(C1*Rzcont1+F1*Rycont1+G1*Rxcont1+R1)
       VctNXj=-2.D0*(A2*Rxcont2+G2*Rzcont2+H2*Rycont2+P2)
       VctNYj=-2.D0*(B2*Rycont2+F2*Rzcont2+H2*Rxcont2+Q2)
       VctNZj=-2.D0*(C2*Rzcont2+F2*Rycont2+G2*Rxcont2+R2)
       VctNXYZi=VctNXi*VctNXi+VctNYi*VctNYi+VctNZi*VctNZi
       VctNXYZj=VctNXj*VctNXj+VctNYj*VctNYj+VctNZj*VctNZj
       VctNXi=VctNXi/sqrt(VctNXYZi+xsmall)
       VctNYi=VctNYi/sqrt(VctNXYZi+xsmall)
       VctNZi=VctNZi/sqrt(VctNXYZi+xsmall)
       VctNXj=VctNXj/sqrt(VctNXYZj+xsmall)
       VctNYj=VctNYj/sqrt(VctNXYZj+xsmall)
       VctNZj=VctNZj/sqrt(VctNXYZj+xsmall)
c**********************************************************
c   Latitude angle for particle i
c        UPi=OrenXA(i)*VctNxi+OrenYA(i)*VctNYi+OrenZA(i)*VctNZi
c        DNi=OrenXA(i)*OrenXA(i)+OrenYA(i)*OrenYA(i)+OrenZA(i)*OrenZA(i)
c        SINangLi=abs(UPi)/sqrt(DNi+xsmall)
c        EccI2=1.D0-RadC(i)*RadC(i)/RadA(i)/RadA(i)
c        RcurveI=RadA(i)*sqrt(1.D0-EccI2)/(1.D0-EccI2*SINangLi*SINangLi)

c        UPj=OrenXA(j)*VctNxj+OrenYA(j)*VctNYj+OrenZA(j)*VctNZj
c        DNj=OrenXA(j)*OrenXA(j)+OrenYA(j)*OrenYA(j)+OrenZA(j)*OrenZA(j)
c        SINangLj=abs(UPj)/sqrt(DNj+xsmall)
c        EccJ2=1.D0-RadC(j)*RadC(j)/RadA(j)/RadA(j)
c        RcurveJ=RadA(j)*sqrt(1.D0-EccJ2)/(1.D0-EccJ2*SINangLj*SINangLj)
c--------------------------------------
c        RSTAR=RcurveI*RcurveJ/(RcurveI+RcurveJ)       
c--------------------------------------------------------------------       
        UPi=OrenXA(i)*VctNxi+OrenYA(i)*VctNYi+OrenZA(i)*VctNZi
        DNi=OrenXA(i)**2+OrenYA(i)**2+OrenZA(i)**2
        SINangLi=abs(UPi)/sqrt(DNi+xsmall)
        EccI2=1-RadC(i)**2/RadA(i)**2
        WWWI=sqrt(1-EccI2*SINangLi**2)
        RcurveIX=RadA(i)*(1-EccI2)/WWWI**3
        RcurveIY=RadA(i)/WWWI
        RcurveI=sqrt(RcurveIX*RcurveIY)

        UPj=OrenXA(j)*VctNxi+OrenYA(j)*VctNYi+OrenZA(j)*VctNZi
        DNj=OrenXA(j)**2+OrenYA(j)**2+OrenZA(j)**2
        SINangLj=abs(UPj)/sqrt(DNj+xsmall)
        EccJ2=1-RadC(j)**2/RadA(j)**2
        WWWJ=sqrt(1-EccJ2*SINangLj**2)
        RcurveJX=RadA(j)*(1-EccJ2)/WWWJ**3
        RcurveJY=RadA(j)/WWWJ
        RcurveJ=sqrt(RcurveJX*RcurveJY)
        RSTAR=RcurveI*RcurveJ/(RcurveI+RcurveJ)
c----------------------------------------------------------------
c        olap_max1=max(olap_max1,dsnorm)
c        olap_min1=min(olap_min1,dsnorm)
c        XmaxRstar=max(XmaxRstar,Rstar)
c        XminRstar=min(XminRstar,Rstar)
c---------------------------------------
        VctNXi=(VctNXi-VctNXj)*0.5D0
        VctNYi=(VctNYi-VctNYj)*0.5D0
        VctNZi=(VctNZi-VctNZj)*0.5D0
        VctNXj=-VctNXi
        VctNYj=-VctNYi
        VctNZj=-VctNZi
      ENDIF
c--------------------------------------------------------------
c   Van Waals forces--Everaers and Ejtehadi. 2003
c--------------------------------------------------------------
       IF(IFvdw.eq.1.or.IFvdw.eq.-1)then
            Hamaker=6.5D-20
            Fvander=0.0D0
         IF(ishape.eq.0)then
         cosgama=RadA(i)*diam*1.0
         ELSE
            RcurveIX=RcurveIX*diam
            RcurveIY=RcurveIY*diam
            RcurveJX=RcurveJX*diam
            RcurveJY=RcurveJY*diam

         UPi=OrenXA(i)*OrenXA(j)+OrenYA(i)*OrenYA(j)+OrenZA(i)*OrenZA(j)
         if(DNi.ne.0.0.and.DNj.ne.0.0)then
         theta=acos(abs(UPi)/sqrt(DNi*DNj+1.0D-30))
          cosgama=2.0/sqrt((1.0/RcurveIX-1.0/RcurveIY)
     &               *(1.0/RcurveJX-1.0/RcurveJY)*(sin(theta))**2.0+
     &  (1.0/RcurveIX+1.0/RcurveJX)*(1.0/RcurveIY+1.0/RcurveJY)+1.0D-30)
         else
         cosgama=2.0/sqrt((1.0/RcurveIX+1.0/RcurveJX)
     &         *(1.0/RcurveIY+1.0/RcurveJY)+1.0D-30)
         endif
         ENDIF
c------------attractive force --------
c----- direction: particle i from point 1 to point 2 --------------
c       IF(IFvdw.eq.1)then
         Hamaker=Hamaker/36.0
         if(ishape.eq.2.or.ishape.eq.0)then   
           fac1=RadA(i)*diam/(RadA(i)*diam+0.5*closedis*diam)
           fac2=RadC(i)*diam/(RadC(i)*diam+0.5*closedis*diam)
           fac3=2.0*(RadA(i)*diam)**4.0/(RadA(i)*diam+
     &          0.5*Closedis*diam)**5.0
           fac4=(RadC(i)*diam)**2.0/(RadC(i)*diam+
     &          0.5*closedis*diam)**3.0
          
           Fvander=Hamaker*(3*cosgama*(fac1**4*fac2**2)
     &          /(closedis*diam+1.0D-30)**2+(1+3.0*cosgama/
     & (closedis*diam+1.0D-30))*(fac3*fac2**2+fac1**4*fac4))/fac
        endif
        if(ishape.eq.1)then   
           fac1=RadC(i)*diam/(RadC(i)*diam+0.5*closedis*diam)
           fac2=RadA(i)*diam/(RadA(i)*diam+0.5*closedis*diam)
           fac3=2.0*(RadC(i)*diam)**4.0/(RadC(i)*diam+
     &         0.5*closedis*diam)**5.0
           fac4=(RadA(i)*diam)**2.0/(RadA(i)*diam+
     &         0.5*closedis*diam)**3.0
          
           Fvander=Hamaker*(3*cosgama*(fac1**4*fac2**2)
     &          /(closedis*diam+1.0D-30)**2+(1+3.0*cosgama/
     & (closedis*diam+1.0D-30))*(fac3*fac2**2+fac1**4*fac4))/fac
        endif
c---------- if collision, the direction of Fvdw is reverse with-----
c-----------normal contact force, so usw negative sign  -------------     
      if(IFvdw.eq.-1)Fvander=-Fvander 
 
          if(mod(it,50000).eq.0)then
           Fvanforcex(i)=Fvanforcex(i)+Fvander*VctNxi
           Fvanforcey(i)=Fvanforcey(i)+Fvander*VctNyi
           Fvanforcez(i)=Fvanforcez(i)+Fvander*VctNzi
           
           Fvanforcex(j)=Fvanforcex(j)-Fvander*VctNxi
           Fvanforcey(j)=Fvanforcey(j)-Fvander*VctNyi
           Fvanforcez(j)=Fvanforcez(j)-Fvander*VctNzi
c           ivancontact(i)=ivancontact(i)+1
c           ivancontact(j)=ivancontact(j)+1
          endif
c----------------------------------------------------          
           FX(I)=FX(I)+Fvander*VctNXi
           FY(I)=FY(I)+Fvander*VctNYi
           FZ(I)=FZ(I)+Fvander*VctNZi
           FX(J)=FX(J)-Fvander*VctNXi
           FY(J)=FY(J)-Fvander*VctNYi
           FZ(J)=FZ(J)-Fvander*VctNZi
c------------------------------------------------------
c  addition torque because the Fvdw is not necessary --
c--   pass the particle center  ----------------------
c-----------------------------------------------------
       if(Icollision.eq.0)then
       if(apratio.eq.1.0) goto 4757
c-------------------------------------------------
c       Normal toques
c-------------------------------------------------               
         RCI=sqrt((RxCont1-Rx(i))*(RxCont1-Rx(i))+
     &            (RyCont1-Ry(i))*(RyCont1-Ry(i))+
     &            (RzCont1-Rz(i))*(RzCont1-Rz(i)))
         XNUI=(Rx(i)-RxCont1)/RCI
         YNUI=(Ry(i)-RyCont1)/RCI
         ZNUI=(Rz(i)-RzCont1)/RCI
         RXCONI=-XNUI*RCI
         RYCONI=-YNUI*RCI
         RZCONI=-ZNUI*RCI
         RCJ=sqrt((RxCont2-Rxxp)*(RxCont2-Rxxp)+
     &            (RyCont2-Ryyp)*(RyCont2-Ryyp)+
     &            (RzCont2-Rz(j))*(RzCont2-Rz(j)))
         XNUJ=(Rxxp-RxCont2)/RCJ
         YNUJ=(Ryyp-RyCont2)/RCJ
         ZNUJ=(Rz(j)-RzCont2)/RCJ
         RXCONJ=-XNUJ*RCJ
         RYCONJ=-YNUJ*RCJ
         RZCONJ=-ZNUJ*RCJ

          FDNX=Fvander*VctNXi
          FDNY=Fvander*VctNYi
          FDNZ=Fvander*VctNZi
          distNi=-RxconI*VctNXi-RyconI*VctNYi-RzconI*VctNZi
          XpositNi=Rxcont1+distNi*VctNXi
          YpositNi=Rycont1+distNi*VctNYi
          ZpositNi=Rzcont1+distNi*VCtNZi
          XNUNi=XpositNi-Rx(i)
          YNUNi=YpositNi-Ry(i)
          ZNUNi=ZpositNi-Rz(i)
c--------------------------
          distNj=-RxconJ*VctNXj-RyconJ*VctNYj-RzconJ*VctNZj
          XpositNj=Rxcont2+distNj*VctNXj
          YpositNj=Rycont2+distNj*VctNYj
          ZpositNj=Rzcont2+distNj*VCtNZj
          XNUNj=XpositNj-Rxxp
          YNUNj=YpositNj-Ryyp
          ZNUNj=ZpositNj-Rz(j)
          CALL VPROD(XNUNi,YNUNi,ZNUNi,FDNX,FDNY,FDNZ,
     &               DMOMiXN,DMOMiYN,DMOMiZN)
          CALL VPROD(XNUNj,YNUNj,ZNUNj,-FDNX,-FDNY,-FDNZ,
     &               DMOMjXN,DMOMjYN,DMOMjZN)

          RMOMX(I)=RMOMX(I)+DMOMiXN
          RMOMY(I)=RMOMY(I)+DMOMiYN
          RMOMZ(I)=RMOMZ(I)+DMOMiZN

          RMOMX(j)=RMOMX(J)+DMOMjXN
          RMOMY(j)=RMOMY(J)+DMOMjYN
          RMOMZ(j)=RMOMZ(J)+DMOMjZN
          goto 4757      
        endif
      ENDIF

c-------------------------------------------------------
      RMAV=Rmass(I)*Rmass(J)/(Rmass(I)+Rmass(J))
C-------------------------------------------------------
C     Normal contact force: Hertz heory
         HERTZC=ESTARP*SQRT(RSTAR)
          If(IforceMode.eq.1)FNORMT=HERTZC*DSNORM**1.5D0
          If(IforceMode.eq.2)FNORMT=stiffness*RSTAR*DSNORM
               UFNORMT=UFRICP*FNORMT
               DSMAXFP=UFRICP*(2.D0-VPOIS)/(2.D0*(1.D0-VPOIS))
               DSMAX=DSNORM*DSMAXFP
               IF(FRICP(K).LT.UFNORMT)THEN
                  DTLIM=DSMAX*(1.D0-(1.D0-FRICP(K)/UFNORMT)**TWOTH)
                  CALL VMAG(DISPTX(K),DISPTY(K),DISPTZ(K),DISPM)
                  IF(DISPM.GT.DTLIM)
     &            CALL VSCALE(DISPTX(K),DISPTY(K),DISPTZ(K),DTLIM,DISPM)
                ENDIF

c-------------------------------------------------------------
                RCI=sqrt((RxCont-Rx(i))*(RxCont-Rx(i))+
     &                   (RyCont-Ry(i))*(RyCont-Ry(i))+
     &                   (RzCont-Rz(i))*(RzCont-Rz(i)))
                XNUI=(Rx(i)-RxCont)/RCI
                YNUI=(Ry(i)-RyCont)/RCI
                ZNUI=(Rz(i)-RzCont)/RCI
                RXCONI=-XNUI*RCI
                RYCONI=-YNUI*RCI
                RZCONI=-ZNUI*RCI
                CALL VPROD(ANGVX(I),ANGVY(I),ANGVZ(I),RXCONI,RYCONI,
     &                     RZCONI,DXROT,DYROT,DZROT)
                DXCONI=DX(I)+DXROT
                DYCONI=DY(I)+DYROT
                DZCONI=DZ(I)+DZROT
c---------------------------------------------------------
                RCJ=sqrt((RxCont-Rxxp )*(RxCont-Rxxp )+
     &                   (RyCont-Ryyp )*(RyCont-Ryyp )+
     &                   (RzCont-Rz(j))*(RzCont-Rz(j)))
                XNUJ=(Rxxp -RxCont)/RCJ
                YNUJ=(Ryyp -RyCont)/RCJ
                ZNUJ=(Rz(j)-RzCont)/RCJ
                RXCONJ=-XNUJ*RCJ
                RYCONJ=-YNUJ*RCJ
                RZCONJ=-ZNUJ*RCJ
                CALL VPROD(ANGVX(J),ANGVY(J),ANGVZ(J),
     &                 RXCONJ,RYCONJ,RZCONJ,DXROT,DYROT,DZROT)
                DXCONJ=DX(J)+DXROT
                DYCONJ=DY(J)+DYROT
                DZCONJ=DZ(J)+DZROT

                DXCONIJ=DXCONI-DXCONJ
                DYCONIJ=DYCONI-DYCONJ
                DZCONIJ=DZCONI-DZCONJ

              CALL VPROJEC(DXCONIJ,DYCONIJ,DZCONIJ,VctNXi,VctNYi,VctNZi,
     &                     DXTIJ,DYTIJ,DZTIJ,0)
              CALL VPROJEC(DISPTX(K),DISPTY(K),DISPTZ(K),VctNXi,VctNYi,
     &                     VctNZi,X,Y,Z,1)
                X=X+DXTIJ
                Y=Y+DYTIJ
                Z=Z+DZTIJ
                DD=SQRT(X*X+Y*Y+Z*Z)
                IF(DD.GE.DSMAX)THEN
                 CALL VSCALE(X,Y,Z,DSMAX,DD)
                  FFRICP=UFNORMT
                  DISPT=DSMAX
                ELSE
                  FFRICP=UFNORMT*(1.D0-(1.D0-DD/DSMAX)**1.5D0)
                  DISPT=DD
                ENDIF
                DISPTX(K)=X
                DISPTY(K)=Y
                DISPTZ(K)=Z
c************************************
                IF(DISPT.LE.0.000001)THEN
                  FFRICX=0.D0
                  FFRICY=0.D0
                  FFRICZ=0.D0
                ELSE
                  FFRICX=-FFRICP*DISPTX(K)/DISPT
                  FFRICY=-FFRICP*DISPTY(K)/DISPT
                  FFRICZ=-FFRICP*DISPTZ(K)/DISPT
                ENDIF
                FRICP(K)=FFRICP
C    ----------------------------------------------------------------------------
c      Normal damping: normal damping coefficient(CN)*(VREL) in normal direction
C    ---------------------------------------------------------------------------
                VREL=(DXCONIJ*VctNXi+DYCONIJ*VctNYi+DZCONIJ*VctNZi)/DT
                PKN=1.5D0*HERTZC*sqrt(DSNORM)
                CN=DAMPNP2*SQRT(PKN*RMAV)
                IF(IforceMode.eq.1)FDAMPN=-CN*VREL
                IF(IforceMode.eq.2)FDAMPN=-RMAV*VREL*DampCoeffN
C    ----------------------------------------------------------
C      Tangential damping in tangential direction
C    ----------------------------------------------------------
                IF(DD.GE.DSMAX)THEN
                  FDAMPX=0.D0
                  FDAMPY=0.D0
                  FDAMPZ=0.D0
                ELSE
                  PKT=1.5D0*UFNORMT*SQRT(1.D0-DD/DSMAX)/DSMAX
                  IF(IforceMode.eq.1)CT=DAMPTP2*SQRT(PKT*RMAV)
                  IF(IforceMode.eq.2)CT=RMAV*DampCoeffT
                  FDAMPX=-CT*DXTIJ/DT
                  FDAMPY=-CT*DYTIJ/DT
                  FDAMPZ=-CT*DZTIJ/DT
                ENDIF
C    ----------------------------------------------------------
C      Sum all forces in each direction (x, y, z)
C    ----------------------------------------------------------
                FDTX=FFRICX+FDAMPX
                FDTY=FFRICY+FDAMPY
                FDTZ=FFRICZ+FDAMPZ
                FDN=FNORMT+FDAMPN
                FDX=FDTX+FDN*VctNXi
                FDY=FDTY+FDN*VctNYi
                FDZ=FDTZ+FDN*VctNZi
                FX(I)=FX(I)+FDX
                FY(I)=FY(I)+FDY
                FZ(I)=FZ(I)+FDZ
                FX(J)=FX(J)-FDX
                FY(J)=FY(J)-FDY
                FZ(J)=FZ(J)-FDZ 
              fcontactx(i)=fcontactx(i)+FDX
              fcontacty(i)=fcontacty(i)+FDY
              fcontactz(i)=fcontactz(i)+FDZ
              
              fcontactx(j)=fcontactx(j)-FDX
              fcontacty(j)=fcontacty(j)-FDY
              fcontactz(j)=fcontactz(j)-FDZ
              incontact(i)=incontact(i)+1
              incontact(j)=incontact(j)+1
c-----------------------------------------------------------------
C                     ROLLING FRICTION Particle-Particle
c-----------------------------------------------------------------
                ANGIJX1=(ANGVX(I)-ANGVX(J))/dt
                ANGIJY1=(ANGVY(I)-ANGVY(J))/dt
                ANGIJZ1=(ANGVZ(I)-ANGVZ(J))/dt
                XNU=XNUi
                YNU=YNUi
                ZNU=ZNUi
                angijx=angijx1*(znu*znu+ynu*ynu)-angijy1*xnu*ynu
     &           -angijz1*xnu*znu
                angijy=angijy1*(znu*znu+xnu*xnu)-angijx1*xnu*ynu
     &           -angijz1*ynu*znu
                angijz=angijz1*(1-znu*znu)-angijx1*xnu*znu
     &           -angijy1*ynu*znu
                TAANG=SQRT(ANGIJX*ANGIJX+ANGIJY*ANGIJY+
     &              ANGIJZ*ANGIJZ)
                if(taang.eq.0.0)then
                  rijx=0.0
                  rijy=0.0
                  rijz=0.0
                else
                  RIJX=-ANGIJX/TAANG
                  RIJY=-ANGIJY/TAANG
                  RIJZ=-ANGIJZ/TAANG
                endif
                UPPO=UPP*FDN
                ROLLPX=UPPO*RIJX
                ROLLPY=UPPO*RIJY
                ROLLPZ=UPPO*RIJZ
c----------------------------------------------------------
c       Normal toques
c-----------------------------
          FDNX=(FDN+Fvander)*VctNXi
          FDNY=(FDN+Fvander)*VctNYi
          FDNZ=(FDN+Fvander)*VctNZi
          distNi=-RxconI*VctNXi-RyconI*VctNYi-RzconI*VctNZi
          XpositNi=Rxcont+distNi*VctNXi
          YpositNi=Rycont+distNi*VctNYi
          ZpositNi=Rzcont+distNi*VCtNZi
          XNUNi=XpositNi-Rx(i)
          YNUNi=YpositNi-Ry(i)
          ZNUNi=ZpositNi-Rz(i)
c--------------------------
          distNj=-RxconJ*VctNXj-RyconJ*VctNYj-RzconJ*VctNZj
          XpositNj=Rxcont+distNj*VctNXj
          YpositNj=Rycont+distNj*VctNYj
          ZpositNj=Rzcont+distNj*VCtNZj
          XNUNj=XpositNj-Rxxp
          YNUNj=YpositNj-Ryyp
          ZNUNj=ZpositNj-Rz(j)
C---------------------------------------------------------
c     Tangential toques
c------------------------
          ffricxyz=ffricx*ffricx+ffricy*ffricy+ffricz*ffricz
          VctTXi=ffricx/sqrt(ffricxyz+xsmall)
          VctTYi=ffricy/sqrt(ffricxyz+xsmall)
          VctTZi=ffricz/sqrt(ffricxyz+xsmall)
          distTi=-RxconI*VctTXi-RyconI*VctTYi-RzconI*VctTZi
          XpositTi=Rxcont+distTi*VctTXi
          YpositTi=Rycont+distTi*VctTYi
          ZpositTi=Rzcont+distTi*VctTZi
          XNUTi=XpositTi-Rx(i)
          YNUTi=YpositTi-Ry(i)
          ZNUTi=ZpositTi-Rz(i)
c------------------------------------------
          VctTXj=-VctTXi
          VctTYj=-VctTYi
          VctTZj=-VctTZi
          distTj=-RxconJ*VctTXj-RyconJ*VctTYj-RzconJ*VctTZj
          XpositTj=Rxcont+distTj*VctTXj
          YpositTj=Rycont+distTj*VctTYj
          ZpositTj=Rzcont+distTj*VctTZj
          XNUTj=XpositTj-Rxxp
          YNUTj=YpositTj-Ryyp
          ZNUTj=ZpositTj-Rz(j)
c-------------------------------
          CALL VPROD(XNUNi,YNUNi,ZNUNi,FDNX,FDNY,FDNZ,
     &               DMOMiXN,DMOMiYN,DMOMiZN)
          CALL VPROD(XNUNj,YNUNj,ZNUNj,-FDNX,-FDNY,-FDNZ,
     &               DMOMjXN,DMOMjYN,DMOMjZN)
          CALL VPROD(XNUTi,YNUTi,ZNUTi,FDTX,FDTY,FDTZ,
     &               DMOMiXT,DMOMiYT,DMOMiZT)
          CALL VPROD(XNUTj,YNUTj,ZNUTj,-FDTX,-FDTY,-FDTZ,
     &               DMOMjXT,DMOMjYT,DMOMjZT)
C-------------------------------------------------------------
c     TOtal torque
c---------------------------
          RMOMX(I)=RMOMX(I)+DMOMiXT+DMOMiXN+ROLLPX
          RMOMY(I)=RMOMY(I)+DMOMiYT+DMOMiYN+ROLLPY
          RMOMZ(I)=RMOMZ(I)+DMOMiZT+DMOMiZN+ROLLPZ

          RMOMX(J)=RMOMX(J)+DMOMjXT+DMOMjXN-ROLLPX
          RMOMY(J)=RMOMY(J)+DMOMjYT+DMOMjYN-ROLLPY
          RMOMZ(J)=RMOMZ(J)+DMOMjZT+DMOMjZN-ROLLPZ
C    ---------------------------------------------------------------
C    ---------------------------------------------------------------
             ELSE
C    -------------------------------------------------------------
C      No contact 
C    -------------------------------------------------------------
               DISPTX(K)=0.
               DISPTY(K)=0.
               DISPTZ(K)=0.
               FRICP(K)=0.
           ENDIF 
4757     CONTINUE
         MF=ML+1
        ENDIF
4756   CONTINUE
c***************************************************************************
c--------------------------------------------------------------------------
c     This part is to treat particle-wall interaction
C     Calculating forces of Particle-Walls & GRAVITY
c--------------------------------------------------------------------------
       COORDA=COORDA/FLOAT(N)
       COORDB=COORDB/FLOAT(N)
       DO 2015 I=1,N
c--------------------------------------------------------------------------
       CALL HOPP_HIT(I,Estarw,RadA(i),RadB(i),RadC(i),RX(I),RY(I),RZ(I),
     &                Q00(i),Q01(i),Q02(i),Q03(i),IWHIT,
     &                RxCont,RyCont,RzCont,CurveRadi)
c----------------------------------
        if(iwhit.eq.5)icolor(i)=5
        IF(IWHIT.NE.0)THEN
c-----------------------------------------
        DSNORM=disth(0)
        VctNXi=disth(1)
        VctNYi=disth(2)
        VctNZi=disth(3)
        Rstar=CurveRadi
        Rmav=Rmass(i)
c------------------------------------------------------------
c------------------------------------------------------------
        olap_max2=max(olap_max2,dsnorm*1.0)
        olap_min2=min(olap_min2,dsnorm*1.0)
C    --------------------------------------------------------------
C      Wall normal compression: Direction (Maximum -> Contact Point 
C    --------------------------------------------------------------
        HERTZC=ESTARW*SQRT(RSTAR)
        IF(IforceMode.eq.1) FWALL=HERTZC*DSNORM**1.5D0
        IF(IforceMode.eq.2) FWALL=Stiffness*RSTAR*DSNORM

        fcontact(i)=fcontact(i)+FWALL*1.0
        ncontact(i)=ncontact(i)+1
C    ----------------------------------------------------------
C       Wall friction: Direction (NORMAL to FWALL)
C    ----------------------------------------------------------
          DSMAXFW=UFRICW*(2.D0-VPOISW)/(2.D0*(1.D0-VPOISW))
          UFWALL=UFRICW*FWALL
          DSMAX=DSNORM*DSMAXFW
          IF(FRICW(I).LT.UFWALL)THEN
            DTLIM=DSMAX*(1.D0-(1.D0-FRICW(I)/UFWALL)**TWOTH)
            CALL VMAG(DISPTWX(I),DISPTWY(I),DISPTWZ(I),DISPM)
            IF(DISPM.GT.DTLIM)
     &      CALL VSCALE(DISPTWX(I),DISPTWY(I),DISPTWZ(I),DTLIM,DISPM)
          ENDIF
c*******************************************************************
          RCI=sqrt((RxCont-Rx(i))*(RxCont-Rx(i))+
     &             (RyCont-Ry(i))*(RyCont-Ry(i))+
     &             (RzCont-Rz(i))*(RzCont-Rz(i)))
          XNU=(Rx(i)-RxCont)/RCI
          YNU=(Ry(i)-RyCont)/RCI
          ZNU=(Rz(i)-RzCont)/RCI
          RXCONI=-XNU*RCI
          RYCONI=-YNU*RCI
          RZCONI=-ZNU*RCI
          CALL VPROD(ANGVX(I),ANGVY(I),ANGVZ(I),RXCONI,RYCONI,RZCONI,
     &               DXROT,DYROT,DZROT)
c     The velocity of contact point Vc(i)
          DXCONI=DX(I)+DXROT
          DYCONI=DY(I)+DYROT
          DZCONI=DZ(I)+DZROT
          CALL VPROJEC(DXCONI,DYCONI,DZCONI,VctNXi,VctNYi,VctNZi,
     &                 DXTIJ,DYTIJ,DZTIJ,0)
          CALL VPROJEC(DISPTWX(I),DISPTWY(I),DISPTWZ(I),VctNXi,VctNYi,
     &                 VctNZi,X,Y,Z,1)
          X=X+DXTIJ
          Y=Y+DYTIJ
          Z=Z+DZTIJ
          DD=SQRT(X*X+Y*Y+Z*Z)
          IF(DD.GE.DSMAX)THEN
            CALL VSCALE(X,Y,Z,DSMAX,DD)
            FFRICW=UFWALL
            DISPT=DSMAX
          ELSE
            FFRICW=UFWALL*(1.D0-(1.D0-DD/DSMAX)**1.5D0)
            DISPT=DD
          ENDIF
            DISPTWX(I)=X
            DISPTWY(I)=Y
            DISPTWZ(I)=Z
            IF(DISPT.LE.0.000001)THEN
              FFRICX=0.D0
              FFRICY=0.D0
              FFRICZ=0.D0
            ELSE
              FFRICX=-FFRICW*DISPTWX(I)/DISPT
              FFRICY=-FFRICW*DISPTWY(I)/DISPT
              FFRICZ=-FFRICW*DISPTWZ(I)/DISPT
            ENDIF
            FRICW(I)=FFRICW
C    ------------------------------------------------------------------
C      Wall  normal damping
C    ------------------------------------------------------------------
            VREL=(DX(I)*VctNXi+DY(I)*VctNYi+DZ(I)*VctNZi)/DT
c           VREL=(DXCONI*VctNXi+DYCONI*VctNYi+DZCONI*VctNZi)/DT
            PKN=1.5D0*HERTZC*sqrt(DSNORM)
            CN=DAMPNW2*SQRT(PKN*RMAV)
            If(IforceMode.eq.1)FDAMPN=-CN*VREL
            if(IforceMode.eq.2)FDAMPN=-RMAV*VREL*DampCoeffN
C    --------------------------------------------------------------------
C      Wall tangential damping
C    --------------------------------------------------------------------
            IF(DD.GE.DSMAX)THEN
               FDAMPX=0.D0
               FDAMPY=0.D0
               FDAMPZ=0.D0
            ELSE
               PKT=1.5D0*UFWALL*SQRT(1.D0-DD/DSMAX)/DSMAX
               If(IforceMode.eq.1)CT=DAMPTW2*SQRT(PKT*RMAV)
               If(IforceMode.eq.2)CT=RMAV*DampCoeffT
               FDAMPX=-CT*DXTIJ/DT
               FDAMPY=-CT*DYTIJ/DT
               FDAMPZ=-CT*DZTIJ/DT
            ENDIF
C----------------------------------------------------------------------
C       Sum all forces acted by WALL
C----------------------------------------------------------------------
            FDTX=FFRICX+FDAMPX
            FDTY=FFRICY+FDAMPY
            FDTZ=FFRICZ+FDAMPZ
            FDN=FWALL+FDAMPN
            FDX=FDTX+FDN*VctNXi
            FDY=FDTY+FDN*VctNYi
            FDZ=FDTZ+FDN*VctNZi
            FX(I)=FX(I)+FDX
            FY(I)=FY(I)+FDY
            FZ(I)=FZ(I)+FDZ
c-----------------------------------------------------------------
c                   rolling friction particle-wall
c-----------------------------------------------------------------
            angijx1=angvx(i)/dt
            angijy1=angvy(i)/dt
            angijz1=angvz(i)/dt
c            XNU=VctNXi
c            YNU=VctNYi
c            ZNU=VctNZi
            angijx=angijx1*(znu*znu+ynu*ynu)-angijy1*xnu*ynu
     $         -angijz1*xnu*znu
            angijy=angijy1*(znu*znu+xnu*xnu)-angijx1*xnu*ynu
     $         -angijz1*ynu*znu
            angijz=angijz1*(1-znu*znu)-angijx1*xnu*znu
     $         -angijy1*ynu*znu
            taangw=sqrt(angijx*angijx+angijy*angijy
     &                   +angijz*angijz)
            if(taangw.eq.0.0)then
               rijx=0.0
               rijy=0.0
               rijz=0.0
            else
               RIJX=-ANGIJX/TAANGW
               RIJY=-ANGIJY/TAANGW
               RIJZ=-ANGIJZ/TAANGW
            endif
            upwo=upw*FDN
            rollpx=upwo*rijx
            rollpy=upwo*rijy
            rollpz=upwo*rijz
c------------------------------
c       Normal toques
c-----------------------------
          FDNX=FDN*VctNXi
          FDNY=FDN*VctNYi
          FDNZ=FDN*VctNZi
          distNi=-RxconI*VctNXi-RyconI*VctNYi-RzconI*VctNZi
          XpositNi=Rxcont+distNi*VctNXi
          YpositNi=Rycont+distNi*VctNYi
          ZpositNi=Rzcont+distNi*VCtNZi
          XNUNi=XpositNi-Rx(i)
          YNUNi=YpositNi-Ry(i)
          ZNUNi=ZpositNi-Rz(i)
c--------------------------------------
c       Tangential toque
c---------------------------------------
          ffricxyz=ffricx*ffricx+ffricy*ffricy+ffricz*ffricz
          vctTXi=ffricx/sqrt(ffricxyz+xsmall)
          vctTYi=ffricy/sqrt(ffricxyz+xsmall)
          vctTZi=ffricz/sqrt(ffricxyz+xsmall)
          distTi=-RxconI*VctTXi-RyconI*VctTYi-RzconI*VctTZi
          XpositTi=Rxcont+distTi*VctTXi
          YpositTi=Rycont+distTi*VctTYi
          ZpositTi=Rzcont+distTi*VctTZi
          XNUTi=XpositTi-Rx(i)
          YNUTi=YpositTi-Ry(i)
          ZNUTi=ZpositTi-Rz(i)
cc--------------------------------------------
          CALL VPROD(XNUNi,YNUNi,ZNUNi,FDNX,FDNY,FDNZ,
     &               DMOMiXN,DMOMiYN,DMOMiZN)
          CALL VPROD(XNUTi,YNUTi,ZNUTi,FDTX,FDTY,FDTZ,
     &               DMOMiXT,DMOMiYT,DMOMiZT)
c-------------------------
c       Total toques
c--------------------------
          RMOMX(I)=RMOMX(I)+DMOMiXT+DMOMiXN+ROLLPX
          RMOMY(I)=RMOMY(I)+DMOMiYT+DMOMiYN+ROLLPY
          RMOMZ(I)=RMOMZ(I)+DMOMiZT+DMOMiZN+ROLLPZ
c---------------------------
          ELSE
            DISPTWX(I)=0.D0
            DISPTWY(I)=0.D0
            DISPTWZ(I)=0.D0
            FRICW(I)=0.D0
        ENDIF
C    -------------------------------------------------------
C    Gravity g=1.0
C    -------------------------------------------------------
          FZ(I)=FZ(I)-RMASS(I)*1.D0
2015   CONTINUE
C ******************************************************************
C *               UPDATE particle POSITIONS                        *
C ******************************************************************
       Energy=0.0
       Energy1=0.0
       Energy2=0.0
       Energy3=0.0
       inumber=0
       hmaxrz=0.0
       pvolume=0.0
c-----------------------------
       DO 2020 I=1,N
       if (rz(i).lt.0.0) goto 2020
c-------------------------------------------
c--------Translational motions -------------
c-------------------------------------------
        dx(i)=adamp*dx(i)+bdamp*fx(i)*DT2/rmass(i)
        dy(i)=adamp*dy(i)+bdamp*fy(i)*DT2/rmass(i)
        dz(i)=adamp*dz(i)+bdamp*fz(i)*DT2/rmass(i)
        velo_x(i)=dx(i)*1.0/dt
        velo_y(i)=dy(i)*1.0/dt
        velo_z(i)=dz(i)*1.0/dt
        RX(I)=RX(I)+DX(I)
        ry(i)=ry(i)+dy(i)
        rz(i)=rz(i)+dz(i)
        QX(I)=QX(I)+DX(I)
        QY(I)=QY(I)+DY(I)
        QZ(I)=QZ(I)+DZ(I)
        RMOV2=QX(I)*QX(I)+QY(I)*QY(I)+QZ(I)*QZ(I)
        IF(RMOV2.GT.ANN2) IALL=1
        hmaxrz=max(hmaxrz,rz(i)*1.0)
        pvolume=pvolume+pi*diamA*diamB*diamC/6.
c---------------------------------------------------------------------
c--------------------------------------------------------------------------
c    Calculation of angular velocity (it should be noted that the 
c    folloing algorithems are only applied to SPHEOID shaped particles
c--------------------------------------------------------------------------
      CALL FunMatrix(1,Q00(i),Q01(i),Q02(i),Q03(i),Ea01,Ea02,Ea03,
     &               Eb01,Eb02,Eb03,Ec01,Ec02,Ec03)
       RMOMXB=RMOMX(i)*Ea01+RMOMY(i)*Ea02+RMOMZ(i)*Ea03
       RMOMYB=RMOMX(i)*Eb01+RMOMY(i)*Eb02+RMOMZ(i)*Eb03
       RMOMZB=RMOMX(i)*Ec01+RMOMY(i)*Ec02+RMOMZ(i)*Ec03
c------------------------------------------------------------
       If(ishape.eq.0.or.ishape.eq.1)then
       angvxB0=angvxB(i)/dt
       angvyB0=angvyB(i)/dt
       angvzB0=angvzB(i)/dt
       YZI=ZZI(i)-XXI(i)
       coeffM=1.D0-cdamp*dt/2.D0
       coeffP=1.D0+cdamp*dt/2.D0
       angvxBt=(angvxB0+RMOMXB*dt/XXI(i)/2.D0)/coeffP
       cbot=ZZI(i)*ZZI(i)*coeffP*coeffP+(YZI*angvxBt*dt/2.D0)**2
       cmd1=RmomYB*dt+ZZI(i)*coeffM*angvyB0+YZI*angvxBt
     &      *dt*angvzB0/2.D0
       cmd2=RmomZB*dt-YZI*angvxBt*dt*angvyB0/2.D0+ZZI(i)*coeffM*angvzB0
       angvxB1=(angvxB0*coeffM+RMOMXB*dt/XXI(i))/coeffP
       angvyB1=( ZZI(i)*coeffP*cmd1+YZI*angvxBt*dt*cmd2/2.D0)/cbot
       angvzB1=(-YZI*angvxBt*dt*cmd1/2.D0+ZZI(i)*coeffP*cmd2)/cbot
       angvxB(i)=angvxB1*dt
       angvyB(i)=angvyB1*dt
       angvzB(i)=angvzB1*dt
      endif
c-----------
      if(ishape.eq.2)then
       angvxB0=angvxB(i)/dt
       angvyB0=angvyB(i)/dt
       angvzB0=angvzB(i)/dt
       YZI=YYI(i)-ZZI(i)
       coeffM=1.D0-cdamp*dt/2.D0
       coeffP=1.D0+cdamp*dt/2.D0
       angvzBt=(angvzB0+RMOMZB*dt/ZZI(i)/2.D0)/coeffP
       cbot=XXI(i)*YYI(i)*coeffP*coeffP+(YZI*angvzBt*dt/2.D0)**2
       cmd1=RmomXB*dt+XXI(i)*coeffM*angvxB0+YZI*angvzBt
     &      *dt*angvyB0/2.D0
       cmd2=RmomYB*dt-YZI*angvzBt*dt*angvzB0/2.D0+YYI(i)*coeffM*angvyB0
       angvxB1=( YYI(i)*coeffP*cmd1+YZI*angvzBt*dt*cmd2/2.D0)/cbot
       angvyB1=(-YZI*angvzBt*dt*cmd1/2.D0+XXI(i)*coeffP*cmd2)/cbot
       angvzB1=(angvzB0*coeffM+RMOMZB*dt/ZZI(i))/coeffP

       angvxB(i)=angvxB1*dt
       angvyB(i)=angvyB1*dt
       angvzB(i)=angvzB1*dt
      endif
c-----------
      if(ishape.eq.3)then
       wxold=angvxB(i)/dt
       wyold=angvyB(i)/dt
       wzold=angvzB(i)/dt
      CAll angwxyz(wxold,wyold,wzold,xxi(i),yyi(i),zzi(i),rmomxB,
     &     rmomyB,rmomzB,dt,cdamp,wxnew,wynew,wznew,ImaxN0)
       angvxB(i)=wxnew*dt
       angvyB(i)=wynew*dt
       angvzB(i)=wznew*dt
      endif
c-------------
      CALL FunMatrix(2,Q00(i),Q01(i),Q02(i),Q03(i),Ea01,Ea02,Ea03,
     &              Eb01,Eb02,Eb03,Ec01,Ec02,Ec03)
       ANGVX(I)=ANGVXB(i)*Ea01+ANGVYB(i)*Ea02+ANGVZB(i)*Ea03
       ANGVY(I)=ANGVXB(i)*Eb01+ANGVYB(i)*Eb02+ANGVZB(i)*Eb03
       ANGVZ(I)=ANGVXB(i)*Ec01+ANGVYB(i)*Ec02+ANGVZB(i)*Ec03
       ANGPX(I)=ANGPX(I)+ANGVX(I)
       ANGPY(I)=ANGPY(I)+ANGVY(I)
       ANGPZ(I)=ANGPZ(I)+ANGVZ(I)
c----------------------------------------------------------
c     Quaternion method
c----------------------------------------------------------
      DTQ00=-0.5D0*( Q01(i)*ANGVXB(i)+Q02(i)*ANGVYB(i)+Q03(i)*ANGVZB(i))
      DTQ01= 0.5D0*( Q00(i)*ANGVXB(i)-Q03(i)*ANGVYB(i)+Q02(i)*ANGVZB(i))
      DTQ02= 0.5D0*( Q03(i)*ANGVXB(i)+Q00(i)*ANGVYB(i)-Q01(i)*ANGVZB(i))
      DTQ03= 0.5D0*(-Q02(i)*ANGVXB(i)+Q01(i)*ANGVYB(i)+Q00(i)*ANGVZB(i))
      Q00(i)=Q00(i)+DTQ00
      Q01(i)=Q01(i)+DTQ01
      Q02(i)=Q02(i)+DTQ02
      Q03(i)=Q03(i)+DTQ03
      QQQ=Q00(I)*Q00(i)+Q01(i)*Q01(i)+Q02(i)*Q02(i)+Q03(i)*Q03(i)
      Q00(i)=Q00(i)/sqrt(QQQ)
      Q01(i)=Q01(i)/sqrt(QQQ)
      Q02(i)=Q02(i)/sqrt(QQQ)
      Q03(i)=Q03(i)/sqrt(QQQ)
c--------------------------------------
      CALL FunMatrix(1,Q00(i),Q01(i),Q02(i),Q03(i),Ea01,Ea02,Ea03,
     &               Eb01,Eb02,Eb03,Ec01,Ec02,Ec03)
       ANGVXB(I)=ANGVX(i)*Ea01+ANGVY(i)*Ea02+ANGVZ(i)*Ea03
       ANGVYB(I)=ANGVX(i)*Eb01+ANGVY(i)*Eb02+ANGVZ(i)*Eb03
       ANGVZB(I)=ANGVX(i)*Ec01+ANGVY(i)*Ec02+ANGVZ(i)*Ec03
c----------------------------------------------------------------
c-----------------------------------------
c------System Engery checking---------
       veloXYZ=dx(i)*dx(i)/dt/dt+dy(i)*dy(i)/dt/dt+dz(i)*dz(i)/dt/dt
       Energy01=0.5*rmass(i)*veloXYZ
       Energy02=0.5*XXI(i)*angvx(i)*angvx(i)/dt/dt+0.5*YYI(i)*angvy(i)*
     &          angvy(i)/dt/dt+0.5*ZZI(i)*angvz(i)*angvz(i)/dt/dt
       Energy03=rmass(i)*rz(i)
       Energy=Energy+Energy01+Energy02+Energy03
       Energy1=Energy1+Energy01
       Energy2=Energy2+Energy02
       Energy3=Energy3+Energy03
c---------------------------------
c---------------------------------
      CALL FunMatrix(2,Q00(i),Q01(i),Q02(i),Q03(i),Ea01,Ea02,Ea03,
     &               Eb01,Eb02,Eb03,Ec01,Ec02,Ec03)
        if(ishape.eq.0.or.ishape.eq.1.or.ishape.eq.3)then
        ORENXA(i)=Ea01*RadA(i)*1.D0
        ORENYA(i)=Eb01*RadA(i)*1.D0
        ORENZA(i)=Ec01*RadA(i)*1.D0
        endif
        if(ishape.eq.2)then
        ORENXA(i)=Ea03*RadC(i)*1.D0
        ORENYA(i)=Eb03*RadC(i)*1.D0
        ORENZA(i)=Ec03*Radc(i)*1.D0
        endif
c------------------------------------
         IF(abs(ry(i)).gt.thick*.5)then
          if(ry(i).gt.0.)then
          ry(i)=ry(i)-thick
          goto 15
          endif
          if(ry(i).lt.0.)then
          ry(i)=ry(i)+thick
          goto 15
          endif
         ENDIF
15       IF(abs(rx(i)).gt.dmt25)then
          if(rx(i).gt.0.)then
          rx(i)=rx(i)-dmt2
          goto 2020
          endif
          if(rx(i).lt.0.)then
          rx(i)=rx(i)+dmt2
          goto 2020
          endif
         ENDIF
2020   CONTINUE
c**********************************************************************
c----------------------------------------------
       bedpor=1.0-pvolume/(Dmt2*thick*(hmaxrz-height0))
       if(mod(it,10000).eq.0.or.it.eq.1)then
        write(6,*)it,'contacts=',ninl,'maxpzon=',maxpzon,n
        write(*,*)it,'contacts=',ninl,'maxpzon=',maxpzon,n
        open(13,file='energy.dat',position='append')
        tsd=it*real_dt
        write(13,989)it,tsd,energy,energy1,energy2,energy3,bedpor,
     &               sphety,height0
989     format((i8),8(g15.6))
       close(13)
      endif
c********************************************************************
C**********************************************************************
C         Store data of particle positons at different time           *
C**********************************************************************
1029  IF(mod(it,50000).eq.0.and.it.ge.1.and.it.le.80000000)then
      open(10,file='particle.dat',position='append')
      write(10,*)'zone'
      do i=1,n
      write(10,991)rx(i)*diam,ry(i)*diam,rz(i)*diam,diaA(i)*diam,
     &      icolor(i),i
991   format(4(g15.6),2(i8))
      enddo
      close(10)
      endif
c*************************************************************************
c     Calculation of Porosity, CN and Orientation Based on the Cut volume
c*************************************************************************
110      IF(mod(it,1).eq.0)THEN
c------------ 
c      DO 559 IJK=1,1
      IPorosity=1
      ICCNN=0
      IOren=0
      IRDF=0
c----------------------
c---------------------
       DO I=1,1000
       NXang1(i)=0
       NXang2(i)=0
       NXang3(i)=0
       ENDDO
c----------------------------------
       CriDis=0.005D0*Ediam/diam
       Pdiswall=0.0

       CNNtot=0.0
       ICNNtot=0.0

       PVolume=0.0
       CNNtot=0.0
       ICNNtot=0.0
       pvolumecut=0.0
       cutvolume=0.0
c---------------------
       DO 189 I=1,N
c************************************
c      Porosity & Packing density
c************************************
      IF(IPorosity.eq.1)then
c------------
       Nin=0
       Ncut=0
       NP05=125
       NP15=2376
       NP20=5245
       NP25=9689
       NP30=15833
       NP40=33975
       NP50=59778
       NNPP=NP15
      open(31,file='location-1-015.dat',status='old')
c    define the domain for porosity calculation
c     pdiswall=2.0
      ph0top=hmaxrz*1.0-20.0
      ph0bot=height0*1.0+0.5
      ph0left=-dmt25+0.5
      ph0right=dmt25-0.5
      ph0rear=thick5-0.5
      ph0front=-thick5+0.5
c      dmtdp=dmt25
c      rmxyz=sqrt(rx(i)**2+ry(i)**2+1.E-30)*1.0
      volumespace=(ph0right-ph0left)*(ph0rear-ph0front)*(ph0top-ph0bot)
c----------
      CALL PartPara(RadA(i),RadB(i),RadC(i),Q00(I),Q01(I),Q02(I),Q03(i),
     &              RX(i),RY(i),RZ(i),A1,B1,C1,F1,G1,H1,P1,Q1,R1,D1)
       do ik=1,NNPP
        read(31,*)rpx(ik),rpy(ik),rpz(ik)
        xsmal=rpx(ik)+rx(i)
        ysmal=rpy(ik)+ry(i)
        zsmal=rpz(ik)+rz(i)
        FUVW=A1*Xsmal*Xsmal+B1*Ysmal*Ysmal+C1*Zsmal*Zsmal+
     &     2.*F1*Ysmal*Zsmal+2.*G1*Zsmal*Xsmal+2.*H1*Xsmal*Ysmal+
     &     2.*P1*Xsmal+2.*Q1*Ysmal+2.*R1*Zsmal+D1
        IF(FUVW.lt.0)Nin=Nin+1
        IF(FUVW.lt.0.and.zsmal.ge.ph0bot.and.zsmal.le.ph0top.and.
     &     xsmal.le.ph0right.and.xsmal.ge.ph0left.and.
     &     ysmal.le.ph0rear.and.ysmal.ge.ph0front)Ncut=Ncut+1
       ENDDO
       Close(31)
       cutvolume=cutvolume+Ncut*1./Nin*pi*diamA*diamB*diamC/6.
      ENDIF
c**************************************
c************************************
c      Coordination Number (CN)
c************************************
      IF(ICCNN.eq.1)Then
c-----------------------
        ICNNtot=ICNNtot+1
        SAI0=RadA(i)+CriDis
        SBI0=RadB(i)+CriDis
        SCI0=RadC(i)+CriDis
         DO 899 J=1,N
           if(I.eq.J)goto 899
           PorDisIJ=sqrt((Rx(i)-Rx(j))**2+(Ry(i)-Ry(j))**2+
     &                   (Rz(i)-Rz(j))**2)
c---------- sphere ----------------------------
           if(ishape.eq.0)then
            if(pordisij.gt.RadA(i)+RadA(j)+Cridis)then
              goto 899
            else
              CNNTOT=CNNtot+1.0
              NCN(i)=NCN(i)+1
              goto 899
            endif
           else
c----------------------------------------------
            if(PorDisIJ.gt.1.5)goto 899
            CALL DETECT(SAI0,SBI0,SCI0,RadA(j),RadB(j),RadC(j),
     &               RX(I),RY(I),RZ(I),RX(J),RY(J),RZ(J),
     &               Q00(I),Q01(I),Q02(I),Q03(I),
     &               Q00(J),Q01(J),Q02(J),Q03(J),
     &               RxCont1,RyCont1,RzCont1,Icontact1)
            CALL DETECT(RadA(j),RadB(j),RadC(j),SAI0,SBI0,SCI0,
     &               RX(J),RY(J),RZ(J),RX(I),RY(I),RZ(I),
     &               Q00(J),Q01(J),Q02(J),Q03(J),
     &               Q00(I),Q01(I),Q02(I),Q03(I),
     &               RxCont2,RyCont2,RzCont2,Icontact2)
            IF(Icontact1.eq.1.and.Icontact2.eq.1)then
            CNNtot=CNNtot+1.0
            NCN(i)=NCN(i)+1
            ENDIF
          endif
899         Continue
      ENDIF
c**********************************
c      Orientation in Planes
c**********************************
       IF(IOren.eq.1)then
c------Equatorial plane is X-O-Y
       deltangle=3.0
       DisRmXY=sqrt(orenXA(i)**2+orenYA(i)**2+1.e-30)
       DISRmxyz=sqrt(orenXA(i)**2+orenYA(i)**2+orenZA(i)**2)
       AngLong1=acos(DISRmXY/DisRmXYZ)*180./pi
       if(orenZA(i).gt.0)AngLong1=AngLong1
c      if(orenYA(i).lt.0)AngLong1=90.0-AngLong1
       Nang1=INT(AngLong1/deltangle)+1
       NXang1(Nang1)=NXang1(Nang1)+1
c------Equatorial plane is X-O-Z
c       DisRmXZ=sqrt(orenXA(i)**2+orenZA(i)**2+1.e-30)
c       AngLong2=acos(abs(orenXA(i))/DisRmXZ)*180./pi
c       if(orenXA(i).gt.0.and.orenZA(i).gt.0)AngLong2=AngLong2
c       if(orenXA(i).lt.0.and.orenZA(i).gt.0)AngLong2=180.0-AngLong2
c       if(orenXA(i).lt.0.and.orenZA(i).lt.0)AngLong2=180.0+AngLong2
c       if(orenXA(i).gt.0.and.orenZA(i).lt.0)AngLong2=360.0-AngLong2
c       if(orenXA(i).eq.0.and.orenZA(i).gt.0)AngLong2=90.0
c       if(orenXA(i).eq.0.and.orenZA(i).lt.0)AngLong2=270.0
c       if(orenZA(i).eq.0.and.orenXA(i).gt.0)AngLong2=0.0
c       if(orenZA(i).eq.0.and.orenXA(i).lt.0)AngLong2=180.0
c       Nang2=INT(AngLong2/deltangle)+1
c       NXang2(Nang2)=NXang2(Nang2)+1
c------Equatorial plane is Y-O-Z
c       DisRmYZ=sqrt(orenYA(i)**2+orenZA(i)**2+1.e-30)
c       AngLong3=acos(abs(orenYA(i))/DisRmYZ)*180./pi
c       if(orenYA(i).gt.0.and.orenZA(i).gt.0)AngLong3=AngLong3
c       if(orenYA(i).lt.0.and.orenZA(i).gt.0)AngLong3=180.0-AngLong3
c       if(orenYA(i).lt.0.and.orenZA(i).lt.0)AngLong3=180.0+AngLong3
c       if(orenYA(i).gt.0.and.orenZA(i).lt.0)AngLong3=360.0-AngLong3
c       if(orenYA(i).eq.0.and.orenZA(i).gt.0)AngLong3=90.0
c       if(orenYA(i).eq.0.and.orenZA(i).lt.0)AngLong3=270.0
c       if(orenZA(i).eq.0.and.orenYA(i).gt.0)AngLong3=0.0
c       if(orenZA(i).eq.0.and.orenYA(i).lt.0)AngLong3=180.0
c       Nang3=INT(AngLong3/deltangle)+1
c       NXang3(Nang3)=NXang3(Nang3)+1
       ENDIF
c--------------
189    CONTINUE
c---------------------------------------------
c-------record the density and CN
       averCN=CNNtot/(ICNNtot*1.0+1.0e-30)
       porosity5=1.-cutvolume/volumespace
       bedaver=1.-porosity5
      open(61,file='ppor.dat',position='append')
      write(61,661)it,apratio,porosity5,bedaver,averCN,sphety,
     & hmaxrz   
661   format((i8),6(g15.6))
      close(61)
	stop
c-----------------
c559   CONTINUE
c*********************************************************
c---------------------------------------------------------
       open(76,file='coordination.dat')
       maxCN=0
       DO I=1,N
       maxCN=max(maxCN,NCN(i))
       write(76,776)Rx(i)*diam,Ry(i)*diam,Rz(i)*diam,diam,
     & Icolor(i),i,NCN(i)
776    format(4(g15.6),3(i8))
       ENDDO
       write(*,*)apratio,maxCN
       CLOSE(76)
c********************************************
c----------------------------------------------
       if(Ioren.eq.1)then      
       open(79,file='orenXY.dat')
       open(78,file='orenXZ.dat')
       open(77,file='orenYZ.dat')
       Totangle=deltangle/2.0
       Iangle=90/Int(deltangle)
       DO 188 I=1,Iangle
       write(79,779)totangle,NXang1(i),NXang1(i)*100./(N*1.0)
c       write(78,779)totangle,NXang2(i)
c       write(77,779)totangle,NXang3(i)
       Totangle=Totangle+deltangle
779    format((g15.6),(i8),(g15.6))
188    CONTINUE
       close(79)
       close(78)
       close(77)
       endif
c----------------------
c      STOP
      ENDIF
c      if(mod(it,500).eq.0)then
c      write(*,*)it,n,energy02,bedaver,height0
c      endif
c----------------------------------------------------------
c-------output all the forces component on the particles --
c----------------------------------------------------------
       if(mod(it,50000).eq.0)then
c       totaldrag=0.0
       totalvanforce=0.0
       totalcontforce=0.0
       tsd=it*real_dt
       DO I=1,N
       Fvanforce(i)=sqrt(Fvanforcex(i)**2+Fvanforcey(i)**2+
     &                   Fvanforcez(i)**2)
       totalvanforce=totalvanforce+Fvanforce(i)
c       contforce=fcontact(i)
       Fcontact(i)=sqrt(Fcontactx(i)**2+Fcontacty(i)**2+
     &                   Fcontactz(i)**2)
       totalcontforce=totalcontforce+fcontact(i)
       ENDDO
       open(23,file='averageforce.dat',position='append')
       write(23,602)tsd,totalvanforce/N,totalcontforce/N,fac
602    format(4(g15.6))
       close(23) 
       endif
c**************************************************************************
c     Write data for ParaView 
c     This is to write a vtu(unstructured grid) data file that can be 
c     applied to the software of ParaView
c*************************************************************************
      IF(mod(it,50000).eq.0.and.N.ne.0)then
c------------------
      Nparticledat=Nparticledat+1
      str2='t'
      str3='.vtu'
      tsd=it*real_dt
      write(str1,"(i3.3)")Nparticledat
c---------------------------
      open(24,file=trim(str2)//trim(str1)//trim(str3))
       write(24,*)'<VTKFile type="UnstructuredGrid" version="0.1"'
       write(24,*) '                   byte_order="LittleEndian">'
       write(24,*)'<UnstructuredGrid>'
       write(24,*)'<Piece NumberOfPoints="',n,'" NumberOfCells="1">'
       write(24,*)'<Points>'
       write(24,*)'<DataArray type="Float32" NumberOfComponents="3">'
       DO I=1,N
       write(24,608)rx(i)*diam,ry(i)*diam,rz(i)*diam
608    format(3(g15.6))
       ENDDO
       write(24,*)'</DataArray>'
       write(24,*)'</Points>'
c----------------------------
       write(24,*)'<Cells>'
       write(24,*)'<DataArray type="Int32" Name="connectivity">'
       DO I=1,N
       write(24,609)I-1
609    format((i8)) 
       ENDDO

       write(24,*)'</DataArray>'
       write(24,*)'<DataArray type="Int32" Name="offsets">'
       write(24,*)N
       write(24,*)'</DataArray>'
       write(24,*)'<DataArray type="UInt8" Name="types"> 1 </DataArray>'
       write(24,*)'</Cells>'
       write(24,*)'<PointData>'
       write(24,*)'<DataArray type="Float32" Name="diameter">'
       DO I=1,N
       write(24,610)diaA(i)*diam
610    format((g15.6))
       ENDDO
       write(24,*)'</DataArray>'

       write(24,*)'<DataArray type="Int32" Name="Color">'
       DO I=1,N
       write(24,611)Icolor(i)
611    format((i8))
       ENDDO
       write(24,*)'</DataArray>'

       write(24,*)'<DataArray type="Float32" Name="bedheight">'
       DO I=1,N
       write(24,606)rz(i)*diam
606    format((g15.6))
       ENDDO
       write(24,*)'</DataArray>'

       write(24,*)'<DataArray type="Float32" Name="vdwforce">'
       DO I=1,N
       Fvanforce(i)=sqrt(Fvanforcex(i)**2+Fvanforcey(i)**2+
     &                   Fvanforcez(i)**2)
       write(24,604)Fvanforce(i)
604    format((g15.6))
       ENDDO
       write(24,*)'</DataArray>'

       write(24,*)'<DataArray type="Float32" Name="contforce">'
       DO I=1,N
       Fcontact(i)=sqrt(Fcontactx(i)**2+Fcontacty(i)**2+
     &                   Fcontactz(i)**2)
       write(24,605)Fcontact(i)
605    format((g15.6))
       ENDDO
       write(24,*)'</DataArray>'

       write(24,*)'<DataArray type="Float32" Name="velvectors"'
       write(24,*)'                 NumberOfComponents="3">'
       DO I=1,N
       write(24,607)velo_x(i),velo_y(i),velo_z(i)
607    format(3(g15.6))
       ENDDO
       write(24,*)'</DataArray>'

       write(24,*)'<DataArray type="Float32" Name="orenvectors"'
       write(24,*)'                 NumberOfComponents="3">'
       DO I=1,N
       write(24,612)ORENXA(i),ORENYA(i),ORENZA(i)
612    format(3(g15.6))
       ENDDO
       write(24,*)'</DataArray>'
       write(24,*)'</PointData>'
       write(24,*)'</Piece>'
       write(24,*)'</UnstructuredGrid>'
       write(24,*)'</VTKFile>'
       Close(24)  
      ENDIF
c----------------------------------------------------
c    Write the geometry of the container
c    It is constructed by 8 points (a slot model with
c    the diameter of dmt25 and thickness of thick
c----------------------------------------------------
      IF(mod(it,50000).eq.0)then
c------------------
      open(22,file='geometry.vtu')
       write(22,*)'<VTKFile type="UnstructuredGrid" version="0.1"'
       write(22,*) '                   byte_order="LittleEndian">'
       write(22,*)'<UnstructuredGrid>'
       write(22,*)'<Piece NumberOfPoints="',8,'" NumberOfCells="1">'
       write(22,*)'<Points>'
       write(22,*)'<DataArray type="Float32" NumberOfComponents="3">'
       write(22,*)-dmt25*diam,-thick5*diam,0.0
       write(22,*) dmt25*diam,-thick5*diam,0.0
       write(22,*)-dmt25*diam, thick5*diam,0.0
       write(22,*) dmt25*diam, thick5*diam,0.0
       write(22,*)-dmt25*diam,-thick5*diam,hopzt*diam/4.
       write(22,*) dmt25*diam,-thick5*diam,hopzt*diam/4.
       write(22,*)-dmt25*diam, thick5*diam,hopzt*diam/4.
       write(22,*) dmt25*diam, thick5*diam,hopzt*diam/4.
       write(22,*)'</DataArray>'
       write(22,*)'</Points>'
c----------------------------
       write(22,*)'<Cells>'
       write(22,*)'<DataArray type="Int32" Name="connectivity">'
       write(22,*)0,1,2,3,4,5,6,7
       write(22,*)'</DataArray>'
       write(22,*)'<DataArray type="Int32" Name="offsets">'
       write(22,*)8
       write(22,*)'</DataArray>'
       write(22,*)'<DataArray type="UInt8" Name="types">11</DataArray>'
       write(22,*)'</Cells>'
       write(22,*)'</Piece>'
       write(22,*)'</UnstructuredGrid>'
       write(22,*)'</VTKFile>'
       Close(22)
       ENDIF
******************************************************************
C    *      write preflow data for re-calculation                     *
C    ******************************************************************
1030  call pclock(tret)
      if(mod(it,1000).eq.0)then
      open(16,file='timesteps.dat',position='append')
      write(16,971)tret,it*1.0
971   format(2(g15.6))
      close(16)
      endif
      if(tret.ge.realtime)then
       itime=0
        open(4,file='preflow.dat')
        write(4,900)it,n,inext,ippp,itime
900     format(5(i8))
        do i=1,n
        write (4,903)rx(i),ry(i),rz(i),dx(i),dy(i),dz(i),
     &        angpx(i),angpy(i),angpz(i),angvx(i),angvy(i),angvz(i),
     &        angvxB(i),angvyB(i),angvzB(i),velo_x(i),velo_z(i),
     &        orenxa(i),orenya(i),orenza(i)
903     format(20(g20.8))
        enddo
        DO 4003 I=1,N
        write(4,905) lfirst(i),llast(i),disptwx(i),disptwy(i),
     &         disptwz(i),fricw(i),Q00(i),Q01(i),Q02(i),Q03(i),icolor(i)
905     format(2(i8),8(g20.8),i8)
4003    CONTINUE
        do k=1,maxil
        write(4,909)llist(k),disptx(k),dispty(k),disptz(k),fricp(k)
909     format(i8,4(g20.8))
        enddo
        close(4)
        realtime=realtime+3600.
        endif
C    ******************************************************************
C    *      stop programme temporarily if running time is 	      *
C    *       beyond 7 hours                                           *             
C    ******************************************************************
       if(tret.gt.5200000.)then
       open(unit=9,file='restart.dat',status='old')
       write(9,*)' 0 ',nparticledat
       CLOSE(9)
       stop
       endif
       close(6)
       if (it.lt.tstop) goto 6
C    ******************************************************************
C    *      Finish off the Run now                                    *
C    ******************************************************************
       close(10)
       close(12)
       close(11)
       Close(19)
       close(25)
       close(63)
       open(unit=8,file='jobfinished')
       write(8,*)'Job finished'
       close(8)
       open(unit=9,file='restart.dat')
       write(9,*)'1',' 1'
       close(9)
       STOP 'program comes to end'
       END
C    ******************************************************************
C    ******************************************************************
C    ******************************************************************
C    *      The following are some subroutines                        *
C    ******************************************************************
       SUBROUTINE LIN_DIS(XP,YP,XA,YA,XB,YB,DIS)
       DOUBLE PRECISION XP,YP,XA,YA,XB,YB,DIS,dp,da,db,canga,cangb,anga 
c      calcs distance from point p to line segment ab using cosine rule
       dp=sqrt((xa-xb)*(xa-xb)+(ya-yb)*(ya-yb))
       da=sqrt((xb-xp)*(xb-xp)+(yb-yp)*(yb-yp))
       db=sqrt((xa-xp)*(xa-xp)+(ya-yp)*(ya-yp))
       if(da.lt..001.or.db.lt..001.or.dp.lt..001) then
       write(6,*) 'error lin_dis'
       dis=0.D0
       return
       stop 'error lin_dis'
       endif

       canga=(db*db+dp*dp-da*da)/(2.D0*db*dp)
       cangb=(da*da+dp*dp-db*db)/(2.D0*da*dp)

       if(canga.le.0.0)then
         dis=db
       elseif(cangb.le.0.0)then
         dis=da
       else
         if(abs(canga).gt.1.01) then
          write(6,*) 'error2 lin_dis'
c         stop 'error2 lin_dis'
         endif
         if(canga.gt.1.0)canga=1.D0
         if(canga.lt.-1.0)canga=-1.D0
         anga=acos(canga)
         dis=db*sin(anga)
       endif
       return
       end
c     **********************************************************************
       SUBROUTINE VMAG(X,Y,Z,RMAG)
       DOUBLE PRECISION X,Y,Z,RMAG
c      vector magnitude 
       rmag=sqrt(x*x+y*y+z*z)
       return
       end
c     **********************************************************************
       SUBROUTINE VPROD(A1,A2,A3,B1,B2,B3,C1,C2,C3)
       DOUBLE PRECISION A1,A2,A3,B1,B2,B3,C1,C2,C3
c      vector product c=axb
       c1=a2*b3-a3*b2
       c2=a3*b1-a1*b3
       c3=a1*b2-a2*b1
       return
       end
c     *******************************************************************
       SUBROUTINE VPROJEC(P1,P2,P3,U1,U2,U3,R1,R2,R3,ITYP)
       DOUBLE PRECISION P1,P2,P3,U1,U2,U3,R1,R2,R3,q1,q2,q3,
     &                  p,r,pr 
       call vprod(u1,u2,u3,p1,p2,p3,q1,q2,q3)
       call vprod(-u1,-u2,-u3,q1,q2,q3,r1,r2,r3)
       if(ityp.eq.0)return
       p=p1*p1+p2*p2+p3*p3
       if(p.lt.0.000001)return
       r=r1*r1+r2*r2+r3*r3
       if(p.gt.r*1000)return
       pr=sqrt(p/r)
       r1=r1*pr
       r2=r2*pr
       r3=r3*pr
       return
       end
c     *******************************************************************
       SUBROUTINE VSCALE(X,Y,Z,P1,P2)
       DOUBLE PRECISION X,Y,Z,P1,P2,ratio
c      vector scaled 
       ratio=p1/p2
       x=x*ratio
       y=y*ratio
       z=z*ratio
       return
       end
c     *******************************************************************
       SUBROUTINE VUNIT(X,Y,Z)
       DOUBLE PRECISION X,Y,Z,rmag
       rmag=sqrt(x*x+y*y+z*z)
       x=x/rmag
       y=y/rmag
       z=z/rmag
       return
       end
c     *******************************************************************
       subroutine pclock(tret) ! returns the CPU time elasped
       implicit none
       real etime, tret
       real :: tarray(2)
       tret = etime(tarray)
       end subroutine pclock
c     ***************************************************************************
       function random()
       double precision new, ia,ib
       data ia,ib /86237,29415/
       new=ia+ib
       new=new-int((ia+ib)/100000)*100000
       ia=ib
       ib=new
       random=new/100000
       return
       end
C*************************************************************************
C*************************************************************************
c-------------------------------------------------------------------------
C         Subroutine   CON_HIT(XI,YI,ZI,RI,IHIT,DISTH,FLAG)  determines  *
C     whether the particle hit the container wall or the medium surface. *
C     If so, it will calculate the hitting distance.
C*************************************************************************
      SUBROUTINE HOPP_HIT(J,Ymodu,SAI,SBI,SCI,RXX,RYY,RZZ,
     &          QG0,QG1,QG2,QG3,IHIT,X00,Y00,Z00,CurveRad)
      include 'common.for'
       DOUBLE PRECISION distF1,distF2,distF3,dis0,dis1,dis2,dis3,
     & dis10,dis11,dis12,dis13,dis20,dis21,dis22,dis23,dis30,dis31,
     & dis32,dis33,Rstar0,Rstar1,Rstar2,Rstar3,CurveRad,Ymodu,
     & SAI,SBI,SCI,RXX,RYY,RZZ,QG0,QG1,QG2,QG3,X00,Y00,Z00,
     & A1,B1,C1,F1,G1,H1,P1,Q1,R1,D1,Xleft,Yleft,Zleft,
     & Xright,Yright,Zright,Xrear,Yrear,Zrear,Xfront,Yfront,Zfront,
     & Xtop,Ytop,Ztop,Xbot,Ybot,Zbot,
     & X01,Y01,Z01,X02,Y02,Z02,X03,Y03,Z03
c      Common/Pcontact/A1,B1,C1,F1,G1,H1,P1,Q1,R1,D1
c      Common/Pcontat1/Xleft,Yleft,Zleft,Xright,Yright,Zright,
c     &                Xfront,Yfront,Zfront,Xrear,Yrear,Zrear,
c     &                Xtop,Ytop,Ztop,Xbot,Ybot,Zbot
c-----------------------------
       ihit=0
       distF1=0.D0
       distF2=0.D0
       distF3=0.D0
       dis0=0.D0
       dis1=0.D0
       dis2=0.D0
       dis3=0.D0
       dis10=0.D0
       dis11=0.D0
       dis12=0.D0
       dis13=0.D0
       dis20=0.D0
       dis21=0.D0
       dis22=0.D0
       dis23=0.D0
       dis30=0.D0
       dis31=0.D0
       dis32=0.D0
       dis33=0.D0
       Rstar0=0.D0
       Rstar1=0.D0
       Rstar2=0.D0
       Rstar3=0.D0
       CurveRad=0.D0
c----------------------------------------------------------------------
c    To determine the ellipsoid is in the wall contacting limit. If the
c    the distance between the center to the wall are larger than the 
c    Longest major semi-radius, then no contacting. 
       if (rzz.lt.0.0) return
       if (rzz.gt.hopzt) return
c       if(abs(rxx).gt.Dmt25.or.abs(ryy).gt.thick5) then
c       write(6,*) 'error rrm > hopx25 ',rxx,ryy,rzz
c       stop 'error rrm > hopx25 '
c       return
c       endif
       if(rzz.ge.height0+SAI)return
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c    By calling the following subroutines, we can determine the maximum 
c    and minimum value of an particle in each direction
c    LEFT and RIGHT: (Xright,Yright,Zright),(Xleft,Yleft,Zleft)
c    Front and REAR: (Xrear,Yrear,Zrear),(Xfront,Yfront,Zfront)
c    Top and Bottom: (Xtop,Ytop,Ztop),(Xbot,Ybot,Zbot)
      CALL PartPara(SAI,SBI,SCI,QG0,QG1,QG2,QG3,RXX,RYY,RZZ,A1,B1,C1,
     &              F1,G1,H1,P1,Q1,R1,D1)
      CALL PartX(A1,B1,C1,F1,G1,H1,P1,Q1,R1,D1,Xtop,Ytop,Ztop,
     &              Xbot,Ybot,Zbot,3)
c***********************************************************************
c----------------------------------------------------------------------
      IF(rzz.lt.height0+SAI.and.rzz.ge.0)THEN
         if(Zbot.le.height0)then
         CALL WDDD(5,Xbot,Ybot,Zbot,dis0,dis1,dis2,dis3,ihit)
         CALL Pxyz(5,height0,A1,B1,C1,F1,G1,H1,P1,Q1,R1,D1,X00,Y00,Z00)
         Call curv(j,dis1,dis2,dis3,Rstar0)
         CurveRad=Rstar0
         goto 100
         endif
      ENDIF
c--------------------------------------
100   disth(0)=dis0
      disth(1)=dis1
      disth(2)=dis2
      disth(3)=dis3
      if(disth(0).le.0.0)ihit=0
      RETURN
      END
C***************************************************************************
C***************************************************************************
C--------------------------------------------------------------------------
C     Subroutine for particle-particle overlap detection
c     This subroutine is to detect the contact of two ellipsoids. The method
c     used here is called Geometric Potential Algorithm
C--------------------------------------------------------------------------
      SUBROUTINE DETECT(SA1,SB1,SC1,SA2,SB2,SC2,X01,Y01,Z01,X02,Y02,
     &                  Z02,QI0,QI1,QI2,QI3,QJ0,QJ1,QJ2,QJ3,RxC,RyC,
     &                  RzC,Kcontact)
c--------------------------
      double precision X01,Y01,Z01,X02,Y02,Z02,QI0,QI1,QI2,QI3,
     &  SA1,SB1,SC1,SA2,SB2,SC2,QJ0,QJ1,QJ2,QJ3,
     &  A1,B1,C1,F1,G1,H1,P1,Q1,R1,D1,A2,B2,C2,F2,G2,H2,P2,Q2,R2,D2,
     &  xxa1,xxa2,xxa12,fxyz1,fxyz2,fxyz12,Xsmall,DDD,DDX,DDY,DDZ,Fdis,
     &  RxC,RyC,RzC
c------------------------------
      Kcontact=0
c-----Euler angles for each particle: matrix A(transposed)
       CALL PartPara(SA1,SB1,SC1,QI0,QI1,QI2,QI3,X01,Y01,Z01,
     &                   A1,B1,C1,F1,G1,H1,P1,Q1,R1,D1)
       CALL PartPara(SA2,SB2,SC2,QJ0,QJ1,QJ2,QJ3,X02,Y02,Z02,
     &                   A2,B2,C2,F2,G2,H2,P2,Q2,R2,D2)
c------------------------------------------------------
c-----The minimum value must be in [-10000,0]
      Xsmall=0.000001D0
      xxa1=100.D0
      xxa2=-0.1D0
      DO I=1,1000
      xxa12=0.5D0*(xxa1+xxa2)
      Call FunFxyz(A1,A2,B1,B2,C1,C2,F1,F2,G1,G2,H1,H2,P1,
     &     P2,Q1,Q2,R1,R2,D1,D2,xxa1,FXYZ1,DDD,DDX,DDY,DDZ)
      Call FunFxyz(A1,A2,B1,B2,C1,C2,F1,F2,G1,G2,H1,H2,P1,
     &     P2,Q1,Q2,R1,R2,D1,D2,xxa2,FXYZ2,DDD,DDX,DDY,DDZ)
      Call FunFxyz(A1,A2,B1,B2,C1,C2,F1,F2,G1,G2,H1,H2,P1,
     &     P2,Q1,Q2,R1,R2,D1,D2,xxa12,FXYZ12,DDD,DDX,DDY,DDZ)
      if(fxyz1*fxyz12.lt.0.0)xxa2=xxa12
      if(fxyz2*fxyz12.lt.0.0)xxa1=xxa12
      if(abs(xxa1-xxa2).le.Xsmall)goto 901
c      if(abs(FXYZ12).le.Small)goto 901
      ENDDO
c----------------------------
901   Fdis=A2*DDX*DDX+B2*DDY*DDY+C2*DDZ*DDZ+
     &       2.D0*F2*DDY*DDZ+2.D0*G2*DDZ*DDX+
     &       2.D0*H2*DDX*DDY+2.D0*P2*DDX*DDD+
     &       2.D0*Q2*DDY*DDD+2.D0*R2*DDZ*DDD+D2*DDD*DDD
      If(Fdis.lt.0.0)THEN
      KContact=1
      RxC=DDX/DDD
      RyC=DDY/DDD
      RzC=DDZ/DDD
      ENDIF
      RETURN
      END
c-----------------------------------------------------------------
       SUBROUTINE FunFxyz(XA1,XA2,XB1,XB2,XC1,XC2,XF1,XF2,XG1,XG2,
     &                    XH1,XH2,XP1,XP2,XQ1,XQ2,XR1,XR2,XD1,XD2,
     &                    xa,FFxyz,XDD,XDX,XDY,XDZ)
       Double Precision XA1,XB1,XC1,XF1,XG1,XH1,XP1,XQ1,XR1,XD1,
     &                  XA2,XB2,XC2,XF2,XG2,XH2,XP2,XQ2,XR2,XD2,
     &                  xa,ffxyz,XDD,XDX,XDY,XDZ,
     &                  A11,A12,A13,A14,A21,A22,A23,A24,A31,A32,A33,A34
c-------------------------
       A11= XA2+xa*XA1
       A12= XH2+xa*XH1
       A13= XG2+xa*XG1
       A14=-XP2-xa*XP1
       A21= XH2+xa*XH1
       A22= XB2+xa*XB1
       A23= XF2+xa*XF1
       A24=-XQ2-xa*XQ1
       A31= XG2+xa*XG1
       A32= XF2+xa*XF1
       A33= XC2+xa*XC1
       A34=-XR2-xa*XR1
       XDD=A11*A22*A33+A12*A23*A31+A13*A21*A32-
     &    A11*A23*A32-A12*A21*A33-A13*A22*A31
       XDX=A14*A22*A33+A12*A23*A34+A13*A24*A32-
     &   A14*A23*A32-A12*A24*A33-A13*A22*A34
       XDY=A11*A24*A33+A14*A23*A31+A13*A21*A34-
     &   A11*A23*A34-A14*A21*A33-A13*A24*A31
       XDZ=A11*A22*A34+A12*A24*A31+A14*A21*A32-
     &   A11*A24*A32-A12*A21*A34-A14*A22*A31
       FFXYZ=XA1*XDX*XDX+XB1*XDY*XDY+XC1*XDZ*XDZ+
     &   2.D0*XF1*XDY*XDZ+2.D0*XG1*XDX*XDZ+2.D0*XH1*XDX*XDY+
     &   2.D0*XP1*XDX*XDD+2.D0*XQ1*XDY*XDD+2.D0*XR1*XDZ*XDD+XD1*XDD*XDD
      RETURN
      END
c-----------------------------------------------------------------
c-----------------------------------------------------------------
       SUBROUTINE PartPara(SAA,SBB,SCC,QQ0,QQ1,QQ2,QQ3,XX01,YY01,ZZ01,
     &                   ZA1,ZB1,ZC1,ZF1,ZG1,ZH1,ZP1,ZQ1,ZR1,ZD1)
       Double Precision SAA,SBB,SCC,QQ0,QQ1,QQ2,QQ3,XX01,YY01,ZZ01,
     &                  ZA1,ZB1,ZC1,ZF1,ZG1,ZH1,ZP1,ZQ1,ZR1,ZD1,
     &                  SA11,SB11,SC11,AA1,BB1,CC1,DD1,EE1,FF1,
     &                  Ea01,Ea02,Ea03,EB01,EB02,EB03,EC01,EC02,EC03
c-------------------------
      SA11=SAA*SAA
      SB11=SBB*SBB
      SC11=SCC*SCC
      CALL FunMatrix(1,QQ0,QQ1,QQ2,QQ3,Ea01,Ea02,Ea03,
     &               Eb01,Eb02,Eb03,Ec01,Ec02,Ec03)
      AA1=EA01*EA01/SA11+EB01*EB01/SB11+EC01*EC01/SC11
      BB1=EA02*EA02/SA11+EB02*EB02/SB11+EC02*EC02/SC11
      CC1=EA03*EA03/SA11+EB03*EB03/SB11+EC03*EC03/SC11
      DD1=EA01*EA02/SA11+EB01*EB02/SB11+EC01*EC02/SC11
      EE1=EA01*EA03/SA11+EB01*EB03/SB11+EC01*EC03/SC11
      FF1=EA03*EA02/SA11+EB03*EB02/SB11+EC03*EC02/SC11
      ZA1=AA1
      ZB1=BB1
      ZC1=CC1
      ZF1=FF1
      ZG1=EE1
      ZH1=DD1
      ZP1=-AA1*XX01-DD1*YY01-EE1*ZZ01
      ZQ1=-DD1*XX01-BB1*YY01-FF1*ZZ01
      ZR1=-EE1*XX01-FF1*YY01-CC1*ZZ01
      ZD1=AA1*XX01*XX01+BB1*YY01*YY01+CC1*ZZ01*ZZ01+
     &   2.D0*DD1*XX01*YY01+2.D0*EE1*XX01*ZZ01+2.D0*FF1*YY01*ZZ01-1.D0
      RETURN
      END
C*****************************************************************
c-----------------------------------------------------------------
c     Subrotine PartX is to calcuate the Maximum and Mininum values 
c     of an ellipsoid position in X, Y, Z direction
c-----------------------------------------------------------------
      SUBROUTINE PartX(A1,B1,C1,F1,G1,H1,P1,Q1,R1,D1,X001,Y001,Z001,
     &                 X002,Y002,Z002,IXYZ)
      Double Precision A1,B1,C1,F1,G1,H1,P1,Q1,R1,D1,X001,Y001,Z001,
     &   X002,Y002,Z002,ev1,ev2,ev3,ev4,TEMP,xsmall,aaa,bbb,ccc
c------------------------------
c     IXYZ=1: Max and Min at X direction
      xsmall=0.0D0
      IF(IXYZ.eq.1)THEN
      TEMP=F1*F1-B1*C1
      if(TEMP.eq.0)xsmall=1.0D-30
      ev1=-(G1*F1-H1*C1)/(F1*F1-B1*C1+xsmall)
      ev2=-(R1*F1-Q1*C1)/(F1*F1-B1*C1+xsmall)
      ev3= (B1*G1-H1*F1)/(F1*F1-B1*C1+xsmall)
      ev4= (B1*R1-Q1*F1)/(F1*F1-B1*C1+xsmall)
      aaa=A1+C1*ev3*ev3+B1*ev1*ev1+2*G1*ev3+2*H1*ev1+2*F1*ev1*ev3
      bbb=2*B1*ev1*ev2+2*C1*ev3*ev4+2*F1*ev2*ev3+2*F1*ev1*ev4+
     *    2*G1*ev4+2*H1*ev2+2*P1+2*Q1*ev1+2*R1*ev3
      ccc=2*R1*ev4+B1*ev2*ev2+2*F1*ev2*ev4+2*Q1*ev2+C1*ev4*ev4+D1
c     X001,Y001,Z001 are the Maximum Values (right, front or top walls)
      X001=(-bbb+sqrt(bbb*bbb-4*aaa*ccc))/2/aaa
      Y001=ev1*X001+ev2
      Z001=ev3*X001+ev4
c     X002,Y002,Z002 are Minimum Values (left, rear, or bottom walls)
      X002=(-bbb-sqrt(bbb*bbb-4*aaa*ccc))/2/aaa
      Y002=ev1*X002+ev2
      Z002=ev3*X002+ev4
      ENDIF
c     IXYZ=2: Max and Min at Y direction
      IF(IXYZ.eq.2)THEN
      TEMP=G1*G1-A1*C1
      if(TEMP.eq.0)xsmall=1.0D-30
      ev1=(H1*C1-F1*G1)/(G1*G1-A1*C1+xsmall)
      ev2=(P1*C1-R1*G1)/(G1*G1-A1*C1+xsmall)
      ev3=(A1*F1-H1*G1)/(G1*G1-A1*C1+xsmall)
      ev4=(R1*A1-P1*G1)/(G1*G1-A1*C1+xsmall)
      aaa=A1*ev1*ev1+B1+C1*ev3*ev3+2*F1*ev3+2*G1*ev1*ev3+2*H1*ev1
      bbb=2*A1*ev1*ev2+2*C1*ev3*ev4+2*F1*ev4+2*G1*ev2*ev3+
     *    2*G1*ev1*ev4+2*H1*ev2+2*P1*ev1+2*R1*ev3+2*Q1
      ccc=A1*ev2*ev2+C1*ev4*ev4+2*P1*ev2+2*G1*ev2*ev4+2*R1*ev4+D1
c     X001,Y001,Z001 are the Maximum Values (right, front or top walls)
      Y001=(-bbb+sqrt(bbb*bbb-4*aaa*ccc))/2/aaa
      X001=ev1*Z001+ev2
      Z001=ev3*Z001+ev4
c     X002,Y002,Z002 are Minimum Values (left, rear, or bottom walls)
      Y002=(-bbb-sqrt(bbb*bbb-4*aaa*ccc))/2/aaa
      X002=ev1*Z002+ev2
      Z002=ev3*Z002+ev4
      ENDIF
c     IXYZ=3: Max and Min at Z direction
      IF(IXYZ.eq.3)THEN
      TEMP=H1*H1-A1*B1
      if(TEMP.eq.0)xsmall=1.0D-30
      ev1=(A1*F1-G1*H1)/(H1*H1-A1*B1+xsmall)
      ev2=(A1*Q1-P1*H1)/(H1*H1-A1*B1+xsmall)
      ev3=(G1*B1-H1*F1)/(H1*H1-A1*B1+xsmall)
      ev4=(P1*B1-Q1*H1)/(H1*H1-A1*B1+xsmall)
      aaa=A1*ev3*ev3+B1*ev1*ev1+C1+2*G1*ev3+2*H1*ev1*ev3+2*F1*ev1
      bbb=2*A1*ev3*ev4+2*B1*ev1*ev2+2*H1*ev1*ev4+2*H1*ev2*ev3+
     *    2*P1*ev3+2*Q1*ev1+2*F1*ev2+2*G1*ev4+2*R1
      ccc=A1*ev4*ev4+B1*ev2*ev2+2*H1*ev2*ev4+2*P1*ev4+2*Q1*ev2+D1
c     X001,Y001,Z001 are the Maximum Values (right, front or top walls)
      Z001=(-bbb+sqrt(bbb*bbb-4*aaa*ccc))/2/aaa
      Y001=ev1*Z001+ev2
      X001=ev3*Z001+ev4
c     X002,Y002,Z002 are Minimum Values (left, rear, or bottom walls)
      Z002=(-bbb-sqrt(bbb*bbb-4*aaa*ccc))/2/aaa
      Y002=ev1*Z002+ev2
      X002=ev3*Z002+ev4
      ENDIF
c-------------------------
      RETURN
      END
c-----------------------------------------------------------------------------
c     This subroutine is to calculate the elements of Matrix of transformation
c     between space-fixed coordinate system AND body-fixed coordinate system
c----------------------------------------------------------------------------
      SUBROUTINE FunMatrix(Jmatrix,Qj00,Qj01,Qj02,Qj03,
     &           Ea001,Ea002,Ea003,Eb001,Eb002,Eb003,
     &           Ec001,Ec002,Ec003)
      Double Precision Qj00,Qj01,Qj02,Qj03,Ea001,Ea002,Ea003,
     &                 Eb001,Eb002,Eb003,Ec001,Ec002,Ec003
c-------------------------
c    A(+1)A(+1)A(+1)A(+1)A(+1)A(+1)A(+1)A(+1)A(+1)A(+1)A(+1)
      If(Jmatrix.eq.1)then
      Ea001=Qj00*Qj00+Qj01*Qj01-Qj02*Qj02-Qj03*Qj03
      Ea002=2.D0*(Qj01*Qj02+Qj00*Qj03)
      Ea003=2.D0*(Qj01*Qj03-Qj00*Qj02)
      Eb001=2.D0*(Qj01*Qj02-Qj00*Qj03)
      Eb002=Qj00*Qj00-Qj01*Qj01+Qj02*Qj02-Qj03*Qj03
      Eb003=2.D0*(Qj02*Qj03+Qj00*Qj01)
      Ec001=2.D0*(Qj01*Qj03+Qj00*Qj02)
      Ec002=2.D0*(Qj02*Qj03-Qj00*Qj01)
      Ec003=Qj00*Qj00-Qj01*Qj01-Qj02*Qj02+Qj03*Qj03
      endif
c    A(-1)A(-1)A(-1)A(-1)A(-1)A(-1)A(-1)A(-1)A(-1)A(-1)A(-1)A(-1)
      if(Jmatrix.eq.2)then
      Ea001=Qj00*Qj00+Qj01*Qj01-Qj02*Qj02-Qj03*Qj03
      Ea002=2.D0*(Qj01*Qj02-Qj00*Qj03)
      Ea003=2.D0*(Qj01*Qj03+Qj00*Qj02)
      Eb001=2.D0*(Qj01*Qj02+Qj00*Qj03)
      Eb002=Qj00*Qj00-Qj01*Qj01+Qj02*Qj02-Qj03*Qj03
      Eb003=2.D0*(Qj02*Qj03-Qj00*Qj01)
      Ec001=2.D0*(Qj01*Qj03-Qj00*Qj02)
      Ec002=2.D0*(Qj02*Qj03+Qj00*Qj01)
      Ec003=Qj00*Qj00-Qj01*Qj01-Qj02*Qj02+Qj03*Qj03
      endif
      RETURN
      END
c--------------------------------------------------------------------
c-------------------------------------------------------------------
c    This subroutine is to calculate the distance between contact
c    point and max or min of point
c    1: rear(+); 2: fornt(-); 3: left(-); 4: right(+); 0, bottom(-)
c    CALL WDDD(5,Xbot,Ybot,Zbot,X00,Y00,Z00,disth(0),disth(1),
c              disth(2),disth(3),ihit)
      SUBROUTINE WDDD(Jwall,walla1,wallb1,wallc1,walld0,walld1,
     *                walld2,walld3,iconwall)
      INCLUDE 'common.for'
      Double Precision walla1,wallb1,wallc1,walld0,walld1,walld2,walld3,
     &            frontrear5,leftright5,topheight
      FRONTREAR5=thick*0.5D0
      LEFTRIGHT5=dmt2*0.5D0
      TOPHEIGHT=hopzt  
c     8.0 is the thickness btween front and rear walls
c     20.0 is the diameter of geometry (distance betwneen left and right walls)
      If(Jwall.eq.1)then
      iconwall=1
      walld0=wallb1-frontrear5
      walld1= 0.D0
      walld2=-1.D0
      walld3= 0.D0
      endif
      If(Jwall.eq.2)then
      iconwall=2
      walld0=abs(wallb1+frontrear5)
      walld1= 0.D0
      walld2= 1.D0
      walld3= 0.D0
      endif
      If(Jwall.eq.3)then
      iconwall=3
      walld0=abs(walla1+leftright5)
      walld1= 1.D0
      walld2= 0.D0
      walld3= 0.D0
      endif
      If(Jwall.eq.4)then
      iconwall=4
      walld0=(walla1-leftright5)
      walld1=-1.D0
      walld2= 0.D0
      walld3= 0.D0
      endif
      If(Jwall.eq.5)then
      iconwall=5
      walld0=abs(height0-wallc1)
      walld1= 0.D0
      walld2= 0.D0
      walld3= 1.D0
      endif
      If(Jwall.eq.6)then
      iconwall=6
      walld0=abs(wallc1-topheight)
      walld1= 0.D0
      walld2= 0.D0
      walld3=-1.D0
      endif
      RETURN
      END
c------------------------------------------------------------------
c-------------------------------------------------------------------
c    This subroutine is to caclulate the position of contact point
c    of an ellipsoid at the wall 
c    IPxyz: 1 (Front) 2(rear) 3(left) 4 (right) 5(top)
c    CALL Pxyz(2,-thick5,A1,B1,C1,F1,G1,H1,P1,Q1,R1,D1,X00,Y00,Z00)
      SUBROUTINE Pxyz(IPxyz,EEE,A1,B1,C1,F1,G1,H1,P1,Q1,R1,D1,
     &                X000,Y000,Z000)
       Double Precision EEE,A1,B1,C1,F1,G1,H1,P1,Q1,R1,D1,
     &                X000,Y000,Z000,XEE,YEE,ZEE
c-----------------------
      If(IPxyz.eq.1.or.IPxyz.eq.2)Then
      YEE=EEE
      X000=((H1*YEE*C1+P1*C1)-(F1*YEE*G1+R1*G1))/(G1*G1-A1*C1)
      Y000=YEE
      Z000=((F1*YEE*A1+R1*A1)-(H1*YEE*G1+P1*G1))/(G1*G1-A1*C1)
      Endif
      If(IPxyz.eq.3.or.IPxyz.eq.4)Then
      XEE=EEE
      X000=XEE
      Y000=((H1*XEE*C1+Q1*C1)-(G1*XEE*F1+R1*F1))/(F1*F1-B1*C1)
      Z000=((G1*XEE*B1+R1*B1)-(H1*XEE*F1+Q1*F1))/(F1*F1-B1*C1)
      Endif
      If(IPxyz.eq.5.or.IPxyz.eq.6)Then
      ZEE=EEE
      X000=((G1*ZEE*B1+P1*B1)-(F1*ZEE*H1+Q1*H1))/(H1*H1-A1*B1)
      Y000=((F1*ZEE*A1+A1*Q1)-(G1*ZEE*H1+P1*H1))/(H1*H1-A1*B1)
      Z000=ZEE
      Endif
      RETURN
      END
c------------------------------------------------------------------
c-------------------------------------------------------------------
      subroutine curv(ipN,VctNXi,VctNYi,VctNZi,Rstar0)
      include 'common.for'
      Double Precision VctNXi,VctNYi,VctNZi,Rstar0,xsmall,
     &   RcurveI,UPi,DNi,SINangLi,EccI2
c-----------------------
       xsmall=1.0D-30
       UPi=OrenXA(iPN)*VctNxi+OrenYA(iPN)*VctNYi+OrenZA(iPN)*VctNZi
       DNi=OrenXA(iPN)**2+OrenYA(iPN)**2+OrenZA(iPN)**2
       SINangLi=abs(UPi)/sqrt(DNi+xsmall)
c--------------------------------------
       EccI2=1.D0-RadC(iPN)*RadC(iPN)/RadA(iPN)/RadA(iPN)
       RcurveI=RadA(iPN)*sqrt(1.D0-EccI2)/(1.D0-EccI2*SINangLi*SINangLi)
       Rstar0=RcurveI
c     Rstar0=0.5D0
      RETURN
      END
c------------------------------------------------------------------
c-------------------------------------------------------------------
      subroutine angwxyz(wx,wy,wz,xi,yi,zi,rmx,rmy,rmz,delt,alfa,
     &                   wx1,wy1,wz1,ImaxN)
      Double precision wx,wy,wz,xi,yi,zi,rmx,rmy,rmz,delt,alfa,
     &                 wx1,wy1,wz1,wx1old,wy1old,wz1old,wx0,wy0,wz0,
     &                 deltwx,deltwy,deltwz,deltmax
c---------------------------------------
        wx1=0.D0
        wy1=0.D0
        wz1=0.D0
        DO 1291 I=1,100000
        wx1old=wx1
        wy1old=wy1
        wz1old=wz1
        wx0=wx+rmx*delt/xi+delt*(yi-zi)*(wy+wy1)*(wz+wz1)/4.D0/xi
     &      -alfa*delt*(wx+wx1)/2.D0
        wy0=wy+rmy*delt/yi+delt*(zi-xi)*(wx+wx1)*(wz+wz1)/4.D0/yi
     &      -alfa*delt*(wy+wy1)/2.D0
        wz0=wz+rmz*delt/zi+delt*(xi-yi)*(wx+wx1)*(wy+wy1)/4.D0/zi
     &      -alfa*delt*(wz+wz1)/2.D0
        wx1=wx0
        wy1=wy0
        wz1=wz0
        deltwx=abs(wx1-wx1old)
        deltwy=abs(wy1-wy1old)
        deltwz=abs(wz1-wz1old)
        deltmax=max(deltwx,deltwy,deltwz)
        if(deltmax.le.1.E-20) goto 1292
1291    CONTINUE
1292    ImaxN=i
        RETURN
        END
c****************
