      PARAMETER (MI=102,MJ=22,imax=21,jmax=101)
      PARAMETER (NN=8100,MAXIL=800000)
      PARAMETER (NZSX=20,NZSY=20,NZSZ=220,NZSP=120)
c---------------------------------
c     Particle parameters (DEM)
c----------------------------------
      COMMON /Constant/ N,IHBASE
      COMMON /PP/RX(NN),RY(NN),RZ(NN),DX(NN),DY(NN),DZ(NN),RAD(NN),
     &  FX(NN),FY(NN),fZ(NN),velo_x(nn),velo_y(nn),velo_z(nn),
     &  ORENXA(NN),ORENYA(NN),ORENZA(NN),FPPz(NN),
     &  DIAA(NN),DIAB(NN),DIAC(NN),RADA(NN),RADB(NN),RADC(NN),
     &  Q00(NN),Q01(NN),Q02(NN),Q03(NN),disth(0:3)
       DOUBLE PRECISION
     &  RX,RY,RZ,DX,DY,DZ,RAD,FX,FY,FZ,ORENXA,ORENYA,ORENZA,Fppz,
     &  DIAA,DIAB,DIAC,RADA,RADB,RADC,Q00,Q01,Q02,Q03,disth

       COMMON /PT/dmt1,dmt2,zht1,zht2,dmt15,dmt25,hgt1,hgt2,
     &            hopzt,DT,PI,diam,thick,thick5,height0
       DOUBLE PRECISION
     & dmt1,dmt2,zht1,zht2,dmt15,dmt25,hgt1,hgt2,
     & hopzt,DT,PI,diam,thick,thick5,height0

      COMMON/ht/TP(NN),TPOLD(NN),heatflux(nn),densp(nn),emods(nn),
     &           conduct_par(nn),cp_par(nn)
      DOUBLE PRECISION
     &  tp,tpold
c----------------------------
c    FLuid parameters(CFD)
c---------------------------
      COMMON/A/AP(MI,MJ),AN(MI,MJ),AS(MI,MJ),AE(MI,MJ),AW(MI,MJ),
     &         B(MI,MJ)
      COMMON/G/DDX(MI),SDDX(MI),DDY(MJ),SDDY(MJ)
      COMMON/CT/RELU,RELV,RELE,RELK,RELT,RELP,RELVS,RELY,PRELU,PRELV,
     &         KEY,CMAX,CTOTAL,IC,JC,CPRICISION,NK,NI,NJ,NII,NJJ,
     &         VIS,isource,NMIST,REAL_DT,ITERT,NITERT,ICFD,PBCF,DTIME,
     &         T0,Tinlet,CFDELX,REAL_CFDELX,PPTOP,PPBOTTOM,REAL_HEIGHT,
     &         REAL_D2,DENF
       COMMON/U/U(MI,MJ),V(MI,MJ),E(MI,MJ),XK(MI,MJ),P(MI,MJ),
     &         CP(MI,MJ),EVIS(MI,MJ),GU(MI,MJ),GV(MI,MJ),
     &    GX(0:MI,0:MI),GZ(0:MI,0:MI),UUU(MJ,MI),VVV(MJ,MI),PPP(MJ,MI),
     &         EVISS(MJ,MI),gxold(mi,mi),gzold(mi,mi),T(mi,mi)
       COMMON/DEN/DEN(MI,MJ),DENU(MI,MJ),DENV(MI,MJ)
       COMMON/FD/FFU(MI,MJ),FFV(MI,MJ),DDU(MI,MJ),DDV(MI,MJ)
       COMMON/UD/UD(MI,MJ),VD(MI,MJ),UU(MI,MJ),VV(MI,MJ)
       COMMON/CO/CD,TC1,TC2,DELTK,DELTE,DELTT,XKIN,EIN,EVISIN,ZL,DL,V10,
     &           SMALL,BIGNO
       COMMON/GGG/GGG(MI,MJ)
       COMMON/SOUR/sourceu(mi,mj),sourcev(mi,mj),apu(mi,mj),apv(mi,mj)
       COMMON/OLD/UOLD(MI,MJ),VOLD(MI,MJ),CPOLD(MI,MJ),DENOLD(MI,MJ),
     &       XKOLD(MI,MJ),EOLD(MI,MJ),POLD(MI,MJ),suu(mi,mi),svv(mi,mi),
     &       sup(mi,mi),svp(mi,mi),icell(mi,mi),xvolume(mi,mi),
     &       por(0:mi,0:mi),ppor(mi,mi),NPCELL(mi,mi,mi),NPCELLN(mi,mi),
     &       pcn(mi,mi),ppcn(mi,mi),TG(mi,mi),tsp(mi,mi),told(mi,mi),
     &       tsource(mi,mi)
       Common/mean/velx(180,180),velz(180,180),icell01(180,180),
     *            icell02(180,180),dmean(180,180),
     *            area01(180,180),area02(180,180),
     *            ffxcell(180,180),ffzcell(180,180)

