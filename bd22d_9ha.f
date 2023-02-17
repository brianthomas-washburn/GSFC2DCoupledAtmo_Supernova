C ======================================================================
C
      BLOCK DATA BD22D
C
C ======================================================================
C

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      include  'com_aers.h'
ccccccc      include 'species.i'

      CHARACTER *24 SPID(NTSP),PLQID(NPLQ)
      COMMON/SID/SPID,PLQID
      CHARACTER *24 QID(NJR),JRID(NJR)
      COMMON/QID/QID,JRID
      COMMON/QJ/QJ(NWLI,3,NJR),QTMP(3,NJR),QYIELD(NJR)
      CHARACTER *24 RE (NKR)
      COMMON/REAC/RE
      COMMON/RATE/RR  1(NRP)
     +,RR  2(NRP),RR  3(NRP),RR  4(NRP),RR  5(NRP),RR  6(NRP),RR  7(NRP)
     +,RR  8(NRP),RR  9(NRP),RR 10(NRP),RR 11(NRP),RR 12(NRP),RR 13(NRP)
     +,RR 14(NRP),RR 15(NRP),RR 16(NRP),RR 17(NRP),RR 18(NRP),RR 19(NRP)
     +,RR 20(NRP),RR 21(NRP),RR 22(NRP),RR 23(NRP),RR 24(NRP),RR 25(NRP)
     +,RR 26(NRP),RR 27(NRP),RR 28(NRP),RR 29(NRP),RR 30(NRP),RR 31(NRP)
     +,RR 32(NRP),RR 33(NRP),RR 34(NRP),RR 35(NRP),RR 36(NRP),RR 37(NRP)
     +,RR 38(NRP),RR 39(NRP),RR 40(NRP),RR 41(NRP),RR 42(NRP),RR 43(NRP)
     +,RR 44(NRP),RR 45(NRP),RR 46(NRP),RR 47(NRP),RR 48(NRP),RR 49(NRP)
     +,RR 50(NRP),RR 51(NRP),RR 52(NRP),RR 53(NRP),RR 54(NRP),RR 55(NRP)
     +,RR 56(NRP),RR 57(NRP),RR 58(NRP),RR 59(NRP),RR 60(NRP),RR 61(NRP)
     +,RR 62(NRP),RR 63(NRP),RR 64(NRP),RR 65(NRP),RR 66(NRP),RR 67(NRP)
     +,RR 68(NRP),RR 69(NRP),RR 70(NRP),RR 71(NRP),RR 72(NRP),RR 73(NRP)
     +,RR 74(NRP),RR 75(NRP),RR 76(NRP),RR 77(NRP),RR 78(NRP),RR 79(NRP)
     +,RR 80(NRP),RR 81(NRP),RR 82(NRP),RR 83(NRP),RR 84(NRP),RR 85(NRP)
     +,RR 86(NRP),RR 87(NRP),RR 88(NRP),RR 89(NRP),RR 90(NRP),RR 91(NRP)
     +,RR 92(NRP),RR 93(NRP),RR 94(NRP),RR 95(NRP),RR 96(NRP),RR 97(NRP)
     +,RR 98(NRP),RR 99(NRP),RR100(NRP),RR101(NRP),RR102(NRP),RR103(NRP)
     +,RR104(NRP),RR105(NRP),RR106(NRP),RR107(NRP),RR108(NRP),RR109(NRP)
     +,RR110(NRP),RR111(NRP),RR112(NRP),RR113(NRP),RR114(NRP),RR115(NRP)
     +,RR116(NRP),RR117(NRP),RR118(NRP),RR119(NRP),RR120(NRP),RR121(NRP)
     +,RR122(NRP),RR123(NRP),RR124(NRP),RR125(NRP),RR126(NRP),RR127(NRP)
     +,RR128(NRP),RR129(NRP),RR130(NRP),RR131(NRP),RR132(NRP),RR133(NRP)
     +,RR134(NRP),RR135(NRP),RR136(NRP),RR137(NRP),RR138(NRP),RR139(NRP)
     +,RR140(NRP),RR141(NRP),RR142(NRP),RR143(NRP),RR144(NRP),RR145(NRP)
     +,RR146(NRP),RR147(NRP),RR148(NRP),RR149(NRP),RR150(NRP),RR151(NRP)
     +,RR152(NRP),RR153(NRP),RR154(NRP),RR155(NRP),RR156(NRP),RR157(NRP)
     +,RR158(NRP),RR159(NRP),RR160(NRP),RR161(NRP),RR162(NRP),RR163(NRP)
     +,RR164(NRP),RR165(NRP),RR166(NRP),RR167(NRP),RR168(NRP),RR169(NRP)
     +,RR170(NRP),RR171(NRP),RR172(NRP),RR173(NRP),RR174(NRP),RR175(NRP)
     +,RR176(NRP),RR177(NRP),RR178(NRP),RR179(NRP),RR180(NRP),RR181(NRP)
     +,RR182(NRP),RR183(NRP),RR184(NRP),RR185(NRP),RR186(NRP),RR187(NRP)
     +,RR188(NRP),RR189(NRP),RR190(NRP),RR191(NRP),RR192(NRP),RR193(NRP)
     +,RR194(NRP),RR195(NRP),RR196(NRP),RR197(NRP),RR198(NRP),RR199(NRP)
     +,RR200(NRP),RR201(NRP),RR202(NRP),RR203(NRP),RR204(NRP),RR205(NRP)
     +,RR206(NRP),RR207(NRP),RR208(NRP),RR209(NRP),RR210(NRP),RR211(NRP)
     +,RR212(NRP),RR213(NRP),RR214(NRP),RR215(NRP),RR216(NRP),RR217(NRP)
     +,RR218(NRP),RR219(NRP),RR220(NRP),RR221(NRP),RR222(NRP),RR223(NRP)
     +,RR224(NRP),RR225(NRP),RR226(NRP),RR227(NRP),RR228(NRP),RR229(NRP)
     +,RR230(NRP),RR231(NRP),RR232(NRP),RR233(NRP),RR234(NRP),RR235(NRP)
     +,RR236(NRP),RR237(NRP),RR238(NRP),RR239(NRP),RR240(NRP),RR241(NRP)
     +,RR242(NRP),RR243(NRP),RR244(NRP),RR245(NRP),RR246(NRP),RR247(NRP)
     +,RR248(NRP),RR249(NRP),RR250(NRP),RR251(NRP),RR252(NRP),RR253(NRP)
     +,RR254(NRP),RR255(NRP),RR256(NRP),RR257(NRP),RR258(NRP),RR259(NRP)
     +,RR260(NRP),RR261(NRP),RR262(NRP),RR263(NRP),RR264(NRP),RR265(NRP)
     +,RR266(NRP),RR267(NRP),RR268(NRP),RR269(NRP),RR270(NRP),RR271(NRP)
     +,RR272(NRP),RR273(NRP),RR274(NRP),RR275(NRP),RR276(NRP),RR277(NRP)
     +,RR278(NRP),RR279(NRP),RR280(NRP),RR281(NRP),RR282(NRP),RR283(NRP)
     +,RR284(NRP),RR285(NRP),RR286(NRP),RR287(NRP),RR288(NRP),RR289(NRP)
     +,RR290(NRP),RR291(NRP),RR292(NRP),RR293(NRP),RR294(NRP),RR295(NRP)
     +,RR296(NRP),RR297(NRP),RR298(NRP),RR299(NRP),RR300(NRP),RR301(NRP)
     +,RR302(NRP),RR303(NRP),RR304(NRP),RR305(NRP),RR306(NRP),RR307(NRP)
     +,RR308(NRP),RR309(NRP),RR310(NRP),RR311(NRP),RR312(NRP),RR313(NRP)
     +,RR314(NRP),RR315(NRP),RR316(NRP),RR317(NRP),RR318(NRP),RR319(NRP)
     +,RR320(NRP),RR321(NRP),RR322(NRP),RR323(NRP),RR324(NRP),RR325(NRP)
     +,RR326(NRP),RR327(NRP),RR328(NRP),RR329(NRP),RR330(NRP),RR331(NRP)
     +,RR332(NRP),RR333(NRP),RR334(NRP),RR335(NRP),RR336(NRP),RR337(NRP)
     +,RR338(NRP),RR339(NRP),RR340(NRP),RR341(NRP),RR342(NRP),RR343(NRP)
     +,RR344(NRP),RR345(NRP),RR346(NRP),RR347(NRP),RR348(NRP),RR349(NRP)
     +,RR350(NRP),RR351(NRP),RR352(NRP),RR353(NRP),RR354(NRP),RR355(NRP)
     +,RR356(NRP),RR357(NRP),RR358(NRP),RR359(NRP),RR360(NRP),RR361(NRP)
     +,RR362(NRP),RR363(NRP),RR364(NRP),RR365(NRP),RR366(NRP),RR367(NRP)
     +,RR368(NRP),RR369(NRP),RR370(NRP),RR371(NRP),RR372(NRP),RR373(NRP)
     +,RR374(NRP),RR375(NRP),RR376(NRP),RR377(NRP),RR378(NRP),RR379(NRP)
     +,RR380(NRP),RR381(NRP),RR382(NRP)
      COMMON/CLOSURE/IFAM(NFAM),ICF(NFAM-1,NFP),ISF(NFSP,NFAM-1)
      DATA IFAM/ 0, 3, 4, 8,11,20,31,40,48,54
     +/
      DATA ICF/ 0, 0, 0, 0, 1, 2, 9, 1, 0, 0, 0, 0, 0,95,92,91,93, 0, 0
     +, 0, 0, 0, 0,48,18, 0, 0, 0, 0, 0, 0, 0, 0,19, 0, 0, 0, 0, 0, 0, 0
     +, 0,20, 0, 0, 0, 0, 0, 0, 0, 0,25, 0, 0, 0, 0, 0, 0, 0, 0,31, 0, 0
     +, 0, 0, 0, 0, 0, 0,45, 0, 0, 0, 0, 0, 0, 0, 0,46, 0, 0, 0, 0, 0, 0
     +, 0, 0,54, 0, 0
     +/
      DATA ISF/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
     +, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
     +, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
     +, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
     +, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
     +, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
     +, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
     +, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
     +, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
     +, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
     +, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 1, 1, 1, 1, 1, 0, 0, 0
     +, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
     +, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
     +, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 0, 0, 0, 0
     +, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0
     +, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1
     +, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0
     +, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
     +, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
     +, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
     +, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
     +, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
     +, 0, 0, 0, 0, 0
     +/
      COMMON/TREAT/LTPW(NRSP),INDBC(NRSP)
      DATA LTPW/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     + 1,1,1,1,1,1,1,2,2,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1
     +/
      DATA SPID(  1)/'O1D     '/
      DATA SPID(  2)/'O       '/
      DATA SPID(  3)/'O2D     '/
      DATA SPID(  4)/'O3      '/
      DATA SPID(  5)/'H       '/
      DATA SPID(  6)/'OH      '/
      DATA SPID(  7)/'HO2     '/
      DATA SPID(  8)/'H2O2    '/
      DATA SPID(  9)/'MEO2    '/
      DATA SPID( 10)/'MEOH    '/
      DATA SPID( 11)/'CH2O    '/
      DATA SPID( 12)/'I       '/
      DATA SPID( 13)/'IO      '/
      DATA SPID( 14)/'I2      '/
      DATA SPID( 15)/'I2O2    '/
      DATA SPID( 16)/'HI      '/
      DATA SPID( 17)/'HOI     '/
      DATA SPID( 18)/'INO     '/
      DATA SPID( 19)/'INO2    '/
      DATA SPID( 20)/'INO3    '/
      DATA SPID( 21)/'CL      '/
      DATA SPID( 22)/'CLO     '/
      DATA SPID( 23)/'HCL     '/
      DATA SPID( 24)/'HOCL    '/
      DATA SPID( 25)/'CLNO3   '/
      DATA SPID( 26)/'CLO3    '/
      DATA SPID( 27)/'OCLO    '/
      DATA SPID( 28)/'CLOO    '/
      DATA SPID( 29)/'CL2O2   '/
      DATA SPID( 30)/'CL2     '/
      DATA SPID( 31)/'CLNO2   '/
      DATA SPID( 32)/'N       '/
      DATA SPID( 33)/'NO      '/
      DATA SPID( 34)/'NO2     '/
      DATA SPID( 35)/'NO3     '/
      DATA SPID( 36)/'N2O5    '/
      DATA SPID( 37)/'HO2NO2  '/
      DATA SPID( 38)/'HONO    '/
      DATA SPID( 39)/'HNO3    '/
      DATA SPID( 40)/'HOONO   '/
      DATA SPID( 41)/'BR      '/
      DATA SPID( 42)/'BRO     '/
      DATA SPID( 43)/'HBR     '/
      DATA SPID( 44)/'HOBR    '/
      DATA SPID( 45)/'BRNO3   '/
      DATA SPID( 46)/'BRONO   '/
      DATA SPID( 47)/'BR2     '/
      DATA SPID( 48)/'BRCL    '/
      DATA SPID( 49)/'ETHO2   '/
      DATA SPID( 50)/'ETHOH   '/
      DATA SPID( 51)/'CH3CHO  '/
      DATA SPID( 52)/'CH3CO3  '/
      DATA SPID( 53)/'CH3CO3H '/
      DATA SPID( 54)/'ETHO2NO2'/
      DATA SPID( 55)/'N2O     '/
      DATA SPID( 56)/'CH4     '/
      DATA SPID( 57)/'C2H6    '/
      DATA SPID( 58)/'H2      '/
      DATA SPID( 59)/'CO      '/
      DATA SPID( 60)/'CH3CL   '/
      DATA SPID( 61)/'CCL4    '/
      DATA SPID( 62)/'CFCL3   '/
      DATA SPID( 63)/'CF2CL2  '/
      DATA SPID( 64)/'CH3CCL3 '/
      DATA SPID( 65)/'CF2O    '/
      DATA SPID( 66)/'CH3BR   '/
      DATA SPID( 67)/'CH2BR2  '/
      DATA SPID( 68)/'CHBR3   '/
      DATA SPID( 69)/'CHCLF2  '/
      DATA SPID( 70)/'C2CL3F3 '/
      DATA SPID( 71)/'C2CL2F4 '/
      DATA SPID( 72)/'C2CLF5  '/
      DATA SPID( 73)/'C2H3FCL2'/
      DATA SPID( 74)/'C2H3F2CL'/
      DATA SPID( 75)/'C2HF3CL2'/
      DATA SPID( 76)/'CBRCLF2 '/
      DATA SPID( 77)/'CBRF3   '/
      DATA SPID( 78)/'CBR2F2  '/
      DATA SPID( 79)/'C2BR2F4 '/
      DATA SPID( 80)/'CH3I    '/
      DATA SPID( 81)/'CF3I    '/
      DATA SPID( 82)/'CH2FCF3 '/
      DATA SPID( 83)/'CF3CH3  '/
      DATA SPID( 84)/'CHF3    '/
      DATA SPID( 85)/'CH2F2   '/
      DATA SPID( 86)/'CHF2CF3 '/
      DATA SPID( 87)/'CH3CHF2 '/
      DATA SPID( 88)/'C3HF7   '/
      DATA SPID( 89)/'C3H3F5  '/
      DATA SPID( 90)/'OX      '/
      DATA SPID( 91)/'NOX     '/
      DATA SPID( 92)/'CLX     '/
      DATA SPID( 93)/'BRX     '/
      DATA SPID( 94)/'FX      '/
      DATA SPID( 95)/'IX      '/
      DATA SPID( 96)/'PAN     '/
      DATA SPID( 97)/'SHNO3   '/
      DATA SPID( 98)/'H2O     '/
      DATA SPID( 99)/'SH2O    '/
      DATA SPID(100)/'ODP     '/
      DATA SPID(101)/'CS2     '/
      DATA SPID(102)/'DMS     '/
      DATA SPID(103)/'H2S     '/
      DATA SPID(104)/'MSA     '/
      DATA SPID(105)/'OCS     '/
      DATA SPID(106)/'SO2     '/
      DATA SPID(107)/'SO3     '/
      DATA SPID(108)/'H2SO4   '/
      DATA SPID(109)/'CH3CN   '/
      DATA SPID(110)/'HCN     '/
      DATA SPID(111)/'CO2     '/
      DATA SPID(112)/'NOY     '/
      DATA SPID(113)/'CLY     '/
      DATA SPID(114)/'O2      '/
      DATA SPID(115)/'N2      '/
      DATA SPID(116)/'M       '/
C
      DATA PLQID(  1)/'N2OPROD                 '/
      DATA PLQID(  2)/'N2OLLOS                 '/
      DATA PLQID(  3)/'CH4PROD                 '/
      DATA PLQID(  4)/'CH4LLOS                 '/
      DATA PLQID(  5)/'C2H6PROD                '/
      DATA PLQID(  6)/'C2H6LLOS                '/
      DATA PLQID(  7)/'H2PROD                  '/
      DATA PLQID(  8)/'H2LLOS                  '/
      DATA PLQID(  9)/'COPROD                  '/
      DATA PLQID( 10)/'COLLOS                  '/
      DATA PLQID( 11)/'CH3CLPROD               '/
      DATA PLQID( 12)/'CH3CLLLOS               '/
      DATA PLQID( 13)/'CCL4PROD                '/
      DATA PLQID( 14)/'CCL4LLOS                '/
      DATA PLQID( 15)/'CFCL3PROD               '/
      DATA PLQID( 16)/'CFCL3LLOS               '/
      DATA PLQID( 17)/'CF2CL2PROD              '/
      DATA PLQID( 18)/'CF2CL2LLOS              '/
      DATA PLQID( 19)/'CH3CCL3PROD             '/
      DATA PLQID( 20)/'CH3CCL3LLOS             '/
      DATA PLQID( 21)/'CF2OPROD                '/
      DATA PLQID( 22)/'CF2OLLOS                '/
      DATA PLQID( 23)/'CH3BRPROD               '/
      DATA PLQID( 24)/'CH3BRLLOS               '/
      DATA PLQID( 25)/'CH2BR2PROD              '/
      DATA PLQID( 26)/'CH2BR2LLOS              '/
      DATA PLQID( 27)/'CHBR3PROD               '/
      DATA PLQID( 28)/'CHBR3LLOS               '/
      DATA PLQID( 29)/'CHCLF2PROD              '/
      DATA PLQID( 30)/'CHCLF2LLOS              '/
      DATA PLQID( 31)/'C2CL3F3PROD             '/
      DATA PLQID( 32)/'C2CL3F3LLOS             '/
      DATA PLQID( 33)/'C2CL2F4PROD             '/
      DATA PLQID( 34)/'C2CL2F4LLOS             '/
      DATA PLQID( 35)/'C2CLF5PROD              '/
      DATA PLQID( 36)/'C2CLF5LLOS              '/
      DATA PLQID( 37)/'C2H3FCL2PROD            '/
      DATA PLQID( 38)/'C2H3FCL2LLOS            '/
      DATA PLQID( 39)/'C2H3F2CLPROD            '/
      DATA PLQID( 40)/'C2H3F2CLLLOS            '/
      DATA PLQID( 41)/'C2HF3CL2PROD            '/
      DATA PLQID( 42)/'C2HF3CL2LLOS            '/
      DATA PLQID( 43)/'CBRCLF2PROD             '/
      DATA PLQID( 44)/'CBRCLF2LLOS             '/
      DATA PLQID( 45)/'CBRF3PROD               '/
      DATA PLQID( 46)/'CBRF3LLOS               '/
      DATA PLQID( 47)/'CBR2F2PROD              '/
      DATA PLQID( 48)/'CBR2F2LLOS              '/
      DATA PLQID( 49)/'C2BR2F4PROD             '/
      DATA PLQID( 50)/'C2BR2F4LLOS             '/
      DATA PLQID( 51)/'CH3IPROD                '/
      DATA PLQID( 52)/'CH3ILLOS                '/
      DATA PLQID( 53)/'CF3IPROD                '/
      DATA PLQID( 54)/'CF3ILLOS                '/
      DATA PLQID( 55)/'CH2FCF3PROD             '/
      DATA PLQID( 56)/'CH2FCF3LLOS             '/
      DATA PLQID( 57)/'CF3CH3PROD              '/
      DATA PLQID( 58)/'CF3CH3LLOS              '/
      DATA PLQID( 59)/'CHF3PROD                '/
      DATA PLQID( 60)/'CHF3LLOS                '/
      DATA PLQID( 61)/'CH2F2PROD               '/
      DATA PLQID( 62)/'CH2F2LLOS               '/
      DATA PLQID( 63)/'CHF2CF3PROD             '/
      DATA PLQID( 64)/'CHF2CF3LLOS             '/
      DATA PLQID( 65)/'CH3CHF2PROD             '/
      DATA PLQID( 66)/'CH3CHF2LLOS             '/
      DATA PLQID( 67)/'C3HF7PROD               '/
      DATA PLQID( 68)/'C3HF7LLOS               '/
      DATA PLQID( 69)/'C3H3F5PROD              '/
      DATA PLQID( 70)/'C3H3F5LLOS              '/
      DATA PLQID( 71)/'OXPROD                  '/
      DATA PLQID( 72)/'OXLLOS                  '/
      DATA PLQID( 73)/'OXQLOS                  '/
      DATA PLQID( 74)/'NOXPROD                 '/
      DATA PLQID( 75)/'NOXLLOS                 '/
      DATA PLQID( 76)/'NOXQLOS                 '/
      DATA PLQID( 77)/'CLXPROD                 '/
      DATA PLQID( 78)/'CLXLLOS                 '/
      DATA PLQID( 79)/'BRXPROD                 '/
      DATA PLQID( 80)/'BRXLLOS                 '/
      DATA PLQID( 81)/'FXPROD                  '/
      DATA PLQID( 82)/'FXLLOS                  '/
      DATA PLQID( 83)/'IXPROD                  '/
      DATA PLQID( 84)/'IXLLOS                  '/
      DATA PLQID( 85)/'PANPROD                 '/
      DATA PLQID( 86)/'PANLLOS                 '/
      DATA PLQID( 87)/'SHNO3PROD               '/
      DATA PLQID( 88)/'SHNO3LLOS               '/
      DATA PLQID( 89)/'H2OPROD                 '/
      DATA PLQID( 90)/'H2OLLOS                 '/
      DATA PLQID( 91)/'H2OQLOS                 '/
      DATA PLQID( 92)/'SH2OPROD                '/
      DATA PLQID( 93)/'SH2OLLOS                '/
      DATA PLQID( 94)/'ODPPROD                 '/
      DATA PLQID( 95)/'ODPLLOS                 '/
      DATA PLQID( 96)/'CS2PROD                 '/
      DATA PLQID( 97)/'CS2LLOS                 '/
      DATA PLQID( 98)/'DMSPROD                 '/
      DATA PLQID( 99)/'DMSLLOS                 '/
      DATA PLQID(100)/'H2SPROD                 '/
      DATA PLQID(101)/'H2SLLOS                 '/
      DATA PLQID(102)/'MSAPROD                 '/
      DATA PLQID(103)/'MSALLOS                 '/
      DATA PLQID(104)/'OCSPROD                 '/
      DATA PLQID(105)/'OCSLLOS                 '/
      DATA PLQID(106)/'SO2PROD                 '/
      DATA PLQID(107)/'SO2LLOS                 '/
      DATA PLQID(108)/'SO3PROD                 '/
      DATA PLQID(109)/'SO3LLOS                 '/
      DATA PLQID(110)/'H2SO4PROD               '/
      DATA PLQID(111)/'H2SO4LLOS               '/
      DATA PLQID(112)/'CH3CNPROD               '/
      DATA PLQID(113)/'CH3CNLLOS               '/
      DATA PLQID(114)/'HCNPROD                 '/
      DATA PLQID(115)/'HCNLLOS                 '/
C
C
* NCL=NUMBER OF CHLORINE IN ODP SPECIE
* NBR=NUMBER OF BROMINE IN ODP SPECIE
* NFL=NUMBER OF FLOURINE IN ODP SPECIE
      COMMON /ODP/NCL,NBR,NFL
      DATA NCL,NBR,NFL/0,0,0/
C
C
      CHARACTER*60 SCCSID
      COMMON/SCCS1/SCCSID
      DATA SCCSID/'$Id: bd22d.f,v 1.69 2007/04/23 16:16:14 dkweis Exp$'/
      END
