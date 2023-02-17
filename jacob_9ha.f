C ****************************************************************************
C
      SUBROUTINE JACOB(IFAM)
C
C ****************************************************************************
C
C
      IMPLICIT DOUBLE PRECISION (A-Z)
      INTEGER NCL,NBR,NFL,IFL,IFAM

      include 'com_aers.h'
ccccccc      include 'species.i'

      COMMON/PROD/PR(NCSP)/LOSS/LO(NCSP)/QLOS/QL(NCSP)
      COMMON/MATRIX/A 1(NFSP)
     +,A 2(NFSP),A 3(NFSP),A 4(NFSP),A 5(NFSP),A 6(NFSP)
     +,A 7(NFSP),A 8(NFSP),A 9(NFSP),A10(NFSP),A11(NFSP)
     +,A12(NFSP),A13(NFSP),A14(NFSP),A15(NFSP),A16(NFSP)
     +,A17(NFSP),A18(NFSP),A19(NFSP),A20(NFSP),A21(NFSP)
     +,A22(NFSP),A23(NFSP),A24(NFSP),A25(NFSP),A26(NFSP)
     +,A27(NFSP),A28(NFSP),A29(NFSP),A30(NFSP),A31(NFSP)
     +,A32(NFSP),A33(NFSP),A34(NFSP),A35(NFSP),A36(NFSP)
     +,A37(NFSP),A38(NFSP),A39(NFSP),A40(NFSP),A41(NFSP)
     +,A42(NFSP),A43(NFSP),A44(NFSP),A45(NFSP),A46(NFSP)
     +,A47(NFSP),A48(NFSP),A49(NFSP),A50(NFSP),A51(NFSP)
     +,A52(NFSP),A53(NFSP),A54(NFSP)
      COMMON/JRATE/J  1
     +,J  2,J  3,J  4,J  5,J  6,J  7,J  8,J  9,J 10,J 11,J 12,J 13
     +,J 14,J 15,J 16,J 17,J 18,J 19,J 20,J 21,J 22,J 23,J 24,J 25
     +,J 26,J 27,J 28,J 29,J 30,J 31,J 32,J 33,J 34,J 35,J 36,J 37
     +,J 38,J 39,J 40,J 41,J 42,J 43,J 44,J 45,J 46,J 47,J 48,J 49
     +,J 50,J 51,J 52,J 53,J 54,J 55,J 56,J 57,J 58,J 59,J 60,J 61
     +,J 62,J 63,J 64,J 65,J 66,J 67,J 68,J 69,J 70,J 71,J 72,J 73
     +,J 74,J 75,J 76,J 77,J 78,J 79,J 80,J 81,J 82,J 83,J 84,J 85
     +,J 86,J 87,J 88,J 89,J 90,J 91,J 92,J 93,J 94,J 95,J 96,J 97
     +,J 98,J 99,J100,J101,J102,J103,J104,J105,J106,J107,J108,J109
     +,J110,J111,J112,J113,J114,J115,J116,J117,J118,J119,J120,J121
      COMMON/KRATE/R  1
     +,R  2,R  3,R  4,R  5,R  6,R  7,R  8,R  9,R 10,R 11,R 12,R 13
     +,R 14,R 15,R 16,R 17,R 18,R 19,R 20,R 21,R 22,R 23,R 24,R 25
     +,R 26,R 27,R 28,R 29,R 30,R 31,R 32,R 33,R 34,R 35,R 36,R 37
     +,R 38,R 39,R 40,R 41,R 42,R 43,R 44,R 45,R 46,R 47,R 48,R 49
     +,R 50,R 51,R 52,R 53,R 54,R 55,R 56,R 57,R 58,R 59,R 60,R 61
     +,R 62,R 63,R 64,R 65,R 66,R 67,R 68,R 69,R 70,R 71,R 72,R 73
     +,R 74,R 75,R 76,R 77,R 78,R 79,R 80,R 81,R 82,R 83,R 84,R 85
     +,R 86,R 87,R 88,R 89,R 90,R 91,R 92,R 93,R 94,R 95,R 96,R 97
     +,R 98,R 99,R100,R101,R102,R103,R104,R105,R106,R107,R108,R109
     +,R110,R111,R112,R113,R114,R115,R116,R117,R118,R119,R120,R121
     +,R122,R123,R124,R125,R126,R127,R128,R129,R130,R131,R132,R133
     +,R134,R135,R136,R137,R138,R139,R140,R141,R142,R143,R144,R145
     +,R146,R147,R148,R149,R150,R151,R152,R153,R154,R155,R156,R157
     +,R158,R159,R160,R161,R162,R163,R164,R165,R166,R167,R168,R169
     +,R170,R171,R172,R173,R174,R175,R176,R177,R178,R179,R180,R181
     +,R182,R183,R184,R185,R186,R187,R188,R189,R190,R191,R192,R193
     +,R194,R195,R196,R197,R198,R199,R200,R201,R202,R203,R204,R205
     +,R206,R207,R208,R209,R210,R211,R212,R213,R214,R215,R216,R217
     +,R218,R219,R220,R221,R222,R223,R224,R225,R226,R227,R228,R229
     +,R230,R231,R232,R233,R234,R235,R236,R237,R238,R239,R240,R241
     +,R242,R243,R244,R245,R246,R247,R248,R249,R250,R251,R252,R253
     +,R254,R255,R256,R257,R258,R259,R260,R261,R262,R263,R264,R265
     +,R266,R267,R268,R269,R270,R271,R272,R273,R274,R275,R276,R277
     +,R278,R279,R280,R281,R282,R283,R284,R285,R286,R287,R288,R289
     +,R290,R291,R292,R293,R294,R295,R296,R297,R298,R299,R300,R301
     +,R302,R303,R304,R305,R306,R307,R308,R309,R310,R311,R312,R313
     +,R314,R315,R316,R317,R318,R319,R320,R321,R322,R323,R324,R325
     +,R326,R327,R328,R329,R330,R331,R332,R333,R334,R335,R336,R337
     +,R338,R339,R340,R341,R342,R343,R344,R345,R346,R347,R348,R349
     +,R350,R351,R352,R353,R354,R355,R356,R357,R358,R359,R360,R361
     +,R362,R363,R364,R365,R366,R367,R368,R369,R370,R371,R372,R373
     +,R374,R375,R376,R377,R378,R379,R380,R381,R382
      COMMON/SP/O1D     
     +,O        ,O2D      ,O3       ,H        ,OH       ,HO2     
     +,H2O2     ,MEO2     ,MEOH     ,CH2O     ,I        ,IO      
     +,I2       ,I2O2     ,HI       ,HOI      ,INO      ,INO2    
     +,INO3     ,CL       ,CLO      ,HCL      ,HOCL     ,CLNO3   
     +,CLO3     ,OCLO     ,CLOO     ,CL2O2    ,CL2      ,CLNO2   
     +,N        ,NO       ,NO2      ,NO3      ,N2O5     ,HO2NO2  
     +,HONO     ,HNO3     ,HOONO    ,BR       ,BRO      ,HBR     
     +,HOBR     ,BRNO3    ,BRONO    ,BR2      ,BRCL     ,ETHO2   
     +,ETHOH    ,CH3CHO   ,CH3CO3   ,CH3CO3H  ,ETHO2NO2 ,N2O     
     +,CH4      ,C2H6     ,H2       ,CO       ,CH3CL    ,CCL4    
     +,CFCL3    ,CF2CL2   ,CH3CCL3  ,CF2O     ,CH3BR    ,CH2BR2  
     +,CHBR3    ,CHCLF2   ,C2CL3F3  ,C2CL2F4  ,C2CLF5   ,C2H3FCL2
     +,C2H3F2CL ,C2HF3CL2 ,CBRCLF2  ,CBRF3    ,CBR2F2   ,C2BR2F4 
     +,CH3I     ,CF3I     ,CH2FCF3  ,CF3CH3   ,CHF3     ,CH2F2  
     +,CHF2CF3  ,CH3CHF2  ,C3HF7    ,C3H3F5   ,OX       ,NOX     
     +,CLX      ,BRX      ,FX       ,IX       ,PAN      ,SHNO3   
     +,H2O      ,SH2O     ,ODP      ,CS2      ,DMS      ,H2S     
     +,MSA      ,OCS      ,SO2      ,SO3      ,H2SO4    ,CH3CN   
     +,HCN      ,CO2      ,NOY      ,CLY      ,O2       ,N2      
     +,M       
      COMMON/ODP/NCL,NBR,NFL
      COMMON/FLAGS/IFL(90)
C
      IF(IFL(4).EQ.0)  THEN
         WRITE(16,*) '$Id: jacob.f,v 1.41 2007/04/23 17:20:37 dkweis $'
         IFL(4)=1
      END IF
C
      GO TO( 1, 2, 3, 4, 5, 6, 7, 8, 9,10),IFAM
      WRITE(16,*) ' FAMILY ',IFAM,' NOT FOUND IN JACOB'
      STOP 'ERROR IN JACOB'
C
    1 CONTINUE
C CALC. OF PRODUCTION TERM FOR O1D      SPECIE NO.   1
      PR(  1)=
     + +   O3      *J  3
     + +   N2O     *J 13
     + +   CO2     *J 81
     + +   HNO3    *J100
     + +   O2      *J119
C CALC. OF LOSS TERM FOR O1D            SPECIE NO.   1
      LO(  1)=
     + +   ODP     *R  2
     + +   H2O     *R  7
     + +   CH4     *R  8
     + +   H2      *R  9
     + +   N2O     *R 46
     + +   N2O     *R 47
     + +   CFCL3   *R 48
     + +   CF2CL2  *R 49
     + +   CH3CL   *R 50
     + +   CCL4    *R 51
     + +   CF2O    *R 55
     + +   CHCLF2  *R106
     + +   O2      *R127
     + +   O3      *R128
     + +   O3      *R129
     + +   CH4     *R130
     + +   CH3BR   *R143
     + +   CBRCLF2 *R144
     + +   CBRF3   *R145
     + +   C2CL3F3 *R146
     + +   HBR     *R173
     + +   C2CL2F4 *R175
     + +   C2CLF5  *R176
     + +   C2H3FCL2*R177
     + +   C2H3F2CL*R179
     + +   C2HF3CL2*R181
     + +   CBR2F2  *R205
     + +   C2BR2F4 *R206
     + +   N2      *R211
     + +   N2      *R212
     + +   CH3CN   *R242
     + +   HCN     *R245
     + +   CH2BR2  *R261
     + +   CHBR3   *R262
     + +   HCL     *R263
     + +   HCL     *R264
     + +   HCL     *R265
     + +   CL2     *R266
     + +R286
     + +   CO2     *R292
     + +   FX      *R293
     + +   HBR     *R294
     + +   HBR     *R295
     + +   CL2     *R296
     + +   CH2F2   *R299
     + +   CH2F2   *R300
     + +   CHF3    *R301
     + +   CHF3    *R302
     + +   CHCLF2  *R303
     + +   CFCL3   *R304
     + +   CF2CL2  *R305
     + +   CBRCLF2 *R306
     + +   CBRF3   *R307
     + +   CBR2F2  *R308
     + +   CH3CHF2 *R309
     + +   CH3CHF2 *R310
     + +   C2H3FCL2*R311
     + +   C2H3F2CL*R312
     + +   CH2FCF3 *R313
     + +   CH2FCF3 *R314
     + +   C2HF3CL2*R315
     + +   CHF2CF3 *R316
     + +   CHF2CF3 *R317
     + +   C2CL3F3 *R318
     + +   C2CL2F4 *R319
     + +   C2CLF5  *R320
     + +   C2BR2F4 *R321
     + +   CF3CH3  *R322
     + +   CF3CH3  *R323
     + +   CCL4    *R324
     + +   CHBR3   *R325
     + +   CH2BR2  *R326
     + +   CH3CCL3 *R327
     + +   CH3CCL3 *R328
     + +   CH3CL   *R329
     + +   CH3BR   *R330
     + +   C3HF7   *R331
     + +   C3HF7   *R332
     + +   C3H3F5  *R333
     + +   C3H3F5  *R334
C CALC. OF PRODUCTION TERM FOR O        SPECIE NO.   2
      PR(  2)=
     + +   O2      *J  1
     + +   O2      *J  1
     + +   O3      *J  2
     + +   NO      *J  5
     + +   NO2     *J  6
     + +   NO3     *J  7
     + +2 *CLO3    *J 21
     + +   BRO     *J 37
     + +   OCLO    *J 40
     + +   CLO     *J 61
     + +   SO3     *J 71
     + +   HO2     *J 78
     + +   H2O     *J 79
     + +   CO2     *J 80
     + +   IO      *J 90
     + +   N2O5    *J 98
     + +   HNO3    *J 99
     + +   CL2O2   *J106
     + +   O2      *J119
     + +   N       *O2      *R  4
     + +   N       *NO      *R  5
     + +   OH      *OH      *R 29
     + +   O2D     *O3      *R 72
     + +   O1D     *O2      *R127
     + +   O1D     *O3      *R129
     + +   O1D     *O3      *R129
     + +   N       *NO2     *R183
     + +   H       *HO2     *R210
     + +   O1D     *N2      *R211
     + +   O1D     *HCL     *R265
     + +   O1D     *R286
     + +   O1D     *CO2     *R292
     + +   O1D     *FX      *R293
     + +   O1D     *HBR     *R294
     + +   O1D     *CL2     *R296
     + +   O1D     *CH2F2   *R300
     + +   O1D     *CHF3    *R302
     + +   O1D     *CHCLF2  *R303
     + +   O1D     *CFCL3   *R304
     + +   O1D     *CF2CL2  *R305
     + +   O1D     *CBRCLF2 *R306
     + +   O1D     *CBRF3   *R307
     + +   O1D     *CBR2F2  *R308
     + +   O1D     *CH3CHF2 *R310
     + +   O1D     *C2H3FCL2*R311
     + +   O1D     *C2H3F2CL*R312
     + +   O1D     *CH2FCF3 *R314
     + +   O1D     *C2HF3CL2*R315
     + +   O1D     *CHF2CF3 *R317
     + +   O1D     *C2CL3F3 *R318
     + +   O1D     *C2CL2F4 *R319
     + +   O1D     *C2CLF5  *R320
     + +   O1D     *C2BR2F4 *R321
     + +   O1D     *CF3CH3  *R323
     + +   O1D     *CCL4    *R324
     + +   O1D     *CHBR3   *R325
     + +   O1D     *CH2BR2  *R326
     + +   O1D     *CH3CCL3 *R328
     + +   O1D     *CH3CL   *R329
     + +   O1D     *CH3BR   *R330
     + +   O1D     *C3HF7   *R332
     + +   O1D     *C3H3F5  *R334
C CALC. OF LOSS TERM FOR O              SPECIE NO.   2
      LO(  2)=
     + +   HO2     *R 15
     + +   OH      *R 24
     + +   CLO     *R 36
     + +   HCL     *R 38
     + +   O2      *R 39
     + +   O3      *R 40
     + +   NO2     *R 45
     + +   HOCL    *R 56
     + +   BRO     *R 61
     + +   HO2NO2  *R 79
     + +   CLNO3   *R 80
     + +   HBR     *R101
     + +   NO3     *R131
     + +   CH2O    *R142
     + +   H2O2    *R147
     + +   OCLO    *R149
     + +   NO      *R156
     + +   NO2     *R157
     + +   HOBR    *R170
     + +   BRNO3   *R172
     + +   CS2     *R185
     + +   CS2     *R186
     + +   OCS     *R188
     + +   H2S     *R190
     + +   DMS     *R196
     + +   O       *M       *R208
     + +   O       *M       *R208
     + +   CO      *M       *R213
     + +   I2      *R214
     + +   IO      *R215
     + +   CH3CN   *R241
     + +   HCN     *R244
     + +   OCLO    *R248
     + +   SO2     *R249
     + +   OCLO    *R268
     + +   HNO3    *R372
C CALC. OF PRODUCTION TERM FOR O2D      SPECIE NO.   3
      PR(  3)=
     + +   O3      *J  3
C CALC. OF LOSS TERM FOR O2D            SPECIE NO.   3
      LO(  3)=
     + +   O2      *R 58
     + +   CLO     *M       *R 59
     + +   O3      *R 72
     + +R287
     + +   NO      *R288
* CALCULATION OF THE JACOBIAN FOR FAMILY NO.  1 INCLUDING SPECIES
*  O1D     ,O       ,O2D     ,
      A 1( 1)=
     + -   ODP     *R  2                -   H2O     *R  7               
     + -   CH4     *R  8                -   H2      *R  9               
     + -   N2O     *R 46                -   N2O     *R 47               
     + -   CFCL3   *R 48                -   CF2CL2  *R 49               
     + -   CH3CL   *R 50                -   CCL4    *R 51               
     + -   CF2O    *R 55                -   CHCLF2  *R106               
     + -   O2      *R127                -   O3      *R128               
     + -   O3      *R129                -   CH4     *R130               
     + -   CH3BR   *R143                -   CBRCLF2 *R144               
     + -   CBRF3   *R145                -   C2CL3F3 *R146               
     + -   HBR     *R173                -   C2CL2F4 *R175               
     + -   C2CLF5  *R176                -   C2H3FCL2*R177               
     + -   C2H3F2CL*R179                -   C2HF3CL2*R181               
     + -   CBR2F2  *R205                -   C2BR2F4 *R206               
     + -   N2      *R211                -   N2      *R212               
     + -   CH3CN   *R242                -   HCN     *R245               
     + -   CH2BR2  *R261                -   CHBR3   *R262               
     + -   HCL     *R263                -   HCL     *R264               
     + -   HCL     *R265                -   CL2     *R266               
     + -   R286                         -   CO2     *R292               
     + -   FX      *R293                -   HBR     *R294               
     + -   HBR     *R295                -   CL2     *R296               
     + -   CH2F2   *R299                -   CH2F2   *R300               
     + -   CHF3    *R301                -   CHF3    *R302               
     + -   CHCLF2  *R303                -   CFCL3   *R304               
     + -   CF2CL2  *R305                -   CBRCLF2 *R306               
     + -   CBRF3   *R307                -   CBR2F2  *R308               
     + -   CH3CHF2 *R309                -   CH3CHF2 *R310               
     + -   C2H3FCL2*R311                -   C2H3F2CL*R312               
     + -   CH2FCF3 *R313                -   CH2FCF3 *R314               
     + -   C2HF3CL2*R315                -   CHF2CF3 *R316               
     + -   CHF2CF3 *R317                -   C2CL3F3 *R318               
     + -   C2CL2F4 *R319                -   C2CLF5  *R320               
     + -   C2BR2F4 *R321                -   CF3CH3  *R322               
     + -   CF3CH3  *R323                -   CCL4    *R324               
     + -   CHBR3   *R325                -   CH2BR2  *R326               
     + -   CH3CCL3 *R327                -   CH3CCL3 *R328               
     + -   CH3CL   *R329                -   CH3BR   *R330               
     + -   C3HF7   *R331                -   C3HF7   *R332               
     + -   C3H3F5  *R333                -   C3H3F5  *R334               
      A 4( 1)=
     + +   J  3                         -   O1D     *R128               
     + -   O1D     *R129               
      A23( 1)=
     + -   O1D     *R263                -   O1D     *R264               
     + -   O1D     *R265               
      A30( 1)=
     + -   O1D     *R266                -   O1D     *R296               
      A39( 1)=
     + +   J100                        
      A43( 1)=
     + -   O1D     *R173                -   O1D     *R294               
     + -   O1D     *R295               
      A 1( 2)=
     + +   O2      *R127                +2 *O3      *R129               
     + +   N2      *R211                +   HCL     *R265               
     + +   R286                         +   CO2     *R292               
     + +   FX      *R293                +   HBR     *R294               
     + +   CL2     *R296                +   CH2F2   *R300               
     + +   CHF3    *R302                +   CHCLF2  *R303               
     + +   CFCL3   *R304                +   CF2CL2  *R305               
     + +   CBRCLF2 *R306                +   CBRF3   *R307               
     + +   CBR2F2  *R308                +   CH3CHF2 *R310               
     + +   C2H3FCL2*R311                +   C2H3F2CL*R312               
     + +   CH2FCF3 *R314                +   C2HF3CL2*R315               
     + +   CHF2CF3 *R317                +   C2CL3F3 *R318               
     + +   C2CL2F4 *R319                +   C2CLF5  *R320               
     + +   C2BR2F4 *R321                +   CF3CH3  *R323               
     + +   CCL4    *R324                +   CHBR3   *R325               
     + +   CH2BR2  *R326                +   CH3CCL3 *R328               
     + +   CH3CL   *R329                +   CH3BR   *R330               
     + +   C3HF7   *R332                +   C3H3F5  *R334               
      A 2( 2)=
     + -   HO2     *R 15                -   OH      *R 24               
     + -   CLO     *R 36                -   HCL     *R 38               
     + -   O2      *R 39                -   O3      *R 40               
     + -   NO2     *R 45                -   HOCL    *R 56               
     + -   BRO     *R 61                -   HO2NO2  *R 79               
     + -   CLNO3   *R 80                -   HBR     *R101               
     + -   NO3     *R131                -   CH2O    *R142               
     + -   H2O2    *R147                -   OCLO    *R149               
     + -   NO      *R156                -   NO2     *R157               
     + -   HOBR    *R170                -   BRNO3   *R172               
     + -   CS2     *R185                -   CS2     *R186               
     + -   OCS     *R188                -   H2S     *R190               
     + -   DMS     *R196                -4 *O       *M       *R208      
     + -   CO      *M       *R213       -   I2      *R214               
     + -   IO      *R215                -   CH3CN   *R241               
     + -   HCN     *R244                -   OCLO    *R248               
     + -   SO2     *R249                -   OCLO    *R268               
     + -   HNO3    *R372               
      A 3( 2)=
     + +   O3      *R 72               
      A 4( 2)=
     + +   J  2                         -   O       *R 40               
     + +   O2D     *R 72                +2 *O1D     *R129               
      A 5( 2)=
     + +   HO2     *R210               
      A 6( 2)=
     + -   O       *R 24                +2 *OH      *R 29               
      A 7( 2)=
     + +   J 78                         -   O       *R 15               
     + +   H       *R210               
      A 8( 2)=
     + -   O       *R147               
      A11( 2)=
     + -   O       *R142               
      A13( 2)=
     + +   J 90                         -   O       *R215               
      A14( 2)=
     + -   O       *R214               
      A22( 2)=
     + +   J 61                         -   O       *R 36               
      A23( 2)=
     + -   O       *R 38                +   O1D     *R265               
      A24( 2)=
     + -   O       *R 56               
      A25( 2)=
     + -   O       *R 80               
      A26( 2)=
     + +2 *J 21                        
      A27( 2)=
     + +   J 40                         -   O       *R149               
     + -   O       *R248                -   O       *R268               
      A29( 2)=
     + +   J106                        
      A30( 2)=
     + +   O1D     *R296               
      A32( 2)=
     + +   O2      *R  4                +   NO      *R  5               
     + +   NO2     *R183               
      A33( 2)=
     + +   J  5                         +   N       *R  5               
     + -   O       *R156               
      A34( 2)=
     + +   J  6                         -   O       *R 45               
     + -   O       *R157                +   N       *R183               
      A35( 2)=
     + +   J  7                         -   O       *R131               
      A36( 2)=
     + +   J 98                        
      A37( 2)=
     + -   O       *R 79               
      A39( 2)=
     + +   J 99                         -   O       *R372               
      A42( 2)=
     + +   J 37                         -   O       *R 61               
      A43( 2)=
     + -   O       *R101                +   O1D     *R294               
      A44( 2)=
     + -   O       *R170               
      A45( 2)=
     + -   O       *R172               
      A 3( 3)=
     + -   O2      *R 58                -   CLO     *M       *R 59      
     + -   O3      *R 72                -   R287                        
     + -   NO      *R288               
      A 4( 3)=
     + +   J  3                         -   O2D     *R 72               
      A22( 3)=
     + -   O2D     *M       *R 59      
      A33( 3)=
     + -   O2D     *R288               
      RETURN
    2 CONTINUE
C CALC. OF PRODUCTION TERM FOR O3       SPECIE NO.   4
      PR(  4)=
     + +            J121
     + +   O       *O2      *R 39
     + +   BRO     *HO2     *R 66
     + +   CLO     *HO2     *R 73
C CALC. OF LOSS TERM FOR O3             SPECIE NO.   4
      LO(  4)=
     + +J  2
     + +J  3
     + +J110
     + +   H       *R 11
     + +   OH      *R 12
     + +   HO2     *R 14
     + +   CL      *R 34
     + +   O       *R 40
     + +   NO      *R 43
     + +   NO2     *R 44
     + +   BR      *R 60
     + +   O2D     *R 72
     + +   BRO     *R 98
     + +   O1D     *R128
     + +   O1D     *R129
     + +   SO2     *R203
     + +   I       *R224
     + +   IO      *R232
     + +   OCLO    *R269
     + +   N       *R376
     + +   HONO    *R378
     + +   CLO     *R380
     + +   CLO     *R381
* CALCULATION OF THE JACOBIAN FOR FAMILY NO.  2 INCLUDING SPECIES
*  O3      ,
      A 1( 4)=
     + -   O3      *R128                -   O3      *R129               
      A 2( 4)=
     + +   O2      *R 39                -   O3      *R 40               
      A 3( 4)=
     + -   O3      *R 72               
      A 4( 4)=
     + -   J  2                         -   J  3                        
     + -   J110                         -   H       *R 11               
     + -   OH      *R 12                -   HO2     *R 14               
     + -   CL      *R 34                -   O       *R 40               
     + -   NO      *R 43                -   NO2     *R 44               
     + -   BR      *R 60                -   O2D     *R 72               
     + -   BRO     *R 98                -   O1D     *R128               
     + -   O1D     *R129                -   SO2     *R203               
     + -   I       *R224                -   IO      *R232               
     + -   OCLO    *R269                -   N       *R376               
     + -   HONO    *R378                -   CLO     *R380               
     + -   CLO     *R381               
      A 5( 4)=
     + -   O3      *R 11               
      A 6( 4)=
     + -   O3      *R 12               
      A 7( 4)=
     + -   O3      *R 14                +   BRO     *R 66               
     + +   CLO     *R 73               
      A12( 4)=
     + -   O3      *R224               
      A13( 4)=
     + -   O3      *R232               
      A21( 4)=
     + -   O3      *R 34               
      A22( 4)=
     + +   HO2     *R 73                -   O3      *R380               
     + -   O3      *R381               
      A27( 4)=
     + -   O3      *R269               
      A32( 4)=
     + -   O3      *R376               
      A33( 4)=
     + -   O3      *R 43               
      A34( 4)=
     + -   O3      *R 44               
      A38( 4)=
     + -   O3      *R378               
      A41( 4)=
     + -   O3      *R 60               
      A42( 4)=
     + +   HO2     *R 66                -   O3      *R 98               
      RETURN
    3 CONTINUE
C CALC. OF PRODUCTION TERM FOR H        SPECIE NO.   5
      PR(  5)=
     + +   H2O     *J 12
     + +   HCL     *J 14
     + +   CH3CL   *J 19
     + +   CH3CCL3 *J 24
     + +   CH2O    *J 30
     + +   CHCLF2  *J 43
     + +   CH3CHO  *J 50
     + +   CH3CHO  *J 52
     + +   C2H3F2CL*J 66
     + +   CHBR3   *J 97
     + +   HBR     *J102
     + +   H2O2    *J105
     + +   ODP     *OH      *R  1
     + +   H2      *O1D     *R  9
     + +   OH      *CO      *R 23
     + +   OH      *O       *R 24
     + +   OH      *H2      *R 30
     + +   CL      *H2      *R 31
     + +   O1D     *CH3CL   *R 50
     + +2 *OH      *CH3BR   *R100
     + +   CHCLF2  *O1D     *R106
     + +   O1D     *C2H3FCL2*R177
     + +   O1D     *C2H3F2CL*R179
     + +   O1D     *C2HF3CL2*R181
     + +   CS2     *OH      *R184
     + +   OCS     *OH      *R187
     + +   H2S     *OH      *R189
     + +   N       *OH      *R209
     + +   OH      *CF3I    *R219
     + +2 *CL      *CH3BR   *R250
     + +   OH      *CH2BR2  *R251
     + +   CL      *CH2BR2  *R253
     + +   CL      *C2H3FCL2*R258
     + +   O1D     *CHBR3   *R262
     + +   O1D     *HCL     *R263
     + +   OH      *CCL4    *R289
     + +   OH      *CF2CL2  *R291
     + +   O1D     *HBR     *R295
     + +   O1D     *CH3CCL3 *R327
     + +   OH      *C2CL2F4 *R337
     + +   OH      *CBR2F2  *R339
     + +   OH      *CBRCLF2 *R341
     + +   OH      *C2BR2F4 *R342
     + +   OH      *CH2F2   *R343
     + +   OH      *CH3CHF2 *R345
     + +   OH      *CH2FCF3 *R347
     + +   OH      *C3H3F5  *R349
     + +   CL      *CH3CHF2 *R353
     + +   CL      *CH2FCF3 *R355
C CALC. OF LOSS TERM FOR H              SPECIE NO.   5
      LO(  5)=
     + +   O2      *R 10
     + +   O3      *R 11
     + +   HO2     *R 21
     + +   HO2     *R 22
     + +   NO2     *R133
     + +   HO2     *R210
C CALC. OF PRODUCTION TERM FOR OH       SPECIE NO.   6
      PR(  6)=
     + +   HNO3    *J  8
     + +   H2O2    *J 11
     + +   H2O2    *J 11
     + +   H2O     *J 12
     + +   HOCL    *J 15
     + +   HOBR    *J 22
     + +   HO2NO2  *J 33
     + +   MEOH    *J 34
     + +   ETHOH   *J 48
     + +   CH3CO3H *J 53
     + +   HONO    *J 58
     + +   H2SO4   *J 70
     + +   HO2     *J 78
     + +   HOI     *J 91
     + +   HOONO   *J107
     + +2 *H2O     *O1D     *R  7
     + +   CH4     *O1D     *R  8
     + +   H2      *O1D     *R  9
     + +   H       *O3      *R 11
     + +   HO2     *O3      *R 14
     + +   HO2     *O       *R 15
     + +   HO2     *NO      *R 16
     + +2 *H       *HO2     *R 22
     + +   O       *HCL     *R 38
     + +   O       *HOCL    *R 56
     + +   CL      *HOCL    *R 57
     + +   HO2NO2  *O       *R 79
     + +   O       *HBR     *R101
     + +   MEOH    *OH      *R104
     + +   ETHOH   *OH      *R111
     + +   CH3CO3H *OH      *R119
     + +   H       *NO2     *R133
     + +   HO2     *CL      *R134
     + +   CH2O    *O       *R142
     + +   O       *H2O2    *R147
     + +   O       *HOBR    *R170
     + +   O1D     *HBR     *R173
     + +   O1D     *HCL     *R264
     + +   NO3     *R285
     + +   HOONO   *R370
     + +   O       *HNO3    *R372
     + +   H2O     *NO2     *R374
C CALC. OF LOSS TERM FOR OH             SPECIE NO.   6
      LO(  6)=
     + +   ODP     *R  1
     + +   CH3CL   *R  6
     + +   O3      *R 12
     + +   H2O2    *R 13
     + +   HO2     *R 19
     + +   CO      *R 23
     + +   O       *R 24
     + +   HCL     *R 25
     + +   CH4     *R 26
     + +   NO2     *R 27
     + +   HNO3    *R 28
     + +   OH      *R 29
     + +   OH      *R 29
     + +   H2      *R 30
     + +   HOCL    *R 54
     + +   HBR     *R 70
     + +   CH3CCL3 *R 71
     + +   CLO     *R 74
     + +   CLO     *R 75
     + +   CH2O    *R 76
     + +   HO2NO2  *R 78
     + +   C2H6    *R 82
     + +   MEOH    *R 86
     + +   CLNO3   *R 87
     + +   CLNO3   *R 88
     + +   CLNO3   *R 89
     + +   BRO     *R 99
     + +   CH3BR   *R100
     + +   MEOH    *R104
     + +   CHCLF2  *R105
     + +   ETHOH   *R110
     + +   ETHOH   *R111
     + +   CH3CHO  *R112
     + +   PAN     *R117
     + +   CH3CO3H *R118
     + +   CH3CO3H *R119
     + +   NO      *R125
     + +   HONO    *R126
     + +   OCLO    *R150
     + +   CLNO2   *R151
     + +   OH      *R155
     + +   OH      *R155
     + +   BRO     *R171
     + +   C2H3FCL2*R178
     + +   C2H3F2CL*R180
     + +   C2HF3CL2*R182
     + +   CS2     *R184
     + +   OCS     *R187
     + +   H2S     *R189
     + +   DMS     *R193
     + +   DMS     *R194
     + +   DMS     *R195
     + +   SO2     *R202
     + +   N       *R209
     + +   I2      *R216
     + +   HI      *R217
     + +   CH3I    *R218
     + +   CF3I    *R219
     + +   CH3CN   *R238
     + +   CH3CN   *R239
     + +   HCN     *R243
     + +   BR2     *R246
     + +   CO      *R247
     + +   CH2BR2  *R251
     + +   CHBR3   *R252
     + +   CL2     *R267
     + +   CCL4    *R289
     + +   CFCL3   *R290
     + +   CF2CL2  *R291
     + +   CL2O2   *R297
     + +   N2O     *R335
     + +   C2CL3F3 *R336
     + +   C2CL2F4 *R337
     + +   C2CLF5  *R338
     + +   CBR2F2  *R339
     + +   CBRF3   *R340
     + +   CBRCLF2 *R341
     + +   C2BR2F4 *R342
     + +   CH2F2   *R343
     + +   CHF3    *R344
     + +   CH3CHF2 *R345
     + +   CF3CH3  *R346
     + +   CH2FCF3 *R347
     + +   CHF2CF3 *R348
     + +   C3H3F5  *R349
     + +   C3HF7   *R350
     + +   NO2     *R369
     + +   HOONO   *R371
     + +   NO3     *R373
C CALC. OF PRODUCTION TERM FOR HO2      SPECIE NO.   7
      PR(  7)=
     + +   CH2O    *J 30
     + +   HO2NO2  *J 32
     + +   MEOH    *J 34
     + +   ETHOH   *J 48
     + +   CH3CO3H *J 53
     + +   H2SO4   *J 70
     + +   H2O2    *J105
     + +   H       *O2      *R 10
     + +   OH      *O3      *R 12
     + +   OH      *H2O2    *R 13
     + +   CL      *H2O2    *R 32
     + +   HO2NO2  *R 53
     + +   BR      *H2O2    *R 69
     + +   CLO     *OH      *R 75
     + +   CH2O    *OH      *R 76
     + +   CL      *CH2O    *R 77
     + +   MEO2    *NO      *R 85
     + +   BR      *CH2O    *R 95
     + +   BRO     *OH      *R 99
     + +   ETHO2   *NO      *R108
     + +   CH2O    *O       *R142
     + +   O1D     *CH3BR   *R143
     + +   O       *H2O2    *R147
     + +   SO2     *OH      *R202
     + +   OH      *CO      *R247
     + +   OH      *N2O     *R335
     + +   OH      *NO3     *R373
C CALC. OF LOSS TERM FOR HO2            SPECIE NO.   7
      LO(  7)=
     + +J 78
     + +J109
     + +   O3      *R 14
     + +   O       *R 15
     + +   NO      *R 16
     + +   CLO     *R 17
     + +   HO2     *R 18
     + +   HO2     *R 18
     + +   OH      *R 19
     + +   CL      *R 20
     + +   H       *R 21
     + +   H       *R 22
     + +   NO2     *R 52
     + +   BR      *R 64
     + +   BRO     *R 65
     + +   BRO     *R 66
     + +   CLO     *R 73
     + +   MEO2    *R 84
     + +   ETHO2   *R109
     + +   CH3CO3  *R114
     + +   CL      *R134
     + +   H       *R210
     + +   I       *R220
     + +   IO      *R221
     + +   NO2     *R375
C CALC. OF PRODUCTION TERM FOR H2O2     SPECIE NO.   8
      PR(  8)=
     + +   HO2     *HO2     *R 18
     + +   OH      *OH      *R155
C CALC. OF LOSS TERM FOR H2O2           SPECIE NO.   8
      LO(  8)=
     + +J 11
     + +J 26
     + +J105
     + +   OH      *R 13
     + +   CL      *R 32
     + +   BR      *R 69
     + +   O       *R147
* CALCULATION OF THE JACOBIAN FOR FAMILY NO.  3 INCLUDING SPECIES
*  H       ,OH      ,HO2     ,H2O2    ,
      A 1( 5)=
     + +   H2      *R  9                +   CH3CL   *R 50               
     + +   CHCLF2  *R106                +   C2H3FCL2*R177               
     + +   C2H3F2CL*R179                +   C2HF3CL2*R181               
     + +   CHBR3   *R262                +   HCL     *R263               
     + +   HBR     *R295                +   CH3CCL3 *R327               
      A 2( 5)=
     + +   OH      *R 24               
      A 4( 5)=
     + -   H       *R 11               
      A 5( 5)=
     + -   O2      *R 10                -   O3      *R 11               
     + -   HO2     *R 21                -   HO2     *R 22               
     + -   NO2     *R133                -   HO2     *R210               
      A 6( 5)=
     + +   ODP     *R  1                +   CO      *R 23               
     + +   O       *R 24                +   H2      *R 30               
     + +2 *CH3BR   *R100                +   CS2     *R184               
     + +   OCS     *R187                +   H2S     *R189               
     + +   N       *R209                +   CF3I    *R219               
     + +   CH2BR2  *R251                +   CCL4    *R289               
     + +   CF2CL2  *R291                +   C2CL2F4 *R337               
     + +   CBR2F2  *R339                +   CBRCLF2 *R341               
     + +   C2BR2F4 *R342                +   CH2F2   *R343               
     + +   CH3CHF2 *R345                +   CH2FCF3 *R347               
     + +   C3H3F5  *R349               
      A 7( 5)=
     + -   H       *R 21                -   H       *R 22               
     + -   H       *R210               
      A 8( 5)=
     + +   J105                        
      A11( 5)=
     + +   J 30                        
      A21( 5)=
     + +   H2      *R 31                +2 *CH3BR   *R250               
     + +   CH2BR2  *R253                +   C2H3FCL2*R258               
     + +   CH3CHF2 *R353                +   CH2FCF3 *R355               
      A23( 5)=
     + +   J 14                         +   O1D     *R263               
      A32( 5)=
     + +   OH      *R209               
      A34( 5)=
     + -   H       *R133               
      A43( 5)=
     + +   J102                         +   O1D     *R295               
      A51( 5)=
     + +   J 50                         +   J 52                        
      A 1( 6)=
     + +2 *H2O     *R  7                +   CH4     *R  8               
     + +   H2      *R  9                +   HBR     *R173               
     + +   HCL     *R264               
      A 2( 6)=
     + +   HO2     *R 15                -   OH      *R 24               
     + +   HCL     *R 38                +   HOCL    *R 56               
     + +   HO2NO2  *R 79                +   HBR     *R101               
     + +   CH2O    *R142                +   H2O2    *R147               
     + +   HOBR    *R170                +   HNO3    *R372               
      A 4( 6)=
     + +   H       *R 11                -   OH      *R 12               
     + +   HO2     *R 14               
      A 5( 6)=
     + +   O3      *R 11                +2 *HO2     *R 22               
     + +   NO2     *R133               
      A 6( 6)=
     + -   ODP     *R  1                -   CH3CL   *R  6               
     + -   O3      *R 12                -   H2O2    *R 13               
     + -   HO2     *R 19                -   CO      *R 23               
     + -   O       *R 24                -   HCL     *R 25               
     + -   CH4     *R 26                -   NO2     *R 27               
     + -   HNO3    *R 28                -4 *OH      *R 29               
     + -   H2      *R 30                -   HOCL    *R 54               
     + -   HBR     *R 70                -   CH3CCL3 *R 71               
     + -   CLO     *R 74                -   CLO     *R 75               
     + -   CH2O    *R 76                -   HO2NO2  *R 78               
     + -   C2H6    *R 82                -   MEOH    *R 86               
     + -   CLNO3   *R 87                -   CLNO3   *R 88               
     + -   CLNO3   *R 89                -   BRO     *R 99               
     + -   CH3BR   *R100                -   MEOH    *R104               
     + +   MEOH    *R104                -   CHCLF2  *R105               
     + -   ETHOH   *R110                -   ETHOH   *R111               
     + +   ETHOH   *R111                -   CH3CHO  *R112               
     + -   PAN     *R117                -   CH3CO3H *R118               
     + -   CH3CO3H *R119                +   CH3CO3H *R119               
     + -   NO      *R125                -   HONO    *R126               
     + -   OCLO    *R150                -   CLNO2   *R151               
     + -4 *OH      *R155                -   BRO     *R171               
     + -   C2H3FCL2*R178                -   C2H3F2CL*R180               
     + -   C2HF3CL2*R182                -   CS2     *R184               
     + -   OCS     *R187                -   H2S     *R189               
     + -   DMS     *R193                -   DMS     *R194               
     + -   DMS     *R195                -   SO2     *R202               
     + -   N       *R209                -   I2      *R216               
     + -   HI      *R217                -   CH3I    *R218               
     + -   CF3I    *R219                -   CH3CN   *R238               
     + -   CH3CN   *R239                -   HCN     *R243               
     + -   BR2     *R246                -   CO      *R247               
     + -   CH2BR2  *R251                -   CHBR3   *R252               
     + -   CL2     *R267                -   CCL4    *R289               
     + -   CFCL3   *R290                -   CF2CL2  *R291               
     + -   CL2O2   *R297                -   N2O     *R335               
     + -   C2CL3F3 *R336                -   C2CL2F4 *R337               
     + -   C2CLF5  *R338                -   CBR2F2  *R339               
     + -   CBRF3   *R340                -   CBRCLF2 *R341               
     + -   C2BR2F4 *R342                -   CH2F2   *R343               
     + -   CHF3    *R344                -   CH3CHF2 *R345               
     + -   CF3CH3  *R346                -   CH2FCF3 *R347               
     + -   CHF2CF3 *R348                -   C3H3F5  *R349               
     + -   C3HF7   *R350                -   NO2     *R369               
     + -   HOONO   *R371                -   NO3     *R373               
      A 7( 6)=
     + +   J 78                         +   O3      *R 14               
     + +   O       *R 15                +   NO      *R 16               
     + -   OH      *R 19                +2 *H       *R 22               
     + +   CL      *R134               
      A 8( 6)=
     + +2 *J 11                         -   OH      *R 13               
     + +   O       *R147               
      A10( 6)=
     + +   J 34                         -   OH      *R 86               
     + -   OH      *R104                +   OH      *R104               
      A11( 6)=
     + -   OH      *R 76                +   O       *R142               
      A14( 6)=
     + -   OH      *R216               
      A16( 6)=
     + -   OH      *R217               
      A17( 6)=
     + +   J 91                        
      A21( 6)=
     + +   HOCL    *R 57                +   HO2     *R134               
      A22( 6)=
     + -   OH      *R 74                -   OH      *R 75               
      A23( 6)=
     + -   OH      *R 25                +   O       *R 38               
     + +   O1D     *R264               
      A24( 6)=
     + +   J 15                         -   OH      *R 54               
     + +   O       *R 56                +   CL      *R 57               
      A25( 6)=
     + -   OH      *R 87                -   OH      *R 88               
     + -   OH      *R 89               
      A27( 6)=
     + -   OH      *R150               
      A29( 6)=
     + -   OH      *R297               
      A30( 6)=
     + -   OH      *R267               
      A31( 6)=
     + -   OH      *R151               
      A32( 6)=
     + -   OH      *R209               
      A33( 6)=
     + +   HO2     *R 16                -   OH      *R125               
      A34( 6)=
     + -   OH      *R 27                +   H       *R133               
     + -   OH      *R369                +   H2O     *R374               
      A35( 6)=
     + +   R285                         -   OH      *R373               
      A37( 6)=
     + +   J 33                         -   OH      *R 78               
     + +   O       *R 79               
      A38( 6)=
     + +   J 58                         -   OH      *R126               
      A39( 6)=
     + +   J  8                         -   OH      *R 28               
     + +   O       *R372               
      A40( 6)=
     + +   J107                         +   R370                        
     + -   OH      *R371               
      A42( 6)=
     + -   OH      *R 99                -   OH      *R171               
      A43( 6)=
     + -   OH      *R 70                +   O       *R101               
     + +   O1D     *R173               
      A44( 6)=
     + +   J 22                         +   O       *R170               
      A47( 6)=
     + -   OH      *R246               
      A50( 6)=
     + +   J 48                         -   OH      *R110               
     + -   OH      *R111                +   OH      *R111               
      A51( 6)=
     + -   OH      *R112               
      A53( 6)=
     + +   J 53                         -   OH      *R118               
     + -   OH      *R119                +   OH      *R119               
      A 1( 7)=
     + +   CH3BR   *R143               
      A 2( 7)=
     + -   HO2     *R 15                +   CH2O    *R142               
     + +   H2O2    *R147               
      A 4( 7)=
     + +   OH      *R 12                -   HO2     *R 14               
      A 5( 7)=
     + +   O2      *R 10                -   HO2     *R 21               
     + -   HO2     *R 22                -   HO2     *R210               
      A 6( 7)=
     + +   O3      *R 12                +   H2O2    *R 13               
     + -   HO2     *R 19                +   CLO     *R 75               
     + +   CH2O    *R 76                +   BRO     *R 99               
     + +   SO2     *R202                +   CO      *R247               
     + +   N2O     *R335                +   NO3     *R373               
      A 7( 7)=
     + -   J 78                         -   J109                        
     + -   O3      *R 14                -   O       *R 15               
     + -   NO      *R 16                -   CLO     *R 17               
     + -4 *HO2     *R 18                -   OH      *R 19               
     + -   CL      *R 20                -   H       *R 21               
     + -   H       *R 22                -   NO2     *R 52               
     + -   BR      *R 64                -   BRO     *R 65               
     + -   BRO     *R 66                -   CLO     *R 73               
     + -   MEO2    *R 84                -   ETHO2   *R109               
     + -   CH3CO3  *R114                -   CL      *R134               
     + -   H       *R210                -   I       *R220               
     + -   IO      *R221                -   NO2     *R375               
      A 8( 7)=
     + +   J105                         +   OH      *R 13               
     + +   CL      *R 32                +   BR      *R 69               
     + +   O       *R147               
      A 9( 7)=
     + -   HO2     *R 84                +   NO      *R 85               
      A10( 7)=
     + +   J 34                        
      A11( 7)=
     + +   J 30                         +   OH      *R 76               
     + +   CL      *R 77                +   BR      *R 95               
     + +   O       *R142               
      A12( 7)=
     + -   HO2     *R220               
      A13( 7)=
     + -   HO2     *R221               
      A21( 7)=
     + -   HO2     *R 20                +   H2O2    *R 32               
     + +   CH2O    *R 77                -   HO2     *R134               
      A22( 7)=
     + -   HO2     *R 17                -   HO2     *R 73               
     + +   OH      *R 75               
      A33( 7)=
     + -   HO2     *R 16                +   MEO2    *R 85               
     + +   ETHO2   *R108               
      A34( 7)=
     + -   HO2     *R 52                -   HO2     *R375               
      A35( 7)=
     + +   OH      *R373               
      A37( 7)=
     + +   J 32                         +   R 53                        
      A41( 7)=
     + -   HO2     *R 64                +   H2O2    *R 69               
     + +   CH2O    *R 95               
      A42( 7)=
     + -   HO2     *R 65                -   HO2     *R 66               
     + +   OH      *R 99               
      A49( 7)=
     + +   NO      *R108                -   HO2     *R109               
      A50( 7)=
     + +   J 48                        
      A52( 7)=
     + -   HO2     *R114               
      A53( 7)=
     + +   J 53                        
      A 2( 8)=
     + -   H2O2    *R147               
      A 6( 8)=
     + -   H2O2    *R 13                +2 *OH      *R155               
      A 7( 8)=
     + +2 *HO2     *R 18               
      A 8( 8)=
     + -   J 11                         -   J 26                        
     + -   J105                         -   OH      *R 13               
     + -   CL      *R 32                -   BR      *R 69               
     + -   O       *R147               
      A21( 8)=
     + -   H2O2    *R 32               
      A41( 8)=
     + -   H2O2    *R 69               
      RETURN
    4 CONTINUE
C CALC. OF PRODUCTION TERM FOR MEO2     SPECIE NO.   9
      PR(  9)=
     + +   CH3BR   *J 42
     + +   CH3CHO  *J 52
     + +   CH3I    *J 84
     + +   CH4     *O1D     *R  8
     + +   OH      *CH4     *R 26
     + +   CL      *CH4     *R 35
     + +   MEOH    *OH      *R 86
     + +   CH3CO3  *NO      *R113
     + +   DMS     *O       *R196
C CALC. OF LOSS TERM FOR MEO2           SPECIE NO.   9
      LO(  9)=
     + +   HO2     *R 84
     + +   NO      *R 85
C CALC. OF PRODUCTION TERM FOR MEOH     SPECIE NO.  10
      PR( 10)=
     + +   MEO2    *HO2     *R 84
     + +   DMS     *OH      *R193
     + +   DMS     *OH      *R194
     + +   OH      *CH3I    *R218
C CALC. OF LOSS TERM FOR MEOH           SPECIE NO.  10
      LO( 10)=
     + +J 34
     + +J 35
     + +   OH      *R 86
     + +   OH      *R104
C CALC. OF PRODUCTION TERM FOR CH2O     SPECIE NO.  11
      PR( 11)=
     + +   MEOH    *J 34
     + +   CH3CO3H *J 53
     + +            J120
     + +   OH      *CH3CL   *R  6
     + +   O1D     *CH3CL   *R 50
     + +   MEO2    *NO      *R 85
     + +   OH      *CH3BR   *R100
     + +   MEOH    *OH      *R104
     + +   PAN     *OH      *R117
     + +   CH3CO3H *OH      *R119
     + +   O1D     *CH4     *R130
     + +   O1D     *CH3BR   *R143
     + +   DMS     *OH      *R193
     + +   DMS     *OH      *R194
     + +   DMS     *OH      *R195
     + +   DMS     *O       *R196
     + +   CL      *CH3CL   *R255
     + +   O1D     *CH2BR2  *R261
C CALC. OF LOSS TERM FOR CH2O           SPECIE NO.  11
      LO( 11)=
     + +J 30
     + +J 31
     + +J 36
     + +   OH      *R 76
     + +   CL      *R 77
     + +   BR      *R 95
     + +   O       *R142
* CALCULATION OF THE JACOBIAN FOR FAMILY NO.  4 INCLUDING SPECIES
*  MEO2    ,MEOH    ,CH2O    ,
      A 1( 9)=
     + +   CH4     *R  8               
      A 2( 9)=
     + +   DMS     *R196               
      A 6( 9)=
     + +   CH4     *R 26                +   MEOH    *R 86               
      A 7( 9)=
     + -   MEO2    *R 84               
      A 9( 9)=
     + -   HO2     *R 84                -   NO      *R 85               
      A10( 9)=
     + +   OH      *R 86               
      A21( 9)=
     + +   CH4     *R 35               
      A33( 9)=
     + -   MEO2    *R 85                +   CH3CO3  *R113               
      A51( 9)=
     + +   J 52                        
      A52( 9)=
     + +   NO      *R113               
      A 6(10)=
     + -   MEOH    *R 86                -   MEOH    *R104               
     + +   DMS     *R193                +   DMS     *R194               
     + +   CH3I    *R218               
      A 7(10)=
     + +   MEO2    *R 84               
      A 9(10)=
     + +   HO2     *R 84               
      A10(10)=
     + -   J 34                         -   J 35                        
     + -   OH      *R 86                -   OH      *R104               
      A 1(11)=
     + +   CH3CL   *R 50                +   CH4     *R130               
     + +   CH3BR   *R143                +   CH2BR2  *R261               
      A 2(11)=
     + -   CH2O    *R142                +   DMS     *R196               
      A 6(11)=
     + +   CH3CL   *R  6                -   CH2O    *R 76               
     + +   CH3BR   *R100                +   MEOH    *R104               
     + +   PAN     *R117                +   CH3CO3H *R119               
     + +   DMS     *R193                +   DMS     *R194               
     + +   DMS     *R195               
      A 9(11)=
     + +   NO      *R 85               
      A10(11)=
     + +   J 34                         +   OH      *R104               
      A11(11)=
     + -   J 30                         -   J 31                        
     + -   J 36                         -   OH      *R 76               
     + -   CL      *R 77                -   BR      *R 95               
     + -   O       *R142               
      A21(11)=
     + -   CH2O    *R 77                +   CH3CL   *R255               
      A33(11)=
     + +   MEO2    *R 85               
      A41(11)=
     + -   CH2O    *R 95               
      A53(11)=
     + +   J 53                         +   OH      *R119               
      RETURN
    5 CONTINUE
C CALC. OF PRODUCTION TERM FOR I        SPECIE NO.  12
      PR( 12)=
     + +2 *I2      *J 85
     + +   INO     *J 86
     + +   INO2    *J 87
     + +   INO3    *J 89
     + +   IO      *J 90
     + +   HOI     *J 91
     + +2 *I2O2    *J 92
     + +   O       *I2      *R214
     + +   O       *IO      *R215
     + +   OH      *I2      *R216
     + +   OH      *HI      *R217
     + +   NO3     *HI      *R222
     + +   IO      *NO      *R228
     + +   IO      *CLO     *R230
     + +   IO      *BRO     *R231
     + +   IO      *O3      *R232
     + +2 *IO      *IO      *R234
C CALC. OF LOSS TERM FOR I              SPECIE NO.  12
      LO( 12)=
     + +   HO2     *R220
     + +   O3      *R224
     + +   NO      *R225
     + +   NO2     *R226
     + +   BRO     *R227
C CALC. OF PRODUCTION TERM FOR IO       SPECIE NO.  13
      PR( 13)=
     + +   INO3    *J 88
     + +   O       *I2      *R214
     + +   I       *O3      *R224
     + +   I       *BRO     *R227
     + +   I2O2    *R235
     + +   I2O2    *R235
C CALC. OF LOSS TERM FOR IO             SPECIE NO.  13
      LO( 13)=
     + +J 90
     + +   O       *R215
     + +   HO2     *R221
     + +   NO      *R228
     + +   NO2     *R229
     + +   CLO     *R230
     + +   BRO     *R231
     + +   O3      *R232
     + +   IO      *R233
     + +   IO      *R233
     + +   IO      *R234
     + +   IO      *R234
C CALC. OF PRODUCTION TERM FOR I2       SPECIE NO.  14
      PR( 14)=
     + +   INO     *INO     *R236
     + +   INO2    *INO2    *R237
C CALC. OF LOSS TERM FOR I2             SPECIE NO.  14
      LO( 14)=
     + +J 85
     + +   O       *R214
     + +   OH      *R216
C CALC. OF PRODUCTION TERM FOR I2O2     SPECIE NO.  15
      PR( 15)=
     + +   IO      *IO      *R233
C CALC. OF LOSS TERM FOR I2O2           SPECIE NO.  15
      LO( 15)=
     + +J 92
     + +R235
C CALC. OF PRODUCTION TERM FOR HI       SPECIE NO.  16
      PR( 16)=
     + +   HO2     *I       *R220
C CALC. OF LOSS TERM FOR HI             SPECIE NO.  16
      LO( 16)=
     + +   OH      *R217
     + +   NO3     *R222
C CALC. OF PRODUCTION TERM FOR HOI      SPECIE NO.  17
      PR( 17)=
     + +   OH      *I2      *R216
     + +   HO2     *IO      *R221
C CALC. OF LOSS TERM FOR HOI            SPECIE NO.  17
      LO( 17)=
     + +J 91
C CALC. OF PRODUCTION TERM FOR INO      SPECIE NO.  18
      PR( 18)=
     + +   I       *NO      *R225
C CALC. OF LOSS TERM FOR INO            SPECIE NO.  18
      LO( 18)=
     + +J 86
     + +   INO     *R236
     + +   INO     *R236
C CALC. OF PRODUCTION TERM FOR INO2     SPECIE NO.  19
      PR( 19)=
     + +   I       *NO2     *R226
C CALC. OF LOSS TERM FOR INO2           SPECIE NO.  19
      LO( 19)=
     + +J 87
     + +   INO2    *R237
     + +   INO2    *R237
C CALC. OF PRODUCTION TERM FOR INO3     SPECIE NO.  20
      PR( 20)=
     + +   IO      *NO2     *R229
C CALC. OF LOSS TERM FOR INO3           SPECIE NO.  20
      LO( 20)=
     + +J 88
     + +J 89
* CALCULATION OF THE JACOBIAN FOR FAMILY NO.  5 INCLUDING SPECIES
*  I       ,IO      ,I2      ,I2O2    ,HI      ,HOI     
*  INO     ,INO2    ,INO3    ,
      A 2(12)=
     + +   I2      *R214                +   IO      *R215               
      A 4(12)=
     + -   I       *R224                +   IO      *R232               
      A 6(12)=
     + +   I2      *R216                +   HI      *R217               
      A 7(12)=
     + -   I       *R220               
      A12(12)=
     + -   HO2     *R220                -   O3      *R224               
     + -   NO      *R225                -   NO2     *R226               
     + -   BRO     *R227               
      A13(12)=
     + +   J 90                         +   O       *R215               
     + +   NO      *R228                +   CLO     *R230               
     + +   BRO     *R231                +   O3      *R232               
     + +4 *IO      *R234               
      A14(12)=
     + +2 *J 85                         +   O       *R214               
     + +   OH      *R216               
      A15(12)=
     + +2 *J 92                        
      A16(12)=
     + +   OH      *R217                +   NO3     *R222               
      A17(12)=
     + +   J 91                        
      A18(12)=
     + +   J 86                        
      A19(12)=
     + +   J 87                        
      A20(12)=
     + +   J 89                        
      A22(12)=
     + +   IO      *R230               
      A33(12)=
     + -   I       *R225                +   IO      *R228               
      A34(12)=
     + -   I       *R226               
      A35(12)=
     + +   HI      *R222               
      A42(12)=
     + -   I       *R227                +   IO      *R231               
      A 2(13)=
     + +   I2      *R214                -   IO      *R215               
      A 4(13)=
     + +   I       *R224                -   IO      *R232               
      A 7(13)=
     + -   IO      *R221               
      A12(13)=
     + +   O3      *R224                +   BRO     *R227               
      A13(13)=
     + -   J 90                         -   O       *R215               
     + -   HO2     *R221                -   NO      *R228               
     + -   NO2     *R229                -   CLO     *R230               
     + -   BRO     *R231                -   O3      *R232               
     + -4 *IO      *R233                -4 *IO      *R234               
      A14(13)=
     + +   O       *R214               
      A15(13)=
     + +2 *R235                        
      A20(13)=
     + +   J 88                        
      A22(13)=
     + -   IO      *R230               
      A33(13)=
     + -   IO      *R228               
      A34(13)=
     + -   IO      *R229               
      A42(13)=
     + +   I       *R227                -   IO      *R231               
      A 2(14)=
     + -   I2      *R214               
      A 6(14)=
     + -   I2      *R216               
      A14(14)=
     + -   J 85                         -   O       *R214               
     + -   OH      *R216               
      A18(14)=
     + +2 *INO     *R236               
      A19(14)=
     + +2 *INO2    *R237               
      A13(15)=
     + +2 *IO      *R233               
      A15(15)=
     + -   J 92                         -   R235                        
      A 6(16)=
     + -   HI      *R217               
      A 7(16)=
     + +   I       *R220               
      A12(16)=
     + +   HO2     *R220               
      A16(16)=
     + -   OH      *R217                -   NO3     *R222               
      A35(16)=
     + -   HI      *R222               
      A 6(17)=
     + +   I2      *R216               
      A 7(17)=
     + +   IO      *R221               
      A13(17)=
     + +   HO2     *R221               
      A14(17)=
     + +   OH      *R216               
      A17(17)=
     + -   J 91                        
      A12(18)=
     + +   NO      *R225               
      A18(18)=
     + -   J 86                         -4 *INO     *R236               
      A33(18)=
     + +   I       *R225               
      A12(19)=
     + +   NO2     *R226               
      A19(19)=
     + -   J 87                         -4 *INO2    *R237               
      A34(19)=
     + +   I       *R226               
      A13(20)=
     + +   NO2     *R229               
      A20(20)=
     + -   J 88                         -   J 89                        
      A34(20)=
     + +   IO      *R229               
      RETURN
    6 CONTINUE
C CALC. OF PRODUCTION TERM FOR CL       SPECIE NO.  21
      PR( 21)=
     + +   CLNO3   *J 10
     + +   HCL     *J 14
     + +   HOCL    *J 15
     + +   CL2     *J 38
     + +   CL2     *J 38
     + +   CL2O2   *J 41
     + +   CLNO2   *J 47
     + +   BRCL    *J 57
     + +   CLO     *J 61
     + +   CLOO    *J104
     + +   CL2O2   *J106
     + +   OH      *HCL     *R 25
     + +   CLO     *O       *R 36
     + +   CLO     *NO      *R 37
     + +   O       *HCL     *R 38
     + +   CLO     *OH      *R 75
     + +   CLO     *CLO     *R139
     + +   CLO     *CLO     *R140
     + +   IO      *CLO     *R230
     + +   O1D     *HCL     *R264
     + +   O1D     *CL2     *R266
     + +   OH      *CL2     *R267
     + +   CLOO    *R272
     + +   CL      *CCL4    *R357
     + +   CL      *CFCL3   *R358
     + +   CL      *CF2CL2  *R359
     + +   CL      *C2CL3F3 *R360
     + +   CL      *C2CL2F4 *R361
     + +   CL      *C2CLF5  *R362
     + +   CL      *CBRCLF2 *R363
     + +   CL      *CBRF3   *R364
     + +   CL      *CBR2F2  *R365
     + +   CL      *C2BR2F4 *R366
     + +   NO3     *HCL     *R379
C CALC. OF LOSS TERM FOR CL             SPECIE NO.  21
      LO( 21)=
     + +   ODP     *R  3
     + +   HO2     *R 20
     + +   H2      *R 31
     + +   H2O2    *R 32
     + +   O3      *R 34
     + +   CH4     *R 35
     + +   HOCL    *R 57
     + +   CH2O    *R 77
     + +   C2H6    *R 81
     + +   HO2     *R134
     + +   NO2     *R135
     + +   NO3     *R136
     + +   CLNO3   *R137
     + +   OCLO    *R153
     + +   H2S     *R191
     + +   DMS     *R198
     + +   CL2O2   *R207
     + +   CH3I    *R223
     + +   CH3CN   *R240
     + +   CH3BR   *R250
     + +   CH2BR2  *R253
     + +   CHBR3   *R254
     + +   CH3CL   *R255
     + +   CHCLF2  *R256
     + +   CH3CCL3 *R257
     + +   C2H3FCL2*R258
     + +   C2H3F2CL*R259
     + +   C2HF3CL2*R260
     + +   O2      *R271
     + +   CLOO    *R273
     + +   CLOO    *R274
     + +   CH2F2   *R351
     + +   CHF3    *R352
     + +   CH3CHF2 *R353
     + +   CF3CH3  *R354
     + +   CH2FCF3 *R355
     + +   CHF2CF3 *R356
     + +   CCL4    *R357
     + +   CFCL3   *R358
     + +   CF2CL2  *R359
     + +   C2CL3F3 *R360
     + +   C2CL2F4 *R361
     + +   C2CLF5  *R362
     + +   CBRCLF2 *R363
     + +   CBRF3   *R364
     + +   CBR2F2  *R365
     + +   C2BR2F4 *R366
     + +   C3HF7   *R367
     + +   C3H3F5  *R368
C CALC. OF PRODUCTION TERM FOR CLO      SPECIE NO.  22
      PR( 22)=
     + +   CLO3    *J 21
     + +   OCLO    *J 40
     + +   CLNO3   *J 60
     + +   CL2O2   *J101
     + +   CL2O2   *J101
     + +   CL2O2   *J106
     + +   CL      *O3      *R 34
     + +   OH      *HOCL    *R 54
     + +   O       *HOCL    *R 56
     + +   CLNO3   *O       *R 80
     + +   OH      *CLNO3   *R 88
     + +   CL2O2   *R103
     + +   CL2O2   *R103
     + +   HO2     *CL      *R134
     + +   CL      *NO3     *R136
     + +   O       *OCLO    *R149
     + +   NO      *OCLO    *R152
     + +   CL      *OCLO    *R153
     + +   CL      *OCLO    *R153
     + +   BR      *OCLO    *R154
     + +   O1D     *HCL     *R263
     + +   O1D     *CL2     *R266
     + +   OCLO    *O       *R268
     + +   O3      *OCLO    *R269
     + +   CL      *CLOO    *R274
     + +   CL      *CLOO    *R274
C CALC. OF LOSS TERM FOR CLO            SPECIE NO.  22
      LO( 22)=
     + +J 61
     + +   HO2     *R 17
     + +   NO2     *R 33
     + +   O       *R 36
     + +   NO      *R 37
     + +   O2D     *M       *R 59
     + +   BRO     *R 68
     + +   HO2     *R 73
     + +   OH      *R 74
     + +   OH      *R 75
     + +   HNO3    *R 94
     + +   BRO     *R 96
     + +   CLO     *R102
     + +   CLO     *R102
     + +   BRO     *R123
     + +   NO3     *R138
     + +   CLO     *R139
     + +   CLO     *R139
     + +   CLO     *R140
     + +   CLO     *R140
     + +   CLO     *R141
     + +   CLO     *R141
     + +   DMS     *R199
     + +   IO      *R230
     + +   O3      *R380
     + +   O3      *R381
C CALC. OF PRODUCTION TERM FOR HCL      SPECIE NO.  23
      PR( 23)=
     + +   HO2     *CL      *R 20
     + +   CL      *H2      *R 31
     + +   CL      *H2O2    *R 32
     + +   CL      *CH4     *R 35
     + +   CLO     *HO2     *R 73
     + +   CLO     *OH      *R 74
     + +   CL      *CH2O    *R 77
     + +   C2H6    *CL      *R 81
     + +   OH      *CLNO3   *R 89
     + +   CLO     *HNO3    *R 94
     + +   H2S     *CL      *R191
     + +   DMS     *CL      *R198
     + +   CL      *CH3I    *R223
     + +   CH3CN   *CL      *R240
     + +   CL      *CH3BR   *R250
     + +   CL      *CH2BR2  *R253
     + +   CL      *CHBR3   *R254
     + +   CL      *CH3CL   *R255
     + +   CL      *CHCLF2  *R256
     + +   CL      *CH3CCL3 *R257
     + +   CL      *C2H3FCL2*R258
     + +   CL      *C2H3F2CL*R259
     + +   CL      *C2HF3CL2*R260
     + +   O1D     *HCL     *R265
     + +   CL2              *R275
     + +   CL      *CH2F2   *R351
     + +   CL      *CHF3    *R352
     + +   CL      *CH3CHF2 *R353
     + +   CL      *CF3CH3  *R354
     + +   CL      *CH2FCF3 *R355
     + +   CL      *CHF2CF3 *R356
     + +   CL      *C3HF7   *R367
     + +   CL      *C3H3F5  *R368
C CALC. OF LOSS TERM FOR HCL            SPECIE NO.  23
      LO( 23)=
     + +J 14
     + +   OH      *R 25
     + +   O       *R 38
     + +   HOCL/(HCL+1.)    *R 90
     + +   CLNO3/(HCL+1.)   *R 91
     + +   N2O5    *R107
     + +   HOBR/(HCL+1.)    *R124
     + +   N2O5/(HCL+1.)    *R159
     + +   CLNO3/(HCL+1.)   *R161
     + +   HOCL/(HCL+1.)    *R162
     + +   N2O5/(HCL+1.)    *R164
     + +   CLNO3/(HCL+1.)   *R166
     + +   HOCL/(HCL+1.)    *R167
     + +   HOBR/(HCL+1.)    *R168
     + +   O1D     *R263
     + +   O1D     *R264
     + +   O1D     *R265
     + +   BRNO3/(HCL+1.)   *R283
     + +   BRNO3/(HCL+1.)   *R284
     + +   NO3     *R379
C CALC. OF PRODUCTION TERM FOR HOCL     SPECIE NO.  24
      PR( 24)=
     + +   HO2     *CLO     *R 17
     + +   OH      *CLNO3   *R 87
     + +   CLNO3   *R 93
     + +   OH      *OCLO    *R150
     + +   OH      *CLNO2   *R151
     + +   CLNO3   *R160
     + +   CLNO3   *R165
     + +   DMS     *CLO     *R199
     + +   OH      *CL2     *R267
     + +   OH      *CL2O2   *R297
C CALC. OF LOSS TERM FOR HOCL           SPECIE NO.  24
      LO( 24)=
     + +J 15
     + +   OH      *R 54
     + +   O       *R 56
     + +   CL      *R 57
     + +            R 90
     + +            R162
     + +            R167
     + +   HBR     *R174
     + +            R278
C CALC. OF PRODUCTION TERM FOR CLNO3    SPECIE NO.  25
      PR( 25)=
     + +   CLO     *NO2     *R 33
C CALC. OF LOSS TERM FOR CLNO3          SPECIE NO.  25
      LO( 25)=
     + +J 10
     + +J 60
     + +   O       *R 80
     + +   OH      *R 87
     + +   OH      *R 88
     + +   OH      *R 89
     + +            R 91
     + +R 93
     + +   CL      *R137
     + +R160
     + +            R161
     + +R165
     + +            R166
     + +            R279
     + +            R280
C CALC. OF PRODUCTION TERM FOR CLO3     SPECIE NO.  26
      PR( 26)=
     + +   CLO     *O2D     *M       *R 59
     + +   OCLO    *O       *R248
C CALC. OF LOSS TERM FOR CLO3           SPECIE NO.  26
      LO( 26)=
     + +J 21
C CALC. OF PRODUCTION TERM FOR OCLO     SPECIE NO.  27
      PR( 27)=
     + +   BRO     *CLO     *R 96
     + +   CLO     *CLO     *R139
     + +   CLO     *O3      *R381
C CALC. OF LOSS TERM FOR OCLO           SPECIE NO.  27
      LO( 27)=
     + +J 40
     + +   O       *R149
     + +   OH      *R150
     + +   NO      *R152
     + +   CL      *R153
     + +   BR      *R154
     + +   O       *R248
     + +   O       *R268
     + +   O3      *R269
C CALC. OF PRODUCTION TERM FOR CLOO     SPECIE NO.  28
      PR( 28)=
     + +   CL2O2   *J 41
     + +   BRO     *CLO     *R 68
     + +   CLO     *NO3     *R138
     + +   CLO     *CLO     *R140
     + +   CL      *CL2O2   *R207
     + +   CL      *O2      *R271
     + +   OH      *CL2O2   *R297
     + +   BR      *CL2O2   *R298
     + +   CLO     *O3      *R380
C CALC. OF LOSS TERM FOR CLOO           SPECIE NO.  28
      LO( 28)=
     + +J104
     + +R272
     + +   CL      *R273
     + +   CL      *R274
C CALC. OF PRODUCTION TERM FOR CL2O2    SPECIE NO.  29
      PR( 29)=
     + +   CLO     *CLO     *R102
C CALC. OF LOSS TERM FOR CL2O2          SPECIE NO.  29
      LO( 29)=
     + +J 41
     + +J101
     + +J106
     + +R103
     + +   CL      *R207
     + +   OH      *R297
     + +   BR      *R298
C CALC. OF PRODUCTION TERM FOR CL2      SPECIE NO.  30
      PR( 30)=
     + +   CL      *HOCL    *R 57
     + +   HOCL             *R 90
     + +   CLNO3            *R 91
     + +   CL      *CLNO3   *R137
     + +   CLO     *CLO     *R141
     + +   CLNO3            *R161
     + +   HOCL             *R162
     + +   CLNO3            *R166
     + +   HOCL             *R167
     + +   CL      *CL2O2   *R207
     + +   CL      *CLOO    *R273
     + +   O1D     *CL2     *R296
C CALC. OF LOSS TERM FOR CL2            SPECIE NO.  30
      LO( 30)=
     + +J 38
     + +   O1D     *R266
     + +   OH      *R267
     + +            R275
     + +   O1D     *R296
C CALC. OF PRODUCTION TERM FOR CLNO2    SPECIE NO.  31
      PR( 31)=
     + +   N2O5    *HCL     *R107
     + +   CL      *NO2     *R135
     + +   N2O5             *R159
     + +   N2O5             *R164
C CALC. OF LOSS TERM FOR CLNO2          SPECIE NO.  31
      LO( 31)=
     + +J 47
     + +   OH      *R151
* CALCULATION OF THE JACOBIAN FOR FAMILY NO.  6 INCLUDING SPECIES
*  CL      ,CLO     ,HCL     ,HOCL    ,CLNO3   ,CLO3    
*  OCLO    ,CLOO    ,CL2O2   ,CL2     ,CLNO2   ,
      A 1(21)=
     + +   HCL     *R264                +   CL2     *R266               
      A 2(21)=
     + +   CLO     *R 36                +   HCL     *R 38               
      A 4(21)=
     + -   CL      *R 34               
      A 6(21)=
     + +   HCL     *R 25                +   CLO     *R 75               
     + +   CL2     *R267               
      A 7(21)=
     + -   CL      *R 20                -   CL      *R134               
      A 8(21)=
     + -   CL      *R 32               
      A11(21)=
     + -   CL      *R 77               
      A13(21)=
     + +   CLO     *R230               
      A21(21)=
     + -   ODP     *R  3                -   HO2     *R 20               
     + -   H2      *R 31                -   H2O2    *R 32               
     + -   O3      *R 34                -   CH4     *R 35               
     + -   HOCL    *R 57                -   CH2O    *R 77               
     + -   C2H6    *R 81                -   HO2     *R134               
     + -   NO2     *R135                -   NO3     *R136               
     + -   CLNO3   *R137                -   OCLO    *R153               
     + -   H2S     *R191                -   DMS     *R198               
     + -   CL2O2   *R207                -   CH3I    *R223               
     + -   CH3CN   *R240                -   CH3BR   *R250               
     + -   CH2BR2  *R253                -   CHBR3   *R254               
     + -   CH3CL   *R255                -   CHCLF2  *R256               
     + -   CH3CCL3 *R257                -   C2H3FCL2*R258               
     + -   C2H3F2CL*R259                -   C2HF3CL2*R260               
     + -   O2      *R271                -   CLOO    *R273               
     + -   CLOO    *R274                -   CH2F2   *R351               
     + -   CHF3    *R352                -   CH3CHF2 *R353               
     + -   CF3CH3  *R354                -   CH2FCF3 *R355               
     + -   CHF2CF3 *R356                -   CCL4    *R357               
     + +   CCL4    *R357                -   CFCL3   *R358               
     + +   CFCL3   *R358                -   CF2CL2  *R359               
     + +   CF2CL2  *R359                -   C2CL3F3 *R360               
     + +   C2CL3F3 *R360                -   C2CL2F4 *R361               
     + +   C2CL2F4 *R361                -   C2CLF5  *R362               
     + +   C2CLF5  *R362                -   CBRCLF2 *R363               
     + +   CBRCLF2 *R363                -   CBRF3   *R364               
     + +   CBRF3   *R364                -   CBR2F2  *R365               
     + +   CBR2F2  *R365                -   C2BR2F4 *R366               
     + +   C2BR2F4 *R366                -   C3HF7   *R367               
     + -   C3H3F5  *R368               
      A22(21)=
     + +   J 61                         +   O       *R 36               
     + +   NO      *R 37                +   OH      *R 75               
     + +2 *CLO     *R139                +2 *CLO     *R140               
     + +   IO      *R230               
      A23(21)=
     + +   J 14                         +   OH      *R 25               
     + +   O       *R 38                +   O1D     *R264               
     + +   NO3     *R379               
      A24(21)=
     + +   J 15                         -   CL      *R 57               
      A25(21)=
     + +   J 10                         -   CL      *R137               
      A27(21)=
     + -   CL      *R153               
      A28(21)=
     + +   J104                         +   R272                        
     + -   CL      *R273                -   CL      *R274               
      A29(21)=
     + +   J 41                         +   J106                        
     + -   CL      *R207               
      A30(21)=
     + +2 *J 38                         +   O1D     *R266               
     + +   OH      *R267               
      A31(21)=
     + +   J 47                        
      A33(21)=
     + +   CLO     *R 37               
      A34(21)=
     + -   CL      *R135               
      A35(21)=
     + -   CL      *R136                +   HCL     *R379               
      A48(21)=
     + +   J 57                        
      A 1(22)=
     + +   HCL     *R263                +   CL2     *R266               
      A 2(22)=
     + -   CLO     *R 36                +   HOCL    *R 56               
     + +   CLNO3   *R 80                +   OCLO    *R149               
     + +   OCLO    *R268               
      A 3(22)=
     + -   CLO     *M       *R 59      
      A 4(22)=
     + +   CL      *R 34                +   OCLO    *R269               
     + -   CLO     *R380                -   CLO     *R381               
      A 6(22)=
     + +   HOCL    *R 54                -   CLO     *R 74               
     + -   CLO     *R 75                +   CLNO3   *R 88               
      A 7(22)=
     + -   CLO     *R 17                -   CLO     *R 73               
     + +   CL      *R134               
      A13(22)=
     + -   CLO     *R230               
      A21(22)=
     + +   O3      *R 34                +   HO2     *R134               
     + +   NO3     *R136                +2 *OCLO    *R153               
     + +2 *CLOO    *R274               
      A22(22)=
     + -   J 61                         -   HO2     *R 17               
     + -   NO2     *R 33                -   O       *R 36               
     + -   NO      *R 37                -   O2D     *M       *R 59      
     + -   BRO     *R 68                -   HO2     *R 73               
     + -   OH      *R 74                -   OH      *R 75               
     + -   HNO3    *R 94                -   BRO     *R 96               
     + -4 *CLO     *R102                -   BRO     *R123               
     + -   NO3     *R138                -4 *CLO     *R139               
     + -4 *CLO     *R140                -4 *CLO     *R141               
     + -   DMS     *R199                -   IO      *R230               
     + -   O3      *R380                -   O3      *R381               
      A23(22)=
     + +   O1D     *R263               
      A24(22)=
     + +   OH      *R 54                +   O       *R 56               
      A25(22)=
     + +   J 60                         +   O       *R 80               
     + +   OH      *R 88               
      A26(22)=
     + +   J 21                        
      A27(22)=
     + +   J 40                         +   O       *R149               
     + +   NO      *R152                +2 *CL      *R153               
     + +   BR      *R154                +   O       *R268               
     + +   O3      *R269               
      A28(22)=
     + +2 *CL      *R274               
      A29(22)=
     + +2 *J101                         +   J106                        
     + +2 *R103                        
      A30(22)=
     + +   O1D     *R266               
      A33(22)=
     + -   CLO     *R 37                +   OCLO    *R152               
      A34(22)=
     + -   CLO     *R 33               
      A35(22)=
     + +   CL      *R136                -   CLO     *R138               
      A39(22)=
     + -   CLO     *R 94               
      A41(22)=
     + +   OCLO    *R154               
      A42(22)=
     + -   CLO     *R 68                -   CLO     *R 96               
     + -   CLO     *R123               
      A 1(23)=
     + -   HCL     *R263                -   HCL     *R264               
     + -   HCL     *R265                +   HCL     *R265               
      A 2(23)=
     + -   HCL     *R 38               
      A 6(23)=
     + -   HCL     *R 25                +   CLO     *R 74               
     + +   CLNO3   *R 89               
      A 7(23)=
     + +   CL      *R 20                +   CLO     *R 73               
      A 8(23)=
     + +   CL      *R 32               
      A11(23)=
     + +   CL      *R 77               
      A21(23)=
     + +   HO2     *R 20                +   H2      *R 31               
     + +   H2O2    *R 32                +   CH4     *R 35               
     + +   CH2O    *R 77                +   C2H6    *R 81               
     + +   H2S     *R191                +   DMS     *R198               
     + +   CH3I    *R223                +   CH3CN   *R240               
     + +   CH3BR   *R250                +   CH2BR2  *R253               
     + +   CHBR3   *R254                +   CH3CL   *R255               
     + +   CHCLF2  *R256                +   CH3CCL3 *R257               
     + +   C2H3FCL2*R258                +   C2H3F2CL*R259               
     + +   C2HF3CL2*R260                +   CH2F2   *R351               
     + +   CHF3    *R352                +   CH3CHF2 *R353               
     + +   CF3CH3  *R354                +   CH2FCF3 *R355               
     + +   CHF2CF3 *R356                +   C3HF7   *R367               
     + +   C3H3F5  *R368               
      A22(23)=
     + +   HO2     *R 73                +   OH      *R 74               
     + +   HNO3    *R 94               
      A23(23)=
     + -   J 14                         -   OH      *R 25               
     + -   O       *R 38                -   HOCL/(HCL+1.)    *R 90               
     + -   CLNO3/(HCL+1.)   *R 91       -   N2O5    *R107               
     + -   HOBR/(HCL+1.)    *R124       -   N2O5/(HCL+1.)    *R159               
     + -   CLNO3/(HCL+1.)   *R161       -   HOCL/(HCL+1.)    *R162               
     + -   N2O5/(HCL+1.)    *R164       -   CLNO3/(HCL+1.)   *R166               
     + -   HOCL/(HCL+1.)    *R167       -   HOBR/(HCL+1.)    *R168               
     + -   O1D     *R263                -   O1D     *R264               
     + -   O1D     *R265                +   O1D     *R265               
     + -   BRNO3/(HCL+1.)   *R283       -   BRNO3/(HCL+1.)   *R284               
     + -   NO3     *R379               
      A24(23)=
     + -            R 90                -            R162               
     + -            R167               
      A25(23)=
     + +   OH      *R 89                -            R 91               
     + -            R161                -            R166               
      A30(23)=
     + +            R275               
      A35(23)=
     + -   HCL     *R379               
      A36(23)=
     + -   HCL     *R107                -            R159               
     + -            R164               
      A39(23)=
     + +   CLO     *R 94               
      A43(23)=
     + +   CL2/(HBR+1.)     *R275               
      A44(23)=
     + -            R124                -            R168               
      A45(23)=
     + -            R283                -            R284               
      A 2(24)=
     + -   HOCL    *R 56               
      A 6(24)=
     + -   HOCL    *R 54                +   CLNO3   *R 87               
     + +   OCLO    *R150                +   CLNO2   *R151               
     + +   CL2     *R267                +   CL2O2   *R297               
      A 7(24)=
     + +   CLO     *R 17               
      A21(24)=
     + -   HOCL    *R 57               
      A22(24)=
     + +   HO2     *R 17                +   DMS     *R199               
      A23(24)=
     + -   HOCL/(HCL+1.)    *R 90       -   HOCL/(HCL+1.)    *R162               
     + -   HOCL/(HCL+1.)    *R167               
      A24(24)=
     + -   J 15                         -   OH      *R 54               
     + -   O       *R 56                -   CL      *R 57               
     + -            R 90                -            R162               
     + -            R167                -   HBR     *R174               
     + -            R278               
      A25(24)=
     + +   OH      *R 87                +   R 93                        
     + +   R160                         +   R165                        
      A27(24)=
     + +   OH      *R150               
      A29(24)=
     + +   OH      *R297               
      A30(24)=
     + +   OH      *R267               
      A31(24)=
     + +   OH      *R151               
      A43(24)=
     + -   HOCL    *R174                -   HOCL/(HBR+1.)    *R278               
      A 2(25)=
     + -   CLNO3   *R 80               
      A 6(25)=
     + -   CLNO3   *R 87                -   CLNO3   *R 88               
     + -   CLNO3   *R 89               
      A21(25)=
     + -   CLNO3   *R137               
      A22(25)=
     + +   NO2     *R 33               
      A23(25)=
     + -   CLNO3/(HCL+1.)   *R 91       -   CLNO3/(HCL+1.)   *R161               
     + -   CLNO3/(HCL+1.)   *R166               
      A25(25)=
     + -   J 10                         -   J 60                        
     + -   O       *R 80                -   OH      *R 87               
     + -   OH      *R 88                -   OH      *R 89               
     + -            R 91                -   R 93                        
     + -   CL      *R137                -   R160                        
     + -            R161                -   R165                        
     + -            R166                -            R279               
     + -            R280               
      A34(25)=
     + +   CLO     *R 33               
      A43(25)=
     + -   CLNO3/(HBR+1.)   *R279       -   CLNO3/(HBR+1.)   *R280               
      A 2(26)=
     + +   OCLO    *R248               
      A 3(26)=
     + +   CLO     *M       *R 59      
      A22(26)=
     + +   O2D     *M       *R 59      
      A26(26)=
     + -   J 21                        
      A27(26)=
     + +   O       *R248               
      A 2(27)=
     + -   OCLO    *R149                -   OCLO    *R248               
     + -   OCLO    *R268               
      A 4(27)=
     + -   OCLO    *R269                +   CLO     *R381               
      A 6(27)=
     + -   OCLO    *R150               
      A21(27)=
     + -   OCLO    *R153               
      A22(27)=
     + +   BRO     *R 96                +2 *CLO     *R139               
     + +   O3      *R381               
      A27(27)=
     + -   J 40                         -   O       *R149               
     + -   OH      *R150                -   NO      *R152               
     + -   CL      *R153                -   BR      *R154               
     + -   O       *R248                -   O       *R268               
     + -   O3      *R269               
      A33(27)=
     + -   OCLO    *R152               
      A41(27)=
     + -   OCLO    *R154               
      A42(27)=
     + +   CLO     *R 96               
      A 4(28)=
     + +   CLO     *R380               
      A 6(28)=
     + +   CL2O2   *R297               
      A21(28)=
     + +   CL2O2   *R207                +   O2      *R271               
     + -   CLOO    *R273                -   CLOO    *R274               
      A22(28)=
     + +   BRO     *R 68                +   NO3     *R138               
     + +2 *CLO     *R140                +   O3      *R380               
      A28(28)=
     + -   J104                         -   R272                        
     + -   CL      *R273                -   CL      *R274               
      A29(28)=
     + +   J 41                         +   CL      *R207               
     + +   OH      *R297                +   BR      *R298               
      A35(28)=
     + +   CLO     *R138               
      A41(28)=
     + +   CL2O2   *R298               
      A42(28)=
     + +   CLO     *R 68               
      A 6(29)=
     + -   CL2O2   *R297               
      A21(29)=
     + -   CL2O2   *R207               
      A22(29)=
     + +2 *CLO     *R102               
      A29(29)=
     + -   J 41                         -   J101                        
     + -   J106                         -   R103                        
     + -   CL      *R207                -   OH      *R297               
     + -   BR      *R298               
      A41(29)=
     + -   CL2O2   *R298               
      A 1(30)=
     + -   CL2     *R266                -   CL2     *R296               
     + +   CL2     *R296               
      A 6(30)=
     + -   CL2     *R267               
      A21(30)=
     + +   HOCL    *R 57                +   CLNO3   *R137               
     + +   CL2O2   *R207                +   CLOO    *R273               
      A22(30)=
     + +2 *CLO     *R141               
      A23(30)=
     + +   HOCL/(HCL+1.)    *R 90       +   CLNO3/(HCL+1.)   *R 91               
     + +   CLNO3/(HCL+1.)   *R161       +   HOCL/(HCL+1.)    *R162               
     + +   CLNO3/(HCL+1.)   *R166       +   HOCL/(HCL+1.)    *R167               
      A24(30)=
     + +   CL      *R 57                +            R 90               
     + +            R162                +            R167               
      A25(30)=
     + +            R 91                +   CL      *R137               
     + +            R161                +            R166               
      A28(30)=
     + +   CL      *R273               
      A29(30)=
     + +   CL      *R207               
      A30(30)=
     + -   J 38                         -   O1D     *R266               
     + -   OH      *R267                -            R275               
     + -   O1D     *R296                +   O1D     *R296               
      A43(30)=
     + -   CL2/(HBR+1.)     *R275               
      A 6(31)=
     + -   CLNO2   *R151               
      A21(31)=
     + +   NO2     *R135               
      A23(31)=
     + +   N2O5    *R107                +   N2O5/(HCL+1.)    *R159               
     + +   N2O5/(HCL+1.)    *R164               
      A31(31)=
     + -   J 47                         -   OH      *R151               
      A34(31)=
     + +   CL      *R135               
      A36(31)=
     + +   HCL     *R107                +            R159               
     + +            R164               
      RETURN
    7 CONTINUE
C CALC. OF PRODUCTION TERM FOR N        SPECIE NO.  32
      PR( 32)=
     + +   NO      *J  5
C CALC. OF LOSS TERM FOR N              SPECIE NO.  32
      LO( 32)=
     + +   O2      *R  4
     + +   NO      *R  5
     + +   NO2     *R183
     + +   OH      *R209
     + +   O3      *R376
C CALC. OF PRODUCTION TERM FOR NO       SPECIE NO.  33
      PR( 33)=
     + +   NO2     *J  6
     + +   HONO    *J 58
     + +   NO3     *J 59
     + +   INO     *J 86
     + +   N2O5    *J 98
     + +   N       *O2      *R  4
     + +   NO2     *O       *R 45
     + +   H       *NO2     *R133
     + +   N       *OH      *R209
     + +2 *INO     *INO     *R236
     + +   O2D     *NO      *R288
     + +   N       *O3      *R376
     + +   NO2     *NO3     *R377
C CALC. OF LOSS TERM FOR NO             SPECIE NO.  33
      LO( 33)=
     + +J  5
     + +   N       *R  5
     + +   HO2     *R 16
     + +   CLO     *R 37
     + +   O3      *R 43
     + +   BRO     *R 63
     + +   MEO2    *R 85
     + +   ETHO2   *R108
     + +   CH3CO3  *R113
     + +   OH      *R125
     + +   NO3     *R132
     + +   OCLO    *R152
     + +   O       *R156
     + +   I       *R225
     + +   IO      *R228
     + +   O2D     *R288
C CALC. OF PRODUCTION TERM FOR NO2      SPECIE NO.  34
      PR( 34)=
     + +   NO3     *J  7
     + +   HNO3    *J  8
     + +   N2O5    *J  9
     + +   HO2NO2  *J 32
     + +   CLNO2   *J 47
     + +   PAN     *J 55
     + +   CLNO3   *J 60
     + +   BRNO3   *J 62
     + +   INO2    *J 87
     + +   INO3    *J 88
     + +   BRONO   *J103
     + +   HOONO   *J107
     + +   HO2     *NO      *R 16
     + +   CLO     *NO      *R 37
     + +   N2O5    *R 42
     + +   NO      *O3      *R 43
     + +   HO2NO2  *R 53
     + +   BRO     *NO      *R 63
     + +   HO2NO2  *OH      *R 78
     + +   HO2NO2  *O       *R 79
     + +   MEO2    *NO      *R 85
     + +   OH      *CLNO3   *R 89
     + +   CLO     *HNO3    *R 94
     + +   ETHO2   *NO      *R108
     + +   CH3CO3  *NO      *R113
     + +   PAN     *R116
     + +   PAN     *OH      *R117
     + +   ETHO2NO2*R121
     + +   OH      *HONO    *R126
     + +   NO3     *O       *R131
     + +   NO3     *NO      *R132
     + +   NO3     *NO      *R132
     + +   CL      *NO3     *R136
     + +   CLO     *NO3     *R138
     + +2 *NO3     *NO3     *R148
     + +   OH      *CLNO2   *R151
     + +   NO      *OCLO    *R152
     + +   O       *NO      *R156
     + +   IO      *NO      *R228
     + +2 *INO2    *INO2    *R237
     + +   HOONO   *R370
     + +   OH      *NO3     *R373
     + +   NO2     *NO3     *R377
     + +   BR      *NO3     *R382
C CALC. OF LOSS TERM FOR NO2            SPECIE NO.  34
      LO( 34)=
     + +J  6
     + +   OH      *R 27
     + +   CLO     *R 33
     + +   NO3     *R 41
     + +   O3      *R 44
     + +   O       *R 45
     + +   HO2     *R 52
     + +   BRO     *R 67
     + +   CH3CO3  *R115
     + +   ETHO2   *R120
     + +   H       *R133
     + +   CL      *R135
     + +   O       *R157
     + +   N       *R183
     + +   I       *R226
     + +   IO      *R229
     + +   BR      *R270
     + +   OH      *R369
     + +   H2O     *R374
     + +   HO2     *R375
     + +   NO3     *R377
C CALC. OF PRODUCTION TERM FOR NO3      SPECIE NO.  35
      PR( 35)=
     + +   N2O5    *J  9
     + +   CLNO3   *J 10
     + +   BRNO3   *J 23
     + +   HO2NO2  *J 33
     + +   INO3    *J 89
     + +   N2O5    *J 98
     + +   OH      *HNO3    *R 28
     + +   N2O5    *R 42
     + +   NO2     *O3      *R 44
     + +   CLNO3   *O       *R 80
     + +   OH      *CLNO3   *R 87
     + +   CL      *CLNO3   *R137
     + +   O       *NO2     *R157
     + +   O       *BRNO3   *R172
     + +   HOONO   *OH      *R371
     + +   O       *HNO3    *R372
C CALC. OF LOSS TERM FOR NO3            SPECIE NO.  35
      LO( 35)=
     + +J  7
     + +J 59
     + +   NO2     *R 41
     + +   O       *R131
     + +   NO      *R132
     + +   CL      *R136
     + +   CLO     *R138
     + +   NO3     *R148
     + +   NO3     *R148
     + +   DMS     *R197
     + +   HI      *R222
     + +R285
     + +   OH      *R373
     + +   NO2     *R377
     + +   HCL     *R379
     + +   BR      *R382
C CALC. OF PRODUCTION TERM FOR N2O5     SPECIE NO.  36
      PR( 36)=
     + +   NO2     *NO3     *R 41
C CALC. OF LOSS TERM FOR N2O5           SPECIE NO.  36
      LO( 36)=
     + +J  9
     + +J 98
     + +R 42
     + +   H2O     *R 83
     + +R 92
     + +   HCL     *R107
     + +R158
     + +            R159
     + +R163
     + +            R164
     + +            R276
     + +            R277
C CALC. OF PRODUCTION TERM FOR HO2NO2   SPECIE NO.  37
      PR( 37)=
     + +   HO2     *NO2     *R 52
C CALC. OF LOSS TERM FOR HO2NO2         SPECIE NO.  37
      LO( 37)=
     + +J 32
     + +J 33
     + +R 53
     + +   OH      *R 78
     + +   O       *R 79
C CALC. OF PRODUCTION TERM FOR HONO     SPECIE NO.  38
      PR( 38)=
     + +   HNO3    *J 99
     + +   HNO3    *J100
     + +   OH      *NO      *R125
     + +   H2O     *NO2     *R374
     + +   HO2     *NO2     *R375
C CALC. OF LOSS TERM FOR HONO           SPECIE NO.  38
      LO( 38)=
     + +J 58
     + +   OH      *R126
     + +   O3      *R378
C CALC. OF PRODUCTION TERM FOR HNO3     SPECIE NO.  39
      PR( 39)=
     + +   OH      *NO2     *R 27
     + +   N2O5    *H2O     *R 83
     + +   N2O5    *H2O     *R 83
     + +   OH      *CLNO3   *R 88
     + +   CLNO3            *R 91
     + +   N2O5    *R 92
     + +   N2O5    *R 92
     + +   CLNO3   *R 93
     + +   N2O5    *HCL     *R107
     + +   BRNO3   *R122
     + +   N2O5    *R158
     + +   N2O5    *R158
     + +   N2O5             *R159
     + +   CLNO3   *R160
     + +   CLNO3            *R161
     + +   N2O5    *R163
     + +   N2O5    *R163
     + +   N2O5             *R164
     + +   CLNO3   *R165
     + +   CLNO3            *R166
     + +   BRNO3   *R169
     + +   NO3     *HI      *R222
     + +   N2O5             *R276
     + +   N2O5             *R277
     + +   CLNO3            *R279
     + +   CLNO3            *R280
     + +   BRNO3            *R283
     + +   BRNO3            *R284
     + +   NO3     *R285
     + +   O3      *HONO    *R378
     + +   NO3     *HCL     *R379
C CALC. OF LOSS TERM FOR HNO3           SPECIE NO.  39
      LO( 39)=
     + +J  8
     + +J 25
     + +J 99
     + +J100
     + +   OH      *R 28
     + +   CLO     *R 94
     + +   O       *R372
C CALC. OF PRODUCTION TERM FOR HOONO    SPECIE NO.  40
      PR( 40)=
     + +   OH      *NO2     *R369
C CALC. OF LOSS TERM FOR HOONO          SPECIE NO.  40
      LO( 40)=
     + +J107
     + +J108
     + +R370
     + +   OH      *R371
* CALCULATION OF THE JACOBIAN FOR FAMILY NO.  7 INCLUDING SPECIES
*  N       ,NO      ,NO2     ,NO3     ,N2O5    ,HO2NO2  
*  HONO    ,HNO3    ,HOONO   ,
      A 4(32)=
     + -   N       *R376               
      A 6(32)=
     + -   N       *R209               
      A32(32)=
     + -   O2      *R  4                -   NO      *R  5               
     + -   NO2     *R183                -   OH      *R209               
     + -   O3      *R376               
      A33(32)=
     + +   J  5                         -   N       *R  5               
      A34(32)=
     + -   N       *R183               
      A 2(33)=
     + +   NO2     *R 45                -   NO      *R156               
      A 3(33)=
     + -   NO      *R288                +   NO      *R288               
      A 4(33)=
     + -   NO      *R 43                +   N       *R376               
      A 5(33)=
     + +   NO2     *R133               
      A 6(33)=
     + -   NO      *R125                +   N       *R209               
      A 7(33)=
     + -   NO      *R 16               
      A 9(33)=
     + -   NO      *R 85               
      A12(33)=
     + -   NO      *R225               
      A13(33)=
     + -   NO      *R228               
      A18(33)=
     + +   J 86                         +4 *INO     *R236               
      A22(33)=
     + -   NO      *R 37               
      A27(33)=
     + -   NO      *R152               
      A32(33)=
     + +   O2      *R  4                -   NO      *R  5               
     + +   OH      *R209                +   O3      *R376               
      A33(33)=
     + -   J  5                         -   N       *R  5               
     + -   HO2     *R 16                -   CLO     *R 37               
     + -   O3      *R 43                -   BRO     *R 63               
     + -   MEO2    *R 85                -   ETHO2   *R108               
     + -   CH3CO3  *R113                -   OH      *R125               
     + -   NO3     *R132                -   OCLO    *R152               
     + -   O       *R156                -   I       *R225               
     + -   IO      *R228                -   O2D     *R288               
     + +   O2D     *R288               
      A34(33)=
     + +   J  6                         +   O       *R 45               
     + +   H       *R133                +   NO3     *R377               
      A35(33)=
     + +   J 59                         -   NO      *R132               
     + +   NO2     *R377               
      A36(33)=
     + +   J 98                        
      A38(33)=
     + +   J 58                        
      A42(33)=
     + -   NO      *R 63               
      A49(33)=
     + -   NO      *R108               
      A52(33)=
     + -   NO      *R113               
      A 2(34)=
     + -   NO2     *R 45                +   HO2NO2  *R 79               
     + +   NO3     *R131                +   NO      *R156               
     + -   NO2     *R157               
      A 4(34)=
     + +   NO      *R 43                -   NO2     *R 44               
      A 5(34)=
     + -   NO2     *R133               
      A 6(34)=
     + -   NO2     *R 27                +   HO2NO2  *R 78               
     + +   CLNO3   *R 89                +   PAN     *R117               
     + +   HONO    *R126                +   CLNO2   *R151               
     + -   NO2     *R369                +   NO3     *R373               
      A 7(34)=
     + +   NO      *R 16                -   NO2     *R 52               
     + -   NO2     *R375               
      A 9(34)=
     + +   NO      *R 85               
      A12(34)=
     + -   NO2     *R226               
      A13(34)=
     + +   NO      *R228                -   NO2     *R229               
      A19(34)=
     + +   J 87                         +4 *INO2    *R237               
      A20(34)=
     + +   J 88                        
      A21(34)=
     + -   NO2     *R135                +   NO3     *R136               
      A22(34)=
     + -   NO2     *R 33                +   NO      *R 37               
     + +   HNO3    *R 94                +   NO3     *R138               
      A25(34)=
     + +   J 60                         +   OH      *R 89               
      A27(34)=
     + +   NO      *R152               
      A31(34)=
     + +   J 47                         +   OH      *R151               
      A32(34)=
     + -   NO2     *R183               
      A33(34)=
     + +   HO2     *R 16                +   CLO     *R 37               
     + +   O3      *R 43                +   BRO     *R 63               
     + +   MEO2    *R 85                +   ETHO2   *R108               
     + +   CH3CO3  *R113                +2 *NO3     *R132               
     + +   OCLO    *R152                +   O       *R156               
     + +   IO      *R228               
      A34(34)=
     + -   J  6                         -   OH      *R 27               
     + -   CLO     *R 33                -   NO3     *R 41               
     + -   O3      *R 44                -   O       *R 45               
     + -   HO2     *R 52                -   BRO     *R 67               
     + -   CH3CO3  *R115                -   ETHO2   *R120               
     + -   H       *R133                -   CL      *R135               
     + -   O       *R157                -   N       *R183               
     + -   I       *R226                -   IO      *R229               
     + -   BR      *R270                -   OH      *R369               
     + -   H2O     *R374                -   HO2     *R375               
     + -   NO3     *R377                +   NO3     *R377               
      A35(34)=
     + +   J  7                         -   NO2     *R 41               
     + +   O       *R131                +2 *NO      *R132               
     + +   CL      *R136                +   CLO     *R138               
     + +4 *NO3     *R148                +   OH      *R373               
     + -   NO2     *R377                +   NO2     *R377               
     + +   BR      *R382               
      A36(34)=
     + +   J  9                         +   R 42                        
      A37(34)=
     + +   J 32                         +   R 53                        
     + +   OH      *R 78                +   O       *R 79               
      A38(34)=
     + +   OH      *R126               
      A39(34)=
     + +   J  8                         +   CLO     *R 94               
      A40(34)=
     + +   J107                         +   R370                        
      A41(34)=
     + -   NO2     *R270                +   NO3     *R382               
      A42(34)=
     + +   NO      *R 63                -   NO2     *R 67               
      A45(34)=
     + +   J 62                        
      A46(34)=
     + +   J103                        
      A49(34)=
     + +   NO      *R108                -   NO2     *R120               
      A52(34)=
     + +   NO      *R113                -   NO2     *R115               
      A54(34)=
     + +   R121                        
      A 2(35)=
     + +   CLNO3   *R 80                -   NO3     *R131               
     + +   NO2     *R157                +   BRNO3   *R172               
     + +   HNO3    *R372               
      A 4(35)=
     + +   NO2     *R 44               
      A 6(35)=
     + +   HNO3    *R 28                +   CLNO3   *R 87               
     + +   HOONO   *R371                -   NO3     *R373               
      A16(35)=
     + -   NO3     *R222               
      A20(35)=
     + +   J 89                        
      A21(35)=
     + -   NO3     *R136                +   CLNO3   *R137               
      A22(35)=
     + -   NO3     *R138               
      A23(35)=
     + -   NO3     *R379               
      A25(35)=
     + +   J 10                         +   O       *R 80               
     + +   OH      *R 87                +   CL      *R137               
      A33(35)=
     + -   NO3     *R132               
      A34(35)=
     + -   NO3     *R 41                +   O3      *R 44               
     + +   O       *R157                -   NO3     *R377               
      A35(35)=
     + -   J  7                         -   J 59                        
     + -   NO2     *R 41                -   O       *R131               
     + -   NO      *R132                -   CL      *R136               
     + -   CLO     *R138                -4 *NO3     *R148               
     + -   DMS     *R197                -   HI      *R222               
     + -   R285                         -   OH      *R373               
     + -   NO2     *R377                -   HCL     *R379               
     + -   BR      *R382               
      A36(35)=
     + +   J  9                         +   J 98                        
     + +   R 42                        
      A37(35)=
     + +   J 33                        
      A39(35)=
     + +   OH      *R 28                +   O       *R372               
      A40(35)=
     + +   OH      *R371               
      A41(35)=
     + -   NO3     *R382               
      A45(35)=
     + +   J 23                         +   O       *R172               
      A23(36)=
     + -   N2O5    *R107                -   N2O5/(HCL+1.)    *R159               
     + -   N2O5/(HCL+1.)    *R164               
      A34(36)=
     + +   NO3     *R 41               
      A35(36)=
     + +   NO2     *R 41               
      A36(36)=
     + -   J  9                         -   J 98                        
     + -   R 42                         -   H2O     *R 83               
     + -   R 92                         -   HCL     *R107               
     + -   R158                         -            R159               
     + -   R163                         -            R164               
     + -            R276                -            R277               
      A43(36)=
     + -   N2O5/(HBR+1.)    *R276       -   N2O5 /(HBR+1.)   *R277               
      A 2(37)=
     + -   HO2NO2  *R 79               
      A 6(37)=
     + -   HO2NO2  *R 78               
      A 7(37)=
     + +   NO2     *R 52               
      A34(37)=
     + +   HO2     *R 52               
      A37(37)=
     + -   J 32                         -   J 33                        
     + -   R 53                         -   OH      *R 78               
     + -   O       *R 79               
      A 4(38)=
     + -   HONO    *R378               
      A 6(38)=
     + +   NO      *R125                -   HONO    *R126               
      A 7(38)=
     + +   NO2     *R375               
      A33(38)=
     + +   OH      *R125               
      A34(38)=
     + +   H2O     *R374                +   HO2     *R375               
      A38(38)=
     + -   J 58                         -   OH      *R126               
     + -   O3      *R378               
      A39(38)=
     + +   J 99                         +   J100                        
      A 2(39)=
     + -   HNO3    *R372               
      A 4(39)=
     + +   HONO    *R378               
      A 6(39)=
     + +   NO2     *R 27                -   HNO3    *R 28               
     + +   CLNO3   *R 88               
      A16(39)=
     + +   NO3     *R222               
      A22(39)=
     + -   HNO3    *R 94               
      A23(39)=
     + +   CLNO3/(HCL+1.)   *R 91                +   N2O5    *R107               
     + +   N2O5/(HCL+1.)    *R159       +   CLNO3/(HCL+1.)   *R161               
     + +   N2O5/(HCL+1.)    *R164       +   CLNO3/(HCL+1.)   *R166               
     + +   BRNO3/(HCL+1.)   *R283       +   BRNO3/(HCL+1.)   *R284               
     + +   NO3     *R379               
      A25(39)=
     + +   OH      *R 88                +            R 91               
     + +   R 93                         +   R160                        
     + +            R161                +   R165                        
     + +            R166                +            R279               
     + +            R280               
      A34(39)=
     + +   OH      *R 27               
      A35(39)=
     + +   HI      *R222                +   R285                        
     + +   HCL     *R379               
      A36(39)=
     + +2 *H2O     *R 83                +2 *R 92                        
     + +   HCL     *R107                +2 *R158                        
     + +            R159                +2 *R163                        
     + +            R164                +            R276               
     + +            R277               
      A38(39)=
     + +   O3      *R378               
      A39(39)=
     + -   J  8                         -   J 25                        
     + -   J 99                         -   J100                        
     + -   OH      *R 28                -   CLO     *R 94               
     + -   O       *R372               
      A43(39)=
     + +   N2O5/(HBR+1.)    *R276       +   N2O5/(HBR+1.)    *R277               
     + +   CLNO3/(HBR+1.)   *R279       +   CLNO3/(HBR+1.)   *R280               
      A45(39)=
     + +   R122                         +   R169                        
     + +            R283                +            R284               
      A 6(40)=
     + +   NO2     *R369                -   HOONO   *R371               
      A34(40)=
     + +   OH      *R369               
      A40(40)=
     + -   J107                         -   J108                        
     + -   R370                         -   OH      *R371               
      RETURN
    8 CONTINUE
C CALC. OF PRODUCTION TERM FOR BR       SPECIE NO.  41
      PR( 41)=
     + +   HOBR    *J 22
     + +   BRNO3   *J 23
     + +   BRO     *J 37
     + +   BR2     *J 39
     + +   BR2     *J 39
     + +   BRCL    *J 57
     + +   HBR     *J102
     + +   BRONO   *J103
     + +   BRO     *O       *R 61
     + +2 *BRO     *BRO     *R 62
     + +   BRO     *NO      *R 63
     + +   BRO     *CLO     *R 68
     + +   OH      *HBR     *R 70
     + +   BRO     *CLO     *R 96
     + +   BRO     *O3      *R 98
     + +   BRO     *OH      *R 99
     + +   O       *HBR     *R101
     + +   O1D     *HBR     *R173
     + +   I       *BRO     *R227
     + +   IO      *BRO     *R231
     + +   OH      *BR2     *R246
C CALC. OF LOSS TERM FOR BR             SPECIE NO.  41
      LO( 41)=
     + +   O3      *R 60
     + +   HO2     *R 64
     + +   H2O2    *R 69
     + +   CH2O    *R 95
     + +   OCLO    *R154
     + +   H2S     *R192
     + +   DMS     *R200
     + +   NO2     *R270
     + +   CL2O2   *R298
     + +   NO3     *R382
C CALC. OF PRODUCTION TERM FOR BRO      SPECIE NO.  42
      PR( 42)=
     + +   BRNO3   *J 62
     + +   BR      *O3      *R 60
     + +   BR      *OCLO    *R154
     + +   O       *HOBR    *R170
     + +   O       *BRNO3   *R172
     + +   O1D     *HBR     *R295
     + +   BR      *NO3     *R382
C CALC. OF LOSS TERM FOR BRO            SPECIE NO.  42
      LO( 42)=
     + +J 37
     + +   O       *R 61
     + +   BRO     *R 62
     + +   BRO     *R 62
     + +   NO      *R 63
     + +   HO2     *R 65
     + +   HO2     *R 66
     + +   NO2     *R 67
     + +   CLO     *R 68
     + +   CLO     *R 96
     + +   BRO     *R 97
     + +   BRO     *R 97
     + +   O3      *R 98
     + +   OH      *R 99
     + +   CLO     *R123
     + +   OH      *R171
     + +   DMS     *R201
     + +   I       *R227
     + +   IO      *R231
C CALC. OF PRODUCTION TERM FOR HBR      SPECIE NO.  43
      PR( 43)=
     + +   BR      *HO2     *R 64
     + +   BRO     *HO2     *R 66
     + +   BR      *H2O2    *R 69
     + +   BR      *CH2O    *R 95
     + +   BRO     *OH      *R171
     + +   H2S     *BR      *R192
     + +   DMS     *BR      *R200
     + +   O1D     *HBR     *R294
C CALC. OF LOSS TERM FOR HBR            SPECIE NO.  43
      LO( 43)=
     + +J102
     + +   OH      *R 70
     + +   O       *R101
     + +   O1D     *R173
     + +   HOCL    *R174
     + +   CL2/(HBR+1.)     *R275
     + +   N2O5/(HBR+1.)    *R276
     + +   N2O5/(HBR+1.)    *R277
     + +   HOCL/(HBR+1.)    *R278
     + +   CLNO3/(HBR+1.)   *R279
     + +   CLNO3/(HBR+1.)   *R280
     + +   HOBR/(HBR+1.)    *R281
     + +   HOBR/(HBR+1.)    *R282
     + +   O1D     *R294
     + +   O1D     *R295
C CALC. OF PRODUCTION TERM FOR HOBR     SPECIE NO.  44
      PR( 44)=
     + +   BRO     *HO2     *R 65
     + +   BRNO3   *R122
     + +   BRNO3   *R169
     + +   DMS     *BRO     *R201
     + +   OH      *BR2     *R246
C CALC. OF LOSS TERM FOR HOBR           SPECIE NO.  44
      LO( 44)=
     + +J 22
     + +            R124
     + +            R168
     + +   O       *R170
     + +            R281
     + +            R282
C CALC. OF PRODUCTION TERM FOR BRNO3    SPECIE NO.  45
      PR( 45)=
     + +   BRO     *NO2     *R 67
C CALC. OF LOSS TERM FOR BRNO3          SPECIE NO.  45
      LO( 45)=
     + +J 23
     + +J 62
     + +R122
     + +R169
     + +   O       *R172
     + +            R283
     + +            R284
C CALC. OF PRODUCTION TERM FOR BRONO    SPECIE NO.  46
      PR( 46)=
     + +   BR      *NO2     *R270
     + +   N2O5             *R276
     + +   N2O5             *R277
C CALC. OF LOSS TERM FOR BRONO          SPECIE NO.  46
      LO( 46)=
     + +J103
C CALC. OF PRODUCTION TERM FOR BR2      SPECIE NO.  47
      PR( 47)=
     + +   BRO     *BRO     *R 97
     + +   HOBR             *R281
     + +   HOBR             *R282
C CALC. OF LOSS TERM FOR BR2            SPECIE NO.  47
      LO( 47)=
     + +J 39
     + +   OH      *R246
C CALC. OF PRODUCTION TERM FOR BRCL     SPECIE NO.  48
      PR( 48)=
     + +   BRO     *CLO     *R123
     + +   HOBR             *R124
     + +   HOBR             *R168
     + +   HOCL    *HBR     *R174
     + +   CL2              *R275
     + +   HOCL             *R278
     + +   CLNO3            *R279
     + +   CLNO3            *R280
     + +   BRNO3            *R283
     + +   BRNO3            *R284
     + +   BR      *CL2O2   *R298
C CALC. OF LOSS TERM FOR BRCL           SPECIE NO.  48
      LO( 48)=
     + +J 57
* CALCULATION OF THE JACOBIAN FOR FAMILY NO.  8 INCLUDING SPECIES
*  BR      ,BRO     ,HBR     ,HOBR    ,BRNO3   ,BRONO   
*  BR2     ,BRCL    ,
      A 1(41)=
     + +   HBR     *R173               
      A 2(41)=
     + +   BRO     *R 61                +   HBR     *R101               
      A 4(41)=
     + -   BR      *R 60                +   BRO     *R 98               
      A 6(41)=
     + +   HBR     *R 70                +   BRO     *R 99               
     + +   BR2     *R246               
      A 7(41)=
     + -   BR      *R 64               
      A 8(41)=
     + -   BR      *R 69               
      A11(41)=
     + -   BR      *R 95               
      A12(41)=
     + +   BRO     *R227               
      A13(41)=
     + +   BRO     *R231               
      A22(41)=
     + +   BRO     *R 68                +   BRO     *R 96               
      A27(41)=
     + -   BR      *R154               
      A29(41)=
     + -   BR      *R298               
      A33(41)=
     + +   BRO     *R 63               
      A34(41)=
     + -   BR      *R270               
      A35(41)=
     + -   BR      *R382               
      A41(41)=
     + -   O3      *R 60                -   HO2     *R 64               
     + -   H2O2    *R 69                -   CH2O    *R 95               
     + -   OCLO    *R154                -   H2S     *R192               
     + -   DMS     *R200                -   NO2     *R270               
     + -   CL2O2   *R298                -   NO3     *R382               
      A42(41)=
     + +   J 37                         +   O       *R 61               
     + +4 *BRO     *R 62                +   NO      *R 63               
     + +   CLO     *R 68                +   CLO     *R 96               
     + +   O3      *R 98                +   OH      *R 99               
     + +   I       *R227                +   IO      *R231               
      A43(41)=
     + +   J102                         +   OH      *R 70               
     + +   O       *R101                +   O1D     *R173               
      A44(41)=
     + +   J 22                        
      A45(41)=
     + +   J 23                        
      A46(41)=
     + +   J103                        
      A47(41)=
     + +2 *J 39                         +   OH      *R246               
      A48(41)=
     + +   J 57                        
      A 1(42)=
     + +   HBR     *R295               
      A 2(42)=
     + -   BRO     *R 61                +   HOBR    *R170               
     + +   BRNO3   *R172               
      A 4(42)=
     + +   BR      *R 60                -   BRO     *R 98               
      A 6(42)=
     + -   BRO     *R 99                -   BRO     *R171               
      A 7(42)=
     + -   BRO     *R 65                -   BRO     *R 66               
      A12(42)=
     + -   BRO     *R227               
      A13(42)=
     + -   BRO     *R231               
      A22(42)=
     + -   BRO     *R 68                -   BRO     *R 96               
     + -   BRO     *R123               
      A27(42)=
     + +   BR      *R154               
      A33(42)=
     + -   BRO     *R 63               
      A34(42)=
     + -   BRO     *R 67               
      A35(42)=
     + +   BR      *R382               
      A41(42)=
     + +   O3      *R 60                +   OCLO    *R154               
     + +   NO3     *R382               
      A42(42)=
     + -   J 37                         -   O       *R 61               
     + -4 *BRO     *R 62                -   NO      *R 63               
     + -   HO2     *R 65                -   HO2     *R 66               
     + -   NO2     *R 67                -   CLO     *R 68               
     + -   CLO     *R 96                -4 *BRO     *R 97               
     + -   O3      *R 98                -   OH      *R 99               
     + -   CLO     *R123                -   OH      *R171               
     + -   DMS     *R201                -   I       *R227               
     + -   IO      *R231               
      A43(42)=
     + +   O1D     *R295               
      A44(42)=
     + +   O       *R170               
      A45(42)=
     + +   J 62                         +   O       *R172               
      A 1(43)=
     + -   HBR     *R173                -   HBR     *R294               
     + +   HBR     *R294                -   HBR     *R295               
      A 2(43)=
     + -   HBR     *R101               
      A 6(43)=
     + -   HBR     *R 70                +   BRO     *R171               
      A 7(43)=
     + +   BR      *R 64                +   BRO     *R 66               
      A 8(43)=
     + +   BR      *R 69               
      A11(43)=
     + +   BR      *R 95               
      A24(43)=
     + -   HBR     *R174                -            R278               
      A25(43)=
     + -            R279                -            R280               
      A30(43)=
     + -            R275               
      A36(43)=
     + -            R276                -            R277               
      A41(43)=
     + +   HO2     *R 64                +   H2O2    *R 69               
     + +   CH2O    *R 95                +   H2S     *R192               
     + +   DMS     *R200               
      A42(43)=
     + +   HO2     *R 66                +   OH      *R171               
      A43(43)=
     + -   J102                         -   OH      *R 70               
     + -   O       *R101                -   O1D     *R173               
     + -   HOCL    *R174                -   CL2/(HBR+1.)     *R275               
     + -   N2O5/(HBR+1.)    *R276       -   N2O5/(HBR+1.)    *R277               
     + -   HOCL/(HBR+1.)    *R278       -   CLNO3/(HBR+1.)   *R279               
     + -   CLNO3/(HBR+1.)   *R280       -   HOBR/(HBR+1.)    *R281               
     + -   HOBR/(HBR+1.)    *R282       -   O1D     *R294               
     + +   O1D     *R294                -   O1D     *R295               
      A44(43)=
     + -            R281                -            R282               
      A 2(44)=
     + -   HOBR    *R170               
      A 6(44)=
     + +   BR2     *R246               
      A 7(44)=
     + +   BRO     *R 65               
      A23(44)=
     + -   HOBR/(HCL+1.)    *R124       -   HOBR/(HCL+1.)    *R168               
      A42(44)=
     + +   HO2     *R 65                +   DMS     *R201               
      A43(44)=
     + -   HOBR/(HBR+1.)    *R281       -   HOBR/(HBR+1.)    *R282               
      A44(44)=
     + -   J 22                         -            R124               
     + -            R168                -   O       *R170               
     + -            R281                -            R282               
      A45(44)=
     + +   R122                         +   R169                        
      A47(44)=
     + +   OH      *R246               
      A 2(45)=
     + -   BRNO3   *R172               
      A23(45)=
     + -   BRNO3/(HCL+1.)   *R283       -   BRNO3/(HCL+1.)   *R284               
      A34(45)=
     + +   BRO     *R 67               
      A42(45)=
     + +   NO2     *R 67               
      A45(45)=
     + -   J 23                         -   J 62                        
     + -   R122                         -   R169                        
     + -   O       *R172                -            R283               
     + -            R284               
      A34(46)=
     + +   BR      *R270               
      A36(46)=
     + +            R276                +            R277               
      A41(46)=
     + +   NO2     *R270               
      A43(46)=
     + +   N2O5/(HBR+1.)    *R276       +   N2O5/(HBR+1.)    *R277               
      A46(46)=
     + -   J103                        
      A 6(47)=
     + -   BR2     *R246               
      A42(47)=
     + +2 *BRO     *R 97               
      A43(47)=
     + +   HOBR/(HBR+1.)    *R281       +   HOBR/(HBR+1.)    *R282               
      A44(47)=
     + +            R281                +            R282               
      A47(47)=
     + -   J 39                         -   OH      *R246               
      A22(48)=
     + +   BRO     *R123               
      A23(48)=
     + +   HOBR/(HCL+1.)    *R124       +   HOBR/(HCL+1.)    *R168               
     + +   BRNO3/(HCL+1.)   *R283       +   BRNO3/(HCL+1.)   *R284               
      A24(48)=
     + +   HBR     *R174                +            R278               
      A25(48)=
     + +            R279                +            R280               
      A29(48)=
     + +   BR      *R298               
      A30(48)=
     + +            R275               
      A41(48)=
     + +   CL2O2   *R298               
      A42(48)=
     + +   CLO     *R123               
      A43(48)=
     + +   HOCL    *R174                +   CL2/(HBR+1.)     *R275               
     + +   HOCL/(HBR+1.)    *R278       +   CLNO3/(HBR+1.)   *R279               
     + +   CLNO3/(HBR+1.)   *R280               
      A44(48)=
     + +            R124                +            R168               
      A45(48)=
     + +            R283                +            R284               
      A48(48)=
     + -   J 57                        
      RETURN
    9 CONTINUE
C CALC. OF PRODUCTION TERM FOR ETHO2    SPECIE NO.  49
      PR( 49)=
     + +   C2H6    *CL      *R 81
     + +   C2H6    *OH      *R 82
     + +   ETHOH   *OH      *R110
     + +   ETHO2NO2*R121
C CALC. OF LOSS TERM FOR ETHO2          SPECIE NO.  49
      LO( 49)=
     + +   NO      *R108
     + +   HO2     *R109
     + +   NO2     *R120
C CALC. OF PRODUCTION TERM FOR ETHOH    SPECIE NO.  50
      PR( 50)=
     + +   ETHO2   *HO2     *R109
C CALC. OF LOSS TERM FOR ETHOH          SPECIE NO.  50
      LO( 50)=
     + +J 48
     + +J 49
     + +   OH      *R110
     + +   OH      *R111
C CALC. OF PRODUCTION TERM FOR CH3CHO   SPECIE NO.  51
      PR( 51)=
     + +   ETHOH   *J 48
     + +   ETHO2   *NO      *R108
     + +   ETHOH   *OH      *R111
C CALC. OF LOSS TERM FOR CH3CHO         SPECIE NO.  51
      LO( 51)=
     + +J 50
     + +J 51
     + +J 52
     + +   OH      *R112
C CALC. OF PRODUCTION TERM FOR CH3CO3   SPECIE NO.  52
      PR( 52)=
     + +   CH3CHO  *J 50
     + +   PAN     *J 55
     + +   CH3CHO  *OH      *R112
     + +   PAN     *R116
     + +   CH3CO3H *OH      *R118
C CALC. OF LOSS TERM FOR CH3CO3         SPECIE NO.  52
      LO( 52)=
     + +   NO      *R113
     + +   HO2     *R114
     + +   NO2     *R115
C CALC. OF PRODUCTION TERM FOR CH3CO3H  SPECIE NO.  53
      PR( 53)=
     + +   CH3CO3  *HO2     *R114
C CALC. OF LOSS TERM FOR CH3CO3H        SPECIE NO.  53
      LO( 53)=
     + +J 53
     + +J 54
     + +   OH      *R118
     + +   OH      *R119
C CALC. OF PRODUCTION TERM FOR ETHO2NO2 SPECIE NO.  54
      PR( 54)=
     + +   ETHO2   *NO2     *R120
C CALC. OF LOSS TERM FOR ETHO2NO2       SPECIE NO.  54
      LO( 54)=
     + +R121
* CALCULATION OF THE JACOBIAN FOR FAMILY NO.  9 INCLUDING SPECIES
*  ETHO2   ,ETHOH   ,CH3CHO  ,CH3CO3  ,CH3CO3H ,ETHO2NO2
      A 6(49)=
     + +   C2H6    *R 82                +   ETHOH   *R110               
      A 7(49)=
     + -   ETHO2   *R109               
      A21(49)=
     + +   C2H6    *R 81               
      A33(49)=
     + -   ETHO2   *R108               
      A34(49)=
     + -   ETHO2   *R120               
      A49(49)=
     + -   NO      *R108                -   HO2     *R109               
     + -   NO2     *R120               
      A50(49)=
     + +   OH      *R110               
      A54(49)=
     + +   R121                        
      A 6(50)=
     + -   ETHOH   *R110                -   ETHOH   *R111               
      A 7(50)=
     + +   ETHO2   *R109               
      A49(50)=
     + +   HO2     *R109               
      A50(50)=
     + -   J 48                         -   J 49                        
     + -   OH      *R110                -   OH      *R111               
      A 6(51)=
     + +   ETHOH   *R111                -   CH3CHO  *R112               
      A33(51)=
     + +   ETHO2   *R108               
      A49(51)=
     + +   NO      *R108               
      A50(51)=
     + +   J 48                         +   OH      *R111               
      A51(51)=
     + -   J 50                         -   J 51                        
     + -   J 52                         -   OH      *R112               
      A 6(52)=
     + +   CH3CHO  *R112                +   CH3CO3H *R118               
      A 7(52)=
     + -   CH3CO3  *R114               
      A33(52)=
     + -   CH3CO3  *R113               
      A34(52)=
     + -   CH3CO3  *R115               
      A51(52)=
     + +   J 50                         +   OH      *R112               
      A52(52)=
     + -   NO      *R113                -   HO2     *R114               
     + -   NO2     *R115               
      A53(52)=
     + +   OH      *R118               
      A 6(53)=
     + -   CH3CO3H *R118                -   CH3CO3H *R119               
      A 7(53)=
     + +   CH3CO3  *R114               
      A52(53)=
     + +   HO2     *R114               
      A53(53)=
     + -   J 53                         -   J 54                        
     + -   OH      *R118                -   OH      *R119               
      A34(54)=
     + +   ETHO2   *R120               
      A49(54)=
     + +   NO2     *R120               
      A54(54)=
     + -   R121                        
      RETURN
   10 CONTINUE
C CALC. OF PRODUCTION TERM FOR N2O      SPECIE NO.  55
      PR( 55)=
     + +   N       *NO2     *R183
     + +   O1D     *N2      *R212
C CALC. OF LOSS TERM FOR N2O            SPECIE NO.  55
      LO( 55)=
     + +J 13
     + +   O1D     *R 46
     + +   O1D     *R 47
     + +   OH      *R335
C CALC. OF PRODUCTION TERM FOR CH4      SPECIE NO.  56
      PR( 56)=
     + +   CH3CHO  *J 51
C CALC. OF LOSS TERM FOR CH4            SPECIE NO.  56
      LO( 56)=
     + +J 82
     + +   O1D     *R  8
     + +   OH      *R 26
     + +   CL      *R 35
     + +   O1D     *R130
C CALC. OF PRODUCTION TERM FOR C2H6     SPECIE NO.  57
      PR( 57)=
     + +   0.
C CALC. OF LOSS TERM FOR C2H6           SPECIE NO.  57
      LO( 57)=
     + +   CL      *R 81
     + +   OH      *R 82
C CALC. OF PRODUCTION TERM FOR H2       SPECIE NO.  58
      PR( 58)=
     + +   CH2O    *J 31
     + +   H2O     *J 79
     + +2 *CH4     *J 82
     + +   H       *HO2     *R 21
     + +   O1D     *CH4     *R130
C CALC. OF LOSS TERM FOR H2             SPECIE NO.  58
      LO( 58)=
     + +   O1D     *R  9
     + +   OH      *R 30
     + +   CL      *R 31
C CALC. OF PRODUCTION TERM FOR CO       SPECIE NO.  59
      PR( 59)=
     + +   CFCL3   *J 16
     + +   CCL4    *J 18
     + +   CH3CL   *J 19
     + +   CF2O    *J 20
     + +2 *CH3CCL3 *J 24
     + +   CH2O    *J 30
     + +   CH2O    *J 31
     + +2 *C2CL3F3 *J 44
     + +   CBRF3   *J 46
     + +   CH3CHO  *J 51
     + +   CH3CHO  *J 52
     + +2 *C2CLF5  *J 64
     + +2 *C2H3FCL2*J 65
     + +   C2H3F2CL*J 66
     + +   C2HF3CL2*J 67
     + +   OCS     *J 69
     + +   CO2     *J 80
     + +   CO2     *J 81
     + +   CH4     *J 82
     + +   CH2FCF3 *J111
     + +2 *CF3CH3  *J112
     + +   CH2F2   *J114
     + +2 *CH3CHF2 *J116
     + +   C3H3F5  *J118
     + +   O1D     *CFCL3   *R 48
     + +   O1D     *CCL4    *R 51
     + +   O1D     *CF2O    *R 55
     + +2 *CH3CCL3 *OH      *R 71
     + +   CH2O    *OH      *R 76
     + +   CL      *CH2O    *R 77
     + +   BR      *CH2O    *R 95
     + +   CH3CO3H *OH      *R119
     + +   CH2O    *O       *R142
     + +   O1D     *CBRF3   *R145
     + +2 *O1D     *C2CL3F3 *R146
     + +2 *O1D     *C2CLF5  *R176
     + +2 *O1D     *C2H3FCL2*R177
     + +2 *OH      *C2H3FCL2*R178
     + +   O1D     *C2H3F2CL*R179
     + +   OH      *C2H3F2CL*R180
     + +2 *O1D     *C2HF3CL2*R181
     + +2 *OH      *C2HF3CL2*R182
     + +   CS2     *O       *R186
     + +   OCS     *OH      *R187
     + +   OCS     *O       *R188
     + +   OH      *CF3I    *R219
     + +   CL      *CH3I    *R223
     + +   OH      *CH2BR2  *R251
     + +   OH      *CHBR3   *R252
     + +2 *CL      *C2H3FCL2*R258
     + +2 *CL      *C2H3F2CL*R259
     + +2 *CL      *C2HF3CL2*R260
     + +   O1D     *CHBR3   *R262
     + +   OH      *CCL4    *R289
     + +   OH      *CFCL3   *R290
     + +   O1D     *CH2F2   *R299
     + +   O1D     *CH2FCF3 *R313
     + +2 *O1D     *CH3CCL3 *R327
     + +   O1D     *C3H3F5  *R333
     + +   OH      *C2CL3F3 *R336
     + +2 *OH      *CH3CHF2 *R345
     + +2 *OH      *CF3CH3  *R346
     + +   OH      *C3H3F5  *R349
     + +   CL      *CH2F2   *R351
     + +2 *CL      *CH3CHF2 *R353
     + +2 *CL      *CF3CH3  *R354
     + +   CL      *CCL4    *R357
     + +   CL      *CFCL3   *R358
     + +   CL      *C2CL3F3 *R360
     + +2 *CL      *C3H3F5  *R368
C CALC. OF LOSS TERM FOR CO             SPECIE NO.  59
      LO( 59)=
     + +   OH      *R 23
     + +   O       *M       *R213
     + +   OH      *R247
C CALC. OF PRODUCTION TERM FOR CH3CL    SPECIE NO.  60
      PR( 60)=
     + +   O1D     *CH3CL   *R329
C CALC. OF LOSS TERM FOR CH3CL          SPECIE NO.  60
      LO( 60)=
     + +J 19
     + +   OH      *R  6
     + +   O1D     *R 50
     + +   CL      *R255
     + +   O1D     *R329
C CALC. OF PRODUCTION TERM FOR CCL4     SPECIE NO.  61
      PR( 61)=
     + +   O1D     *CCL4    *R324
C CALC. OF LOSS TERM FOR CCL4           SPECIE NO.  61
      LO( 61)=
     + +J 18
     + +   O1D     *R 51
     + +   OH      *R289
     + +   O1D     *R324
     + +   CL      *R357
C CALC. OF PRODUCTION TERM FOR CFCL3    SPECIE NO.  62
      PR( 62)=
     + +   O1D     *CFCL3   *R304
C CALC. OF LOSS TERM FOR CFCL3          SPECIE NO.  62
      LO( 62)=
     + +J 16
     + +   O1D     *R 48
     + +   OH      *R290
     + +   O1D     *R304
     + +   CL      *R358
C CALC. OF PRODUCTION TERM FOR CF2CL2   SPECIE NO.  63
      PR( 63)=
     + +   O1D     *CF2CL2  *R305
C CALC. OF LOSS TERM FOR CF2CL2         SPECIE NO.  63
      LO( 63)=
     + +J 17
     + +   O1D     *R 49
     + +   OH      *R291
     + +   O1D     *R305
     + +   CL      *R359
C CALC. OF PRODUCTION TERM FOR CH3CCL3  SPECIE NO.  64
      PR( 64)=
     + +   O1D     *CH3CCL3 *R328
C CALC. OF LOSS TERM FOR CH3CCL3        SPECIE NO.  64
      LO( 64)=
     + +J 24
     + +   OH      *R 71
     + +   CL      *R257
     + +   O1D     *R327
     + +   O1D     *R328
C CALC. OF PRODUCTION TERM FOR CF2O     SPECIE NO.  65
      PR( 65)=
     + +   CF2CL2  *J 17
     + +   CHCLF2  *J 43
     + +   CBRCLF2 *J 45
     + +2 *C2CL2F4 *J 63
     + +   C2H3F2CL*J 66
     + +   C2HF3CL2*J 67
     + +   CBR2F2  *J 76
     + +2 *C2BR2F4 *J 77
     + +   CH2FCF3 *J111
     + +   CHF3    *J113
     + +2 *CHF2CF3 *J115
     + +3 *C3HF7   *J117
     + +2 *C3H3F5  *J118
     + +   O1D     *CF2CL2  *R 49
     + +   CHCLF2  *OH      *R105
     + +   CHCLF2  *O1D     *R106
     + +   O1D     *CBRCLF2 *R144
     + +2 *O1D     *C2CL2F4 *R175
     + +   O1D     *C2H3F2CL*R179
     + +   OH      *C2H3F2CL*R180
     + +   O1D     *CBR2F2  *R205
     + +2 *O1D     *C2BR2F4 *R206
     + +   CL      *CHCLF2  *R256
     + +   OH      *CF2CL2  *R291
     + +   O1D     *CHF3    *R301
     + +   O1D     *CH2FCF3 *R313
     + +3 *O1D     *C3HF7   *R331
     + +2 *O1D     *C3H3F5  *R333
     + +   OH      *C2CL3F3 *R336
     + +2 *OH      *C2CL2F4 *R337
     + +2 *OH      *C2CLF5  *R338
     + +   OH      *CBR2F2  *R339
     + +   OH      *CBRF3   *R340
     + +   OH      *CBRCLF2 *R341
     + +2 *OH      *C2BR2F4 *R342
     + +   OH      *CH2F2   *R343
     + +   OH      *CHF3    *R344
     + +2 *OH      *CH2FCF3 *R347
     + +2 *OH      *CHF2CF3 *R348
     + +2 *OH      *C3H3F5  *R349
     + +3 *OH      *C3HF7   *R350
     + +   CL      *CHF3    *R352
     + +2 *CL      *CH2FCF3 *R355
     + +2 *CL      *CHF2CF3 *R356
     + +   CL      *CF2CL2  *R359
     + +   CL      *C2CL3F3 *R360
     + +2 *CL      *C2CL2F4 *R361
     + +2 *CL      *C2CLF5  *R362
     + +   CL      *CBRCLF2 *R363
     + +   CL      *CBRF3   *R364
     + +   CL      *CBR2F2  *R365
     + +2 *CL      *C2BR2F4 *R366
     + +3 *CL      *C3HF7   *R367
     + +   CL      *C3H3F5  *R368
C CALC. OF LOSS TERM FOR CF2O           SPECIE NO.  65
      LO( 65)=
     + +J 20
     + +J 56
     + +   O1D     *R 55
C CALC. OF PRODUCTION TERM FOR CH3BR    SPECIE NO.  66
      PR( 66)=
     + +   O1D     *CH3BR   *R330
C CALC. OF LOSS TERM FOR CH3BR          SPECIE NO.  66
      LO( 66)=
     + +J 42
     + +   OH      *R100
     + +   O1D     *R143
     + +   CL      *R250
     + +   O1D     *R330
C CALC. OF PRODUCTION TERM FOR CH2BR2   SPECIE NO.  67
      PR( 67)=
     + +   O1D     *CH2BR2  *R326
C CALC. OF LOSS TERM FOR CH2BR2         SPECIE NO.  67
      LO( 67)=
     + +J 96
     + +   OH      *R251
     + +   CL      *R253
     + +   O1D     *R261
     + +   O1D     *R326
C CALC. OF PRODUCTION TERM FOR CHBR3    SPECIE NO.  68
      PR( 68)=
     + +   O1D     *CHBR3   *R325
C CALC. OF LOSS TERM FOR CHBR3          SPECIE NO.  68
      LO( 68)=
     + +J 97
     + +   OH      *R252
     + +   CL      *R254
     + +   O1D     *R262
     + +   O1D     *R325
C CALC. OF PRODUCTION TERM FOR CHCLF2   SPECIE NO.  69
      PR( 69)=
     + +   O1D     *CHCLF2  *R303
C CALC. OF LOSS TERM FOR CHCLF2         SPECIE NO.  69
      LO( 69)=
     + +J 43
     + +   OH      *R105
     + +   O1D     *R106
     + +   CL      *R256
     + +   O1D     *R303
C CALC. OF PRODUCTION TERM FOR C2CL3F3  SPECIE NO.  70
      PR( 70)=
     + +   O1D     *C2CL3F3 *R318
C CALC. OF LOSS TERM FOR C2CL3F3        SPECIE NO.  70
      LO( 70)=
     + +J 44
     + +   O1D     *R146
     + +   O1D     *R318
     + +   OH      *R336
     + +   CL      *R360
C CALC. OF PRODUCTION TERM FOR C2CL2F4  SPECIE NO.  71
      PR( 71)=
     + +   O1D     *C2CL2F4 *R319
C CALC. OF LOSS TERM FOR C2CL2F4        SPECIE NO.  71
      LO( 71)=
     + +J 63
     + +   O1D     *R175
     + +   O1D     *R319
     + +   OH      *R337
     + +   CL      *R361
C CALC. OF PRODUCTION TERM FOR C2CLF5   SPECIE NO.  72
      PR( 72)=
     + +   O1D     *C2CLF5  *R320
C CALC. OF LOSS TERM FOR C2CLF5         SPECIE NO.  72
      LO( 72)=
     + +J 64
     + +   O1D     *R176
     + +   O1D     *R320
     + +   OH      *R338
     + +   CL      *R362
C CALC. OF PRODUCTION TERM FOR C2H3FCL2 SPECIE NO.  73
      PR( 73)=
     + +   O1D     *C2H3FCL2*R311
C CALC. OF LOSS TERM FOR C2H3FCL2       SPECIE NO.  73
      LO( 73)=
     + +J 65
     + +   O1D     *R177
     + +   OH      *R178
     + +   CL      *R258
     + +   O1D     *R311
C CALC. OF PRODUCTION TERM FOR C2H3F2CL SPECIE NO.  74
      PR( 74)=
     + +   O1D     *C2H3F2CL*R312
C CALC. OF LOSS TERM FOR C2H3F2CL       SPECIE NO.  74
      LO( 74)=
     + +J 66
     + +   O1D     *R179
     + +   OH      *R180
     + +   CL      *R259
     + +   O1D     *R312
C CALC. OF PRODUCTION TERM FOR C2HF3CL2 SPECIE NO.  75
      PR( 75)=
     + +   O1D     *C2HF3CL2*R315
C CALC. OF LOSS TERM FOR C2HF3CL2       SPECIE NO.  75
      LO( 75)=
     + +J 67
     + +   O1D     *R181
     + +   OH      *R182
     + +   CL      *R260
     + +   O1D     *R315
C CALC. OF PRODUCTION TERM FOR CBRCLF2  SPECIE NO.  76
      PR( 76)=
     + +   O1D     *CBRCLF2 *R306
C CALC. OF LOSS TERM FOR CBRCLF2        SPECIE NO.  76
      LO( 76)=
     + +J 45
     + +   O1D     *R144
     + +   O1D     *R306
     + +   OH      *R341
     + +   CL      *R363
C CALC. OF PRODUCTION TERM FOR CBRF3    SPECIE NO.  77
      PR( 77)=
     + +   O1D     *CBRF3   *R307
C CALC. OF LOSS TERM FOR CBRF3          SPECIE NO.  77
      LO( 77)=
     + +J 46
     + +   O1D     *R145
     + +   O1D     *R307
     + +   OH      *R340
     + +   CL      *R364
C CALC. OF PRODUCTION TERM FOR CBR2F2   SPECIE NO.  78
      PR( 78)=
     + +   O1D     *CBR2F2  *R308
C CALC. OF LOSS TERM FOR CBR2F2         SPECIE NO.  78
      LO( 78)=
     + +J 76
     + +   O1D     *R205
     + +   O1D     *R308
     + +   OH      *R339
     + +   CL      *R365
C CALC. OF PRODUCTION TERM FOR C2BR2F4  SPECIE NO.  79
      PR( 79)=
     + +   O1D     *C2BR2F4 *R321
C CALC. OF LOSS TERM FOR C2BR2F4        SPECIE NO.  79
      LO( 79)=
     + +J 77
     + +   O1D     *R206
     + +   O1D     *R321
     + +   OH      *R342
     + +   CL      *R366
C CALC. OF PRODUCTION TERM FOR CH3I     SPECIE NO.  80
      PR( 80)=
     + +   0.
C CALC. OF LOSS TERM FOR CH3I           SPECIE NO.  80
      LO( 80)=
     + +J 84
     + +   OH      *R218
     + +   CL      *R223
C CALC. OF PRODUCTION TERM FOR CF3I     SPECIE NO.  81
      PR( 81)=
     + +   0.
C CALC. OF LOSS TERM FOR CF3I           SPECIE NO.  81
      LO( 81)=
     + +J 83
     + +   OH      *R219
C CALC. OF PRODUCTION TERM FOR CH2FCF3  SPECIE NO.  82
      PR( 82)=
     + +   O1D     *CH2FCF3 *R314
C CALC. OF LOSS TERM FOR CH2FCF3        SPECIE NO.  82
      LO( 82)=
     + +J111
     + +   O1D     *R313
     + +   O1D     *R314
     + +   OH      *R347
     + +   CL      *R355
C CALC. OF PRODUCTION TERM FOR CF3CH3   SPECIE NO.  83
      PR( 83)=
     + +   O1D     *CF3CH3  *R323
C CALC. OF LOSS TERM FOR CF3CH3         SPECIE NO.  83
      LO( 83)=
     + +J112
     + +   O1D     *R322
     + +   O1D     *R323
     + +   OH      *R346
     + +   CL      *R354
C CALC. OF PRODUCTION TERM FOR CHF3     SPECIE NO.  84
      PR( 84)=
     + +   O1D     *CHF3    *R302
C CALC. OF LOSS TERM FOR CHF3           SPECIE NO.  84
      LO( 84)=
     + +J113
     + +   O1D     *R301
     + +   O1D     *R302
     + +   OH      *R344
     + +   CL      *R352
C CALC. OF PRODUCTION TERM FOR CH2F2   SPECIE NO.  85
      PR( 85)=
     + +   O1D     *CH2F2  *R300
C CALC. OF LOSS TERM FOR CH2F2         SPECIE NO.  85
      LO( 85)=
     + +J114
     + +   O1D     *R299
     + +   O1D     *R300
     + +   OH      *R343
     + +   CL      *R351
C CALC. OF PRODUCTION TERM FOR CHF2CF3  SPECIE NO.  86
      PR( 86)=
     + +   O1D     *CHF2CF3 *R317
C CALC. OF LOSS TERM FOR CHF2CF3        SPECIE NO.  86
      LO( 86)=
     + +J115
     + +   O1D     *R316
     + +   O1D     *R317
     + +   OH      *R348
     + +   CL      *R356
C CALC. OF PRODUCTION TERM FOR CH3CHF2  SPECIE NO.  87
      PR( 87)=
     + +   O1D     *CH3CHF2 *R310
C CALC. OF LOSS TERM FOR CH3CHF2        SPECIE NO.  87
      LO( 87)=
     + +J116
     + +   O1D     *R309
     + +   O1D     *R310
     + +   OH      *R345
     + +   CL      *R353
C CALC. OF PRODUCTION TERM FOR C3HF7    SPECIE NO.  88
      PR( 88)=
     + +   O1D     *C3HF7   *R332
C CALC. OF LOSS TERM FOR C3HF7          SPECIE NO.  88
      LO( 88)=
     + +J117
     + +   O1D     *R331
     + +   O1D     *R332
     + +   OH      *R350
     + +   CL      *R367
C CALC. OF PRODUCTION TERM FOR C3H3F5   SPECIE NO.  89
      PR( 89)=
     + +   O1D     *C3H3F5  *R334
C CALC. OF LOSS TERM FOR C3H3F5         SPECIE NO.  89
      LO( 89)=
     + +J118
     + +   O1D     *R333
     + +   O1D     *R334
     + +   OH      *R349
     + +   CL      *R368
C CALC. OF PRODUCTION TERM FOR OX       SPECIE NO.  90
      PR( 90)=
     + +   0.
C CALC. OF LOSS TERM FOR OX             SPECIE NO.  90
      LO( 90)=
     + +   0.
C CALC. OF PRODUCTION TERM FOR NOX      SPECIE NO.  91
      PR( 91)=
     + +   CH3CN   *J 94
     + +   HCN     *J 95
     + +2 *N2O     *O1D     *R 47
     + +   CH3CN   *OH      *R239
     + +   CH3CN   *CL      *R240
     + +   CH3CN   *O       *R241
     + +   CH3CN   *O1D     *R242
     + +   HCN     *OH      *R243
     + +   HCN     *O       *R244
     + +   HCN     *O1D     *R245
C CALC. OF LOSS TERM FOR NOX            SPECIE NO.  91
      LO( 91)=
     + +   0.
C CALC. OF PRODUCTION TERM FOR CLX      SPECIE NO.  92
      PR( 92)=
     + +   ODP     *J  4
     + +3 *CFCL3   *J 16
     + +2 *CF2CL2  *J 17
     + +4 *CCL4    *J 18
     + +   CH3CL   *J 19
     + +3 *CH3CCL3 *J 24
     + +   CHCLF2  *J 43
     + +3 *C2CL3F3 *J 44
     + +   CBRCLF2 *J 45
     + +2 *C2CL2F4 *J 63
     + +   C2CLF5  *J 64
     + +2 *C2H3FCL2*J 65
     + +   C2H3F2CL*J 66
     + +2 *C2HF3CL2*J 67
     + +   ODP     *OH      *R  1
     + +   ODP     *O1D     *R  2
     + +   ODP     *CL      *R  3
     + +   OH      *CH3CL   *R  6
     + +3 *O1D     *CFCL3   *R 48
     + +2 *O1D     *CF2CL2  *R 49
     + +   O1D     *CH3CL   *R 50
     + +4 *O1D     *CCL4    *R 51
     + +3 *CH3CCL3 *OH      *R 71
     + +   CHCLF2  *OH      *R105
     + +   CHCLF2  *O1D     *R106
     + +   O1D     *CBRCLF2 *R144
     + +3 *O1D     *C2CL3F3 *R146
     + +2 *O1D     *C2CL2F4 *R175
     + +   O1D     *C2CLF5  *R176
     + +2 *O1D     *C2H3FCL2*R177
     + +2 *OH      *C2H3FCL2*R178
     + +   O1D     *C2H3F2CL*R179
     + +   OH      *C2H3F2CL*R180
     + +2 *O1D     *C2HF3CL2*R181
     + +2 *OH      *C2HF3CL2*R182
     + +   CL      *CH3CL   *R255
     + +   CL      *CHCLF2  *R256
     + +3 *CL      *CH3CCL3 *R257
     + +2 *CL      *C2H3FCL2*R258
     + +   CL      *C2H3F2CL*R259
     + +2 *CL      *C2HF3CL2*R260
     + +4 *OH      *CCL4    *R289
     + +3 *OH      *CFCL3   *R290
     + +2 *OH      *CF2CL2  *R291
     + +3 *O1D     *CH3CCL3 *R327
     + +3 *OH      *C2CL3F3 *R336
     + +2 *OH      *C2CL2F4 *R337
     + +   OH      *C2CLF5  *R338
     + +   OH      *CBRCLF2 *R341
     + +4 *CL      *CCL4    *R357
     + +3 *CL      *CFCL3   *R358
     + +2 *CL      *CF2CL2  *R359
     + +3 *CL      *C2CL3F3 *R360
     + +2 *CL      *C2CL2F4 *R361
     + +   CL      *C2CLF5  *R362
     + +   CL      *CBRCLF2 *R363
C CALC. OF LOSS TERM FOR CLX            SPECIE NO.  92
      LO( 92)=
     + +J 27
C CALC. OF PRODUCTION TERM FOR BRX      SPECIE NO.  93
      PR( 93)=
     + +   ODP     *J  4
     + +   CH3BR   *J 42
     + +   CBRCLF2 *J 45
     + +   CBRF3   *J 46
     + +2 *CBR2F2  *J 76
     + +2 *C2BR2F4 *J 77
     + +2 *CH2BR2  *J 96
     + +3 *CHBR3   *J 97
     + +   ODP     *OH      *R  1
     + +   ODP     *O1D     *R  2
     + +   ODP     *CL      *R  3
     + +   OH      *CH3BR   *R100
     + +   O1D     *CH3BR   *R143
     + +   O1D     *CBRCLF2 *R144
     + +   O1D     *CBRF3   *R145
     + +2 *O1D     *CBR2F2  *R205
     + +2 *O1D     *C2BR2F4 *R206
     + +   CL      *CH3BR   *R250
     + +2 *OH      *CH2BR2  *R251
     + +3 *OH      *CHBR3   *R252
     + +2 *CL      *CH2BR2  *R253
     + +3 *CL      *CHBR3   *R254
     + +2 *O1D     *CH2BR2  *R261
     + +3 *O1D     *CHBR3   *R262
     + +2 *OH      *CBR2F2  *R339
     + +   OH      *CBRF3   *R340
     + +   OH      *CBRCLF2 *R341
     + +2 *OH      *C2BR2F4 *R342
     + +   CL      *CBRCLF2 *R363
     + +   CL      *CBRF3   *R364
     + +2 *CL      *CBR2F2  *R365
     + +2 *CL      *C2BR2F4 *R366
C CALC. OF LOSS TERM FOR BRX            SPECIE NO.  93
      LO( 93)=
     + +J 28
C CALC. OF PRODUCTION TERM FOR FX       SPECIE NO.  94
      PR( 94)=
     + +   ODP     *J  4
     + +   CFCL3   *J 16
     + +2 *CF2O    *J 20
     + +3 *C2CL3F3 *J 44
     + +3 *CBRF3   *J 46
     + +5 *C2CLF5  *J 64
     + +   C2H3FCL2*J 65
     + +   C2HF3CL2*J 67
     + +3 *CF3I    *J 83
     + +2 *CH2FCF3 *J111
     + +3 *CF3CH3  *J112
     + +   CHF3    *J113
     + +2 *CH2F2   *J114
     + +   CHF2CF3 *J115
     + +2 *CH3CHF2 *J116
     + +   C3HF7   *J117
     + +   C3H3F5  *J118
     + +   ODP     *OH      *R  1
     + +   ODP     *O1D     *R  2
     + +   ODP     *CL      *R  3
     + +   O1D     *CFCL3   *R 48
     + +2 *O1D     *CF2O    *R 55
     + +3 *O1D     *CBRF3   *R145
     + +3 *O1D     *C2CL3F3 *R146
     + +5 *O1D     *C2CLF5  *R176
     + +   O1D     *C2H3FCL2*R177
     + +   OH      *C2H3FCL2*R178
     + +3 *O1D     *C2HF3CL2*R181
     + +3 *OH      *C2HF3CL2*R182
     + +3 *OH      *CF3I    *R219
     + +   CL      *C2H3FCL2*R258
     + +2 *CL      *C2H3F2CL*R259
     + +3 *CL      *C2HF3CL2*R260
     + +   OH      *CFCL3   *R290
     + +   O1D     *FX      *R293
     + +2 *O1D     *CH2F2   *R299
     + +   O1D     *CHF3    *R301
     + +2 *O1D     *CH3CHF2 *R309
     + +2 *O1D     *CH2FCF3 *R313
     + +   O1D     *CHF2CF3 *R316
     + +3 *O1D     *CF3CH3  *R322
     + +   O1D     *C3HF7   *R331
     + +   O1D     *C3H3F5  *R333
     + +   OH      *C2CL3F3 *R336
     + +   OH      *C2CLF5  *R338
     + +   OH      *CBRF3   *R340
     + +   OH      *CHF3    *R344
     + +2 *OH      *CH3CHF2 *R345
     + +3 *OH      *CF3CH3  *R346
     + +   OH      *CHF2CF3 *R348
     + +   OH      *C3H3F5  *R349
     + +   OH      *C3HF7   *R350
     + +2 *CL      *CH2F2   *R351
     + +   CL      *CHF3    *R352
     + +2 *CL      *CH3CHF2 *R353
     + +3 *CL      *CF3CH3  *R354
     + +   CL      *CHF2CF3 *R356
     + +   CL      *CFCL3   *R358
     + +   CL      *C2CL3F3 *R360
     + +   CL      *C2CLF5  *R362
     + +   CL      *CBRF3   *R364
     + +   CL      *C3HF7   *R367
     + +3 *CL      *C3H3F5  *R368
C CALC. OF LOSS TERM FOR FX             SPECIE NO.  94
      LO( 94)=
     + +J 29
     + +   O1D     *R293
C CALC. OF PRODUCTION TERM FOR IX       SPECIE NO.  95
      PR( 95)=
     + +   CF3I    *J 83
     + +   CH3I    *J 84
     + +   OH      *CH3I    *R218
     + +   OH      *CF3I    *R219
     + +   CL      *CH3I    *R223
C CALC. OF LOSS TERM FOR IX             SPECIE NO.  95
      LO( 95)=
     + +J 93
C CALC. OF PRODUCTION TERM FOR PAN      SPECIE NO.  96
      PR( 96)=
     + +   CH3CO3  *NO2     *R115
     + +   DMS     *NO3     *R197
C CALC. OF LOSS TERM FOR PAN            SPECIE NO.  96
      LO( 96)=
     + +J 55
     + +R116
     + +   OH      *R117
C CALC. OF PRODUCTION TERM FOR SHNO3    SPECIE NO.  97
      PR( 97)=
     + +   0.
C CALC. OF LOSS TERM FOR SHNO3          SPECIE NO.  97
      LO( 97)=
     + +   0.
C CALC. OF PRODUCTION TERM FOR H2O      SPECIE NO.  98
      PR( 98)=
     + +   CH3CL   *J 19
     + +   CH3CCL3 *J 24
     + +   C2H3FCL2*J 65
     + +   C2H3F2CL*J 66
     + +   CH2BR2  *J 96
     + +   CH3CHF2 *J116
     + +   C3H3F5  *J118
     + +   OH      *CH3CL   *R  6
     + +   OH      *H2O2    *R 13
     + +   HO2     *OH      *R 19
     + +   OH      *HCL     *R 25
     + +   OH      *CH4     *R 26
     + +   OH      *HNO3    *R 28
     + +   OH      *OH      *R 29
     + +   OH      *H2      *R 30
     + +   OH      *HOCL    *R 54
     + +   OH      *HBR     *R 70
     + +2 *CH3CCL3 *OH      *R 71
     + +   CH2O    *OH      *R 76
     + +   HO2NO2  *OH      *R 78
     + +   C2H6    *OH      *R 82
     + +   MEOH    *OH      *R 86
     + +   HOCL             *R 90
     + +   MEOH    *OH      *R104
     + +   CHCLF2  *OH      *R105
     + +   ETHOH   *OH      *R110
     + +   ETHOH   *OH      *R111
     + +   CH3CHO  *OH      *R112
     + +   PAN     *OH      *R117
     + +   CH3CO3H *OH      *R118
     + +   CH3CO3H *OH      *R119
     + +   HOBR             *R124
     + +   OH      *HONO    *R126
     + +   HOCL             *R162
     + +   HOCL             *R167
     + +   HOBR             *R168
     + +   HOCL    *HBR     *R174
     + +   O1D     *C2H3FCL2*R177
     + +2 *OH      *C2H3FCL2*R178
     + +   O1D     *C2H3F2CL*R179
     + +2 *OH      *C2H3F2CL*R180
     + +   OH      *C2HF3CL2*R182
     + +   H2S     *OH      *R189
     + +   H2S     *O       *R190
     + +   DMS     *NO3     *R197
     + +   SO3     *H2O     *H2O     *R204
     + +   H       *HO2     *R210
     + +   OH      *HI      *R217
     + +   CL      *CH3I    *R223
     + +   OH      *CH2BR2  *R251
     + +   OH      *CHBR3   *R252
     + +   CL      *CH3CCL3 *R257
     + +   HOCL             *R278
     + +   HOBR             *R281
     + +   HOBR             *R282
     + +   O1D     *CH3CHF2 *R309
     + +   O1D     *CH3CCL3 *R327
     + +   O1D     *C3H3F5  *R333
     + +   OH      *CH2F2   *R343
     + +   OH      *CHF3    *R344
     + +   OH      *CH3CHF2 *R345
     + +   OH      *CF3CH3  *R346
     + +   OH      *CH2FCF3 *R347
     + +   OH      *CHF2CF3 *R348
     + +   OH      *C3H3F5  *R349
     + +   OH      *C3HF7   *R350
     + +   HOONO   *OH      *R371
C CALC. OF LOSS TERM FOR H2O            SPECIE NO.  98
      LO( 98)=
     + +J 12
     + +J 79
     + +   O1D     *R  7
     + +   N2O5    *R 83
     + +   SO3     *H2O     *R204
     + +   SO3     *H2O     *R204
     + +   NO2     *R374
C CALC. OF PRODUCTION TERM FOR SH2O     SPECIE NO.  99
      PR( 99)=
     + +   0.
C CALC. OF LOSS TERM FOR SH2O           SPECIE NO.  99
      LO( 99)=
     + +   0.
C CALC. OF PRODUCTION TERM FOR ODP      SPECIE NO. 100
      PR(100)=
     + +   0.
C CALC. OF LOSS TERM FOR ODP            SPECIE NO. 100
      LO(100)=
     + +J  4
     + +   OH      *R  1
     + +   O1D     *R  2
     + +   CL      *R  3
C CALC. OF PRODUCTION TERM FOR CS2      SPECIE NO. 101
      PR(101)=
     + +   0.
C CALC. OF LOSS TERM FOR CS2            SPECIE NO. 101
      LO(101)=
     + +J 68
     + +   OH      *R184
     + +   O       *R185
     + +   O       *R186
C CALC. OF PRODUCTION TERM FOR DMS      SPECIE NO. 102
      PR(102)=
     + +   0.
C CALC. OF LOSS TERM FOR DMS            SPECIE NO. 102
      LO(102)=
     + +   OH      *R193
     + +   OH      *R194
     + +   OH      *R195
     + +   O       *R196
     + +   NO3     *R197
     + +   CL      *R198
     + +   CLO     *R199
     + +   BR      *R200
     + +   BRO     *R201
C CALC. OF PRODUCTION TERM FOR H2S      SPECIE NO. 103
      PR(103)=
     + +   0.
C CALC. OF LOSS TERM FOR H2S            SPECIE NO. 103
      LO(103)=
     + +   OH      *R189
     + +   O       *R190
     + +   CL      *R191
     + +   BR      *R192
C CALC. OF PRODUCTION TERM FOR MSA      SPECIE NO. 104
      PR(104)=
     + +   DMS     *OH      *R195
C CALC. OF LOSS TERM FOR MSA            SPECIE NO. 104
      LO(104)=
     + +J 74
C CALC. OF PRODUCTION TERM FOR OCS      SPECIE NO. 105
      PR(105)=
     + +   CS2     *J 68
     + +   CS2     *OH      *R184
     + +   CS2     *O       *R185
C CALC. OF LOSS TERM FOR OCS            SPECIE NO. 105
      LO(105)=
     + +J 69
     + +   OH      *R187
     + +   O       *R188
C CALC. OF PRODUCTION TERM FOR SO2      SPECIE NO. 106
      PR(106)=
     + +   CS2     *J 68
     + +   OCS     *J 69
     + +   SO3     *J 71
     + +   CS2     *OH      *R184
     + +   CS2     *O       *R185
     + +2 *CS2     *O       *R186
     + +   OCS     *OH      *R187
     + +   OCS     *O       *R188
     + +   H2S     *OH      *R189
     + +   H2S     *O       *R190
     + +   H2S     *CL      *R191
     + +   H2S     *BR      *R192
     + +   DMS     *OH      *R193
     + +   DMS     *OH      *R194
     + +   DMS     *O       *R196
     + +   DMS     *NO3     *R197
     + +   DMS     *CL      *R198
     + +   DMS     *CLO     *R199
     + +   DMS     *BR      *R200
     + +   DMS     *BRO     *R201
C CALC. OF LOSS TERM FOR SO2            SPECIE NO. 106
      LO(106)=
     + +J 72
     + +J 75
     + +   OH      *R202
     + +   O3      *R203
     + +   O       *R249
C CALC. OF PRODUCTION TERM FOR SO3      SPECIE NO. 107
      PR(107)=
     + +   H2SO4   *J 70
     + +   SO2     *OH      *R202
     + +   SO2     *O3      *R203
     + +   O       *SO2     *R249
C CALC. OF LOSS TERM FOR SO3            SPECIE NO. 107
      LO(107)=
     + +J 71
     + +   H2O     *H2O     *R204
C CALC. OF PRODUCTION TERM FOR H2SO4    SPECIE NO. 108
      PR(108)=
     + +   SO3     *H2O     *H2O     *R204
C CALC. OF LOSS TERM FOR H2SO4          SPECIE NO. 108
      LO(108)=
     + +J 70
     + +J 73
C CALC. OF PRODUCTION TERM FOR CH3CN    SPECIE NO. 109
      PR(109)=
     + +   0.
C CALC. OF LOSS TERM FOR CH3CN          SPECIE NO. 109
      LO(109)=
     + +J 94
     + +   OH      *R238
     + +   OH      *R239
     + +   CL      *R240
     + +   O       *R241
     + +   O1D     *R242
C CALC. OF PRODUCTION TERM FOR HCN      SPECIE NO. 110
      PR(110)=
     + +   CH3CN   *OH      *R238
C CALC. OF LOSS TERM FOR HCN            SPECIE NO. 110
      LO(110)=
     + +J 95
     + +   OH      *R243
     + +   O       *R244
     + +   O1D     *R245
      RETURN
      END
