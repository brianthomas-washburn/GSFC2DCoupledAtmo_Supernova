C 
	SUBROUTINE KMAP(KR,NKR,IJC,IKC)


        include "com2d.h"


        INTEGER NKR
        DOUBLE PRECISION KR(NKR)

        SAVE


C  GAS PHASE REACTION RATES from GSFC array K(RB$=284,L$,Z$) 
C    mapped into AER KR(382) array -use IAGMAP(382) in COMMON
C    if GSFC does not have an AER reaction  (eg, ODP, Ix, sulfur, 
C         hydrocarbons, etc), IAGMAP=0 and set reaction=0.
C
C  for HET reactions, IAGMAP=0 and reaction=0 here - it gets set in HETCHEM
C

C                                   first initialize            
        do 100 iik=1,NKR
            KR(iik) = 0.0E0
            if (IAGMAP(iik) .ne. 0) KR(iik) = K(IAGMAP(iik),IJC,IKC)
 100	CONTINUE


	RETURN
	END
