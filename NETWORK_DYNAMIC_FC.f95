!==================================================================================================================================
!
!       VERSÃO 8 
!
!       IMPLEMENTAÇÃO DO NEURONIO GL
!       GANHO NEURONAL DINAMICO --> GAMA
!       DEPRESSÃO SINAPTICA DINAMICA --> W
!
!       MODELO DE NEURONIO USADO NO PAPER Sci. Report --> Kinouchi 2019.
!
!       VERSÃO 8 CALCULA A ATIVIDADE NA CRITICALIDADE W = GAMA = 1
!       CALCULA O TAMANHO E A DURAÇÃO DAS AVALANCHES, E A ATIVIDADE DE CADA MÓDULO 
!       DA REDE DO SISTEMA NA VERSÃO CRITÍCA W = GAMA = 1.0 E TBM NA VERSÃO COM
!       GANHO DINÂMICO "COMPLETO" E GANHO DINÂMICO SIMPLIFICADO
!
!       CALCULA A RELAÇÃO DE SETHNA PARA DIFERENTES NIVEIS HIERARQUICOS H  
!
!       MAIO DE 2022
!
!==================================================================================================================================
PROGRAM dinamicaGL_TCT_HMN
  IMPLICIT NONE
  INTEGER, PARAMETER :: tipoInteiro = 4, steps = 5000, transiente = 1.d4
  DOUBLE PRECISION, PARAMETER :: uTheta = 0.1d0, uW = 0.1d0 
  DOUBLE PRECISION, PARAMETER :: TauTheta = 1.d3, TauW = 1.d3, A = 1.d0 
  DOUBLE PRECISION, PARAMETER :: mu = 0.d0, curr = 0.d0, g = 0.d0, W = 1.d0, gma = 1.d0
  CHARACTER(LEN=40) :: nomeLista = 'LISTA_DE_CONEXOES_H1_N16384.dat'
  INTEGER(KIND=tipoInteiro), DIMENSION(:,:), ALLOCATABLE :: X
  INTEGER, DIMENSION(steps) :: tempoSilencio
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: indiceNeuronio
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: listaNeuronioNeuronio, listaNeuronioModulo
  DOUBLE PRECISION, DIMENSION(0:steps) :: rho 
  DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: Gama 
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: potencial, PesoSyn
  INTEGER :: modulos, verticesModulo, ii, jj, kk, contaSilencio
  LOGICAL :: avalanche
!
  CALL SD_CarregaRede
  ALLOCATE(X(0:steps, verticesModulo))
  ALLOCATE(Gama(modulos, verticesModulo))
  ALLOCATE(potencial(0:steps, modulos, verticesModulo))
  ALLOCATE(PesoSyn(modulos, verticesModulo, MAXVAL(indiceNeuronio) + verticesModulo))
!  
! ESTABELECIMENTO DAS CONDIÇÕES INICIAIS
  DO ii = 1, modulos
    DO jj = 1, verticesModulo
      Gama(ii, jj) = gma  !F_SorteioReal(0.1d0, 1.d0) 
      potencial(0, ii, jj) = F_SorteioReal(0.1d0, 1.d0) 
      DO kk = 1, verticesModulo + indiceNeuronio(ii, jj)
        PesoSyn(ii, jj, kk) = W !F_SorteioReal(0.1d0, 1.d0) !1.d2*
      ENDDO    
    ENDDO  
  ENDDO
!
  X(:, :) = 0
  !X(0, 256) = IBSET(X(0, 256), 0)
!
  !X(0, 4) = IBSET(X(0, 4), 0)
  ii = 0
  kk = 0
  contaSilencio = 0
  avalanche = .FALSE.
  DO WHILE(ii < steps)
    IF(ALL(X(ii, :) == 0)) THEN 
      X(ii, 4) = IBSET(X(ii, 4), 0) 
!
      IF(avalanche) THEN
        contaSilencio = contaSilencio + 1
        tempoSilencio(contaSilencio) = ii 
      ENDIF
!
      avalanche = .FALSE.
    ELSEIF(.NOT. avalanche) THEN
      avalanche = .TRUE.
      contaSilencio = contaSilencio + 1
      tempoSilencio(contaSilencio) = ii 
    ENDIF
!
    CALL SE_EvolucaoDoSistema(ii, potencial(ii, :, :), potencial(ii+1, :, :)) 
!
    IF(kk < transiente) THEN
      kk = kk + 1
      potencial(ii, :, :) = potencial(ii+1, :, :)
      X(ii, :) = X(ii+1, :)
    ELSE
      ii = ii + 1
      rho(ii) = DFLOAT(SUM(POPCNT(X(ii, :))))/DFLOAT(modulos*verticesModulo) 
    ENDIF
!
  ENDDO
  CALL SD_SalvaResultados
!
CONTAINS
  DOUBLE PRECISION FUNCTION F_Gama(GamaLocal, ALocal, uLocal, tauLocal, estadoLocal)!----------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: GamaLocal, ALocal, uLocal, tauLocal
    LOGICAL, INTENT(IN) :: estadoLocal
!
    F_Gama = gma        !GamaLocal + (ALocal - GamaLocal)/tauLocal 
    IF(estadoLocal) F_Gama = gma        !F_Gama - uLocal*GamaLocal
!
  END FUNCTION F_Gama!-------------------------------------------------------------------------------------------------------------
!
!  DOUBLE PRECISION FUNCTION F_Gama(GamaLocal, tauLocal, estadoLocal)!-------------------------------------------------------------
!    IMPLICIT NONE
!    DOUBLE PRECISION, INTENT(IN) :: GamaLocal, tauLocal
!    LOGICAL, INTENT(IN) :: estadoLocal
!
!    F_Gama = GamaLocal*(1.d0 + 1.d0/tauLocal) 
!    IF(estadoLocal) F_Gama = F_Gama - GamaLocal
!
!  END FUNCTION F_Gama!------------------------------------------------------------------------------------------------------------
!
  DOUBLE PRECISION FUNCTION F_Phi(potencialLocal, gamaLocal)!----------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: potencialLocal, gamaLocal
!
      F_Phi = 1.d0 - 1.d0/(1.d0 + gamaLocal*potencialLocal)
!
  END FUNCTION F_Phi!--------------------------------------------------------------------------------------------------------------
!
  DOUBLE PRECISION FUNCTION F_EvolucaoPotencial(PotenciaLocal, muLocal, currLocal, acoplamentoExcitLocal, estadoLocal)!------------
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: PotenciaLocal, muLocal, currLocal, acoplamentoExcitLocal
    LOGICAL, INTENT(IN) :: estadoLocal
!
    IF(estadoLocal) THEN
      F_EvolucaoPotencial = 0.d0 
    ELSE
      F_EvolucaoPotencial = muLocal*potenciaLocal + currLocal + acoplamentoExcitLocal  
    ENDIF
!
  END FUNCTION F_EvolucaoPotencial!------------------------------------------------------------------------------------------------
!
  DOUBLE PRECISION FUNCTION F_SorteioReal(minimo, maximo)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: minimo, maximo
!
    F_SorteioReal = minimo + (maximo - minimo)*RAND()
  END FUNCTION F_SorteioReal
!
  DOUBLE PRECISION FUNCTION F_AcoplamentoExcitatorio(verticeA, moduloA, indice, listaNN, listaNM, estado, pesosW)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: verticeA, moduloA, indice
    INTEGER, INTENT(IN), DIMENSION(indice) :: listaNN, listaNM
    DOUBLE PRECISION, INTENT(IN), DIMENSION(verticesModulo + indice) :: pesosW
    INTEGER(KIND=tipoInteiro), DIMENSION(verticesModulo) :: estado 
    INTEGER :: i, j, k, l
!
    F_AcoplamentoExcitatorio = 0.d0
    DO i = 1, verticesModulo
      IF(i == verticeA) CYCLE
      IF(BTEST(estado(i), moduloA - 1)) F_AcoplamentoExcitatorio = F_AcoplamentoExcitatorio + W!F_AcoplamentoExcitatorio+pesosW(i) 
    ENDDO  
!
    DO i = 1, indice 
      j = i + verticesModulo
      k = listaNM(i)
      l = listaNN(i)
!
      IF(BTEST(estado(l), k - 1)) F_AcoplamentoExcitatorio = F_AcoplamentoExcitatorio + W  !F_AcoplamentoExcitatorio + pesosW(j) 
    ENDDO 
!
    F_AcoplamentoExcitatorio = F_AcoplamentoExcitatorio/DFLOAT(verticesModulo + indice - 1)    
!  
  END FUNCTION F_AcoplamentoExcitatorio

  SUBROUTINE SE_EvolucaoDoSistema(i, potencialEntrada, potencialSaida)!----------------------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i
    DOUBLE PRECISION, INTENT(IN), DIMENSION(modulos, verticesModulo) :: potencialEntrada
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(modulos, verticesModulo) :: potencialSaida
    DOUBLE PRECISION :: LimiarSaturacao, Phi, acoplamentoExcitLocal
    INTEGER :: j, k
!   
    DO j = 1, modulos
      DO k = 1, verticesModulo
        acoplamentoExcitLocal = F_AcoplamentoExcitatorio(k, j, indiceNeuronio(j, k), &
            listaNeuronioNeuronio(j, k, :), listaNeuronioModulo(j, k, :), X(i, :), &
            PesoSyn(j, k, 1:verticesModulo + indiceNeuronio(j, k))/g)
!
        potencialSaida(j, k) = F_EvolucaoPotencial(potencialEntrada(j, k), mu, curr, & 
                                                    acoplamentoExcitLocal, BTEST(X(i, k), j-1)) 
!       
        !Gama(j, k) = F_Gama(Gama(j, k), TauTheta, BTEST(X(i, k), j-1))
        !Gama(j, k) = F_Gama(Gama(j, k), A, uTheta, TauTheta, BTEST(X(i, k), j-1))
        Gama(j, k) = 1.d0
        Phi = F_Phi(potencialSaida(j, k), Gama(j, k))
!
        IF(RAND() < Phi) X(i+1, k) = IBSET(X(i, k), j-1)
!   
        !CALL SE_AtualizaW(indiceNeuronio(j, k), PesoSyn(j, k, 1:verticesModulo + indiceNeuronio(j, k)), & 
        !        k, X(i+1, :), j) 
!
        !WRITE(22, '(3I15, 3F15.3, L15, I15, 2F15.3)') i, j, k, potencialEntrada(j, k), Gama(j, k), &     
        !                   Phi, BTEST(X(i, k), j-1), indiceNeuronio(j, k), RAND(), acoplamentoExcitLocal 
!
      ENDDO    
    ENDDO      
  END SUBROUTINE SE_EvolucaoDoSistema!---------------------------------------------------------------------------------------------     
!
  SUBROUTINE SE_AtualizaW(indiceLocal, WLocal, verticeA, estado, moduloA)!---------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: indiceLocal, verticeA, moduloA
    INTEGER, INTENT(IN), DIMENSION(verticesModulo) :: estado 
    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(verticesModulo + indiceLocal) :: WLocal
    INTEGER :: i, j, k, l
!
    DO i = 1, verticesModulo
      IF(i == verticeA) CYCLE
!
      WLocal(i) = W 
    ENDDO
!
    DO i = 1, indiceLocal
      j = i + verticesModulo
      k = listaNeuronioModulo(moduloA, verticeA, i)
      l = listaNeuronioNeuronio(moduloA, verticeA, i)
!
      WLocal(j) = W  
    ENDDO  
!
  END SUBROUTINE SE_AtualizaW!-----------------------------------------------------------------------------------------------------
!
  SUBROUTINE SD_CarregaRede!-------------------------------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: i, j, conexoesTotal, moduloLido, verticeLido
!
    OPEN(12, FILE=TRIM(nomeLista), STATUS='old')
    READ(12, *) modulos, verticesModulo, conexoesTotal
    IF(modulos > 1) THEN
      ALLOCATE(indiceNeuronio(modulos, verticesModulo))
      ALLOCATE(listaNeuronioModulo(modulos, verticesModulo, conexoesTotal))
      ALLOCATE(listaNeuronioNeuronio(modulos, verticesModulo, conexoesTotal))
!
      DO i = 1, modulos
        DO j = 1, verticesModulo
          READ(12, *) moduloLido, verticeLido, indiceNeuronio(moduloLido, verticeLido), &
              listaNeuronioModulo(moduloLido, verticeLido, 1:indiceNeuronio(moduloLido, verticeLido))
          READ(12, *) moduloLido, verticeLido, indiceNeuronio(moduloLido, verticeLido), &
              listaNeuronioNeuronio(moduloLido, verticeLido, 1:indiceNeuronio(moduloLido, verticeLido))        
        ENDDO      
      ENDDO
    ENDIF
    CLOSE(12)
  END SUBROUTINE SD_CarregaRede!---------------------------------------------------------------------------------------------------
!
  SUBROUTINE SD_SalvaResultados!--------------------------------------------------------------------------------------------------- 
    IMPLICIT NONE
    INTEGER, DIMENSION(contaSilencio/2) :: s, t
    INTEGER :: i, j, k, sMax, sMin, sSum, tMax, tMin, tSum, contagemT
    DOUBLE PRECISION :: mediaS
!   
    !OPEN(15, FILE='RASTER_PLOT_DINAMICA_GL_l3_N8192_CRITICO.dat', STATUS='unknown')
    OPEN(16, FILE='ATIVIDADE_REDE_H1_N16384_CG.dat', STATUS='unknown')
    OPEN(17, FILE='AVALANCHES_H1_N16384_CG.dat', STATUS = 'unknown')
    OPEN(18, FILE='DURACAO_AVALANCHES_H1_N16384_CG.dat', STATUS = 'unknown')
    OPEN(19, FILE='RELAÇÃO_SETHNA_H1_N16384_CG.dat', STATUS = 'unknown')
    
    !DO i = 1, modulos
    !  DO j = 1, verticesModulo
    !    DO k = 0, steps
    !      WRITE(15, '(B2.1)', ADVANCE='no') IBITS(X(k, j), i-1, 1) 
    !    ENDDO   
    !    WRITE(15, *)
    !  ENDDO    
    !ENDDO  
    
    DO i = 1, steps
      WRITE(16, *) i, DFLOAT(SUM(POPCNT(X(i, :))))/DFLOAT(modulos*verticesModulo)
      !DO j = 1, modulos
      !  WRITE(19, '(I5, 2F9.6)', ADVANCE = 'no') i, DFLOAT(COUNT(BTEST(X(i, :), j-1)))/DFLOAT(verticesModulo), &
      !                                           SUM(Gama(i, j, :))/DFLOAT(verticesModulo)
      !ENDDO
      !WRITE(19, *)
    ENDDO
!   !  
    DO i = 1, contaSilencio/2
      s(i) = SUM(POPCNT(X(tempoSilencio(2*i-1) : tempoSilencio(2*i) - 1, :)))
      t(i) = tempoSilencio(2*i) - tempoSilencio(2*i-1)
    ENDDO
!   
    sMin = MINVAL(s)
    sMax = MAXVAL(s)
    sSum = SUM(s)
!
    tMin = MINVAL(t)
    tMax = MAXVAL(t)
    tSum = SUM(t)
!
    DO i = sMin, sMax
      IF( COUNT(s == i) == 0) CYCLE
      WRITE(17, *) i, COUNT(s == i), DFLOAT(COUNT(s == i))/DFLOAT(sSum)    
    ENDDO
!
    DO i = tMin, tMax
      IF( COUNT(t == i) == 0) CYCLE
      WRITE(18, *) i, COUNT(t == i), DFLOAT(COUNT(t == i))/DFLOAT(tSum)    
    ENDDO
!
    DO i = tMin, tMax
      contagemT = COUNT(t == i)
      IF(contagemT == 0) CYCLE
      mediaS = DFLOAT(SUM(s, MASK = (t == i)))/DFLOAT(contagemT)    
      WRITE(19, *) i, mediaS
    ENDDO
!
    !CLOSE(15)
    CLOSE(16)
    CLOSE(17)
    CLOSE(18)  
    CLOSE(19)
  END SUBROUTINE SD_SalvaResultados!--------------------------------------------------------------------------------------------------------

END PROGRAM dinamicaGL_TCT_HMN
