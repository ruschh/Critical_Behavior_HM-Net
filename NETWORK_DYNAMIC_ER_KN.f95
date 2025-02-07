!PROGRAMA dinamica_gl.f95 IMPLEMENTA A DINAMICA ESTOCÁSTICA DOS NEURONIOS GL 
! 
!       ESSA VERSÃO CALCULA O TAMANHO DAS AVALANCHES, AS DURAÇÕES E A RELAÇÃO DE SETHNA DE
!       UMA REDE HM COM TOPOLOGIA DEPENDENTE DO ARQUIVO DE ENTRADA QUE REPRESENTA A REDE DE
!       GRAFOS 
!
!       VERSAO OFICIAL -- RESULTADOS DESTE CÓDIGO SERÃO APRESENTADOS NO RPB2023 EM PARIS
!
!Flavio Rusch
!
!NEUROMAT - USP
!FEV/2023
!
!==================================================================================================================================
PROGRAM dinamicaGL
IMPLICIT NONE
INTEGER, PARAMETER :: transiente = 1.d4, iteracoes = 1.d4
DOUBLE PRECISION, PARAMETER :: A = 1.d0, tau = 1.d4, u = 1d-1, mu = 0.d0, curr = 0.d0, g = 0.d0
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: raster, rhoModulo 
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: rhoMedio, Var, ISImedio, varISI
DOUBLE PRECISION, DIMENSION(iteracoes) :: rho, meanW, GamaMean
DOUBLE PRECISION :: somaW = 0.d0, gamaMax = 1.d0, wMax = 1.d0, r
INTEGER :: nNeuron, nModulo, idum = -936782451 

TYPE caractNeuron
  INTEGER :: indice, nSpike                                                                !Número de ligações de cada neuronio
  INTEGER(KIND = 1), DIMENSION(iteracoes) :: disparos, tempoSpike = 0
  INTEGER, DIMENSION(:), ALLOCATABLE :: listaConexoes, tipoConexoes
  DOUBLE PRECISION :: V, phi, gama
  DOUBLE PRECISION, DIMENSION(iteracoes) :: GuardaW, GuardaGama, X
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: W
END TYPE caractNeuron

TYPE(caractNeuron), DIMENSION(:), ALLOCATABLE :: neuronio
!
!CALL SYSTEM ('clear')

neuronio%nSpike = 0
CALL CarregaLista
CALL condicoesIniciais
!
ALLOCATE(raster(nNeuron, iteracoes), rhoModulo(nModulo, iteracoes))
ALLOCATE(rhoMedio(nModulo), Var(nModulo), ISImedio(nModulo), varISI(nModulo)) 

CALL evoluiSistema
!CALL estatisticas
!CALL determinaISI
!CALL salvaEstatistica
!CALL salvaAtividade
CALL salvaResultados
!CALL gif

!WRITE(*, '(A20)') 'Programa concluído!'
  
CONTAINS
  SUBROUTINE CarregaLista!=========================================================================================================
    IMPLICIT NONE
    INTEGER :: i, j
    
    OPEN(1, FILE='lista_de_listas_8192_H0.dat', STATUS='old')
    READ(1, *) nNeuron, nModulo 
    ALLOCATE(neuronio(nNeuron))
        
    DO i = 1, nNeuron
      READ(1, '(2I7)', ADVANCE = 'no') j, neuronio(i)%indice
      ALLOCATE (neuronio(i)%listaConexoes(neuronio(i)%indice))
      ALLOCATE (neuronio(i)%tipoConexoes(neuronio(i)%indice))
      ALLOCATE (neuronio(i)%W(neuronio(i)%indice))
      READ(1, *) neuronio(i)%listaConexoes, neuronio(i)%tipoConexoes
    ENDDO
    CLOSE(1) 
!
    !ALOCA A DIMENSÃO DE caractNeuron, X DEPENDENTE DO TEMPO E DOS NEURONIOS
    !DO i = 1, nNeuron
    !  ALLOCATE(neuronio(i)%X(iteracoes))
    !ENDDO
!     
  !  WRITE(*, '(A17)') 'Rede carregada...'
  END SUBROUTINE CarregaLista!=====================================================================================================
!  
  SUBROUTINE condicoesIniciais!====================================================================================================
    IMPLICIT NONE
    INTEGER :: i, jj
!
    DO i = 1, nNeuron
      neuronio(i)%X = 0.d0
      neuronio(i)%V = RAND()
      neuronio(i)%W = 1.d0
      neuronio(i)%gama = RAND()*gamaMax           
    ENDDO  
    !WRITE(*, '(A35)') 'Condições iniciais estabelecidas...'
  ENDSUBROUTINE condicoesIniciais!=================================================================================================
!  
  SUBROUTINE evoluiSistema!========================================================================================================
    IMPLICIT NONE
    INTEGER :: i, j, jj, k, l, ll, neuroniosModulo
    INTEGER :: sMax, sMin, sSum, tMax, tMin, tSum, contagemT, contaSilencio
    INTEGER, DIMENSION(nNeuron) :: tempArray 
    INTEGER, DIMENSION(iteracoes) :: s, t, tempoZeroRho
    DOUBLE PRECISION :: mediaS
    LOGICAL :: avalanche

   ! WRITE(*, '(A22)') 'Evoluindo o sistema...'
   !OPEN(43, FILE = 'VOLTAGEM_ESPECTRO_POTENCIA_ER_H2_N32768_STATIC.dat', STATUS = 'unknown')
   !OPEN(43, FILE = 'SPIKE_TRAIN_ESPECTRO_POTENCIA_ER_H2_N32768_STATIC.dat', STATUS = 'unknown')
!    
    i = 0
    j = 1
    k = 0
    contaSilencio = 0
    avalanche = .FALSE.
!
    DO WHILE(i < iteracoes)
!   
      !CALL determina_W(i)               !P/ O CÁLCULO C/ MEC. HOMEOSTÁTICO
      CALL determinaGamma(i)
   
      rho(i) = SUM(neuronio(:)%X(i))/DFLOAT(nNeuron)
      raster(:, i) = SUM(neuronio(:)%X(i))             
      IF(rho(i) == 0.d0) THEN
        r = ran1(idum)                                        !CONDIÇÃO INICIAL PARA SE TER UM NEURONIO ALEATÓRIO DISPARANDO NA REDE
        l = INT((nNeuron - 1.d0)*r + 1.d0) 
        neuronio(l)%X(i) = 1.d0
        rho(i) = 1.d0/DFLOAT(nNeuron)
        IF(avalanche) THEN
          contaSilencio = contaSilencio + 1
          tempoZeroRho(contaSilencio) = i 
        ENDIF
        
        avalanche = .FALSE.
      ELSEIF(.NOT. avalanche) THEN
        avalanche = .TRUE.
        contaSilencio = contaSilencio + 1
        tempoZeroRho(contaSilencio) = i
      ENDIF
      
      CALL calculaV(i)
      CALL determinaPhi
      !CALL determina_X
!
      !WRITE(43, *) SUM(neuronio(:)%V)
      !WRITE(44, *) raster(:, i)
!     
      s(contaSilencio) = s(contaSilencio) + SUM(neuronio(:)%X(i))
      i = i + 1 
!
      CALL determina_X(i)
!
    ENDDO
!    
    OPEN(66, FILE='ATIVIDADE_ESPECTRO_POTENCIA_ER_H2_N32768_STATIC.dat', STATUS='unknown') 
    OPEN(77, FILE='RHO_REDE_ESTATICA_H2_N32768_STATIC.dat', STATUS='unknown') 
    DO i = 1, iteracoes
      WRITE(66, *) rho(i)
      WRITE(77, *) i, rho(i)
    ENDDO
    CLOSE(77)
!
    OPEN(99, FILE='AVALANCHES_N8192_H0_HOMEOSTATIC.dat', STATUS='unknown') 
    OPEN(100, FILE='DURACAO_AVALANCHES_N8192_H0_HOMEOSTATIC.dat', STATUS='unknown') 
    OPEN(103, FILE='SETHNA_RELATION_N8192_H0_HOMEOSTATIC.dat', STATUS='unknown') 
!
    ll = 0
    jj = MAXVAL(INT(s(1:contaSilencio)))
    DO i = 0, jj
      l = COUNT(INT(s(1:contaSilencio)) == ll)
      IF(l == 0) THEN
        ll = ll + 1
        CYCLE
      ENDIF
      WRITE(99, *) ll, DFLOAT(l)/SUM(s(1:contaSilencio))   
      ll = ll + 1 
    ENDDO
!
    ll = 0
    jj = MAXVAL(INT(tempoZeroRho(1:contaSilencio)))
    DO i = 1, jj
      l = COUNT(tempoZeroRho(2:contaSilencio) - tempoZeroRho(1:contaSilencio - 1) == ll)
      IF(l == 0) THEN
        ll = ll + 1
        CYCLE
      ENDIF
      !WRITE(100, *) i, COUNT(INT(tempoZeroRho(2:kk) - tempoZeroRho(1:kk-1)) == i)/SUM(tempoZeroRho(:))    
      WRITE(100, *) ll, DFLOAT(l)/DFLOAT(contaSilencio - 1)
      ll = ll + 1
    ENDDO
!
    s = 0.d0
    DO i = 1, contaSilencio/2
      !s(i) = SUM(neuronio(:)%X(tempoZeroRho(2*i-1) : tempoZeroRho(2*i)-1))
      FORALL(j = 1: nNeuron) tempArray(j) = COUNT(neuronio(j)%X(tempoZeroRho(2*i-1) : tempoZeroRho(2*i)-1) == 1)
      s(i) = SUM(tempArray)
      t(i) = tempoZeroRho(2*i) - tempoZeroRho(2*i-1)
    ENDDO
!     
    sMin = MINVAL(s)
    sMax = MAXVAL(s)
    sSum = SUM(s)
!
    tMin = MINVAL(t)    !tempoZeroRho
    tMax = MAXVAL(t) !tempoZeroRho
    tSum = SUM(t)       !tempoZeroRho 
!
    !DO i = sMin, sMax
    !  IF(COUNT(s == i) == 0) CYCLE
    !  WRITE(99, *) i, DFLOAT(COUNT(s == i))/DFLOAT(sSum)
    !ENDDO
!
    !DO i = tMin, tMax
    !  IF(COUNT(t == i) == 0) CYCLE
    !  WRITE(100, *) i, DFLOAT(COUNT(t == i))/DFLOAT(sSum)
    !ENDDO
!
    DO i = tMin, tMax
      contagemT = COUNT(t == i)
      IF(contagemT == 0) CYCLE
      !mediaS = DFLOAT(SUM(s, MASK = (tempoZeroRho == i)))/DFLOAT(contagemT)    
      mediaS = DFLOAT(SUM(s, MASK = (t == i)))/DFLOAT(contagemT)    
      WRITE(103, *) i, mediaS
    ENDDO
!
    CLOSE(99)
    CLOSE(100)
    CLOSE(103)
!        
  ENDSUBROUTINE evoluiSistema!=====================================================================================================
!  
  SUBROUTINE estatisticas
    IMPLICIT NONE
    INTEGER :: i
    
    !WRITE(*, '(A27)') 'Calculando a estatística...'
    
    DO i = 1, nModulo
      rhoMedio(i) = SUM(rhoModulo(i,:))/DFLOAT(iteracoes)
      Var(i) = SUM((rhoModulo(i,:) - rhoMedio(i))**2.d0)/DFLOAT(iteracoes - 1)
    ENDDO    
  ENDSUBROUTINE estatisticas
  
!   SUBROUTINE salvaAtividade
!     IMPLICIT NONE
!     INTEGER :: i, j
!     
!     OPEN(5, FILE='ATIVIDADE_MODULOS_H0_A1.10_I00_GAMA_MAX_1.0.dat', STATUS='UNKNOWN')
!     
!     WRITE(*, '(A35)') 'Salvando a atividade dos módulos...'
!     
!     DO i = 1, iteracoes
!       WRITE(5, '(I6, F10.6, F10.6)', ADVANCE='no') i, rho(i), GamaMean(i)
!       DO j = 1, nModulo      
!         WRITE(5, '(F10.6)', ADVANCE='no') rhoModulo(j,i)
!       ENDDO
!       WRITE(5, *)
!     ENDDO
!     
!     CLOSE(5)
!   ENDSUBROUTINE
  
!   SUBROUTINE determinaISI
!     IMPLICIT NONE
!     INTEGER, DIMENSION(nNeuron, iteracoes) :: intervalo 
!     DOUBLE PRECISION :: soma
!     INTEGER :: i, j, neuroniosModulo, denominador  
!     
!     WRITE(*, '(A36)') 'Calculando o Interval Insterspike...'
!     
!     intervalo = 0
!     DO i = 1, nNeuron
!       DO j = 1, neuronio(i)%nSpike - 1
!         intervalo(i, j) = neuronio(i)%tempoSpike(j + 1) - neuronio(i)%tempoSpike(j)
!       ENDDO      
!     ENDDO    
!     
!     neuroniosModulo = nNeuron/nModulo
!     DO i = 1, nModulo
!       soma = 0.d0
!       denominador = 0
!       DO j = (i - 1)*neuroniosModulo + 1, i*neuroniosModulo
!         denominador = denominador + neuronio(j)%nSpike - 1
!         soma = soma + DFLOAT(SUM(intervalo(j, 1:neuronio(j)%nSpike - 1)))
!       ENDDO
!       ISImedio(i) = soma/DFLOAT(denominador)
!       
!       soma = 0.d0
!       DO j = (i - 1)*neuroniosModulo + 1, i*neuroniosModulo
!         soma = soma + SUM((DFLOAT(intervalo(j, 1:neuronio(j)%nSpike - 1)) - ISImedio(i))**2.d0) 
!       ENDDO
!       varISI(i) = soma/DFLOAT(denominador - 1)
!     ENDDO    
!   ENDSUBROUTINE determinaISI
  
!   SUBROUTINE salvaEstatistica
!     IMPLICIT NONE
!     INTEGER :: i
!     
!     OPEN(4, FILE='ESTATISTICA_MODULOS_H0_A1.70_I025.dat', STATUS='UNKNOWN')
!     
!     WRITE(4, '(4A10, A8, A13, A6)') '# Modulo', 'Rho_Medio', 'Variancia', 'Sincronia', 'ISI', 'DP do ISI', 'CV'
!     WRITE(4, *)
!     DO i = 1, nModulo
!       WRITE(4, '(I10, 6F10.6)') i, rhoMedio(i), Var(i), Var(i)/rhoMedio(i), ISImedio(i), DSQRT(varISI(i)) &
!             , DSQRT(varISI(i))/ISImedio(i)
!     ENDDO  
!     
!     CLOSE(4)
!   ENDSUBROUTINE salvaEstatistica
  
  SUBROUTINE salvaResultados
    IMPLICIT NONE
    INTEGER :: i, j, k
    !WRITE(*, '(A22)') 'Salvando resultados...'
    
    OPEN (1, FILE='RASTER_PLOT_H2_N32768_STATIC.dat', STATUS='UNKNOWN')
    OPEN(2, FILE='GRANDEZAS_H2_N32768_STATIC.dat', STATUS='UNKNOWN')
    OPEN(3, FILE='MATRIZ_W_MEDIO_H0_A1.70_I00_MU00.dat  ', STATUS='UNKNOWN')
    
    DO i = 1, nNeuron
      WRITE(1, '(5000F3.0)') raster(i, :)
    ENDDO
    
    DO i = 1, iteracoes
      WRITE(2, '(I5, F9.4, F9.4)', ADVANCE ='no') i, rho(i), GamaMean(i)
      DO j = 1, nNeuron
        WRITE(2, '(2F9.4)', ADVANCE='no') neuronio(j)%GuardaGama(i), neuronio(j)%GuardaW(i)           
        WRITE(3, '(F9.4)', ADVANCE='no') neuronio(j)%GuardaW(i)                                  !neuronio(j)%GuardaGama(i)
      ENDDO
      WRITE(2, *)
      WRITE(3, *)
    ENDDO
    
    DO i = 1, 3
      CLOSE(i)
    ENDDO
  ENDSUBROUTINE salvaResultados

  SUBROUTINE determinaGamma(tAtual)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: tAtual
    INTEGER :: i
    
    DO i = 1, nNeuron
      !neuronio(i)%gama = gamaMax 
      neuronio(i)%gama = neuronio(i)%gama + (A - neuronio(i)%gama)/tau - u*neuronio(i)%gama*neuronio(i)%X(tAtual)
    ENDDO  
  END SUBROUTINE determinaGamma
  
  SUBROUTINE determinaPhi
    IMPLICIT NONE
    INTEGER :: i
    
    DO i = 1, nNeuron
      neuronio(i)%phi = 1.d0 - 1.d0/(1.d0 + neuronio(i)%gama*neuronio(i)%V)          
    ENDDO  
  END SUBROUTINE determinaPhi
  
  SUBROUTINE determina_W(tAtual)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: tAtual
    INTEGER :: i, j
    
    DO i = 1, nNeuron
      DO j = 1, neuronio(i)%indice
        !neuronio(i)%W(j) = wMax
        neuronio(i)%W(j) = neuronio(i)%W(j) + (A - neuronio(i)%W(j))/tau - &
        u*neuronio(i)%W(j)*neuronio(neuronio(i)%listaConexoes(j))%X(tAtual)
        !IF(neuronio(i)%tipoConexoes(j) == -1) neuronio(i)%W(j) = -g*neuronio(i)%W(j)             !INIBIU
      ENDDO    
    ENDDO     
  ENDSUBROUTINE determina_W
  
  SUBROUTINE determina_X(tAtual)
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: tAtual
    INTEGER :: i
    
    DO i = 1, nNeuron
      IF (RAND() < neuronio(i)%phi) THEN
        neuronio(i)%X(tAtual) = 1
      ELSE
        neuronio(i)%X(tAtual) = 0
      ENDIF
    ENDDO
    
  ENDSUBROUTINE determina_X
  
  SUBROUTINE guardaHistorico(tempo, tAtual)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: tempo, tAtual
    INTEGER :: i
    
    neuronio%disparos(tempo) = INT(neuronio%X(tAtual)) 
    DO i = 1, nNeuron
      IF(neuronio(i)%disparos(tempo) == 1) THEN
        neuronio(i)%nSpike = neuronio(i)%nSpike + 1
        neuronio(i)%tempoSpike(neuronio(i)%nSpike) = tempo
      ENDIF
    ENDDO
    
  ENDSUBROUTINE guardaHistorico
  
!   SUBROUTINE gif
!     IMPLICIT NONE
!     DOUBLE PRECISION, DIMENSION(nNeuron, 2) :: posicao 
!     INTEGER :: i, j
!     
!     OPEN(4, FILE='POSICAO.dat', STATUS='old')
!     OPEN(5, FILE='GIF_PONTOS_HMN3_A1.10_I010.dat', STATUS='UNKNOWN')
!     OPEN(7, FILE='GIF_CONEXOES_HMN3_A1.10_I010.dat', STATUS='UNKNOWN')
!     
!     DO i = 1, nNeuron
!       READ(4, *) posicao(i, :)
!     ENDDO
!     CLOSE(4)
!     
!     DO i = 1, nNeuron
!       WRITE(5, *) posicao(i, :), neuronio(i)%disparos    
!     ENDDO
!     
!     DO i = 1, nNeuron
!       DO j = 1, neuronio(i)%indice
!         WRITE(7, *) posicao(i, :), neuronio(i)%disparos*neuronio(i)%tipoConexoes(j)        
!         WRITE(7, *) posicao(neuronio(i)%listaConexoes(j), :), neuronio(i)%disparos*neuronio(i)%tipoConexoes(j)
!         WRITE(7, *)
!       ENDDO
!     ENDDO
!     
!     CLOSE(5)
!     CLOSE(7)
!   ENDSUBROUTINE gif
  
  SUBROUTINE calculaV(tAtual)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: tAtual
    INTEGER :: i
    
    DO i = 1, nNeuron
      IF (neuronio(i)%X(tAtual) == 1) THEN
        neuronio(i)%V = 0.d0
      ELSE
        neuronio(i)%V = mu*neuronio(i)%V + curr + &
                        SUM(neuronio(i)%W(:)*neuronio(neuronio(i)%listaConexoes(:))%X(tAtual))/(DFLOAT(neuronio(i)%indice))
      ENDIF
!       WRITE(*,*) i, neuronio(i)%V, &
!       SUM(neuronio(i)%W(:)*neuronio(neuronio(i)%listaConexoes(:))%X)/(DFLOAT(neuronio(i)%indice)), &
!       ANY(neuronio(neuronio(i)%listaConexoes(:))%X == 1)
    ENDDO  
  ENDSUBROUTINE calculaV

  DOUBLE PRECISION FUNCTION ran1(input)
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: input
  INTEGER, PARAMETER :: ia = 16807, im = 2147483647, iq = 127773, ir = 2836
  INTEGER, PARAMETER :: ntab = 32, ndiv = 1 + (im-1)/ntab
  DOUBLE PRECISION, PARAMETER :: eps = 1.2d-7, rnmx = 1.d0 - eps, am = 1.d0/FLOAT(im)
  INTEGER :: j, k, iv(ntab), iy
  SAVE :: iv, iy
  DATA iv/ntab*0/, iy/0/
!
  IF (input <= 0 .OR. iy == 0) THEN
    input = MAX(-input, 1)
    DO j = ntab + 8, 1, -1
      k = input/iq
      input = ia*(input - k*iq) - ir*k
      IF (input < 0) input = input + im
      IF (j < ntab) iv(j) = input
    END DO
    iy = iv(1)
  END IF
  k = input/iq
  input = ia*(input - k*iq) - ir*k
  IF (input < 0) input = input + im
  j = 1 + iy/ndiv
  iy = iv(j)
  iv(j) = input
  ran1 = MIN(am*iy,rnmx)
  RETURN
  END FUNCTION ran1 
  
END PROGRAM dinamicaGL




