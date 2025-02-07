! PROGRAMA GERA UM GRAFO ALEATÓRIO E SUA MATRIZ DE ADJACÊNCIA n x n 
! eps É O PARÂMETRO DE CONECTIVIDADE ENTRE OS VÉRTICES
! 
! 
!
! Neuromat - USP  janeiro de 2021
! 
! Flavio Rusch
!==================================================================================================================================
PROGRAM random_graph
IMPLICIT NONE
INTEGER, PARAMETER :: nNeuron = 16384, hMax = 1
DOUBLE PRECISION, PARAMETER :: pi = dacos(-1.d0)
INTEGER, DIMENSION(nNeuron, nNeuron)  :: conecIntern, matriz, Pot 
INTEGER, DIMENSION(nNeuron, nNeuron/10) :: listConect, tipoConex          !Com quem esta conectado | nNeuron/10 max de conexões que
INTEGER, DIMENSION(nNeuron) :: indice = 0                                 !1 neuro pode ter
DOUBLE PRECISION, DIMENSION(nNeuron, 2) :: posicao                        !Qtas (número de) conexoes o n-esimo neuronio tem   
DOUBLE PRECISION, DIMENSION(2**hMax, 2) :: centroModulo
DOUBLE PRECISION :: r, V, eps = 0.01d0, probInib = 0.2d0, Rex = 0.5d0, distancia 
INTEGER :: i, j, d, f, k, l, h, elementosModulo, inicio, fim, m, lateralizacao, verticalizacao 

OPEN(1, FILE='MATRIZ_ADJ_N2048_eps=0.01.dat', STATUS='UNKNOWN')
OPEN(2, FILE='GRAFO.dat', STATUS='UNKNOWN')
OPEN(3, FILE='RANDOM_GRAPH_N2048_eps=0.01.dat', STATUS='UNKNOWN')
OPEN(4, FILE='lista_simples_16384_H1.dat', STATUS='UNKNOWN')
OPEN(5, FILE='lista_de_listas_16384_H1.dat', STATUS='UNKNOWN')
OPEN(7, FILE='POSICAO.dat', STATUS='UNKNOWN')


DO i = 1, nNeuron - 1
  DO j = i + 1, nNeuron
    IF (RAND() <= eps) THEN
      indice(i) = indice(i) + 1
      indice(j) = indice(j) + 1
      listConect(i, indice(i)) = j 
      listConect(j, indice(j)) = i      
      IF (RAND() < probInib) THEN
        tipoConex(i, indice(i)) = -1
        tipoConex(j, indice(j)) = -1
        !WRITE(4, '(2I5, I3)') i, j, -1 
        !WRITE(4, '(2I5, I3)') j, i, -1
        conecIntern(i, j) = -1   !INIBITORIO
      ELSE
        tipoConex(i, indice(i)) = 1
        tipoConex(j, indice(j)) = 1
        !WRITE(4, '(2I5, I3)') i, j, 1 
        !WRITE(4, '(2I5, I3)') j, i, 1
        conecIntern(i, j) = 1    !EXCITATORIO   
      END IF
    ELSE
      conecIntern(i, j) = 0      !S/ CONECT
    ENDIF
  ENDDO
ENDDO

DO i = 1, nNeuron
  IF (indice(i) == 0) THEN
    k = i
    DO WHILE (k == i)
      k = 1 + FLOOR(DFLOAT(nNeuron)*RAND())
    ENDDO
    indice(i) = 1
    listConect(i, 1) = k
    indice(k) = indice(k) + 1
    listConect(k, indice(k)) = i 
  END IF
ENDDO

elementosModulo = nNeuron
DO h = 1, hMax
  elementosModulo = elementosModulo/2
  fim = -elementosModulo
  DO m = 1, 2**(h - 1)
    inicio = fim + elementosModulo + 1
    fim = inicio + elementosModulo - 1  
    DO i = inicio, fim
      !WRITE(*, *) h, i, inicio, fim
      IF (indice(i) == 0) CYCLE
      DO j = 1, indice(i)
        IF (listConect(i, j) > fim .AND. listConect(i, j) <= fim + elementosModulo) THEN !Garante conexões somente dentro do modulo
          IF (tipoConex(i, j) == -1 .OR. RAND() <= Rex) THEN
            k = i 
            DO WHILE (k == i .OR. ANY(k == listConect(i, 1:indice(i))))
              k = inicio + FLOOR(DFLOAT(elementosModulo)*RAND())            
            ENDDO
!
            l = listConect(i, j)
            listConect(i, j) = k
            indice(k) = indice(k) + 1
            listConect(k, indice(k)) = i                                                        !Ligação entre o neurônio i e k
            tipoConex(k, indice(k)) = tipoConex(i, j)
            
            k = l
            DO WHILE (k == l .OR. ANY(k == listConect(l, 1:indice(l))))
              k = ((l-1)/elementosModulo)*elementosModulo + 1 + FLOOR(DFLOAT(elementosModulo)*RAND())            
            ENDDO
            listConect(l, MINLOC(listConect(l, 1:indice(l)), DIM = 1, MASK = (listConect(l, 1:indice(l)) == i))) = k
            indice(k) = indice(k) + 1
            listConect(k, indice(k)) = l                                                        !Ligação entre o neurônio j e k
            tipoConex(k, indice(k)) = tipoConex(i, j)
          ENDIF
        ENDIF
      ENDDO
    ENDDO  
  ENDDO
ENDDO


DO i = 1, nNeuron
  DO j = 1, indice(i)
    WRITE(4, '(3I5)') i, listConect(i, j), tipoConex(i, j) 
  ENDDO  
ENDDO

WRITE(5, '(2I7)') nNeuron, 2**hMax
DO i = 1, nNeuron
  WRITE(5, '(2I7)', ADVANCE='no') i, indice(i)
  DO j = 1, indice(i)
    WRITE(5, '(I7)', ADVANCE='no') listConect(i, j)
  ENDDO  
  DO j = 1, indice(i)
    WRITE(5, '(I7)', ADVANCE='no') tipoConex(i, j) 
  ENDDO  
  WRITE(5, *)
ENDDO

DO i = 1, nNeuron
  posicao(i, 1) = RAND()
  posicao(i, 2) = 2*pi*RAND()
ENDDO

centroModulo(1, :) = 0.d0
distancia = 2.d0/(2**hmax)

IF (hmax > 0) THEN
  lateralizacao = ((hmax + 1)/2)*2 
  verticalizacao = (hmax/2)*2  
  
  k = 0
  DO i = 1, lateralizacao
    
    DO j = 1, verticalizacao
      k = k + 1
      centroModulo(k, 1) = DFLOAT(i)
      centroModulo(k, 2) = DFLOAT(j)    
    ENDDO   
  ENDDO
  
ENDIF

posicao(:, 1) = posicao(:, 1)*distancia*0.9d0

DO i = 1, nNeuron
  IF (indice(i) == 0) CYCLE
  DO j = 1, indice(i)
    k = (i - 1)*2**hMax/nNeuron + 1
    WRITE(2, *) posicao(i, 1)*dcos(posicao(i, 2)) + centroModulo(k, 1), posicao(i, 1)*dsin(posicao(i, 2)) + &
                centroModulo(k, 2), tipoConex(i, j)
    k = (listConect(i, j) - 1)*2**hmax/nNeuron + 1
    WRITE(2, *) posicao(listConect(i,j), 1)*dcos(posicao(listConect(i,j), 2)) + centroModulo(k, 1), & 
                posicao(listConect(i,j), 1)*dsin(posicao(listConect(i,j), 2)) + centroModulo(k, 2), tipoConex(i, j)
    WRITE(2,*)
  ENDDO
    WRITE(7, *) posicao(i, 1)*dcos(posicao(i, 2)) + centroModulo(k, 1), posicao(i, 1)*dsin(posicao(i, 2)) + &
                centroModulo(k, 2)
ENDDO

CLOSE(1)
CLOSE(2)
CLOSE(3)
CLOSE(4)
CLOSE(5)
CLOSE(7)

END PROGRAM random_graph
