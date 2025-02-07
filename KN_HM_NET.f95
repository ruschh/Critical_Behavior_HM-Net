PROGRAM GeraRede_AlgumNome
  IMPLICIT NONE
  INTEGER, PARAMETER :: neuronios = 1024, vizinhos = 2
  DOUBLE PRECISION, PARAMETER :: probabilidadeInibitorio = 2.d-1
  INTEGER, DIMENSION(neuronios) :: indiceConexoes
  INTEGER, DIMENSION(neuronios,2*vizinhos) :: listaConexoes, listaTipoConexao
  INTEGER :: ii, jj
!
  indiceConexoes = 0
  CALL ConectaVerticesAnel
  CALL NovasConexoes  
  CALL CorrigeDuplicatas
  CALL EstabeleceTipoConexao
  CALL SalvaLista('lista_N1024_K2')
  CALL SalvaMatrizAdjacencia
CONTAINS
  SUBROUTINE ConectaVerticesAnel!--------------------------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: i, j, k
!
    DO i = 1, neuronios
      DO j = 1, vizinhos
        k = i + j
        IF (k > neuronios) k = k - neuronios
        indiceConexoes(i) = indiceConexoes(i) + 1
        indiceConexoes(k) = indiceConexoes(k) + 1
        listaConexoes(i,indiceConexoes(i)) = k
        listaConexoes(k,indiceConexoes(k)) = i
      END DO
    END DO
  END SUBROUTINE ConectaVerticesAnel
!
  SUBROUTINE NovasConexoes!--------------------------------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, DIMENSION(neuronios) :: indiceNovasConexoes
    INTEGER :: conexaoNova, conexaoVelha, quebraConexao, i, j, k, l
!
    DO i = 1, neuronios
      DO j = 1, 2*vizinhos
        CALL SorteiaVertice(i, conexaoNova)
        IF (ANY(listaConexoes(i,:) == conexaoNova)) THEN
          k = MINLOC(listaConexoes(i,:), MASK = (listaConexoes(i,:) == conexaoNova), DIM = 1)
          IF (k /= j) THEN
            listaConexoes(i,k) = listaConexoes(i,j)
            listaConexoes(i,j) = conexaoNova
          END IF
          CYCLE
        END IF
!     
        conexaoVelha = listaConexoes(i,j)
        listaConexoes(i,j) = conexaoNova
!
        conexaoNova = conexaoVelha
        DO WHILE (conexaoNova == conexaoVelha)
          quebraConexao = FLOOR(2.d0*FLOAT(vizinhos)*RAND()) + 1
          conexaoNova = listaConexoes(listaConexoes(i,j),quebraConexao)
        END DO  
        listaConexoes(listaConexoes(i,j),quebraConexao) = i
!      
        k = MINLOC(listaConexoes(conexaoVelha,:), MASK = (listaConexoes(conexaoVelha,:) == i), DIM = 1)
        listaConexoes(conexaoVelha,k) = conexaoNova

        k = MINLOC(listaConexoes(conexaoNova,:), MASK = (listaConexoes(conexaoNova,:) == listaConexoes(i,j)), DIM = 1)
        listaConexoes(conexaoNova,k) = conexaoVelha
      END DO
    END DO
  END SUBROUTINE NovasConexoes
!
  SUBROUTINE CorrigeDuplicatas!-----------------------------------------------------------------------------------------------------
    IMPLICIT NONE  
    INTEGER :: posicaoVertice, verticeAuxiliar_1, verticeAuxiliar_2, i, j, k
!
    DO i = 1, neuronios - 1
      DO j = 1, 2*vizinhos - 1
        IF (COUNT(listaConexoes(i,:) == listaConexoes(i,j)) > 1) THEN
          !write(*, *) i, listaConexoes(i,j), listaConexoes(i,:)
          posicaoVertice = MINLOC(listaConexoes(listaConexoes(i,j),:), MASK = (listaConexoes(listaConexoes(i,j),:) == i), DIM = 1)
!
          DO
            CALL SorteiaVertice(i, verticeAuxiliar_1)
            IF (ANY(listaConexoes(i,:) == verticeAuxiliar_1)) CYCLE
!          
            DO k = 1, 2*vizinhos
              IF (ALL(listaConexoes(listaConexoes(i,j),:) /= listaConexoes(verticeAuxiliar_1,k)) .AND.&
                listaConexoes(i,j) /= listaConexoes(verticeAuxiliar_1,k)) EXIT
            END DO
            IF (k <= 2*vizinhos) EXIT
          END DO
!
          listaConexoes(listaConexoes(i,j),posicaoVertice) = listaConexoes(verticeAuxiliar_1,k)
!          
          posicaoVertice = MINLOC(listaConexoes(listaConexoes(verticeAuxiliar_1,k),:),&
            MASK = (listaConexoes(listaConexoes(verticeAuxiliar_1,k),:) == verticeAuxiliar_1), DIM = 1)
          listaConexoes(listaConexoes(verticeAuxiliar_1,k),posicaoVertice) = listaConexoes(i,j)
!
          listaConexoes(i,j) = verticeAuxiliar_1
          listaConexoes(verticeAuxiliar_1,k) = i
        END IF
      END DO
    END DO
  END SUBROUTINE CorrigeDuplicatas
!
  SUBROUTINE SorteiaVertice(vertice_1, vertice_2)!----------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: vertice_1
    INTEGER, INTENT(OUT) :: vertice_2
!
    vertice_2 = vertice_1
    DO WHILE (vertice_1 == vertice_2)
      vertice_2 = FLOOR(FLOAT(neuronios)*RAND()) + 1
    END DO
  END SUBROUTINE SorteiaVertice
!
  SUBROUTINE EstabeleceTipoConexao!-------------------------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: i, j, k, aux, l
    
    k = 0
    listaTipoConexao(:, 1:2*vizinhos) = 1
!
    DO i = 1, neuronios
      IF(COUNT(listaTipoConexao(i, 1:2*vizinhos) == -1) == 2) CYCLE
      DO j = 1, 2*vizinhos
        IF(COUNT(listaTipoConexao(listaConexoes(i, j), 1:2*vizinhos) == -1) == 2) CYCLE
        listaTipoConexao(i, j) = -1
        k = MINLOC(listaConexoes(listaConexoes(i, j), 1:2*vizinhos), DIM = 1, &  
        MASK = listaConexoes(listaConexoes(i, j), 1:2*vizinhos) == i)
        listaTipoConexao(listaConexoes(i, j), k) = -1
        IF(COUNT(listaTipoConexao(i, 1:2*vizinhos) == -1) == 2) EXIT       
      END DO
    END DO
!   
  END SUBROUTINE EstabeleceTipoConexao
!
  SUBROUTINE SalvaLista(nome)!-----------------------------------------------------------------------------------------------------
    IMPLICIT NONE
    CHARACTER(len=14), INTENT(IN) :: nome
    INTEGER :: i, j
!
    OPEN (101, FILE = nome//'.dat', STATUS = 'unknown')
    WRITE(101, '(2I6)') neuronios, MAXVAL(indiceConexoes)
    DO i = 1, neuronios
      WRITE(101, '(2I6)', ADVANCE = 'no') i, indiceConexoes(i)
      DO j = 1, 2*vizinhos
        WRITE(101, '(I6)', ADVANCE = 'no') listaConexoes(i,j)
      END DO
      DO j = 1, 2*vizinhos
        WRITE(101, '(I6)', ADVANCE = 'no') listaTipoConexao(i,j)
      END DO
      WRITE(101, *)
    END DO
    CLOSE(101)
  END SUBROUTINE SalvaLista
!
  SUBROUTINE SalvaMatrizAdjacencia!------------------------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: i, j, k
!    
    OPEN (101, FILE = 'MatrizAdjacencia_N1024', STATUS = 'unknown')
    DO i = 1, neuronios
      DO j = 1, neuronios
        IF (ANY(listaConexoes(i,:) == j)) THEN
          k = MINLOC(listaConexoes(i,:), MASK = (listaConexoes(i,:) == j), DIM = 1)
          WRITE(101, '(I3)', ADVANCE = 'no') listaTipoConexao(i,k)
        ELSE
          WRITE(101, '(I3)', ADVANCE = 'no') 0
        END IF
      END DO
      WRITE(101, *)
    END DO
    CLOSE(101)
  END SUBROUTINE SalvaMatrizAdjacencia
END PROGRAM GeraRede_AlgumNome
