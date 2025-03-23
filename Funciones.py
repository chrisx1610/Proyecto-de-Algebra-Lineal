def mostrar_menu():
    print("MENU:")
    print("1. FACTORIZACIÓN LU")
    print("2. MÉTODO DE JACOBI")
    print("3. SOLUCIÓN DIRECTA POR GAUSS-JORDAN")
    print("4. SALIR")

def ingresar_matriz_y_vector():
    print("Ingresa la matriz 4x4 (fila por fila, separando los números con espacios):")
    matriz = []
    for i in range(4):
        fila = list(map(float, input(f"Fila {i + 1}: ").split()))
        if len(fila) != 4:
            print("Error: Debes ingresar exactamente 4 números por fila.")
            return None, None
        matriz.append(fila)

    print("Ingresa el vector de términos independientes (4 números separados por espacios):")
    vector = list(map(float, input().split()))
    if len(vector) != 4:
        print("Error: Debes ingresar exactamente 4 números para el vector.")
        return None, None

    return matriz, vector

def factorizacion_lu():
    print("\nHas seleccionado Factorización LU.")
    matriz, vector = ingresar_matriz_y_vector()
    if matriz is not None and vector is not None:
        n = 4
        L = [[0.0] * n for _ in range(n)]
        U = [[0.0] * n for _ in range(n)]
        P = [[float(i == j) for j in range(n)] for i in range(n)]  # Matriz de permutación

        # Copiar la matriz original para no modificarla
        A = [fila.copy() for fila in matriz]

        # Factorización LU con pivoteo parcial
        for i in range(n):
            # Pivoteo parcial: encontrar la fila con el máximo elemento en la columna i
            max_row = i
            for k in range(i + 1, n):
                if abs(A[k][i]) > abs(A[max_row][i]):
                    max_row = k

            if max_row != i:
                A[i], A[max_row] = A[max_row], A[i]
                P[i], P[max_row] = P[max_row], P[i]
                if i > 0:
                    L[i][:i], L[max_row][:i] = L[max_row][:i], L[i][:i]

            
            if A[i][i] == 0:
                print("Error: La matriz es singular y no se puede factorizar.")
                return

            
            for j in range(i, n):
                U[i][j] = A[i][j] - sum(L[i][k] * U[k][j] for k in range(i))

            
            for j in range(i + 1, n):
                L[j][i] = (A[j][i] - sum(L[j][k] * U[k][i] for k in range(i))) / U[i][i]

            
            L[i][i] = 1.0

      
        print("\nMatriz P (permutación):")
        for fila in P:
            print(fila)

        print("\nMatriz L (triangular inferior):")
        for fila in L:
            print(fila)

        print("\nMatriz U (triangular superior):")
        for fila in U:
            print(fila)

        # Resuelve el sistema Ly = Pb
        Pb = [sum(P[i][j] * vector[j] for j in range(n)) for i in range(n)]
        y = [0.0] * n
        for i in range(n):
            y[i] = Pb[i] - sum(L[i][j] * y[j] for j in range(i))

        # Resuelve el sistema Ux = y
        x = [0.0] * n
        for i in range(n - 1, -1, -1):
            x[i] = (y[i] - sum(U[i][j] * x[j] for j in range(i + 1, n))) / U[i][i]

        print("\nSolución del sistema (x):")
        print(x)

def es_diagonalmente_dominante(matriz):
    n = len(matriz)
    for i in range(n):
        suma = sum(abs(matriz[i][j]) for j in range(n) if j != i)
        if abs(matriz[i][i]) <= suma:
            return False
    return True

def redefinir_sistema(matriz, vector):
    n = len(matriz)
    for i in range(n):
        suma = sum(abs(matriz[i][j]) for j in range(n) if j != i)
        if abs(matriz[i][i]) <= suma:
            # Ajustar el elemento diagonal para que sea mayor que la suma de los otros elementos
            matriz[i][i] = suma + 1  # Se suma 1 para asegurar que sea estrictamente mayor
    return matriz, vector

def jacobi():
    print("\nHas seleccionado el Método de Jacobi.")
    matriz, vector = ingresar_matriz_y_vector()
    if matriz is None or vector is None:
        return

    if not es_diagonalmente_dominante(matriz):
        print("El sistema no está bien condicionado. Redefiniendo el sistema...")
        matriz, vector = redefinir_sistema(matriz, vector)
        if not es_diagonalmente_dominante(matriz):
            print("No se pudo redefinir el sistema para que sea diagonalmente dominante.")
            return

    n = len(matriz)
    iteraciones = int(input("Ingresa la cantidad de iteraciones: "))
    x = [0.0] * n

    for _ in range(iteraciones):
        x_nuevo = [0.0] * n
        for i in range(n):
            suma = sum(matriz[i][j] * x[j] for j in range(n) if j != i)
            x_nuevo[i] = (vector[i] - suma) / matriz[i][i]
        x = x_nuevo

    print("\nSolución aproximada del sistema (x):")
    print(x)