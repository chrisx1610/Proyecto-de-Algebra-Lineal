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
        P = [[float(i == j) for j in range(n)] for i in range(n)]  

        A = [fila.copy() for fila in matriz]

        # Factorización LU con pivoteo parcial
        for i in range(n):
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
    for col in range(n):
        max_row = col
        for k in range(col + 1, n):
            if abs(matriz[k][col]) > abs(matriz[max_row][col]):
                max_row = k
        
        if max_row != col:
            matriz[col], matriz[max_row] = matriz[max_row], matriz[col]
            vector[col], vector[max_row] = vector[max_row], vector[col]
    return matriz, vector

def jacobi():
    print("\nHas seleccionado el Método de Jacobi.")
    matriz, vector = ingresar_matriz_y_vector()
    if matriz is None or vector is None:
        return

    n = len(matriz)
    matriz, vector = redefinir_sistema(matriz, vector)  # Reordenar filas primero

    if not es_diagonalmente_dominante(matriz):
        print("\nADVERTENCIA: El sistema no es diagonalmente dominante incluso después de reordenar.")
        print("El método de Jacobi puede no converger o dar resultados incorrectos.")
        continuar = input("¿Deseas continuar de todos modos? (s/n): ").lower()
        if continuar != 's':
            return

    iteraciones = int(input("Ingresa la cantidad de iteraciones: "))
    x = [0.0] * n  # Solución inicial :)

    print("\nIteración\tValores de x")
    print(f"0\t\t{x}")

    for it in range(1, iteraciones + 1):
        x_nuevo = [0.0] * n
        for i in range(n):
            suma = sum(matriz[i][j] * x[j] for j in range(n) if j != i)
            x_nuevo[i] = (vector[i] - suma) / matriz[i][i]
        x = x_nuevo
        print(f"{it}\t\t{x}")

    print("\nSolución aproximada después de", iteraciones, "iteraciones:")
    print([round(val, 6) for val in x])  # Redondeo 

def gauss_jordan():
    print("\nHas seleccionado Solución Directa por Gauss-Jordan.")
    matriz, vector = ingresar_matriz_y_vector()
    if matriz is None or vector is None:
        return
    
    n = 4

    aumentada = [fila.copy() for fila in matriz]
    for i in range(n):
        aumentada[i].append(vector[i])
    
    for col in range(n):
        max_row = col
        for k in range(col + 1, n):
            if abs(aumentada[k][col]) > abs(aumentada[max_row][col]):
                max_row = k
        
        if max_row != col:
            aumentada[col], aumentada[max_row] = aumentada[max_row], aumentada[col]
        
        if abs(aumentada[col][col]) < 1e-10:  
            print("Error: La matriz es singular y no tiene solución única.")
            return
        
        pivot = aumentada[col][col]
        for j in range(col, n + 1):
            aumentada[col][j] /= pivot
        
        for i in range(n):
            if i != col:
                factor = aumentada[i][col]
                for j in range(col, n + 1):
                    aumentada[i][j] -= factor * aumentada[col][j]
    
    solucion = [round(aumentada[i][n], 10) for i in range(n)] 
    
    print("\nMatriz aumentada final en forma reducida:")
    for fila in aumentada:
        print([f"{x:.6f}".rstrip('0').rstrip('.') if '.' in f"{x:.6f}" else f"{x:.6f}" for x in fila])
    
    print("\nSolución del sistema (x):")
    print([x if abs(x) > 1e-10 else 0.0 for x in solucion])  
