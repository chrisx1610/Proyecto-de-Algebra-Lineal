
import numpy
from numpy import array, zeros

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

def ingresar_matriz_y_vectorGJ():
    print("Ingresa la matriz 3x3 (fila por fila, separando los números con espacios):")
    matriz = []
    for i in range(3):
        fila = list(map(float, input(f"Fila {i + 1}: ").split()))
        if len(fila) != 3:
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
        for i in range(n): #empieza a iterarse por numero de columna
            # Pivoteo parcial: encontrar la fila con el máximo elemento en la columna i
            max_row = i #aquí empieza desde la primera columna
            for k in range(i + 1, n): #aquí se va a ver cada fila de cada columna
                if abs(A[k][i]) > abs(A[max_row][i]): #aquí se ve examina cual es el mayor elemento de todas las filas
                    max_row = k

            if max_row != i: #si el max_row no es igual al numero de columna se intercambian las filas en A, P y L
                A[i], A[max_row] = A[max_row], A[i] #se definen las dimensiones de la matriz A y P
                P[i], P[max_row] = P[max_row], P[i]
                if i > 0: 
                    L[i][:i], L[max_row][:i] = L[max_row][:i], L[i][:i] #Realiza un intercambio de elementos entre dos filas de la matriz L pero solo para las primeras i columnas< va desde columna 0 hasta columna i-1
            
            
            if A[i][i] == 0: #Si no se encuentra el pivote
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

def Gauss_Jordan():
    matriz,vector=ingresar_matriz_y_vector()
    if matriz is None:
        return
    a=array(matriz, float)
    b=array(vector, float)
    n=len(b)
    # Main loop
    for k in range(n):
        # Partial Pivoting
        if numpy.fabs(a[k, k]) < 1.0e-12:
            for i in range(k + 1, n):
                if numpy.fabs(a[i, k]) > numpy.fabs(a[k, k]):
                    for j in range(k, n):
                        a[k, j], a[i, j] = a[i, j], a[k, j]
                    b[k], b[i] = b[i], b[k]
                    break

        # Division of the pivot row
        pivot = a[k, k]
        for j in range(k, n):
            a[k, j] /= pivot
        b[k] /= pivot

        # Elimination loop
        for i in range(n):
            if i == k or a[i, k] == 0: continue
            factor = a[i, k]
            for j in range(k, n):
                a[i, j] -= factor * a[k, j]
            b[i] -= factor * b[k]

    return b, a