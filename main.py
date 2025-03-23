import Funciones
while True:
        Funciones.mostrar_menu()
        opcion = input("Selecciona una opción (1-4): ")

        if opcion == '1':
            Funciones.factorizacion_lu()
        elif opcion == '2':
            Funciones.jacobi()
        elif opcion == '3':
            A, B=Funciones.Gauss_Jordan()
            print("La matriz transformada:")
            print(B)
            print("La solución es:")
            print(A)
        elif opcion == '4':
            print("Saliendo del programa...")
            break
        else:
            print("Opción no válida. Por favor, selecciona una opción del 1 al 4.")