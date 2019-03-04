Archivo =  SecuencialFinal.c
	Esta fue la primera versión secuencial que se propuso, con todos los tipos de datos inttypes y con el Group Builder funcionando con una estructura de datos llamada ID utilizada para identificar cada grupo que se iba generando.
Ejecución = gcc -o seq SecuencialFinal.c 
			./seq size(tamaño_cadena)


Archivo = SecuencialVersion2
	En esta versión ya no aparece la estructura de datos ID y el arreglo de frontera tiene un funcionamiento diferente, ya se rellena con 1's y 0's, cada que se marca un 1 es un limite superior. Nos quedamos con este método ya que en la primera versión cuando los archivos son muy grandes el GroupBuilder toma mucho tiempo.
Ejecución = gcc -o seq SecuencialVersion2.c
		  ./seq size(tamaño_cadena)

Archivo = intragrupo.c

	En esta versión aparecen los hilos que se despliegan dentro de cada grupo al momento de sacar los caracteres válidos y de la ejecución del Radixsort.

Ejecución = gcc -o intra intragupo.c -fopenmp
	 		./intra size nthreads

Archivo = Embebido.c

	En esta versión cada hilo se encarga de un subgrupo de ordenamiento y a partir de un valor del tamaño de grupo la ejecución de ese grupo se realiza secuencial y no se divide en más.

Ejecución = gcc -o emb Embebido.c -fopenmp
	 		./emb size nthreads

Archivo = intergrupo.c

	En esta versión cada hilo se encarga de un subgrupo de ordenamiento, la ejecución del RadixSort y del SelectValidChar es secuencial
Ejecución = gcc -o inter intergrupo.c -fopenmp
	 		./inter size nthreads

Archivo = intra_seq(threshold).c

	En esta versión aparecen los hilos que se despliegan dentro de cada grupo al momento de sacar los caracteres válidos y de la ejecución del Radixsort. A partir del valor de SEQ_THRESHOLD en un subgrupo la ejecución de ese grupo es secuencial.

Ejecución = gcc -o intra intra_seq(threshold).c -fopenmp
	 		./intra size nthreads


Archivo = 	Sort.py
		Verifica correctitud de todos los archivos que terminen de la forma Out_SACSeq512_size.txt, este es el archivo que se genera al ejecutar el programa. Es necesario descomentar la línea #define archivos del archivo .c que se esté ejecutando.
	Ejecución= python Sort.py
	