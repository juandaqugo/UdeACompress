Archivo =  SecuencialFinal.c
	Esta fue la primera versi�n secuencial que se propuso, con todos los tipos de datos inttypes y con el Group Builder funcionando con una estructura de datos llamada ID utilizada para identificar cada grupo que se iba generando.
Ejecuci�n = gcc -o seq SecuencialFinal.c 
			./seq size(tama�o_cadena)


Archivo = SecuencialVersion2
	En esta versi�n ya no aparece la estructura de datos ID y el arreglo de frontera tiene un funcionamiento diferente, ya se rellena con 1's y 0's, cada que se marca un 1 es un limite superior. Nos quedamos con este m�todo ya que en la primera versi�n cuando los archivos son muy grandes el GroupBuilder toma mucho tiempo.
Ejecuci�n = gcc -o seq SecuencialVersion2.c
		  ./seq size(tama�o_cadena)

Archivo = intragrupo.c

	En esta versi�n aparecen los hilos que se despliegan dentro de cada grupo al momento de sacar los caracteres v�lidos y de la ejecuci�n del Radixsort.

Ejecuci�n = gcc -o intra intragupo.c -fopenmp
	 		./intra size nthreads

Archivo = Embebido.c

	En esta versi�n cada hilo se encarga de un subgrupo de ordenamiento y a partir de un valor del tama�o de grupo la ejecuci�n de ese grupo se realiza secuencial y no se divide en m�s.

Ejecuci�n = gcc -o emb Embebido.c -fopenmp
	 		./emb size nthreads

Archivo = intergrupo.c

	En esta versi�n cada hilo se encarga de un subgrupo de ordenamiento, la ejecuci�n del RadixSort y del SelectValidChar es secuencial
Ejecuci�n = gcc -o inter intergrupo.c -fopenmp
	 		./inter size nthreads

Archivo = intra_seq(threshold).c

	En esta versi�n aparecen los hilos que se despliegan dentro de cada grupo al momento de sacar los caracteres v�lidos y de la ejecuci�n del Radixsort. A partir del valor de SEQ_THRESHOLD en un subgrupo la ejecuci�n de ese grupo es secuencial.

Ejecuci�n = gcc -o intra intra_seq(threshold).c -fopenmp
	 		./intra size nthreads


Archivo = 	Sort.py
		Verifica correctitud de todos los archivos que terminen de la forma Out_SACSeq512_size.txt, este es el archivo que se genera al ejecutar el programa. Es necesario descomentar la l�nea #define archivos del archivo .c que se est� ejecutando.
	Ejecuci�n= python Sort.py
	