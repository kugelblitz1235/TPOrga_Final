# Implementación de algoritmos de alineamiento de secuencias de ADN utilizando instrucciones SIMD

Como opción equivalente al final de la materia Organización del Computador II realizamos la implemetación de una versión paralelizada empleando instrucciones SIMD de los algoritmos de alineamiento de secuencias biológicas Smith-Waterman y Needleman-Wunsch. Dichos algoritmos son fundamentales en el estudio de secuencias de ADN durante el proceso de secuenciación que permite conocer el el genoma de cualquier especie. Adicionalmente, se evaluó su desempeño realizando distintos experimentos que se detallan en el informe. El trabajo fue realizado por Bruno Gomez, Tomás Tropea y Octavio Gianatiempo.

## Instrucciones de compilación

Para compilar ejecutar:

```
cd src
make
make test
make random_test
````
Esto creará los ejecutables `cli`, `test` y `random_test` en el directorio raíz del repositorio.

- `cli`: Interfaz de línea de comandos para ejecutar los algoritmos.
- `test`: Ejecuta los tests estáticos de las implementaciones lineales.
- `random_test`: Ejecuta test aleatorios del resto de las implementaciones y compara los resultados con los obtenidos en las implementaciones lineales.

Adicionalmente, al compilar `test` y `random_test` se ejecutarán una vez dichos tests.
Si se desea, se pueden correr los tests aleatorios hasta encontrar un error (esperemos que no suceda, lo dejamos un buen rato y no explotó 🤞). Para ello, ejecutar el script `random_test_until_fail.sh`.

## Modo de uso

La interfaz de línea de comandos permite alinear dos secuencias tomadas de dos archivos en formato [FASTA](https://es.wikipedia.org/wiki/Formato_FASTA), y tiene las siguientes opciones:

```
Especificar algoritmo con -a
Especificar dos archivos FASTA con las secuencias a alinear con -s y -t
Especificar puntaje de match con -p
Especificar puntaje de mismatch con -q
Especificar puntaje de gap con -r
Los algoritmos existentes son: 
- NW_ASM_AVX	Needlman-Wunsch implementado en ASM usando instrucciones AVX.
- NW_ASM_AVX512	Needlman-Wunsch implementado en ASM usando instrucciones AVX512.
- NW_ASM_LIN	Needlman-Wunsch implementado en ASM de forma lineal.
- NW_ASM_SSE	Needlman-Wunsch implementado en ASM usando instrucciones SSE.
- NW_C_AVX	Needlman-Wunsch implementado en C usando instrucciones AVX.
- NW_C_AVX512	Needlman-Wunsch implementado en C usando instrucciones AVX512.
- NW_C_LIN	Needlman-Wunsch implementado en C de forma lineal.
- NW_C_SSE	Needlman-Wunsch implementado en C usando instrucciones SSE.
- SW_ASM_AVX	Smith-Waterman implementado en ASM usando instrucciones AVX.
- SW_ASM_AVX512	Smith-Waterman implementado en ASM usando instrucciones AVX512.
- SW_ASM_LIN	Smith-Waterman implementado en ASM de forma lineal.
- SW_ASM_SSE	Smith-Waterman implementado en ASM usando instrucciones SSE.
- SW_C_AVX	Smith-Waterman implementado en C usando instrucciones AVX.
- SW_C_AVX512	Smith-Waterman implementado en C usando instrucciones AVX512.
- SW_C_LIN	Smith-Waterman implementado en C de forma lineal.
- SW_C_SSE	Smith-Waterman implementado en C usando instrucciones SSE.
```

A continación, un ejemplo de uso alineando dos genomas virales disponibles en el repositorio:
```
❯ ./cli -a NW_ASM_SSE -s sequences/genomes/Marburg_RefSeq.fasta -t sequences/genomes/Marburg_JN408064.fasta
456111 365631760 11822 19115 19114
```
Los resultados que se imprimen en la salida estándar corresponden a el tiempo de ejecución, la cantidad de celdas en la matriz de programación dinámica, el puntaje del alineamiento, y la longitud de cada una de las secuencias alineadas.

Para correr los algoritmos que emplean instrucciones SIMD necesario tener un procesador compatible con el set de instrucciones necesario.
Por ejemplo, este es un intento de uso del algoritmo que emplea instrucciones AVX-512 en un procesador que no tiene dichas instrucciones:
```
❯ ./cli -a NW_ASM_AVX512 -s sequences/genomes/Marburg_RefSeq.fasta -t sequences/genomes/Marburg_JN408064.fasta
[1]    6487 illegal hardware instruction (core dumped)  ./cli -a NW_ASM_AVX512 -s sequences/genomes/Marburg_RefSeq.fasta -t 
```

Sin embargo, el [_Intel Software Development Emulator_](https://software.intel.com/content/www/us/en/develop/articles/intel-software-development-emulator.html) permite emular la ejecución de instrucciones no compatibles con un procesador para que los desarrolladores puedan ganar familiaridad con sets de instrucciones que todavía están ampliamente disponibles. Dicho emulador es el que permite correr los tests aleatorios sobre todas las implementaciones y puede ser utilizado para realizar un alineamiento particular de la siguiente forma:
```
❯ sde-kit/sde -skx -- ./cli -a NW_ASM_AVX512 -s sequences/genomes/Marburg_RefSeq.fasta -t sequences/genomes/Marburg_JN408064.fasta
8821971 366396992 11822 19115 19114
```
