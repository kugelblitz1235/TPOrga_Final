# TPOrga_Final

* Hacer parser de formato FASTA para leer secuencias de archivos
* Leer las secuencias de archivos para no tener que ingresar una secuencia extensa por consola
* Guardar la respuesta del alineamiento local/global en una struct que se pasa por referencia a los distintos algoritmos con el siguiente formato: 
~~~~
struct Alignment{
  unsigned int length;
  char* sequence_1;
  char* sequence_2;
};
~~~~

## TODO:

- [x] Ver si no se puede shiftear bytes y luego desempaquetar para evitar usar las mascaras raras.
- [x] Adaptar NW withLogicSSE para incluir las nuevas mascaras.
- [x] Dividir version SSE de NW en subrutinas para facilitar el pasaje a ASM.
- [x] Dividir version SSE de SW en subrutinas para facilitar el pasaje a ASM.
- [x] NW_ASM_SSE
- [x] Asegurarse de pushear todo lo que se podria romper antes de llamar a malloc u otras funciones de C.
- [x] Ver bug de score 0 en SW (Ver si quedo arreglado inicializando posicion maxima en 0).
- [x] SW_ASM_SSE
- [x] Hacer versiones AVX-256.
- [x] Hacer versiones AVX-512.
- [x] Documentar los pasos de NW y SW para LIN, SSE, AVX y AVX512 en C
- [x] Dividir el .c  de NW y SW para las distintas versiones porque quedo enorme
- [x] Documentar los pasos de NW y SW para LIN, SSE, AVX y AVX512 en ASM
- [x] Reducir push y pop alrededor de los malloc en NW_ASM_SSE
- [x] Emprolijar mascaras que no se usan y demas basura que quedo.
- [x] Pensar en remover casos base y utilizar diag1 y diag2
- [x] Decidir si agregar o no los cambios a NW y SW de C en AVX512 a los respectivos asm.
- [x] Agregar a SSE diag1 y diag2 para evitar levantar de memoria en cada vuelta del loop principal.
- [x] Opcional: dejar parametros como variables globales en SSE.
- [x] Opcional: agregar tipo de datos union para simular registro en memoria en AVX.
- [x] Cli o las funciones que manejan la selección de algoritmos deben chequear que se cumplan las precondiciones (i.e. secuencias de al menos 8 letras para los SSE).
- [x] Pensar en usar instancia de AWS para usar AVX-512.
- [ ] Chequear que el tamaño de los strings sea manejable en cli y evitar que se use en ese caso.
- [ ] Medir el tiempo en cli.
- [ ] Medir GCUPS en cli.
- [ ] Pensar si habría que cambiar la representación del score a int.