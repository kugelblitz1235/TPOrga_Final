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
- [ ] NW_ASM_SSE
  - [ ] Asegurarse de pushear todo lo que se podria romper antes de llamar a malloc u otras funciones de C.
- [ ] SW_ASM_SSE
- [ ] Hacer versiones AVX-256.
- [ ] Pensar en usar instancia de AWS para usar AVX-512.
- [ ] Chequear que el tamaño de los strings sea manejable en cli y evitar que se use en ese caso.
- [ ] Ver bug de score 0 en SW (Ver si quedo arreglado inicializando posicion maxima en 0).
- [ ] Emprolijar mascaras que no se usan y demas basura que quedo.
- [ ] Cli o las funciones que manejan la selección de algoritmos deben chequear que se cumplan las precondiciones (i.e. secuencias de al menos 8 letras para los SSE).
- [ ] Reducir push y pop alrededor de los malloc en NW_ASM_SSE
