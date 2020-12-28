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

- [ ] Ver si no se puede shiftear bytes y luego desempaquetar para evitar usar las mascaras raras.
- [ ] Macro para reemplazar for loops en las versiones "withlogic".
- [x] Dividir version SSE de NW en subrutinas para facilitar el pasaje a ASM.
- [ ] Dividir version SSE de SW en subrutinas para facilitar el pasaje a ASM.
- [ ] Hacer versiones AVX-256.
- [ ] Pensar en usar instancia de AWS para usar AVX-512.
