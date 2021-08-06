# Dinamica-Molecular-Fortran90

Código modular escrito en Fortran90, para realizar mediciones y simulaciones 3D de Dinámica Molecular. 
----------------------------------------------------------------------------------------------------



Contiene dependencias para editar y compilar con Code::Blocks. El programa principal es main.f90.

Contiene una variedad de subrutinas que cumplen diferentes funciones: 

-ConfigIniAle.f90: Genera una configuración inicial aleatoria, es decir, con posiciones aleatorias. Ideal para sistemas con baja concentración.

-configinifccalt.f90: Genera una configuración inicial cristalina FCC. Ideal para sistemas con muy alta concentración que naturalmente presentan estructura FCC.

-ConfigIniReg.f90: Genera una configuración regular cúbica. Ideal para sistemas con alta concentración. Requiere DistLineaML.f90.

-correctefrain.f90: Corrige la configuración predicha de acuerdo al método de Gear predictor-corrector. Utiliza unidades de posición reducida.

-dev.f90: Lee la variable de una columna específica de un archivo, saca su promedio y desviación.

-DistLineaML.f90: Coloca partículas de manera equidistante dentro de una línea con márgenes. Dependencia de ConfigIniReg.f90.

-fuerzas.f90: Calcula la fuerza ejercida sobre cada partícula. Utiliza el modelo de potencial de Gupta que es ideal para metales de transición que muestran estructura FCC.

-gdr.f90: Calcula la función de distribución radial promedio de las configuraciones guardadas.

-instgdr.f90: Calcula la función de distribución radial de una sola configuración, incluye el cálculo de entropía de exceso S2.

-predictefrain.f90: Predice la siguiente configuración, de acuerdo al método de Gear predictor-corrector. Utiliza unidades de posición reducida.

-veliniale.f90: Genera velocidades iniciales aleatorias y cuida que se respete el teorema de equipartición de la energía, y que el momento lineal total del sistema sea 0.

-wvdt.f90: Calcula el desplazamiento cuadrático medio (MSD), coeficiente de difusión (relación de Einstein) y la función de autocorrelación de velocidades promedio utilizando las configuraciones guardadas.

-----------------------------------

Carpeta CNA y desviaciones contiene archivos que hacen cálculos con lectura de archivos, no dependen de main.f90:

-cnaport.f90: Calcula la microestructura de una configuración.

-desviación.f90: Calcula la desviación media de una columna de datos.


-----------------------------------

Se incluye ejemplo de archivos de resultados .txt generados por el programa.
