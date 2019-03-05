## Análisis de máxima verosimilitud con RAxML








**i. Conexión al servidor**



`ssh { nombre de usuario}  @132.248.13.5 -p 2022`

**Copia el archivo de trabajo a un directorio de trabajo**



```
```

**A. Un script para realizar unos análisis para un solo gen**


```
ls
```

En la menú de ayuda la primera línea indica la información obligatoria para un análisis: 1) `-s` seguido por el nombre del archivo de secuencias, 2) `-n` seguido por el nombre para el archivo para salida y 3) `-m` seguido para el modelo de sustitución. Sin embargo, `RAxML` también devuelva un error al intentar realizar una búsqueda que inicia a partir de un árbol de parsimonia sin proporcionar un número semilla. Para indicar el número, necesitas un cuarto parámetro `-p`. Para realizar un análisis muy simple puedes arreglar los parámetros así:

`raxmlHPC -s 0-615.phy -p 12345 -m GTRGAMMA -n 0-615.phy.tre`

Este análisis realizará una búsqueda heurística en el alineamiento `0-615.phy` y genera cinco archivos de salida,entre ellos el mejor árbol. El modelo `GTRGAMMA` especifica el modelo “general-time reversible” con el parámetro gamma para heterogeneidad de tasas. Revisa las condiciones de la búsqueda y posteriormente el mejor árbol en formato `phylip`.


`cat  RAxML_info.0-615.phy.tre`



Debes poder ver el árbol en formato Newick en el archivo `...bestTree...`. Puedes seleccionar y copiar el archivo Newick con el ratón y pegar el texto en una ventanita afuera del servidor (p. ej., en tu Mac, PC, o máquina de Linux) en [Dendroscope 3](http://dendroscope.org/) (Huson & Scornavacca 2015). Este programa se puede bajar de la red e instalar en tu computadora.

`raxmlHPC -s 0-615.phy -p 12345 -m GTRGAMMA -o devoDSG721 -n 0-615.phy.OG.tre`
```

También puedes usar el método de bootstrap (“remuestreo”) no-paramétrico **1)** como estrategia heurística en la búsqueda del mejor árbol y **2)** estimar el apoyo de las ramas:

`

```




```

**Para analizar particiones con modelos diferentes, pero con las mismas longitudes de ramas para todas las particiones:**

```



nohup raxmlHPC -f a -x 12345 -p12345 -# 300 -m GTRGAMMA -q Plastid_partitions_Gblocks.txt -s aaAustrales_GBlocks_mask_11jan2018.phy -n Plastid_Gblocks_part.bs300 &
```
`RAxML` tal como está instalada en nuestro servidor no permite aplicar longitudes de ramas diferentes para más de 128 particiones. Pero el comando para hacerlo sería:

```
nohup raxmlHPC -M -m GTRGAMMA Australes_n339_partitions.txt testpartitions &
```

 
```
RAxML_bipartitions.part.bs300 

RAxML_bootstrap.part.bs300
```

```
```
___

**B. Un script para analizar una tanda de alineamientos en una carpeta**



```
perl ../applyRAxML2AllFilesInDirectory.pl ./ "/usr/local/bin/raxmlHPC-PTHREADS-AVX -T 12 -m GTRGAMMA -p 685387 -o devoDSG721"

```




**Concatenar todos los árboles en un sólo archivo:**

```
```

**Otra posibilidad**




**Referencias**








**Apéndice**


Instalación de `RAxML` en Mac y Linux (bajo construcción). Para la versión reciente de `RAxML`, consulta la página del [autor]

La manera de instalar varía entre tipos de computadoras. Para Mac (y PC?), puedes dar clic en “Clone or Download” y posteriormente “Download ZIP”. Una vez que termina de bajar, puedes extraer el archivo. Sugiero que colocas la carpeta extraida en un archivo “Phylogenetics” que colocas cerca del directorio root de la computadora. Si estás trabajando en Linux sin acceso a un interfaz gráfica se puede bajar el archivo con todo el código necesario en terminal:

Hay que desempacar el archivo:





Esto debe compilar la versión secuencial `raxmlHPC`, que tal vez conviene si tienes una computadora personal. Para un servidor, puede ser más útil instalar la versión PPThreads. El servidor que usamos en este ejercicio cuenta con `raxmlHPC-THREADS-AVX`.


Probablemente será necesario agregar el programa a su `PATH`. Si eres el administrador en una computadora LINUX, simplemente coloca el programa en un lugar que ya está en el path. Por ejemplo:



Para Mac, si el archivo está ubicado el la carpeta



Agregar la siguiente linea al archivo y posteriormente salir (CTRL-X):


Una vez que saliste del archivo puedes teclear el siguiente para activar los cambios



Para confirmar que RAxML quedó en el `PATH`:

`