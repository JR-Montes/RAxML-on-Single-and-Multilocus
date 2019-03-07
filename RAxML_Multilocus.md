## Análisis de máxima verosimilitud con RAxML

By David Gernandt and José Rubén Montes
Ver. Marzo 2019

**Introducción**

Máxima verosimilitud es un método poderoso y popular de análisis filogenética con datos moleculares. Uno de sus primeras aplicaciones en la filogenética fue por Cavalli-Sforza & Edwards (1967) quienes lo propusieron para la inferencia de un árbol para expresar las relaciones entre poblaciones humanas con base en frecuencia génicas de grupos sanguíneas. Fue desarrollado para secuencias de DNA por Felsenstein (1981). En el contexto de inferencia filogenética, la estimación de verosimilitud representa la probabilidad de los datos (los caracteres) dado un hipótesis (el árbol).


Uno de los limitantes prácticos más importante para la inferencia filogenética de máxima verosimilitud es el gran esfuerzo computacional que requiere para calcular la función de verosimilitud comparado con distancia, parsimonia, o de un menor grado, inferencia Bayesiana. Como resultado, por mucho tiempo fue poco factible realizar análisis filogenéticas de máxima verosimilitud con conjuntos de datos con mucho taxa (p.ej., más de 100) o muchos caracteres. Stamatakis et al. (2005) reportaron avances importantes en este tema en su programa RAxML-III (Randomized A(x)ccelerated Maximum Likelihood). El programa se ha desarrollado activamente desde su forma inicial, permitiendo el uso de más modelos de evolución para DNA, aminoácidos y morfología, más opciones para medir apoyo mediante el remuestreo (“bootstrap) y más maneras de probar hipótesis. También se ha mejorado su habilidad de aprovechar el procesamiento en paralelo con núcleos múltiples que ofrecen las computadoras recientes (Stamatakis 2014).


En este tutorial vamos a realizar análisis de máxima verosimilitud para secuencias de genes de bajo número de copia provenientes de  Pinus  subsección  Australes , un grupo natural de unas 30 especies de pinos distribuidos en el Nuevo Mundo. El servidor que usaremos maneja procesos múltiples (cuenta con dos procesadores, cada uno con 12 núcleos, lo cual permite manejar 24 subprocesos). La versión raxmlHPC-THREADS-AVX fue compilado recientemente a partir del source code y el ejecutable fue transferido a su `PATH` en `/usr/local/bin/`.


**i. Conexión al servidor**
Usaremos un  secure shell p  ara entrar al servidor de BioLinux. Los usuarios de Windows pueden emplear PuTTY, mientros los usuarios de Mac y Linux pueden usar Terminal. Así vas a poder interactuar con el servidor desde la línea de comandos. Es necesario contar con tu propio nombre de usuario y contraseña.


**ii. Usuarios de PuTTY**  hay que indicar el nombre del servidor “pinaceae.ib.unam.mx” y oprimir “Open”. La primera vez que se realiza la conección desde una computadora recibes una advertencia al respecto. Elige “Yes” para confirmar la conexión. Debe salir una ventana que pide “login as:”. Hay que proporcionar tu nombre de usuario y después tu contraseña.
Usuarios de Terminal  deben teclear:

`ssh { nombre de usuario}  @132.248.13.5 -p 2022`

**Copia el archivo de trabajo a un directorio de trabajo**

`RAxML` requiere archivos en formato `PHYLIP`. Hay un archivo comprimido disponible en el servidor que contiene 20 alineamientos en este formato. Estos fueron generados mediante la opción “batch export” en Geneious. También realicé una reducción en el muestreo con la ayuda del script [Beforephylo.pl](https://github.com/qiyunzhu/BeforePhylo) y posteriormente los empaqué en un tarball (`tar -zcvf Australes_t20_n20.tar.gz *.phy`). Dentro de tu directorio home, hay que crear una carpeta `raxml/`, copiar el archivo de ejemplo y descomprimirlo:


```
mkdir  raxml
cd  raxml
ls  /mnt/datos/
cp  /mnt/datos/Australes_t19_n20.tar.gz ./ tar -xvzf  Australes_t19_n20.tar.gz
ls -1 | wc -l
```

**A. Un script para realizar unos análisis para un solo gen**

Para empezar, vamos a inspeccionar uno de los alineamientos en formato `PHYLIP` (`0-615.phy`) y las opciones disponibles en `RAxML`:

```
ls
mkdir  1gen
cp /mnt/datos/0-615.phy ./ cp  0-615.phy 1gen/
cd  1gen/
cat  0-615.phy
raxmlHPC -h
```

En la menú de ayuda la primera línea indica la información obligatoria para un análisis: 1) `-s` seguido por el nombre del archivo de secuencias, 2) `-n` seguido por el nombre para el archivo para salida y 3) `-m` seguido para el modelo de sustitución. Sin embargo, `RAxML` también devuelva un error al intentar realizar una búsqueda que inicia a partir de un árbol de parsimonia sin proporcionar un número semilla. Para indicar el número, necesitas un cuarto parámetro `-p`. Para realizar un análisis muy simple puedes arreglar los parámetros así:

`raxmlHPC -s 0-615.phy -p 12345 -m GTRGAMMA -n 0-615.phy.tre`

Este análisis realizará una búsqueda heurística en el alineamiento `0-615.phy` y genera cinco archivos de salida,entre ellos el mejor árbol. El modelo `GTRGAMMA` especifica el modelo “general-time reversible” con el parámetro gamma para heterogeneidad de tasas. Revisa las condiciones de la búsqueda y posteriormente el mejor árbol en formato `phylip`.


`cat  RAxML_info.0-615.phy.tre`

`cat  RAxML_bestTree.0-615.phy.tre`


Debes poder ver el árbol en formato Newick en el archivo `...bestTree...`. Puedes seleccionar y copiar el archivo Newick con el ratón y pegar el texto en una ventanita afuera del servidor (p. ej., en tu Mac, PC, o máquina de Linux) en [Dendroscope 3](http://dendroscope.org/) (Huson & Scornavacca 2015). Este programa se puede bajar de la red e instalar en tu computadora.
Por defecto, `RAxML` produce árboles no enraizados. Si quieres realizar un análisis y generar un árbol enraizado, puedes emplear la opción `-o`:

`raxmlHPC -s 0-615.phy -p 12345 -m GTRGAMMA -o devoDSG721 -n 0-615.phy.OG.tre`

```
ls
cat  RAxML_bestTree.0-615.phy.OG.tre
```

También puedes usar el método de bootstrap (“remuestreo”) no-paramétrico **1)** como estrategia heurística en la búsqueda del mejor árbol y **2)** estimar el apoyo de las ramas:

`
raxmlHPC -s 0-615.phy -f a -p 12345 -x 12345 -# 100 -m GTRGAMMA -o devoDSG721 -n 0-615.phy.OG.bo.tre`

```
ls

cat  RAxML_bipartitions.0-615.phy.OG.bo.tre


nohup raxmlHPC -x 12345 -p 12345 -# 100 -m GTRGAMMA -s Australes_RAxML_t74n339.phy -n TEST &

nohup raxmlHPC -f a -x 12345 -p12345 -# 300 -m GTRGAMMA -s Australes_RAxML_t74n339.phy -n testbs300 &
```

**Para analizar particiones con modelos diferentes, pero con las mismas longitudes de ramas para todas las particiones:**

```
nohup raxmlHPC -m GTRGAMMA -p 12345 -q Australes_n339_partitions.txt -s Australes_n339.phy -n testpartitions &


nohup raxmlHPC -f a -x 12345 -p12345 -# 300 -m GTRGAMMA -q Australes_n339_partitions.txt -s Australes_n339.phy -n test.part.bs300 &

nohup raxmlHPC -f a -x 12345 -p12345 -# 300 -m GTRGAMMA -q Plastid_partitions_Gblocks.txt -s Australes_GBlocks_mask_11jan2018.phy -n Plastid_Gblocks_part.bs300 &
```
`RAxML` tal como está instalada en nuestro servidor no permite aplicar longitudes de ramas diferentes para más de 128 particiones. Pero el comando para hacerlo sería:

```
nohup raxmlHPC -M -m GTRGAMMA Australes_n339_partitions.txt testpartitions &
```

**Para identificar “rogue taxa”**
 
```
RAxML_bipartitions.part.bs300 

RAxML_bootstrap.part.bs300
```

```
mkdir detect_rogues
cp RAxML_bootstrap.part.bs300
cd detect_rogues
raxmlHPC -J MR_DROP -z RAxML_bootstrap.part.bs300 -m GTRCAT -n TEST
```
___

**B. Un script para analizar una tanda de alineamientos en una carpeta**


Para algunos de los análisis multilocus, hay que compilar árboles para cada gen individualmente. Esto puede ser muy tardado. Afortunadamente existe un script de RAxML (applyRAxML2AllFilesInDirectory.pl) para facilitar el trabajo. Si tienes el script en la carpeta “raxml” y tus alineamientos están en un directorio de raxml, puedes ir al directorio con los alineamientos y teclear:

```
perl ../applyRAxML2AllFilesInDirectory.pl ./ "/usr/local/bin/raxmlHPC-PTHREADS-AVX -T 12 -m GTRGAMMA -p 685387 -o devoDSG721"

perl ../applyRAxML2AllFilesInDirectory.pl ./ "/usr/local/bin/raxmlHPC-PTHREADS-AVX -T 12 -f a -m GTRGAMMA -p 685387 -o devoDSG721 -x 574684 -# autoMRE"
```


Tratar esto para Astral-III (because RAxML manual discourages designating outgroup):

```
perl ../applyRAxML2AllFilesInDirectory.pl ./ "/usr/local/bin/raxmlHPC-PTHREADS-AVX -T 12 -f a -m GTRGAMMA -p 685387 -x 574684 -# autoMRE"
```

**Concatenar todos los árboles en un sólo archivo:**

```
mkdir  besttrees
mv  RAxML_bestTree.* besttrees
cd  besttrees
ls
cat  RAxML_best* > alltrees.tre
nl  -ba -w1 -s ' ' alltrees.tre > alltrees2.tre
awk  '{ printf "Tree "; print }' alltrees2.tre > alltrees3.tre
nano  alltrees3.tre 
```

**Otra posibilidad**

`sed -i [adding tree name = ] & cat [for header and end]`



**Referencias**

Cavalli-Sforza, L.L., and A.W.F. Edwards. 1967. Phylogenetic analysis: models and estimation procedures. Evolution 21: 550–570.

Felsenstein, J. 1981. Evolutionary trees from DNA sequences: a maximum likelihood approach. Journal of Molecular Evolution 17: 368–376.

Huson, D.H., and C. Scornavacca. 2012. Dendroscope 3: an interactive tool for rooted phylogenetic trees and networks. Systematic Biology 61: 1061–1067.

Stamatakis, A. 2014. RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics 30: 1312–1313.

Stamatakis, A., T. Ludwig, and H. Meier. 2005. RAxML-III: a fast program for maximum likelihood-based inference of large phylogenetic trees. Bioinformatics 21: 456–463.

Than, C., D. Ruths, and L. Nakhleh. 2008. PhyloNet: a software package for analyzing and reconstructing reticulate evolutionary relationships. BMC Bioinformatics 9: 322.


**Apéndice**

Apéndice 1.

Instalación de `RAxML` en Mac y Linux (bajo construcción). Para la versión reciente de `RAxML`, consulta la página del [autor]
(http://sco.h-its.org/exelixis/web/software/raxml/). Actualmente está disponible en [GibHub](https://github.com/stamatak/standard-RAxML). Como siempre, es altamente recomendable leer el manual.

La manera de instalar varía entre tipos de computadoras. Para Mac (y PC?), puedes dar clic en “Clone or Download” y posteriormente “Download ZIP”. Una vez que termina de bajar, puedes extraer el archivo. Sugiero que colocas la carpeta extraida en un archivo “Phylogenetics” que colocas cerca del directorio root de la computadora. Si estás trabajando en Linux sin acceso a un interfaz gráfica se puede bajar el archivo con todo el código necesario en terminal:


`wget https://github.com/stamatak/standard-RAxML/archive/master.zip`

Hay que desempacar el archivo:

`unzip master.zip`

En Mac o Linux, hay que usar un terminal para cambiar al directorio y compilar el código. Por ejemplo:

`cd standard-RAxML-master make -f Makefile.gcc`


Esto debe compilar la versión secuencial `raxmlHPC`, que tal vez conviene si tienes una computadora personal. Para un servidor, puede ser más útil instalar la versión PPThreads. El servidor que usamos en este ejercicio cuenta con `raxmlHPC-THREADS-AVX`.

`make -f Makefile.AVX.PTHREADS.gcc`

Probablemente será necesario agregar el programa a su `PATH`. Si eres el administrador en una computadora LINUX, simplemente coloca el programa en un lugar que ya está en el path. Por ejemplo:


`mv raxmlHPC-PTHREADS-AVX /usr/local/bin/`

Para Mac, si el archivo está ubicado el la carpeta
`/users/nombrecompu/Phylogenetics/standard-RAXML-master/`, hay que modificar su archivo `.bash_profile`:

`nano ~/.bash_profile`


Agregar la siguiente linea al archivo y posteriormente salir (CTRL-X):

`PATH="$PATH:Users/ nombrecompu/Phylogenetics/standard-RAXML-master/"`

Una vez que saliste del archivo puedes teclear el siguiente para activar los cambios

`source ~/.bash_profile`


Para confirmar que RAxML quedó en el `PATH`:

`echo $PATH=
`
