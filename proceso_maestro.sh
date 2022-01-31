#Proyecto: Interacción Parásito-Microbiota_Hospedero en modelo murino en infección por Trypanosoma cruzi (Cepa Tulahuen)
#Datos obtenidos por shotgun metagenomics (4Gb por muestra)
#Ratones BALBc y BL6
#Para cada modelo hay tomas a los 0-3-7-10-13-16 días 
#Los días 0 - 10 - 16 tienen un NI
#NI = No Infectado
#DPI: Días Post Infección

#Acceso al clúster. Asegurarse de tener la llave .pem para acceder al clúster

ssh -i sergio.castaneda.pem -p 53841 sergio.castaneda@labcompu-hpc.urosario.edu.co

#Crear o usar una sesión en tmux para trabajar

#Crear:
tmux new -s nombresesion

#acceder
tmux a -t nombresesion

# Ruta de Proyecto: /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com

#Cargar ambiente de trabajo en el clúster
module load conda
source activate /datacnmat01/biologia/gimur/Anaconda/Anaconda3/envs/gimur-tools


#Raw reads

#pwd de los reads: /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com
#Allí se crean varias carpetas:
#En análisis manual quedan cada una de las subcarpetas de análisis
#En herramientas quedan las herramientas que han sido instaladas fuera de ambientes conda


#Preprocesamiento

# Ruta de Proyecto: /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com

#Cargar ambiente de trabajo en el clúster
module load conda
source activate /datacnmat01/biologia/gimur/Anaconda/Anaconda3/envs/gimur-tools

#Se realiza un control de calidad confirmatorio por medio de fastqc.
#Se realizó la installación de fastqc con el siguiente comando:

conda install -c bioconda/label/broken fastqc

#ejecutar análisis. Lecturas se encuentran en la ruta: /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/raw_data_merge/

fastqc *.gz

# Uso de trimmomatic
# La herramienta se instaló con el comando:

conda install -c bioconda/label/broken trimmomatic

# Se efectúa el trimming de las lecturas para confirmar eliminación de adaptadores,
# para mantener secuencias de al menos 150 pb
# para realizar filtrado por calidad (Phred Score)
# Repetir el siguiente comando para cada muestra

trimmomatic PE -threads 16 /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/raw_data_merge/BL6_16_NI_1.fq.gz /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/raw_data_merge/BL6_16_NI_2.fq.gz /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/trimmomatic_out/BL6_16_NI_1_paired_trim.fq.gz /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/trimmomatic_out/BL6_16_NI_1_unpaired_trim.fq.gz /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/trimmomatic_out/BL6_16_NI_2_paired_trim.fq.gz /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/trimmomatic_out/BL6_16_NI_2_unpaired_trim.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads MINLEN:150 AVGQUAL:20 TRAILING:20

#Mapeo con genoma del hospedero con Bowtie2. Se usó script previamente escrito que fue ajustado para el presente análisis.
#Se puede hacer manual. Pare este caso se usa el script adaptado bowtie2_mice:

----------------------------------------------------------------------------------------------------------------------------------------

#!/bin/bash
# Bowtie2 wrapper script to filter out contaminations from metagenomic samples.
# Adapted by: Sergio Castañeda
# Last updated on: 2021-09-11

set -e

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
cat <<HELP_USAGE
Remove contaminant sequences with Bowtie2.
Usage: metabiome bowtie2 [options] -i <in_dir> -o <out_dir> -opts bowtie2_options

Required:
  -i in_dir         Input directory containing FASTQ files.
  -o out_dir        Directory in which results will be saved. This directory will
                    be created if it doesn't exist.

Options:
  -ho host          Host reference genome in FASTA format.
  -ph PhiX          PhiX-174 phage reference genome in FASTA format.
  -m  mus           Mus musculus reference genome in FASTA format.
  -t  NUM           Number of threads to use. (default=1)
  -opts OPTIONS     Bowtie2's options.
  -h, --help        Show this help.
HELP_USAGE
}

# Exit if command is called with no arguments
validate_arguments "$#"

# Parse command line arguments
while (("$#")); do
    case "$1" in
        -h|--help ) usage; exit 0 ;;
        -i )        input_dir=$(readlink -f "$2"); shift 2 ;;
        -o )        out_dir=$(readlink -m "$2"); shift 2 ;;
        -ho )       host=$(readlink -f "$2"); shift 2 ;;
        -t )        threads="$2"; shift 2 ;;
        -ph )       PhiX=$(readlink -f "$2"); shift 2 ;;
        -m  )       Mus=$(readlink -f "$2"); shift 2 ;;
        -opts )     shift; bowtie2_opts="$@"; break ;;
        * )         echo "Option '$1' not recognized"; exit 1 ;;
    esac
done

# Activate Conda environment
#activate_env metabiome-preprocessing

# Download Mice and PhiX reference genomes
# Highlight: Next code downloads PhiX and Mus musculus genome and checks if the
# downloads were successfull
cd "$out_dir"
if [[ ! -e "$Mus" ]]; then
    URL="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/latest_assembly_versions/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz"
    for i in {1..10}; do
        echo "Downloading Mus musculus reference genome"
        wget "$URL" -P "$out_dir"
        if [[ -s GCF_000001635.27_GRCm39_genomic.fna.gz ]]; then
            gunzip GCF_000001635.27_GRCm39_genomic.fna.gz
            echo "Mus musculus genome was downloaded"
            Mus="$(readlink -f GCF_000001635.27_GRCm39_genomic.fna)"
            break
        else
            echo "Try downloading again the Mus musculus reference genome"
        fi
    done
fi

if [[ ! -e "$PhiX" ]]; then
    for i in {1..10}; do
        echo "Downloading PhiX reference genome"
        esearch -db nucleotide -query "NC_001422.1" | efetch -format fasta > PhiX_NC_001422.1.fasta
        if [[ -s PhiX_NC_001422.1.fasta ]]; then
            echo "PhiX genome was downloaded"
            PhiX=$(readlink -f PhiX_NC_001422.1.fasta)
            break
        else
            echo "Try downloading again the PhiX reference genome"
        fi
    done
fi

# Verify that input directory is set and exists
validate_input_dir

# Create output directory if it doesn't exists
validate_output_dir

# Output info
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Input directory: $input_dir"
echo "Output directory: $out_dir"
echo "Number of threads: ${threads:=1}"
echo "Phix Genome: ${PhiX:?'PhiX genome not set'}"
echo "Mus musculus Genome: ${Mus:?'Mus musculus genome not set'}"
echo "Bowtie2 called with options: $bowtie2_opts"

# Generate the mixed fasta and its index
# Small index if the genome is small
if [[ ! -f "$host" ]] && [[ ! -f Mix.rev.2.bt2 ]]; then
    # Pool together the Mus and host genome
    cat "$PhiX" "$Mus" > Mixed.fasta

    # Generate a small index
    echo "Small index needs to be generated"
    bowtie2-build Mixed.fasta Mix --threads "$threads"

# Large index if the genome is large due to the presence of the host genome
elif [[ -f "$host" ]] && [[ ! -f Mix.rev.2.bt2l ]]; then
    # Pool together the host, phage and Mus musculus genome
    cat "$host" "$PhiX" "$Mus" > Mixed.fasta

    # Generate a large index
    echo "Large index needs to be generated"
    bowtie2-build --large-index Mixed.fasta Mix --threads "$threads"
fi

# Alignment
for file in "$input_dir"/*; do
    # Paired end reads
    if [[ "$file" == @(*_R1_*|*_1).@(fq|fastq|fq.gz|fastq.gz) ]]; then
        forward_file="$file"
        core_name=$(get_core_name "$forward_file")
        bowtie2 -x Mix -1 "$forward_file" -2 $(forward_to_reverse "$forward_file") -p "$threads" \
            -q --un-conc-gz "$out_dir"/$(echo "$core_name" | sed 's/_trim/_bt2/').fq.gz \
            2> "$out_dir"/$(echo "$core_name" | sed 's/_trim/_bt2/')_summary.txt \
            $bowtie2_opts > /dev/null # Bowtie2 output to terminal is excesive and we do not need it in this case

    # Unpaired reads
    elif [[ ! "$file" ==  *_@(R1_|1.|R2_|2.)* ]] && [[ "$file" == *.@(fq|fastq|fq.gz|fastq.gz) ]]; then
        unpaired_file="$file"
        core_name=$(get_core_name "$unpaired_file")
        bowtie2 -x Mix -U "$unpaired_file" -q -p "$threads" \
            --un-gz "$out_dir"/$(basename -- "$unpaired_file" | sed 's/_trim/_bt2/') \
            2> "$out_dir"/$(get_core_name "$unpaired_file" | sed 's/_trim/_bt2/')_summary.txt \
            $bowtie2_opts > /dev/null

    # Files that do not match the required extension
    elif [[ ! "$file" == *.@(fq|fastq|fq.gz|fastq.gz) ]]; then
        echo -e "$(basename -- "$file") will not be processed as is not a .fastq or .fq.gz file."
    fi
done

# Rename bowtie2 output files
rename "s/.fq.1.gz/_1.fq.gz/" *.fq.1.gz 2> /dev/null
rename "s/.fq.2.gz/_2.fq.gz/" *.fq.2.gz 2> /dev/null

# Set the correct file extension (.gz) for unpaired output files
rename -f "s/fq/fq.gz/" *.fq 2> /dev/null
rename -f "s/fastq/fastq.gz/"  *.fastq 2> /dev/null

echo "Done."
echo "You can now use clean reads to:"
echo "- Assemble genomes using metabiome metaspades or megahit"
echo "- Perform taxonomic binning with metabiome kraken2, kaiju or metaphlan3"
echo "- Perform functional profiling using metabiome humann"
echo "- Extract 16S sequences with metabiome bbduk"

------------------------------------------------------------------------------------------------------------------------------------------------------------------


#Reads límpios en la carpeta Clean_reads: /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/Clean_reads

#Se realiza la asignación taxonómica exploratoria en la plataforma OneCodex https://app.onecodex.com/
#Se evaluan aspectos generales pero se observa clasificación solo a nivel de familia de Lachnospiraceae bacterium
#Se realiza asignación taxonómica con Kaiju (Web Server)
#Se explora la herramienta Centrifuge que mostró mejor resolución a nivel de especie.

#Clasificación con centrifuge

#Se creó un ambiente conda con la herramienta siguiendo las instrucciones de https://genomics.sschmeier.com/ngs-taxonomic-investigation/index.html#centrifuge:

conda create --yes -n centrifuge centrifuge
conda activate centrifuge

#Se generó la base de datos de centrifuge.

curl -O ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/p_compressed+h+v.tar.gz

# alternatively we can use wget

wget ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/p_compressed+h+v.tar.gz

# once the download is finished, we need to extract the archive content
# It will extract a few files from the archive and may take a moment to finish.

tar -xvzf p_compressed+h+v.tar.gz

#Se corre la Herramienta. En el argumento -x se pone la ruta donde quedó la base previamente generada. EL paso previo genera p_compressed+h+v con números.
#Simplemente ponerlo sin ningún número. -1 y -2 corresponde a los reads PE. Cambiar los nombres de report.txt y results.txt
#Se hace esto para cada muestra

centrifuge -x p_compressed+h+v -1 example.1.fq -2 example.2.fq -U single.fq --report-file report.txt -S results.txt

#Para poder construir gráficas con Pavian en RStudio se necesita un formato tipo Kraken. Por ello se debe transformar el archivo results.txt
#a kreport.txt
#Se hace esto para cada muestra

centrifuge-kreport -x p_compressed+h+v results.txt > kreport.txt

#Con este kreport se puede realizar una visualización de estadísticas y asignación con Pavian en RStudio.
#Se hace esto para cada muestra
#Se abre RStudio y se ejecuta:

if (!require(remotes)) { install.packages("remotes") }
remotes::install_github("fbreitwieser/pavian")

pavian::runApp(port=5000)


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rsamtools")


#Esto genera una ventana interactiva donde se pueden ver varios análisis.
#De allí vale la pena descargar los CSVs con las asignaciones para Bacterias-Virus-Eucariotas que se consideren.

#Se genera una figura exploratoria con los 10 taxones más abundantes.

----------------------------------------------------------------------------------

# library
library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)
library(viridis)
library(hrbrthemes)
library(forcats)
library(cowplot)


#Se pasa de formato ancho a largo para graficar

cruzi_pivot <- BL6_virus %>% pivot_longer(

  cols = 2:11,

  names_to = "Species",

  values_to = "Abundance",

  values_drop_na = TRUE)

# Stacked + percent
ggplot(cruzi_pivot_BL6_virus, aes(fill=Species, y=Abundance, x=Name)) +
  geom_bar(position="fill", stat="identity") +
  scale_fill_brewer(palette = "Spectral")

--------------------------------------------------------------------------------

#Se procede a hacer de manera paralela los ensamblajes.
#Se usaron Spades y Megahit.
# https://github.com/voutcn/megahit
# https://github.com/ablab/spades#sec2

#Se hace esto para cada muestra
#Para el uso de Spades se usa como se muestra en el siguiente ejemplo:

/home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/Herramientas/SPAdes-3.15.3-Linux/bin/metaspades.py -t 24 -1 /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/Clean_reads/BL6/BL610NI_R1.fq.gz -2 /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/Clean_reads/BL6/BL610NI_R2.fq.gz -o BL6_10NI

#Se hace esto para cada muestra
#Para el uso de Megahit como en el siguiente ejemplo:

/home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/Herramientas/megahit/build/megahit -t 24 -1 /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/Clean_reads/BL6/BL60NI_R1.fq.gz -2 /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/Clean_reads/BL6/BL60NI_R2.fq.gz -o BL60NI


#Después de evalua la calidad de los ensamblajes para comparar las herramientas y resultados por mediode metaquast http://cab.cc.spbu.ru/quast/manual.html#sec2.4
#Se hace esto para cada muestra
#La instrucción ejemplo es:

/home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/Herramientas/quast-5.0.2/metaquast.py --threads 18 /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/spades/Balbc/BALBc0NI/contigs.fasta /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/megahit/BALBc/BALBc0NI/final.contigs.fa -o 0_NI
/home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/Herramientas/quast-5.0.2/metaquast.py --threads 18 /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/spades/BL6/BL6_0NI/contigs.fasta /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/megahit/BL6/BL60NI/final.contigs.fa -o 0_NI


#Donde contigs.fasta corresponde al ensamblaje de Metaspades y el final.contigs al de Megahit.
#Cuando se presentan inconvenientes en el análisis que arrojan dificultades de conexión a red o similares, se puede aplicar el parámetro: --max-ref-num <int>
#Indicando --max-ref-num 0 
#De esta manera se corre el análisis de estadísticos básicos y se obtiene el reporte habitual.

#Paralelamente se realiza el análisis funcional a partir de los reads con Humann3 (https://github.com/biobakery/biobakery/wiki/humann3) (https://github.com/biobakery/humann)
#Se sigue el proceso indicado para su instalación y uso https://huttenhower.sph.harvard.edu/humann/
#Importante trabajar la instalación y uso dentro de conda biobakery3 para que las dependencias y demás funcionen correctamente

conda create --name biobakery3 python=3.7
conda activate biobakery3
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels biobakery
conda install humann -c biobakery
humann test

#Base de Datos Pangenoma
humann_databases --download chocophlan full  /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/Herramientas/databases_humann/ --update-config yes


#Base de datos proteínas
humann_databases --download uniref uniref90_diamond /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/Herramientas/databases_humann/ --update-config yes

#Base de datos anotaciones

humann_databases --download utility_mapping full /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/Herramientas/databases_humann/ --update-config yes

# Se instaló dentro de este ambiente Bowtie2 para su correcto funcionamiento

conda install bowtie2

#Se ejecuta el siguiente comando para uso de Humann3. El input en este caso corresponde al mejor ensamblaje obtenido, ya sea por spades o megahit. Se debe trabajar en la carpeta
#donde se encuentra las bases previamente creadas para que funcione.
#Así mismo asegurarse de seguir dentro del ambiente biobakery

humann --threads 18 -i /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/spades/BL6/BL6_0NI/contigs.fasta -o /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/humann/BL6_0NI

humann --threads 18 -i /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/spades/BL6/BL6_7DPI/contigs.fasta -o /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/humann/BL6_7DPI


#Con los siguientes comandos se puede visualizar e inspeccionar los output arrojados que son 3 archivos más una carpeta con temp.


column -t -s $'\t' contigs_genefamilies.tsv | less -S

column -t -s $'\t' contigs_pathabundance.tsv | less -S


# Para facilitar las comparaciones entre muestras con diferentes profundidades de secuenciación, es importante normalizar los valores RPK predeterminados de HUMAnN antes de realizar análisis estadísticos. 
#La suma-normalización a abundancia relativa o unidades de "copias por millón" (CPM) son enfoques comunes para esto; estos métodos se implementan en el humann_renorm_tablescript 
#(con consideraciones especiales para los formatos de salida de HUMAnN, por ejemplo, estratificación de características). Normalicemos nuestra abundancia de genes a unidades CPM

humann_renorm_table --input contigs_genefamilies.tsv --output contigs_genefamilies-cpm.tsv --units cpm --update-snames

#Con los siguientes comandos se puede visualizar e inspeccionar los output arrojados.

column -t -s $'\t' contigs_genefamilies-cpm.tsv | less -S

# Las "unidades" predeterminadas de la función microbiana de HUMAnN son las familias de genes UniRef (que usamos internamente para calcular las abundancias de reacción y, 
# a partir de ahí, las abundancias de la vía). Sin embargo, a partir de la abundancia de familias de genes, es posible reconstruir la abundancia de otras categorías funcionales 
# en un microbioma utilizando el humann_regroup_tablescript. Inspeccione las opciones de reagrupación usando el humann_regroup_table -hcomando.

# Reagrupamos nuestros valores de abundancia de la familia de genes normalizados por CPM a las abundancias de la reacción MetaCyc (RXN), que se incluyen con la instalación predeterminada de HUMAnN


humann_regroup_table --input contigs_genefamilies-cpm.tsv  --output BL616NI_genefamilies_uniref90rxn-cpm.tsv --groups uniref90_rxn


#generalmente es útil adjuntar algunas descripciones legibles por humanos de estos ID para facilitar la interpretación biológica

humann_rename_table --input BL616NI_genefamilies_uniref90rxn-cpm.tsv --output BL616NI_genefamilies_uniref90rxn-cpm-named.tsv --names metacyc-rxn


# LLevar todos los archivos finales a una misma carpeta para unir todos los output -named.tsv. Se ejecuta el siguiente comando

humann_join_tables -i . -o BL6_Final_Humann.tsv --file_name uniref90rxn-cpm

#Para graficar correctamente usando Humann, se debe obtener en Humann3 una tabla estratificada con el siguiente comando:

humann_split_stratified_table --input /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/humann/BL6_Humann_Final/BL6_Final_Humann.tsv --output /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/humann/BL6_Humann_Final/BL6_Final_Humann_stritified


# Se descarga este archivo estratificado para trabajar en Maaslin2 en R Studio (https://github.com/biobakery/biobakery/wiki/maaslin2#1-installing-r)
# Importante revisar la tabla descargada, organizarla de acuerdo a la estructura que requiere la herramienta. No recibe nombres de columnas tan largos
# Es necesario ajustar manualmente 
# En RStudio se trabaja de la siguiente manera:

--------------------------------------------------------------------------------------------------------------------------------

library(dplyr)

# Se puede ralizar de manera manual el cargue de los archivos tsv resultantes

BALBc0_NI = read.table(file = "/Users/scast/Desktop/SERGIO/Doctorado/Tesis/Modelo_Tcruzi/Clean_Reads/humann/BALBc0NI_pathabundance.tsv", header = TRUE, sep = "\t")
BALBc3_DPI = read.table(file = "/Users/scast/Desktop/SERGIO/Doctorado/Tesis/Modelo_Tcruzi/Clean_Reads/humann/BALBc3DPI_pathabundance.tsv", header = TRUE, sep = "\t")
BALBc7_DPI = read.table(file = "/Users/scast/Desktop/SERGIO/Doctorado/Tesis/Modelo_Tcruzi/Clean_Reads/humann/BALBc7DPI_pathabundance.tsv", header = TRUE, sep = "\t")
BALBc10_NI = read.table(file = "/Users/scast/Desktop/SERGIO/Doctorado/Tesis/Modelo_Tcruzi/Clean_Reads/humann/BALBc10NI_pathabundance.tsv", header = TRUE, sep = "\t")
BALBc10_DPI = read.table(file = "/Users/scast/Desktop/SERGIO/Doctorado/Tesis/Modelo_Tcruzi/Clean_Reads/humann/BALBc10DPI_pathabundance.tsv", header = TRUE, sep = "\t")
BALBc13_DPI = read.table(file = "/Users/scast/Desktop/SERGIO/Doctorado/Tesis/Modelo_Tcruzi/Clean_Reads/humann/BALBc13DPI_pathabundance.tsv", header = TRUE, sep = "\t")
BALBc16_NI = read.table(file = "/Users/scast/Desktop/SERGIO/Doctorado/Tesis/Modelo_Tcruzi/Clean_Reads/humann/BALBc16NI_pathabundance.tsv", header = TRUE, sep = "\t")
BALBc16_DPI = read.table(file = "/Users/scast/Desktop/SERGIO/Doctorado/Tesis/Modelo_Tcruzi/Clean_Reads/humann/BALBc16DPI_pathabundance.tsv", header = TRUE, sep = "\t")

# Y posteriormente realizar un Inner Join para posterior análisis
# Sin embargo es mejor hacer el paso a paso de Maaslin2 para normalizar, reagrupar 
# y consolidar la tabla, directamente en línea de comando.

PathTotal= BALBc0_NI %>% inner_join(BALBc3_DPI,by="Pathway") %>% inner_join(BALBc7_DPI,by="Pathway") %>% inner_join(BALBc10_DPI,by="Pathway") %>% inner_join(BALBc10_NI,by="Pathway") %>% inner_join(BALBc13_DPI,by="Pathway") %>% inner_join(BALBc16_DPI,by="Pathway") %>% inner_join(BALBc16_NI,by="Pathway")
write.csv(PathTotal, "PathTotal.csv")


# Ejecución Maaslin2

if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Maaslin2")

library(Maaslin2)
library(readr)

# Como input data se puede poner el resultado de asignación taxonómica 
# así como también los resultados del análisis por Humann3
# Aquí se usa el .tsv consolidado obtenido por Humann3 en línea de comando como input_data

df_input_data = read.table(file = "/Users/scast/Desktop/SERGIO/Doctorado/Tesis/Modelo_Tcruzi/Clean_Reads/humann/BALBc_Final_Humann.tsv", header = FALSE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
df_input_data[1:5, 1:5]
df_input_data=as.data.frame(df_input_data)


# Input metadata correspondiente a las muestras evaluadas

df_input_metadata = read.table(file = "/Users/scast/Desktop/SERGIO/Doctorado/Tesis/Modelo_Tcruzi/Clean_Reads/humann/metadata.tsv", header = TRUE, sep = "\t",
                               row.names = 1,
                               stringsAsFactors = FALSE)
df_input_metadata[1:5, ]


# Ejecución del análisis

fit_func = Maaslin2(
  input_data = df_input_data, 
  input_metadata = df_input_metadata, 
  output = "demo_functional4", 
  fixed_effects = c("Condition"),
  min_abundance = 0.0001,
  min_prevalence = 0.1,
  normalization = "NONE",
  plot_heatmap = TRUE,
  plot_scatter = TRUE)

--------------------------------------------------------------------------------------------------------------------------------

# Ya conociendo que vías tienen significancia estadística se busca directamente con el término correspondiente en la tabla estratificada. A esta tabla se debe agregar la metadata por filas de aquello que queramos considerar 
# para cada muestra en función de la figura a construir

grep 'ENOYL-COA-HYDRAT-RXN' /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/humann/BALBc_stratified.tsv | less -S

#Se genera la figura

humann_barplot --input /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/humann/BALBc_Humann_Final/BALBc_Final_Humann_stritified.tsv/BALBc_stratified_Final.tsv --focal-metadata Time --last-metadata Condition --output BALBc_ENTDB-RXN.png --focal-feature 'ENTDB-RXN' --sort sum metadata --scaling logstack


#Binning

# Primero necesitamos mapear las lecturas contra el ensamblaje para obtener información de cobertura
# Mapear las lecturas para cada muestra con el ensamblaje nos brinda información de "cobertura" 
# para cada contig en cada muestra, lo que nos ayudará con nuestros esfuerzos para recuperar los MAGs.
# Se ejecuta el mapeo dentro del env biobakery que contiene Bowtie2 y samtools.

bowtie2-build /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/spades/Balbc/BALBc10NI/contigs.fasta final.contigs

bowtie2 -x final.contigs -1 /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/Clean_reads/BALBc/BALBc10NI_R1.fq.gz -2 /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/Clean_reads/BALBc/BALBc10NI_R2.fq.gz | samtools view -bS -o BALBc10NI_to_sort.bam

samtools sort BALBc10NI_to_sort.bam -o BALBc10NI.bam

samtools index BALBc10NI.bam


# Instalar dentro de biobakery metabat2 con 

conda install metabat2

# Ejecutar el siguiente comando que incluye el parámetro -m (--minContig) que se refiere al tamaño mínimo del binning, el ensamblaje por spades en este caso
# y el archivo .bam obtenido en pasos previos

runMetaBat.sh -m 1500 -t 18 /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/spades/Balbc/BALBc10NI/contigs.fasta /datacnmat01/biologia/gimur/proyectos/microbiomas/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/binning/BALBc10NI/BALBc10NI.bam


# Instalar dentro de biobakery Maxbin2 

conda install maxbin2

# En el folder correspondiente a cada mx correr el siguiente comando
# Los argumentos son el ensamblaje, las lecturas Fw y Rv de cada mx y el prefijo del
# https://carpentries-incubator.github.io/metagenomics/05-binning/index.html


run_MaxBin.pl -thread 16 -contig /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/spades/Balbc/BALBc10NI/contigs.fasta -reads /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/Clean_reads/BALBc/BALBc10NI_R1.fq.gz -reads2 /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/Clean_reads/BALBc/BALBc10NI_R2.fq.gz -out BALBc10NI


# Ahora se realiza el binning con concoct
# Se realiza la instalación en la carpeta herramientas siguiendo las intrucciones tomadas de https://www.giters.com/antagomir/CONCOCT


git clone https://github.com/BinPro/CONCOCT.git
cd CONCOCT
pip install -r requirements.txt
python setup.py install


# Seguido se ejecuta en la carpeta donde quedarán los resultados por cada mx (https://github.com/BinPro/CONCOCT):
# Cortar contigs en partes más pequeñas


cut_up_fasta.py /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/spades/Balbc/BALBc16NI/contigs.fasta -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa

# Genere una tabla con información de profundidad de cobertura por muestra y subcontig. 
# Este paso asume que el directorio 'mapeo' contiene archivos bam ordenados e indexados donde cada muestra ha sido mapeada contra los contigs originales.

concoct_coverage_table.py contigs_10K.bed /datacnmat01/biologia/gimur/proyectos/microbiomas/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/binning/BALBc16NI/BALBc16NI.bam > coverage_table.tsv

# Ejecutar binning

concoct --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -b concoct_BALBc16NI/

# Fusionar la agrupación en clústeres de subcontigos en la agrupación en clústeres de contig original

merge_cutup_clustering.py concoct_BALBc16NI/clustering_gt1000.csv > concoct_BALBc16NI/clustering_merged.csv

# Extraer bins como FASTA individual

mkdir concoct_BALBc16NI/fasta_bins

extract_fasta_bins.py /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/spades/Balbc/BALBc10NI/contigs.fasta concoct_BALBc10NI/clustering_merged.csv --output_path concoct_BALBc10NI/fasta_bins


#CheckM
# La primera vez que se usa se debe realizar el siguiente procedimiento
# De aquí en adelante esta herramienta se trabaja dentro del env checkm

conda create -n checkm python=3.9
conda activate checkm
conda install numpy matplotlib pysam
conda install hmmer prodigal pplacer
pip3 install checkm-genome

wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar -xvzf 

#Donde se descomprimió el archivo anterior ejecutar: 
checkm data setRoot .

# Comprobar instalación

checkm

# Comprobar funcionamiento

checkm test ~/checkm_test_results

# Ejecutar esto para cada muestra para cada resultado dependiendo de la herramienta usada (metabat, maxbin, concoct)

checkm lineage_wf -t 18 -x fasta /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/binning/Maxbin/BALBc16NI /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/checkM/Muestras_BALBc_Maxbin/BALBc16NI


# Cuando se obtiene el resultado previo, en cada carpeta resultante ingresar y ejectutar el siguiente comando para poder ver los estadísticos de los bins de cada unas de las herramientas usadas.
# qa es la función usada de checkm, lineage.ms es el archivo resultante de la ejecución anterior, el . es para indicar que los output de checkm están allí
# el -o 2 es para indicarle que es en el formato 2 que incluye todos los estadísticos de calidad y el -f es para nombrar el archivo de salida

checkm qa lineage.ms . -o 2 -f qa


# Usar DAS Tool

# Cargar ambiente virtual:

source activate /home/sergio.castaneda/data/conda/conda_envs/das_tool


# cargar modulo usearch

module load usearch

# Desde la carpeta de DAS Tool creada en herramientas inicialmente se deben obtener los archivos necesarios como input para DAS Tool
# Por lo anterior se debe ejecutar la siguiente línea para cada proceso de binning por las 3 herramientas usadas para cada una de las muestras

src/Fasta_to_Scaffolds2Bin.sh -i /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/binning/Maxbin/BALBc16NI -e fasta > /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/binning/DASTool/BALBc16NI/maxbin.scaffolds2bin.tsv
src/Fasta_to_Scaffolds2Bin.sh -i /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/binning/Metabat/BALBc16NI/contigs.fasta.metabat-bins18-20211222_134206 -e fa > /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/binning/DASTool/BALBc16NI/metabat.scaffolds2bin.tsv

# Algunas herramientas de agrupamiento (como CONCOCT) proporcionan una salida tabular separada por comas. 
# Para convertir una coma archivo separado en un archivo .tsv se puede utilizar:

perl -pe "s/,/\tconcoct./g;" /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/binning/concoct/BALBc16NI/concoct_BALBc16NI/clustering_gt1000.csv > /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/binning/DASTool/BALBc16NI/concoct.scaffolds2bin.tsv


# Eejecutar DAS Tool

DAS_Tool -i /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/binning/DASTool/BALBc10NI/maxbin.scaffolds2bin.tsv,/home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/binning/DASTool/BALBc10NI/metabat.scaffolds2bin.tsv,/home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/binning/DASTool/BALBc10NI/concoct.scaffolds2bin.tsv -l maxbin,metabat,concoct -c /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/spades/Balbc/BALBc10NI/contigs.fasta -t 16  --write_bins 1 -o /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/binning/DASTool/BALBc10NI/BALBc10NI_Result

DAS_Tool -i /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/binning/DASTool/BALBc16NI/maxbin.scaffolds2bin.tsv,/home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/binning/DASTool/BALBc16NI/metabat.scaffolds2bin.tsv -l maxbin,metabat -c /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/spades/Balbc/BALBc16NI/contigs.fasta -t 16  --write_bins 1 -o /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/binning/DASTool/BALBc16NI/BALBc16NI_Result


# CheckM para bins refinados


conda activate checkm

checkm lineage_wf -t 18 -x fa /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/binning/DASTool/BALBc0NI/BALBc0NI_Result_DASTool_bins /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/checkM/Muestras_BALBc_DASTool/BALBc10NI -f BALBc10NI.txt
checkm qa lineage.ms . -o 2 -f qa

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Se realiza clasificación taxonómica de los bins resultantes de DASTool por cada una de las muestras.
#Para tal propósito se utiliza la herramienta GTDBTK se acuerdo con referencias previas. (https://ecogenomics.github.io/GTDBTk/index.html)
#Inicialmente se efectúa la instalación de la herramienta usando conda

conda create -n gtdbtk -c conda-forge -c bioconda gtdbtk

#Se debe ejecutar el siguiende comando para que se instale adecuadamente la base de datos. Este paso piede tomar varios días

download-db.sh

#Se ejecuta el Workflow gtdbtk classify wf

#El flujo de trabajo classify wf consta de tres pasos: identify, aligny y classify.

# identify llama a los genes usando Prodigal y usa modelos HMM y el paquete HMMER para identificar los 120 genes marcadores bacterianos y 122 arqueales usados para la inferencia filogenética ( Parks et al., 2018 ). 
# Los alineamientos de secuencias múltiples (MSA) se obtienen alineando los genes marcadores con su modelo HMM respectivo.
# align concatena los genes marcadores alineados y filtra el MSA concatenado a aproximadamente 5000 aminoácidos.
# Finalmente, classify utiliza pplacer para encontrar la ubicación de máxima probabilidad de cada genoma en el árbol de referencia GTDB-Tk. 
# GTDB-Tk clasifica cada genoma en función de su ubicación en el árbol de referencia, su divergencia evolutiva relativa y/o la identidad de nucleótidos promedio (ANI) con respecto a los genomas de referencia.

# se debe especificar el --e para indicar que extensión tiene los archivos (fasta, fa, fna)
# --genome_dir indica donde están los bins a clasificar
# --out_dir donde quedarán los archivos
# --pplacer_cpus 8 --scratch_dir son muy importantes para evitar inconvenientes del uso de memoria del servidor, principalmente lo relacionado
# con la ejecución de pplacer. Se recomienda no usar muchos ya que el consumo que requiere por thread es muy alto.
# se realizaron varios test con diferentes num de cpus y con 8 funcionó apropiadamente


gtdbtk classify_wf --e fa --genome_dir /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/binning/DASTool/BALBc13DPI/BALBc13DPI_Result_DASTool_bins --out_dir /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/GTDBTK/BALBc13DPI --pplacer_cpus 4 --scratch_dir .





conda install -c conda-forge -c bioconda -c defaults prokka
conda install -c bioconda perl-bioperl
prokka --version
prokka --outdir /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/prokka --prefix BALBc0NI /home/sergio.castaneda/7_19_Metagenomics_SAC/usftp21.novogene.com/analisis_manual/binning/DASTool/BALBc0NI/BALBc0NI_Result_DASTool_bins


Pasos con los bins:

Extraer 16s Barnap
Genoma candidato genoma referencia
ANI 
Anotación Prokka 
Filogenia



