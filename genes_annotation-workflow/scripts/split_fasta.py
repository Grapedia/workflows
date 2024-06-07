# Importar los módulos SeqIO y sys de BioPython
from Bio import SeqIO
import sys

# Definir el path al fichero fasta de entrada como el primer argumento
input_file = sys.argv[1]

# Definir el path con el prefijo de los ficheros fasta de salida como el segundo argumento
output_prefix = sys.argv[2]

# Definir el número de partes a dividir
num_parts = 10

# Leer el fichero fasta de entrada y almacenar las secuencias en una lista
sequences = list(SeqIO.parse(input_file, "fasta"))

# Contar el número total de secuencias
num_sequences = len(sequences)

# Calcular el número de secuencias por cada parte
num_sequences_per_part = num_sequences // num_parts

# Inicializar el índice de la parte actual
part_index = 1

# Inicializar el índice de la secuencia actual
sequence_index = 0

# Iterar sobre las partes
while part_index <= num_parts:
    # Definir el nombre del fichero fasta de salida para la parte actual con un formato de dos dígitos y restando uno al índice
    output_file = f"{output_prefix}{part_index - 1:02d}.fasta"

    # Escribir las secuencias correspondientes a la parte actual en el fichero fasta de salida
    SeqIO.write(sequences[sequence_index:sequence_index + num_sequences_per_part], output_file, "fasta")

    # Incrementar el índice de la parte actual
    part_index += 1

    # Incrementar el índice de la secuencia actual
    sequence_index += num_sequences_per_part

# Si hay secuencias sobrantes, añadirlas a la última parte
if sequence_index < num_sequences:
    # Definir el nombre del fichero fasta de salida para la última parte con un formato de dos dígitos y restando uno al índice
    output_file = f"{output_prefix}{num_parts - 1:02d}.fasta"

    # Escribir las secuencias sobrantes en el fichero fasta de salida
    SeqIO.write(sequences[sequence_index:], output_file, "fasta")
