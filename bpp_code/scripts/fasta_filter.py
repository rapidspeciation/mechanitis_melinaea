import sys
import re
import os

def read_fasta(filename):
    sequences = {}
    with open(filename, 'r') as file:
        sequence_id = None
        sequence = ''
        for line in file:
            if line.startswith('>'):
                if sequence_id:
                    sequences[sequence_id] = sequence
                sequence_id = line.strip()[1:]
                sequence = ''
            else:
                sequence += line.strip()
        if sequence_id:
            sequences[sequence_id] = sequence
    return sequences


def read_phylip(filename):
    sequences = {}
    with open(filename, 'r') as file:
        header = file.readline().strip().split()
        num_sequences = int(header[0])
        alignment_length = int(header[1])
        for _ in range(num_sequences):
            line = file.readline().strip().split()
            seq_id = line[0]
            sequence = ''.join(line[1:])
            sequences[seq_id] = sequence
    return sequences


def remove_elements_from_values(dictionary, indices_to_remove):
    for key, value in dictionary.items():
        if isinstance(value, str):
            modified_value = ''.join([value[i] for i in range(len(value)) if i not in indices_to_remove])
            dictionary[key] = modified_value
    return(dictionary)


def read_alignment(filename):
    file_extension = os.path.splitext(filename)[1].lower()
    if file_extension == '.fasta' or file_extension == '.fa':
        return read_fasta(filename)
    elif file_extension == '.phylip':
        return read_phylip(filename)
    else:
        raise ValueError("Unsupported file format. Please provide a FASTA or PHYLIP file.")


def read_population_assignment(filename):
    population_assignment = {}
    with open(filename, 'r') as file:
        for line in file:
            individual, population = line.strip().split()
            population_assignment[individual] = population
    return population_assignment


def read_min_individuals_per_population(filename):
    min_individuals_per_population = {}
    with open(filename, 'r') as file:
        for line in file:
            population, min_individuals = line.strip().split()
            min_individuals_per_population[population] = int(min_individuals)
    return min_individuals_per_population


def filter_alignment(sequences, max_missing_per_sequence, max_missing_per_position, max_iupac_per_position):
    filtered_sequences = {}
    alignment_length = len(next(iter(sequences.values())))
    for seq_id, seq in sequences.items():
        missing_count = seq.count('-') + seq.count('N')  # Counting gaps and ambiguous characters
        if missing_count / alignment_length <= max_missing_per_sequence:
            filtered_sequences[seq_id] = seq

    bad_positions = []
    for i in range(alignment_length):
        column = ''.join(seq[i] for seq in filtered_sequences.values())
        missing_count = column.count('-') + column.count('N')  # Counting gaps and ambiguous characters
        iupac_count = len(re.findall(r'[RYSWKMBDHVN]', column, flags=re.IGNORECASE))  # Counting IUPAC characters
        if missing_count / len(filtered_sequences) > max_missing_per_position or iupac_count / len(filtered_sequences) > max_iupac_per_position or (missing_count+iupac_count)/len(filtered_sequences) > (max_missing_per_position+max_iupac_per_position):
            bad_positions.append(i)

    filtered_alignment = remove_elements_from_values(filtered_sequences, bad_positions)
    return(filtered_alignment)


def check_min_individuals_per_population(population_assignment, filtered_sequences, min_individuals_per_population):
    alnPass = 1
    population_counts = {}
    for individual, population in population_assignment.items():
        if individual in filtered_sequences:
            population_counts[population] = population_counts.get(population, 0) + 1
    for population, min_individuals in min_individuals_per_population.items():
        if population not in population_counts or population_counts[population] < min_individuals:
            print(population, population_counts[population])
            return False
    return True
# def check_min_individuals_per_population(population_assignment, filtered_sequences, min_individuals_per_population):
#     population_counts = {}
#     for individual, population in population_assignment.items():
#         if population in filtered_sequences:
#             population_counts[population] = population_counts.get(population, 0) + 1

#     for population, min_individuals in min_individuals_per_population.items():
#         if population in population_counts and population_counts[population] < min_individuals:
#             return False
#     return True


def write_fasta(sequences, output_filename):
    with open(output_filename, 'w') as file:
        for seq_id, seq in sequences.items():
            file.write('>' + seq_id + '\n')
            file.write(seq + '\n')


if __name__ == "__main__":
    if len(sys.argv) != 8:
        print("Usage: python filter_alignment.py input_file output_file max_missing_per_sequence max_missing_per_position max_iupac_per_position population_assignment_file min_individuals_file")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    max_missing_per_sequence = float(sys.argv[3])
    max_missing_per_position = float(sys.argv[4])
    max_iupac_per_position = float(sys.argv[5])
    population_assignment_file = sys.argv[6]
    min_individuals_file = sys.argv[7]


    sequences = read_alignment(input_file)
    population_assignment = read_population_assignment(population_assignment_file)
    min_individuals_per_population = read_min_individuals_per_population(min_individuals_file)


    filtered_sequences = filter_alignment(sequences, max_missing_per_sequence, max_missing_per_position, max_iupac_per_position)
    
    if check_min_individuals_per_population(population_assignment, filtered_sequences, min_individuals_per_population):
        write_fasta(filtered_sequences, output_file)
        print("Filtered alignment saved to", output_file)
    else:
        print("Minimum number of individuals per population requirement not met.")
