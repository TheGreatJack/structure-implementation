import argparse
import numpy as np


def parse_arguments():
  """Parses command line arguments using argparse.

  Returns:
    A namespace object containing parsed arguments.
  """
  parser = argparse.ArgumentParser(description="Process STRUCTURE data file")
  parser.add_argument("input_file", type=str, help="Path to the STRUCTURE data file")
  return parser.parse_args()


def load_structure_data(filename):
  """Loads genotype data from a STRUCTURE format file.

  Args:
    filename: Path to the STRUCTURE data file.

  Returns:
    A NumPy array containing the genotype data.
  """
  with open(filename, 'r') as f:
    # Read the first line containing meta information
    header_line = f.readline().strip().split()
    num_loci = int(header_line[0])
    num_individuals = int(header_line[1])

    # Initialize NumPy array to store genotypes
    genotype_data = np.zeros((num_individuals, num_loci), dtype=int)

    # Read remaining lines containing genotypes
    for i, line in enumerate(f):
      genotypes = line.strip().split()
      # Convert genotypes to integers (assuming integer encoding)
      genotype_data[i] = np.array([int(g) for g in genotypes])

  return genotype_data


def main():
  # Parse command line arguments
  args = parse_arguments()

  # Load data from STRUCTURE file
  genotype_data = load_structure_data(args.input_file)

  # Print information about loaded data (optional)
  print(f"Loaded genotype data with dimensions: {genotype_data.shape}")

  # You can now use the genotype_data NumPy array for further analysis

if __name__ == "__main__":
  main()
