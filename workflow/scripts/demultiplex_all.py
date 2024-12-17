import os
import sys


def create_token_files(bcl_run_indices, fastq_dirs):
    # Convert the inputs to a dictionary
    bcl_run_index_dir_dict = {
        bcl_run_index: dirs.split(",")
        for bcl_run_index, dirs in zip(
            bcl_run_indices.split(";"), fastq_dirs.split(";")
        )
    }

    # Generate tokens and create files
    tokens = []
    for bcl_run_index, fastq_dirs in bcl_run_index_dir_dict.items():
        for fastq_dir in fastq_dirs:
            # Ensure the directory exists
            os.makedirs(fastq_dir, exist_ok=True)

            # Define the token file path
            token_file = os.path.join(fastq_dir, "sample_demultiplexed.txt")

            # Create the token file
            with open(token_file, "w") as f:
                f.write("success")

            tokens.append(token_file)

    return tokens


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python demultiplex.py <bcl_run_indices> <fastq_dirs>")
        sys.exit(1)

    # Get arguments
    bcl_run_indices = sys.argv[1]
    fastq_dirs = sys.argv[2]

    # Call the function and get the result
    tokens = create_token_files(bcl_run_indices, fastq_dirs)
