import shutil
import os, sys
from pathlib import Path

# Define the path to the zip file
zip_file_path = 'my_archive.zip'

# Define the directory where the contents should be extracted
extract_directory = 'extracted_files'

def unzip_file(zip_file_path, extract_directory):
	zip_file_path = Path(zip_file_path)
	extract_directory = Path(extract_directory)
	if not zip_file_path.exists() or not zip_file_path.is_file():
		sys.exit("The provided path to the zipped file does not exists. Unzipping aborted. Files: ",zip_file_path.as_posix())
	if not extract_directory.exists() or not extract_directory.is_dir():
		sys.exit("The provided output directory to store the unzipped file does not exists. Unzipping aborted. Directory:",
		         extract_directory.as_posix())
	print("* Unzipping file ", zip_file_path.name, " *")
	# Unzip the archive
	try:
		shutil.unpack_archive(zip_file_path, extract_directory)
		print(f"Successfully unzipped to '{extract_directory.name}'")
	except shutil.ReadError as e:
		print(f"Error unzipping file: {e}")
	except FileNotFoundError:
		print(f"Error: The file '{zip_file_path}' was not found.")
	except Exception as e:
		print(f"An unexpected error occurred: {e}")



if __name__ == "__main__":
	if len(sys.argv) > 2:
		zip_file_path = sys.argv[1]
		extract_directory = sys.argv[2]
	else:
		print("Error: Two arguments must be supplied to unzip a file to a given directory:\n"
			"  1 - zip_file_path: A compressed zip file (.zip);\n"
			"  2 - extract_directory: Path to the desired directory to unzip and extract the zip file.\n")
		sys.exit(1)
	# call final report creation
	unzip_file(zip_file_path, extract_directory)