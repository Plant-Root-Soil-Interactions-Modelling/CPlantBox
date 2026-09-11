from pathlib import Path

# Folder where this script is located
root_folder = Path(__file__).resolve().parent

# Go through all files in the folder and subfolders
for file_path in root_folder.rglob("*"):
    if file_path.is_file() and "Zone.Identifier" in file_path.name:
        try:
            file_path.unlink()
            print(f"Deleted: {file_path}")
        except Exception as e:
            print(f"Could not delete {file_path}: {e}")

print("Done.")