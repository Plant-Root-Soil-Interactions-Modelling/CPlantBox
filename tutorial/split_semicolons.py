#!/usr/bin/env python3
""" fixes E702 for subfolder (python3 split_semicolons.py path/to/folder) """

import os
import re
import sys

# Folder to process (can be current folder)
folder = sys.argv[1] if len(sys.argv) > 1 else "."

# Regex to match semicolon-separated statements
semicolon_regex = re.compile(r";\s*")

# Walk through all Python files
for root, dirs, files in os.walk(folder):
    for file in files:
        if file.endswith(".py"):
            file_path = os.path.join(root, file)

            with open(file_path, "r", encoding="utf-8") as f:
                lines = f.readlines()

            new_lines = []
            changed = False

            for line in lines:
                # Skip lines that are only whitespace
                if line.strip() == "":
                    new_lines.append(line)
                    continue

                # Split by semicolon, but preserve indentation
                parts = semicolon_regex.split(line)
                if len(parts) > 1:
                    changed = True
                    indent_match = re.match(r"(\s*)", line)
                    indent = indent_match.group(1) if indent_match else ""
                    for i, part in enumerate(parts):
                        stripped = part.strip()
                        if stripped:
                            new_lines.append(f"{indent}{stripped}\n")
                else:
                    new_lines.append(line)

            # Overwrite file if changes were made
            if changed:
                with open(file_path, "w", encoding="utf-8") as f:
                    f.writelines(new_lines)
                print(f"Updated {file_path}")
