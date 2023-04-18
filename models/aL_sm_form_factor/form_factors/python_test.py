import re

# Specify the path to your header file
header_file_path = "form_factors.h"

# Create an empty list to store function names
function_names = []

# Open the header file in read mode
with open(header_file_path, "r") as header_file:
    # Read the contents of the file
    file_contents = header_file.read()
    # Use regular expressions to find function definitions
    pattern = r"(?<=\n)(?:[a-zA-Z0-9_]+\s)+([a-zA-Z0-9_]+)\(.*?\);"
    matches = re.findall(pattern, file_contents)
    # Add the function names to the list
    function_names.extend(matches)

# Print the list of function names
print(function_names)
