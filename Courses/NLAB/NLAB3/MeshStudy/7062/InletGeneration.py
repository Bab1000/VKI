import numpy as np
 
def velocity_profile(y, U_m=0.3, H=0.41):
    return 4 * U_m * y * (H - y) / H**2
 
# Read the inlet.dat file and process it
file_path = 'example_inlet.dat'
with open(file_path, 'r') as file:
    lines = file.readlines()
 
# Process lines and substitute the velocity
# Process lines and substitute the velocity
new_lines = []
for line in lines:
    # Skip comment lines or empty lines
    if line.startswith('#') or not line.strip():
        new_lines.append(line)
        continue
   
    # Split the line by whitespace (this will handle both spaces and tabs)
    values = line.split()
   
    # Ensure there are at least 6 values (since the file has 6 columns)
    if len(values) >= 6:
        # Get the y value from the second column
        y = float(values[1])
       
        # Calculate the velocity
        velocity = velocity_profile(y)
       
        # Replace the velocity value in the 4th column (index 3)
        values[3] = f"{velocity:.15e}"  # Use scientific notation like the input
       
        # Reassemble the line
        new_line = "\t".join(values) + "\n"
        new_lines.append(new_line)
    else:
        # If a line doesn't have the expected format, just keep it as is
        new_lines.append(line)
 
# Write the modified content back to the file or a new file
output_file_path = 'inlet.dat'
with open(output_file_path, 'w') as file:
    file.writelines(new_lines)
 
print(f"Velocity values updated and saved to {output_file_path}")