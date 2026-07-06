import math

# Read the current clashing coordinate file
with open("ionized.gro", "r") as f:
    lines = f.readlines()

header = lines[:2]
footer = [lines[-1]]
atom_lines = lines[2:-1]

new_atom_lines = []
nitrogens = {}

# Pass 1: Map out the correct position of every Nitrogen
for line in atom_lines:
    if len(line) < 44: continue
    resnum = line[:5].strip()
    atomname = line[10:15].strip()
    try:
        x = float(line[20:28])
        y = float(line[28:36])
        z = float(line[36:44])
        if atomname == "N":
            nitrogens[resnum] = (x, y, z)
    except ValueError:
        pass

# Pass 2: Find runaway amide hydrogens and snap them back to 0.1 nm
fixed_count = 0
for line in atom_lines:
    if len(line) < 44:
        new_atom_lines.append(line)
        continue
    resnum = line[:5].strip()
    atomname = line[10:15].strip()
    
    if atomname == "H" and resnum in nitrogens:
        try:
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
            nx, ny, nz = nitrogens[resnum]
            
            # Calculate distance between N and H
            dx, dy, dz = x - nx, y - ny, z - nz
            dist = math.sqrt(dx**2 + dy**2 + dz**2)
            
            # If the hydrogen bond is unnaturally long, rebuild it cleanly
            if dist > 0.15: 
                if dist > 0:
                    new_x = nx + (dx / dist) * 0.100
                    new_y = ny + (dy / dist) * 0.100
                    new_z = nz + (dz / dist) * 0.100
                else:
                    new_x, new_y, new_z = nx + 0.1, ny, nz
                    
                # Reconstruct the text line keeping strict .gro formatting columns
                line = line[:20] + f"{new_x:8.3f}{new_y:8.3f}{new_z:8.3f}" + line[44:]
                fixed_count += 1
        except ValueError:
            pass
    new_atom_lines.append(line)

# Write the fixed structure out
with open("ionized_fixed.gro", "w") as f:
    f.writelines(header + new_atom_lines + footer)

print(f"Success! Cleaned up {fixed_count} runaway hydrogen atoms across the protein backbone.")