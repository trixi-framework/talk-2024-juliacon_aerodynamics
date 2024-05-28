path = "../mesh_data/"
prefix = "NACA4412_2"

### Assign the geometrical ID to the elements in the mesh ###
# The conversion program "p3d2gmsh" does this "wrong"

file = open(path * prefix * ".msh", "r")
lines = readlines(file)
close(file)

# Open the output file in write mode
out = open(path * prefix * "_boundary_assoc.msh", "w")

# Flag to indicate if "$Elements" line has been read
elements_found = false

# Iterate over the lines
for line in lines
    if line == "\$Elements"
        elements_found = true
    end
    # If "$Elements" line has been read, start replacing
    if elements_found
        # Split the line into words
        words = split(line)
        # Check if there are at least 5 words (sanity check if this is really an element)
        if length(words) >= 5
            # Replace the fifth word with the fourth one to pass on the physical ID
            words[5] = words[4]
        end
        # Write the modified line to the file
        write(out, join(words, ' ') * "\n")
    else
        # If "$Elements" line has not been read, write the original line
        write(out, line * "\n")
    end
end

# Close the file
close(out)

mv(path * prefix * "_boundary_assoc.msh", path * prefix * ".msh", force=true)

### Swap y and z coordinates ###

# Read the file
lines = readlines(path * prefix * ".msh")

### Swap y and z coordinates ###

# Find the start and end of the node list
start = findfirst(x -> x == "\$Nodes", lines) + 2
_end = findfirst(x -> x == "\$EndNodes", lines)

# Swap the third and fourth entries in each line
for i in start:_end-1
    entries = split(lines[i])
    if length(entries) >= 4
        entries[3], entries[4] = entries[4], entries[3]
        lines[i] = join(entries, ' ')
    end
end

# Write the file
open(path * prefix * "_swapped.msh", "w") do f
  for line in lines
      write(f, line * "\n")
  end
end

# Overwrite original file
mv(path * prefix * "_swapped.msh", path * prefix * ".msh", force=true)

### Set z coordinates to zero ###

for i in start:_end-1
  entries = split(lines[i])
  if length(entries) >= 4
      # Truncate to 2D
      entries[4] = "0.0"
      lines[i] = join(entries, ' ')
  end
end

# Write the file
open(path * prefix * "_2D.msh", "w") do f
  for line in lines
      write(f, line * "\n")
  end
end