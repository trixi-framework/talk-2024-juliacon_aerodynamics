# Aim to remove duplicate nodes (due to truncation to 2D)
# Also: Remove some stuff irrelevant for Trixi (e.g. element sets)

path = "../mesh_data/"
prefix = "NACA4412_2"

# Read the file
lines = readlines(path * prefix * "_2D.inp")

### Remove duplicate nodes, i.e, nodes with the same x and y coordinates ###

# Find the start and end of the node list
start = findfirst(x -> x == "*NODE", lines) + 1
_end = findfirst(x -> x == "******* E L E M E N T S *************", lines)

nnodes = _end - start
new_lines = lines[1:start-1]
removed_nodes = []

# Stick with first half of nodes
# CARE: Might not work for all NASA meshes!
for i in start:Int(nnodes/2)+start-1
    push!(new_lines, lines[i])
end
# Remove second half of nodes
for i in Int(nnodes/2)+start:_end-1
    entries = split(lines[i])
    push!(removed_nodes, entries[1])
end

# Convert removed_nodes to an array of integers
removed_nodes = [parse(Int, replace(node, "," => "")) for node in removed_nodes]
#println(removed_nodes)
println("Number removed nodes: ", length(removed_nodes))

# Append the rest of the lines
append!(new_lines, lines[_end:end])

### Remove 2D elements ###

# Initialize an empty array to hold the lines to keep
lines_to_keep = []

# Initialize a flag to indicate whether to start checking the lines
check_lines = false

# Iterate over the lines in the file
for line in new_lines
    # Check if the line starts with "*ELEMENT, type=CPS4,"
    if startswith(line, "*ELEMENT, type=CPS4,")
        # If it does, set the flag to start checking the following lines
        check_lines = true
    elseif check_lines
        # If the flag is set, split the line by comma and check if there are 5 entries
        if length(split(line, ',')) != 5
            # If there are not 5 entries, unset the flag and keep the line
            check_lines = false
            push!(lines_to_keep, line)
        end
    else
        # If the line does not start with "*ELEMENT, type=CPS4," and the flag is not set, keep the line
        push!(lines_to_keep, line)
    end
end

### Turn 3D elements into 2D ###

new_lines = lines_to_keep

# Find the start of the element list
start = findfirst(x -> startswith(x, "*ELEMENT, type=C3D8"), new_lines) + 1

# Modify the element list
for i in start:length(new_lines)
    entries = split(new_lines[i], ',')
    if length(entries) == 9
        new_lines[i] = join(entries[1:5], ",")
    end
end

# Replace "*ELEMENT, type=C3D8," with "*ELEMENT, type=CPS4,"
new_lines = [replace(line, "*ELEMENT, type=C3D8," => "*ELEMENT, type=CPS4,") for line in new_lines]

### Replace underscores in node set names (dashes are not permitted in symbol names in Julia) ###
# Iterate over the lines in new_lines
for i in 1:length(new_lines)
    # Check if the line starts with "*NSET,NSET"
    if startswith(new_lines[i], "*NSET,NSET")
        # If it does, replace every dash "-" with an underscore "_"
        new_lines[i] = replace(new_lines[i], "-" => "_")
    end
end

### Delete irrelevant (for Trixi) elsets ###

# Initialize a boolean variable to track whether we're between "*ELSET" lines
in_elset = false

# Create a new list of lines
new_lines_2 = []

# Iterate over the lines
for line in new_lines
    if startswith(line, "*ELSET")
        in_elset = true
    elseif startswith(line, "*ELSET") || startswith(line, "*NSET")
        in_elset = false
    end

    # If we're not between "*ELSET" lines, add the line to the new list
    if !in_elset
        push!(new_lines_2, line)
    end
end

# Replace new_lines with new_lines_2
new_lines = new_lines_2

# IDEA: Could only pass on the specified nodesets

# Write the new list of lines back to the file, skipping empty lines
open(path * prefix * "_2D_unique.inp", "w") do f
  for line in new_lines
      if !isempty(strip(line))
          write(f, line * "\n")
      end
  end
end