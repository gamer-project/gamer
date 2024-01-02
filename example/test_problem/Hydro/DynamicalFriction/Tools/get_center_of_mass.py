import yt

# Define a range of file numbers
file_numbers = range(0, 17)  # This assumes your file names follow the pattern Data_000000, Data_000001, ..., Data_000016

# Create a list to store the center of mass values
center_of_mass_list = []

# Loop through the file numbers
for file_number in file_numbers:
    # Generate the filename for the current dataset
    file = f'../Data_{file_number:06d}'
    
    # Load data
    ds = yt.load(file)
    ad = ds.all_data()

    # Filter data if required
    #dense = ad.cut_region(['obj["Dens"] > 1.0e2'])

    # Calculate center-of-mass
    cm = ad.quantities.weighted_average_quantity(['particle_position_x', 'particle_position_y', 'particle_position_z'], weight='particle_mass')

    # Append the center of mass to the list
    center_of_mass_list.append(cm)

# Define the output file name
output_file = 'center_of_mass.txt'

# Open the file in append mode and write the results
with open(output_file, 'a') as f:
    f.write("Data_ID\tCoM_X\tCoM_Y\tCoM_Z\n")
    for file_number, cm in zip(file_numbers, center_of_mass_list):
        f.write(f'Data_{file_number:06d}\t {cm[0].in_units("code_length").value}\t {cm[1].in_units("code_length").value} \t {cm[2].in_units("code_length").value} \n')

# Print a message to confirm that the data has been logged
print(f'Center of mass data has been logged to {output_file}')

