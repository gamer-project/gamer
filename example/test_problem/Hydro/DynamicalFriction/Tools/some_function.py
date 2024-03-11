def read_selected_parameters(filename):
    '''Read selected parameters from Input__Parameter and return a dictionary'''
    parameters = {}
    with open(filename, 'r') as file:
        for line in file:
            # Remove comments and strip whitespace
            line = line.split('#')[0].strip()
            
            # Skip empty lines
            if not line:
                continue
            
            words = line.split()
            if len(words) != 2:
                continue

            key, value = words

            # Check if the key is one of the desired parameters
            if key in ["BOX_SIZE_X", "NX0_TOT[0]", "END_T"]:
                # Convert the value to float if possible
                try:
                    value = float(value)
                except ValueError:
                    pass
                
                parameters[key] = value

    return parameters


def read_parameters(filename):
    '''Read the parameters from Input__TestProb and return a dictionary'''
    parameters = {}
    with open(filename, 'r') as file:
        for line in file:
            # Remove comments and strip whitespace
            line = line.split('#')[0].strip()
            
            # Skip empty lines
            if not line:
                continue
            
            words = line.split()
            if len(words) == 1:
                # Skip lines with only one word or assign a default value
                print(f"Warning: Skipping line with single word in {filename}: {line}")
                continue
            elif len(words) != 2:
                print(f"Warning: Unexpected line format in {filename}: {line}")
                continue

            key, value = words
            
            # Convert the value to float if possible
            try:
                value = float(value)
            except ValueError:
                pass
            
            parameters[key] = value
    return parameters