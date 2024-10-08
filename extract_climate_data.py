import json
from scipy.spatial import cKDTree


def find_nearest_coordinates(field_coordinates, json_data):
    """
    Find the nearest coordinates in the JSON data to a given field.

    This function uses a KDTree to find the nearest neighbor to the given field.

    Parameters:
        field_coordinates (tuple): a tuple of (latitude, longitude) for the field
        json_data (list): a list of entries, where each entry is a list with the first
            element being a tuple of (latitude, longitude)

    Returns:
        tuple: A tuple of (latitude, longitude) for the nearest coordinates in the JSON file.
    """
    # Extract coordinates from json_data
    coordinates = [entry[0] for entry in json_data]

    # Build a KDTree from the coordinates
    tree = cKDTree(coordinates)

    # Query the tree to find the nearest neighbor
    dist, idx = tree.query(field_coordinates, k=1)

    # Retrieve the nearest coordinates
    nearest_coordinates = coordinates[idx]

    return nearest_coordinates


# Load JSON file containing mapping of coordinates to row and column numbers
# Note: Make sure that the JSON file corresponds to the climate data you want to extract
json_file = 'latlon-to-rowcol.json'
with open(json_file, 'r') as f:
    latlon_to_rowcol = json.load(f)

import pandas as pd

df = pd.read_csv('2018_climaestation.csv')
X = df['POINT_X']
Y = df['POINT_Y']
output = []
for i,j in zip(X,Y):
    # Field coordinates
    field_coordinates = (i,j)  # Latitude and longitude of field
    print(field_coordinates)
    # Multiple fields coordinates
    # field_coordinates = [
    #     (52.959811, 9.465853),
    #     (52.5200, 13.4050)
    # ]

    # Find nearest coordinates in the JSON file
    nearest_coordinates = find_nearest_coordinates(field_coordinates, latlon_to_rowcol)
    if nearest_coordinates:
        for entry in latlon_to_rowcol:
            coordinates = entry[0]
            if coordinates == nearest_coordinates:
                row_col_numbers = entry[1]
                break
        print(f'Field Coordinates: {field_coordinates}')
        print(f'Nearest coordinates: {nearest_coordinates}')
        print(f'Row and column numbers: {row_col_numbers}')
    else:
        print('No coordinates found.')
    output.append(row_col_numbers)
df['station'] = output
df.to_csv('soil_station.csv',index=False)
# # Field coordinates
# field_coordinates = (53.06243, 11.75081)  # Latitude and longitude of field
#
#
# # Multiple fields coordinates
# # field_coordinates = [
# #     (52.959811, 9.465853),
# #     (52.5200, 13.4050)
# # ]
#
# # Find nearest coordinates in the JSON file
# nearest_coordinates = find_nearest_coordinates(field_coordinates, latlon_to_rowcol)
# if nearest_coordinates:
#     for entry in latlon_to_rowcol:
#         coordinates = entry[0]
#         if coordinates == nearest_coordinates:
#             row_col_numbers = entry[1]
#             break
#     print(f'Field Coordinates: {field_coordinates}')
#     print(f'Nearest coordinates: {nearest_coordinates}')
#     print(f'Row and column numbers: {row_col_numbers}')
# else:
#     print('No coordinates found.')

# For multiple fields coordinates
# for field_coord in field_coordinates:
#     nearest_coordinates = find_nearest_coordinates(field_coord, latlon_to_rowcol)
#     if nearest_coordinates:
#         for entry in latlon_to_rowcol:
#             coordinates = entry[0]
#             if coordinates == nearest_coordinates:
#                 row_col_numbers = entry[1]
#                 break
#         print(f'Field Coordinates: {field_coord}')
#         print(f'Nearest coordinates: {nearest_coordinates}')
#         print(f'Row and column numbers: {row_col_numbers}')
#         print("\n")
#     else:
#         print('No coordinates found.')
