import sys
import os
import pandas as pd
from shapely.geometry import MultiPoint
import csv

def read_points_data(data_file_path):
    data = pd.read_csv(data_file_path)
    return data

def get_convex_hull_coordinates(points_data):
    cell_groups = points_data.groupby('cell')
    convex_hulls = {}

    for cell_name, group in cell_groups:
        if cell_name != '0':  # Skip points not assigned to any cell
            points = group[['x', 'y']].values
            multi_point = MultiPoint(points)
            convex_hull = multi_point.convex_hull
            if convex_hull.geom_type == 'Polygon':
                coords = list(convex_hull.exterior.coords)
                convex_hulls[cell_name] = coords

    return convex_hulls

def write_roi_coordinates_to_csv(roi_coordinates, output_csv_path):
    with open(output_csv_path, 'w', newline='') as csvfile:
        fieldnames = ['ROI_Name', 'X', 'Y']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for roi_name, coordinates in roi_coordinates.items():
            for x, y in coordinates:
                writer.writerow({'ROI_Name': roi_name, 'X': x, 'Y': y})

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Usage: python points_to_roi.py <points_data_csv> <output_folder>')
        sys.exit(1)

    points_data_csv = sys.argv[1]
    output_folder = sys.argv[2]

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    points_data = read_points_data(points_data_csv)
    convex_hulls = get_convex_hull_coordinates(points_data)

    output_csv_path = os.path.join(output_folder, 'convex_hulls.csv')
    write_roi_coordinates_to_csv(convex_hulls, output_csv_path)
    print(f'Wrote convex hull ROI coordinates to {output_csv_path}')