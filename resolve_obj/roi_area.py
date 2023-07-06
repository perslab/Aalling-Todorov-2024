import sys
import os
import glob
import csv
from collections import defaultdict
from shapely.geometry import Polygon

def compute_area_from_csv(csv_file_path):
    roi_coordinates = defaultdict(list)

    with open(csv_file_path, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            roi_name = row['ROI_Name']
            x, y = float(row['X']), float(row['Y'])
            roi_coordinates[roi_name].append((x, y))

    roi_areas = {}
    for roi_name, coordinates in roi_coordinates.items():
        if len(coordinates) >= 4:
            polygon = Polygon(coordinates)
            area = polygon.area
            roi_areas[roi_name] = area
        else:
            print(f'Warning: Not enough coordinates for {roi_name} to form a valid polygon.')

    return roi_areas

def write_area_to_csv(roi_areas, output_csv_path):
    with open(output_csv_path, 'w', newline='') as csvfile:
        fieldnames = ['ROI_Name', 'Area']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for roi_name, area in roi_areas.items():
            writer.writerow({'ROI_Name': roi_name, 'Area': area})

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Usage: python compute_roi_area.py <csv_folder> <output_folder>')
        sys.exit(1)

    csv_folder = sys.argv[1]
    output_folder = sys.argv[2]

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    csv_files = glob.glob(os.path.join(csv_folder, '*.csv'))

    for csv_file in csv_files:
        roi_areas = compute_area_from_csv(csv_file)
        output_csv_path = os.path.join(output_folder, os.path.splitext(os.path.basename(csv_file))[0] + '_areas.csv')
        write_area_to_csv(roi_areas, output_csv_path)
        print(f'Wrote areas for {os.path.basename(csv_file)} to {output_csv_path}')

