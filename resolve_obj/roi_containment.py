import sys
import os
import csv
from collections import defaultdict
from shapely.geometry import Polygon
from rtree import index

def read_roi_coordinates(csv_file_path):
    roi_coordinates = defaultdict(list)

    with open(csv_file_path, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            roi_name = row['ROI_Name']
            x, y = float(row['X']), float(row['Y'])
            roi_coordinates[roi_name].append((x, y))

    return roi_coordinates


def find_roi_containment(original_rois, new_rois):
    containment = defaultdict(list)

    idx = index.Index()
    original_polygons = {}
    original_roi_names = []

    for original_roi_name, original_coords in original_rois.items():
        if len(original_coords) < 3:
            continue  # Skip if there are less than 3 coordinates
        if original_coords[0] != original_coords[-1]:
            original_coords.append(original_coords[0])  # Add the first coordinate to the end if necessary
        original_polygon = Polygon(original_coords)
        if original_polygon.is_valid:
            original_polygons[original_roi_name] = original_polygon
            original_roi_names.append(original_roi_name)
            idx.insert(len(original_roi_names) - 1, original_polygon.bounds)

    for new_roi_name, new_coords in new_rois.items():
        if len(new_coords) < 3:
            continue  # Skip if there are less than 3 coordinates
        if new_coords[0] != new_coords[-1]:
            new_coords.append(new_coords[0])  # Add the first coordinate to the end if necessary
        new_polygon = Polygon(new_coords)
        if new_polygon.is_valid:
            new_centroid = new_polygon.centroid
            for i in idx.intersection(new_polygon.bounds):
                original_roi_name = original_roi_names[i]
                original_polygon = original_polygons[original_roi_name]
                if original_polygon.contains(new_centroid):
                    containment[original_roi_name].append(new_roi_name)

    return containment


def write_containment_to_csv(containment, output_file_path):
    with open(output_file_path, 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Original_ROI', 'New_ROI'])

        for original_roi, new_rois in containment.items():
            for new_roi in new_rois:
                writer.writerow([original_roi, new_roi])

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print('Usage: python roi_containment.py <original_roi_csv> <new_roi_csv> <output_csv>')
        sys.exit(1)

    original_roi_csv = sys.argv[1]
    new_roi_csv = sys.argv[2]
    output_csv = sys.argv[3]

    original_rois = read_roi_coordinates(original_roi_csv)
    new_rois = read_roi_coordinates(new_roi_csv)

    containment = find_roi_containment(original_rois, new_rois)
    write_containment_to_csv(containment, output_csv)