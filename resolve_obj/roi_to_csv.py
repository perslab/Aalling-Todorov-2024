import sys
import os
import glob
import read_roi
import csv

def get_roi_info(roi):
    if roi['type'] == 'rectangle':
        x, y, w, h = roi['left'], roi['top'], roi['width'], roi['height']
    elif roi['type'] == 'freehand':
        x, y = roi['x'], roi['y']
        w, h = max(roi['x']) - min(roi['x']), max(roi['y']) - min(roi['y'])
    else:
        raise ValueError(f'Unsupported ROI type: {roi["type"]}')

    return x, y, w, h

def convert_roi_zip_to_csv(zip_file_path, output_csv_path):
    rois = read_roi.read_roi_zip(zip_file_path)

    with open(output_csv_path, 'w', newline='') as csvfile:
        fieldnames = ['ROI_Name', 'X', 'Y']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for roi_name, roi in rois.items():
            if roi['type'] == 'rectangle':
                x, y, w, h = roi['left'], roi['top'], roi['width'], roi['height']
                writer.writerow({'ROI_Name': roi_name, 'X': x, 'Y': y})
            elif roi['type'] == 'freehand':
                x, y = roi['x'], roi['y']
                for i in range(len(x)):
                    writer.writerow({'ROI_Name': roi_name, 'X': x[i], 'Y': y[i]})
            else:
                raise ValueError(f'Unsupported ROI type: {roi["type"]}')


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Usage: python convert_roi_to_csv.py <input_folder> <output_folder>')
        sys.exit(1)

    input_folder = sys.argv[1]
    output_folder = sys.argv[2]

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    roi_zip_files = glob.glob(os.path.join(input_folder, '*.zip'))

    for roi_zip_file in roi_zip_files:
        output_csv_path = os.path.join(output_folder, os.path.splitext(os.path.basename(roi_zip_file))[0] + '.csv')
        convert_roi_zip_to_csv(roi_zip_file, output_csv_path)
        print(f'Converted {roi_zip_file} to {output_csv_path}')
