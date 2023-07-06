import argparse
import pandas as pd
from shapely.geometry import Polygon, Point
from rtree import index
from tqdm import tqdm


def load_nucleus_boundaries(file_path):
    """Load nucleus boundaries from a CSV file and create Shapely Polygon objects for each nucleus."""
    df = pd.read_csv(file_path)
    nuclei = {}
    for name, group in df.groupby('ROI_Name'):
        vertices = [(x, y) for x, y in zip(group['X'], group['Y'])]
        # Add the last vertex identical to the first to make the polygon closed
        vertices.append(vertices[0])
        nuclei[name] = Polygon(vertices)
    return nuclei


def load_transcripts(file_path):
    """Load transcript coordinates from a CSV file and create Shapely Point objects for each transcript."""
    df = pd.read_csv(file_path)
    transcripts = [Point(x, y) for x, y in zip(df['x'], df['y'])]
    return transcripts


def create_rtree(nuclei):
    """Create an R-tree spatial index for the nucleus boundaries."""
    idx = index.Index()
    for i, nucleus in enumerate(nuclei.values()):
        idx.insert(i, nucleus.bounds, obj=nucleus)
    return idx

#check if points are in any polygons
# ~15min per run
def main(convex_hulls_file, baysor_results_file, output_file):
    """Main script function."""
    print('Loading nucleus boundaries...')
    nuclei = load_nucleus_boundaries(convex_hulls_file)
    print(f'Loaded {len(nuclei)} nuclei.')
    print('Loading transcripts...')
    transcripts = load_transcripts(baysor_results_file)
    print(f'Loaded {len(transcripts)} transcripts.')
    # Create R-tree spatial index for nuclei
    print('Creating R-tree index...')
    idx = create_rtree(nuclei)
    # Check if each transcript is within a nucleus
    print('Checking transcripts...')
    df = pd.read_csv(baysor_results_file)
    results = []
    for i, transcript in tqdm(enumerate(transcripts), total=len(transcripts)):
        # Check if transcript is within any nucleus
        within_nucleus = False
        for j in idx.intersection((transcript.x, transcript.y)):
            if transcript.within(nuclei[list(nuclei.keys())[j]]):
                within_nucleus = True
                break
        results.append(within_nucleus)
    # Write results to output file
    df['within_nucleus'] = results
    df.to_csv(output_file, index=False)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compute whether or not a given transcript is within a given nucleus.')
    parser.add_argument('convex_hulls_file', help='Path to the convex_hulls file')
    parser.add_argument('baysor_results_file', help='Path to the baysor_results file')
    parser.add_argument('output_file', help='Path to the output file')
    args = parser.parse_args()
    main(args.convex_hulls_file, args.baysor_results_file, args.output_file)
