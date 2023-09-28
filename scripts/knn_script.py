import sys
import numpy as np
from sklearn.neighbors import NearestNeighbors
import math

def load_points(file_path):
    points = []
    with open(file_path, 'r') as file:
        for line in file:
            point = list(map(float, line.strip().split()))
            points.append(point)
    return np.array(points)

def calculate_mbr(neighbors):
    min_x = np.min(neighbors[:, 0])
    max_x = np.max(neighbors[:, 0])
    min_y = np.min(neighbors[:, 1])
    max_y = np.max(neighbors[:, 1])
    return min_x, max_x, min_y, max_y

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Generates rectangles containing a specific number of points using a knn algorithm.")
        print("The points themselves are obtained from an input file.")
        print("")
        print("Usage: python knn_script.py <input_file_path> <k> <n> <output_file_path>")
        print("<input_file_path> - Path to input file containing points to knn")
        print("<k> - How many points to collect in a rectangle")
        print("<n> - The total number of rectangles containing k points to generate")
        print("<output_file_path> - The output file to save rectangles to")
        sys.exit(1)

    input_file_path = sys.argv[1]
    k = int(sys.argv[2])
    n = int(sys.argv[3])
    output_file_path = sys.argv[4]

    output_file = open(output_file_path, 'w')

    # Load points from the input file
    print("Loading points from the input file...")
    points = load_points(input_file_path)

    # Create a NearestNeighbors model
    print("Creating a NearestNeighbors model...")
    nbrs = NearestNeighbors(n_neighbors=k, algorithm='ball_tree').fit(points)

    # Load query points from the input file
    query_points = []
    print(f"Loading {n} query points from the input file...")
    skip = math.ceil(len(points)/n)

    for i in range(0, len(points), skip):
        query_points.append(points[i])

    for idx, query_point in enumerate(query_points):
        print(f"\nQuery Point {idx+1}: {query_point}")
        # Find the k Nearest Neighbors of the query point
        distances, indices = nbrs.kneighbors([query_point])

        # Get the k nearest neighbors
        neighbors = points[indices[0]]

        # Calculate the minimum bounding rectangle
        min_x, max_x, min_y, max_y = calculate_mbr(neighbors)

        print(f"The {k} Nearest Neighbors of the point {query_point} are:")
        for i, index in enumerate(indices[0]):
            print(f"Neighbor {i+1}: {points[index]} (Distance: {distances[0][i]})")

        print(f"MBR for query point {idx+1} (min X, min Y, max X, max Y) | {min_x} {min_y} {max_x} {max_y}")
        print()  # Add an empty line between query points

        output_file.write(f"{min_x} {min_y} {max_x} {max_y}\n")

    print(f"{n} MBRs with {k} points written to {output_file_path}")
    output_file.close()
