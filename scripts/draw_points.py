import sys
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection

def parse_input_file(file_path):
    points = []
    rectangles = []

    count = 0
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()

            # Parse points
            x, y = line.split('\t')
            points.append((float(x), float(y)))

            count += 1
            if count % 1000000 == 0:
                print(count)

    return points, rectangles

def draw_points_and_rectangles(points, rectangles):
    fig, ax = plt.subplots()

    # Draw points using PathCollection
    x_points = [p[0] for p in points]
    y_points = [p[1] for p in points]
    ax.scatter(x_points, y_points, color='red', s=0.1)

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Points and Rectangles')
    plt.grid(False)
    plt.show()

# Check if an input file path is provided as a command line argument
if len(sys.argv) > 1:
    file_path = sys.argv[1]
    points, rectangles = parse_input_file(file_path)
    draw_points_and_rectangles(points, rectangles)
else:
    print("Please provide the path to the input file as a command line argument.")
