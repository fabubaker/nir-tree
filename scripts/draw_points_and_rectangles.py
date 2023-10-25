import sys
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection

def parse_input_file(file_path):
    points = []
    rectangles = []

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()

            if line.startswith('('):
                # Parse points
                x, y = line.strip('()').split(',')
                points.append((float(x), float(y)))

            elif line.startswith('[('):
                # Parse rectangles
                coordinates = line.strip()[1:-1].split('), ')
                x1, y1 = map(float, coordinates[0][1:].split(', '))
                x2, y2 = map(float, coordinates[1][1:-1].split(', '))
                rectangles.append([(x1, y1), (x2, y2)])

    return points, rectangles

def draw_points_and_rectangles(points, rectangles):
    fig, ax = plt.subplots()

    # Draw points using PathCollection
    x_points = [p[0] for p in points]
    y_points = [p[1] for p in points]
    ax.scatter(x_points, y_points, color='red', s=0.1)

    # Draw rectangles
    for rectangle in rectangles:
        x_min, y_min = rectangle[0]
        x_max, y_max = rectangle[1]
        width = x_max - x_min
        height = y_max - y_min
        rect = patches.Rectangle(
            (x_min, y_min), width, height,
            linewidth=0.25, edgecolor='black', facecolor='blue', alpha=0.2
        )
        ax.add_patch(rect)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title('Points and Rectangles')
    ax.grid(False)

    # Save the plot as SVG
    plt.show()

# Check if an input file path and output file path are provided as command line arguments
if len(sys.argv) > 1:
    file_path = sys.argv[1]
    points, rectangles = parse_input_file(file_path)
    draw_points_and_rectangles(points, rectangles)
else:
    print("Invalid arguments")
