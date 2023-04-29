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

            if line.startswith('('):
                # Parse points
                x, y = line.strip('()').split(',')
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
    ax.scatter(x_points, y_points, color='red')

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Points and Rectangles')
    plt.grid(True)
    plt.show()

# Input file path
file_path = '/home/fabubaker/Desktop/printed_tree.txt'

# Parse points and rectangles from input file
points, rectangles = parse_input_file(file_path)

# Draw points and rectangles
draw_points_and_rectangles(points, rectangles)
