import matplotlib.pyplot as plt
import matplotlib.patches as patches

def parse_input_file(file_path):
    rectangles = []

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()

            if line.startswith('[('):
                coordinates = line.strip()[1:-1].split('), ')
                x1, y1 = map(float, coordinates[0][1:].split(', '))
                x2, y2 = map(float, coordinates[1][1:-1].split(', '))
                rectangles.append(([(x1, y1), (x2, y2)]))

    return rectangles

def draw_rectangles(rectangles):
    fig, ax = plt.subplots()

    depth_color_map = {}  # Map depth to color

    for rectangle in rectangles:
        x_min, y_min = rectangle[0]
        x_max, y_max = rectangle[1]
        width = x_max - x_min
        height = y_max - y_min
        rect = patches.Rectangle((x_min, y_min), width, height, linewidth=1, edgecolor='black', facecolor='blue', alpha=0.2)
        ax.add_patch(rect)

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Rectangles')
    plt.grid(visible=False)
    plt.show()

input_file = '/home/fabubaker/Desktop/printed_tree.txt'
rectangles = parse_input_file(input_file)
draw_rectangles(rectangles)