import sys
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
                rectangles.append([(x1, y1), (x2, y2)])

    return rectangles

def check_overlap(rect1, rect2):
    # Check if two rectangles overlap
    x1_min, y1_min = rect1[0]
    x1_max, y1_max = rect1[1]
    x2_min, y2_min = rect2[0]
    x2_max, y2_max = rect2[1]

    if (x1_min <= x2_max and x1_max >= x2_min) and (y1_min <= y2_max and y1_max >= y2_min):
        return True
    else:
        return False

def find_overlapping_rectangles(rectangles):
    overlapping_rectangles = []
    for i in range(len(rectangles)):
        for j in range(i + 1, len(rectangles)):
            if check_overlap(rectangles[i], rectangles[j]):
                overlapping_rectangles.append(rectangles[i])
                overlapping_rectangles.append(rectangles[j])

    return overlapping_rectangles

# Check if an input file path is provided as a command line argument
if len(sys.argv) > 1:
    input_file = sys.argv[1]
    rectangles = parse_input_file(input_file)
    print(find_overlapping_rectangles(rectangles))
else:
    print("Please provide the path to the input file as a command line argument.")
