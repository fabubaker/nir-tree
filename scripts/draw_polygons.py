import matplotlib.pyplot as plt
import matplotlib.patches as patches
import random

# Draws plygons of children on a branch where each polygon can be assigned with a colour 


# Function to draw rectangles
def draw_rectangles(polygon_list):
    fig, ax = plt.subplots()

    for polygon in polygon_list:
        polygon_points, color = polygon
        for rectangle in polygon_points:
            bottom_left = rectangle[0]
            upper_right = rectangle[1]
            width = upper_right[0] - bottom_left[0]
            height = upper_right[1] - bottom_left[1]
            rect = patches.Rectangle(bottom_left, width, height, linewidth=2, edgecolor='black', facecolor=color)
            ax.add_patch(rect)
            ax.text(bottom_left[0], bottom_left[1], f'({bottom_left[0]:.3f},{bottom_left[1]:.3f})', fontsize=10, ha='right', va='bottom',color=color)
            ax.text(upper_right[0], upper_right[1], f'({upper_right[0]:.3f},{upper_right[1]:.3f})', fontsize=10, ha='right', va='bottom',color=color)

    ax.set_aspect('equal')
    ax.autoscale_view()

    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.title('Rectangles with Different Colors')
    plt.grid(True)
    plt.show()

# Example: List of rectangles with polygons and colours 
polygon_list = [
([[(0.84292042681106882362, 0.0017563449269028725587), (0.9164790701528610084, 0.9128197253640669695)],
[(0.9164790701528610084, 0.0017563449269028725587), (0.91685910160056716123, 0.44807536366036226916)],
[(0.91685910160056716123, 0.0017563449269028725587), (0.93422081068905993284, 0.9128197253640669695)]], 'red'),
([[(0.43247836698335284655, 0.020140065567065799229), (0.84273791585997970266, 0.99331711050535242968)]], 'blue'),

([[(0.012892368997317941309, 0.00022219234585404366291), (0.42599086762421850549, 0.42564185382859143214)]], 'green'),

([[(0.0075342274300282570387, 0.44807536366036226916), (0.43247836698335284655, 0.99177730796751606412)],
[(0.84273791585997970266, 0.44807536366036226916), (0.84292042681106882362, 0.99177730796751606412)],
[(0.9164790701528610084, 0.44807536366036226916), (0.91685910160056716123, 0.99177730796751606412)],
[(0.84292042681106882362, 0.9128197253640669695), (0.9164790701528610084, 0.99177730796751606412)]], 'orange'),
]
polygon_list1 = [
([[(0.84292042681106882362, 0.0017563449269028725587), (0.9164790701528610084, 0.9128197253640669695)]], 'red'),
([[(0.43247836698335284655, 0.020140065567065799229), (0.84273791585997970266, 0.99331711050535242968)]],'blue'),

([[(0.012892368997317941309, 0.00022219234585404366291), (0.42599086762421850549, 0.42564185382859143214)]], 'green'),
([[(0.0075342274300282570387, 0.44807536366036226916), (0.43247836698335284655, 0.99177730796751606412)],
[(0.84273791585997970266, 0.44807536366036226916), (0.84292042681106882362, 0.99177730796751606412)],
[(0.9164790701528610084, 0.44807536366036226916), (0.91685910160056716123, 0.99177730796751606412)],
[(0.84292042681106882362, 0.9128197253640669695), (0.9164790701528610084, 0.99177730796751606412)]],'orange'),
]
draw_rectangles(polygon_list)
