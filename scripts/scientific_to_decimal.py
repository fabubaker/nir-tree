# Use this script to convert floats in scientific notation to decimal.

USAGE=\
'''USAGE: python3 scientific_to_decimal.py <filename.csv>'''

import sys
import csv

if len(sys.argv) != 2:
    print(USAGE)

filename = sys.argv[1]
input_csv_file = open(filename)
output_csv_file = open(filename + "_output.csv", "w")

csv_reader = csv.reader(input_csv_file, delimiter=',')
csv_writer = csv.writer(output_csv_file, delimiter=',')

for input_row in csv_reader:
    output_row = []

    if '' in input_row:
        continue

    print("Converting row:", input_row)
    for column in input_row:
        column_float = "{:f}".format(float(column))
        output_row.append(column_float)

    csv_writer.writerow(output_row)

input_csv_file.close()
output_csv_file.close()