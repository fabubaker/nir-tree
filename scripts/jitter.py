# Use this script to add jitter to points in a csv file.
# Might need to run `sort -u` afterwards.

USAGE=\
'USAGE: python3 jitter.py <filename.csv>'

import os
import sys
import csv
import random

if len(sys.argv) != 2:
    print(USAGE)

filename = sys.argv[1]
input_csv_file = open(filename)
file_name_without_extension = os.path.splitext(filename)[0]
output_csv_file = open(file_name_without_extension + "_jitter.csv", "w")

csv_reader = csv.reader(input_csv_file, delimiter=' ')
csv_writer = csv.writer(output_csv_file, delimiter=' ')

for input_row in csv_reader:
    # print("Add jitter to row:", input_row)

    x = float(input_row[0])
    y = float(input_row[1])
    x += float(random.randint(0,99))/(10**4)
    y += float(random.randint(0,99))/(10**4)
    input_row[0] = "{:f}".format(x)
    input_row[1] = "{:f}".format(y)

    csv_writer.writerow(input_row)

input_csv_file.close()
output_csv_file.close()
