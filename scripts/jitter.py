# Use this script to add jitter to points in a csv file.
# By default, both the original and the jittered points are kept in the output file,
# effectively doubling the input.
# Comment out the 'NO_RETAIN' variable below to not retain the original points.
# Might need to run `sort -u` afterwards.

NO_RETAIN = False

USAGE='USAGE: python3 jitter.py <filename.csv>'

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
    jittered_row = [None, None]

    x = float(input_row[0])
    y = float(input_row[1])
    x += float(random.randint(1,99))/(10**4)
    y += float(random.randint(1,99))/(10**4)
    jittered_row[0] = "{:f}".format(x)
    jittered_row[1] = "{:f}".format(y)

    if not NO_RETAIN:
        csv_writer.writerow(input_row)

    csv_writer.writerow(jittered_row)

input_csv_file.close()
output_csv_file.close()
