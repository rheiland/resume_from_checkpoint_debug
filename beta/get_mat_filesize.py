# Validate(?) *_cells*.mat file sizes 
#
# python beta/get_mat_filesize.py output


import os
import sys
from pathlib import Path
import glob

# print(len(sys.argv))
if (len(sys.argv) < 2):
  usage_str = "Usage: %s <outputdir> " % (sys.argv[0])
  print(usage_str)
#  print("e.g.:  Needs <outputdir>")
  exit(1)
else:
   dir1 = sys.argv[1]


# svg_files = Path(dir1).glob(f'{dir1}/out_heterog/snap*.svg')
#print("svg_files=",svg_files)


mat_files = glob.glob(f'{dir1}/output*_cells.mat')
# f"{number_two:03d}"
mat_files.sort()
# print(mat_files)

# cells_mat_size= ((num_cells * num_data_entries) + 25 ) * 8
hdr_size = 25.0
num_cells = 2200
num_data_entries = 125

num_data_entries = 119  # Raquel's 3D model
# (base) M1P~/git/resume_from_checkpoint$ python beta/get_mat_filesize.py .
# The size of './output00000856_cells.mat' is: 16714289 bytes
# num_cells =  17557.0
# file_size_diff =  0.0
# The size of './output00000857_cells.mat' is: 16796161 bytes
# num_cells =  17643.0
# file_size_diff =  0.0

for filename in mat_files:
    f = os.path.basename(filename)
    file_path = os.path.join(dir1,f)
    try:
        cells_mat_size = os.path.getsize(file_path)
        print(f"The size of '{file_path}' is: {cells_mat_size} bytes")
    except OSError as e:
        print(f"Error accessing file '{file_path}': {e}")

#    num_data_entries = (float(cells_mat_size) - 25.0) / / 8
    # num_cells = (float(cells_mat_size) - 25.0) / num_data_entries / 8
    num_cells = (float(cells_mat_size) - 25.0) / num_data_entries / 8
    print("num_cells = ",num_cells)
    file_size_diff = num_cells - int(num_cells)
    print("file_size_diff = ",file_size_diff)
    if abs(file_size_diff) > 1.e-4:
        print("Warning: possible corrupt file")
        break

print("\n --- using num_data_entries = ", num_data_entries)

# from /dev/Raquel/PhysiCell-1.14.2/main.cpp
#                         float num_cells = (float(cells_mat_size) - 25.0) / num_data_entries / 8;
#                         // num_cells += 0.42; // do fake test to fail

#                         float file_size_diff = num_cells - int(num_cells);