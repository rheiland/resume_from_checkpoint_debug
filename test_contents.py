import scipy.io
info_dict = {}
full_fname = "output00000108_cells.mat"
full_fname = "output00000109_cells.mat"
full_fname = "fixed.mat"
full_fname = "output00000856_cells.mat"
full_fname = "output00000857_cells.mat"
print("--- ",full_fname)
scipy.io.loadmat(full_fname, info_dict)
keys = info_dict.keys()
# print("keys= ",keys)
num_keys = len(keys)
print("# keys in info_dict = ",num_keys)

keys_list = list(info_dict.keys())
#print(keys_list)
print(keys_list[0])
# print(keys_list[1])
#print(keys_list[2])


#print(info_dict.keys())
cells = info_dict['cells']

print("cells.shape = ",cells.shape)
print("\n-----len(cells) = ",len(cells))
print(cells)
#Out[10]: (119, 20757)
