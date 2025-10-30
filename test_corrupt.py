import scipy.io
info_dict = {}
full_fname = "output00000108_cells.mat"
full_fname = "fixed.mat"
full_fname = "output00000109_cells.mat"
print("--- ",full_fname)
scipy.io.loadmat(full_fname, info_dict)
keys = info_dict.keys()
# print("keys= ",keys)
num_keys = len(keys)
print("# keys in info_dict = ",num_keys)

keys_list = list(info_dict.keys())

for idx in range(1,num_keys):
    print(" ---> deleting corrupt index ",idx)
    key_to_remove = keys_list[idx]
    del info_dict[key_to_remove]
#    removed_value = info_dict.pop(key_to_remove)

#    del info_dict[keys[0]]
#    del info_dict[keys[1]]


if num_keys > 1:
    scipy.io.savemat('fixed.mat', info_dict, format='4')


keys = info_dict.keys()
num_keys = len(keys)
print("# keys in info_dict = ",num_keys)
#print(info_dict)

#print(info_dict.keys())
cells = info_dict['cells']

print(cells.shape)
#Out[10]: (119, 20757)
