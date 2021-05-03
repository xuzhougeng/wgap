import json
from sys import argv
from wgap.scripts.maker_ctl_io import ctl_reader
from wgap.scripts.maker_ctl_io import ctl_writer

input_fn = argv[1]
out_fn = argv[2]

params_dict = ctl_reader(input_fn)
#print("Reading" +  input_fn ,file=stderr)

json_file = argv[3]
with open(json_file, "r") as f:
    params = json.load(f)

# modify the params
if len(params) > 1:
    for key,value in params.items():
        params_dict[key] = value

# write
ctl_writer(params_dict, out_fn)