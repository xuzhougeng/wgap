from collections import OrderedDict

# read the control file
def ctl_reader(file_name):
    ctl_dict = OrderedDict()
    with open(file_name, "r") as f:
        for line in f.readlines():
            if (line.startswith("#")) or line.startswith("\n"):
                continue
            line = line[0:line.find("#")]
            k,v = line.split("=")
            ctl_dict[k] = v
    return ctl_dict



# write the control file
def ctl_writer(ctl_dict, file_name):
    #ctl_dict = OrderedDict()
    with open(file_name, "w") as f:
        for k,v in ctl_dict.items():
            line = "{}={}\n".format(str(k), str(v))
            f.writelines(line)

if __name__ == "__main__":
    pass