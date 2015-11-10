import os
import subprocess
import re


def mv_rename(old, new):
    try:
        subprocess.check_call(["mv",
                               old,
                               new])
    except:
        log(old + " name not changed")

    print(old + " renamed to " + new)


# error logging
def log(message):
    with open(os.path.join('./', 'renaming_log.txt'), 'a') as datafile:
        datafile.write(message + "\n")
        datafile.close()


# Takes in tab-separated list of:
# 0: old sample number
# 1: new sample number
# 2: old sample name
# 3: new sample name
def read_dicts():
    sample_names = dict()
    sample_nums = dict()

    with open(os.path.join('./', 'sample_data.txt'), 'r') as datafile:
        for line in datafile:
            line_list = line.split("\t")
            sample_nums[line_list[0]] = line_list[1]
            sample_names[line_list[2]] = line_list[3].strip("\n")
        datafile.close()

    return (sample_names, sample_nums)


def main():

    (sample_names, sample_nums) = read_dicts()

    if sample_names and sample_nums:
        print("dictionary populated.")

    root = os.getcwd() + "/BaseSpace/cranio1-26049053/"

    for dirname in os.listdir(root):
        # To check it is dir, not file
        if os.path.isdir(root + dirname):
            bcpath = root + dirname + "/Data/Intensities/BaseCalls/"
            for filename in os.listdir(bcpath):
                # get sample number between S1-96
                num_re = re.match(r"\S*_S(\d+)\S*",
                                  filename)

                if num_re:
                    old_num = num_re.group(1)
                    try:
                        new_num = sample_nums.get(old_num)
                        new_filename = re.sub(r"_S(\d+)",
                                              "_S" + str(new_num),
                                              filename)
                        mv_rename(bcpath + filename,
                                  bcpath + new_filename)
                    except KeyError:
                        log(old_num + " not found in dict!")
                else:
                    if (".sam" in filename or
                            ".bam" in filename or
                            ".fastq" in filename):
                        mv_rename(bcpath + filename,
                                  bcpath + filename + "_old_kept")

                if filename is ".DS_Store":
                    subprocess.check_call(["rm", filename])

            name_re = re.match(r"([\S*\s*]*)-300\S*",
                               dirname)

            if name_re:
                old_name = name_re.group(1)
                try:
                    new_name = sample_names.get(old_name)
                    new_dirname = re.sub(r"([\S*\s*]*)(-300\S*)",
                                         new_name + r"\g<2>",
                                         dirname)
                    mv_rename(root + dirname,
                              root + new_dirname)
                except KeyError:
                    log(old_name + " not found in dict!")
            else:
                mv_rename(root + dirname,
                          root + dirname + "_old_kept")

        if dirname is ".DS_Store":
            subprocess.check_call(["rm", dirname])

if __name__ == "__main__":
    main()
