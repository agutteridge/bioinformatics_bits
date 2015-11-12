###############################################################################
# For when your barcodes correspond to the wrong samples!
#
# Make sure this is executed from the run directory.
#
# Ensure the input file is named "sample_data.txt" and is in the same
# directory, with 4 tab-separated columns corresponding to:
# 0: old sample number
# 1: new sample number
# 2: old sample name
# 3: new sample name

import os
import subprocess
import re
import sys

import dir_tools


def mv_rename(old, new):
    try:
        subprocess.check_call(["mv",
                               old,
                               new])
    except:
        log(old + " name not changed")

    log(old + " renamed to " + new)


# error logging
def log(message):
    with open(os.path.join('./', 'renaming_log.txt'), 'a') as datafile:
        datafile.write(message + "\n")
        datafile.close()


# Reads the input sample information and creates a dictionary for both sample
# names and sample numbers separately.
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

    paths = dir_tools.get_all_paths(sys.argv)

    for dirname in os.listdir(paths["run_path"]):
        # To check it is dir, not file
        if os.path.isdir(paths["run_path"] + dirname):
            bcpath = paths["run_path"] + dirname + "/Data/Intensities/BaseCalls/"
            for filename in os.listdir(bcpath):
                # Renaming files in directory
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
                        if new_filename != filename:
                            mv_rename(bcpath + filename,
                                      bcpath + new_filename)
                    except KeyError:
                        log(old_num + " not found in dict!")
                elif (".sam" in filename or
                      ".bam" in filename or
                      ".fastq" in filename):
                    log(filename + " name not changed")

            # Renaming directory
            # Assumes Illumina numeric ID is at least 8 digits long
            name_re = re.match(r"([\S*\s*]*[A|B])-\d\d\d\d\d\d\d\d+\Z",
                               dirname)

            if name_re:
                old_name = name_re.group(1)
                try:
                    new_name = sample_names.get(old_name)
                    new_name = new_name.replace(" ", "_")
                    new_dirname = re.sub(r"([\S*\s*]*[A|B])(-\d\d\d\d\d\d\d\d+)\Z",
                                         new_name + r"\g<2>",
                                         dirname)
                    if dirname != new_dirname:
                        mv_rename(paths["run_path"] + dirname,
                                  paths["run_path"] + new_dirname)
                except KeyError:
                    log(old_name + " not found in dict!")
            else:
                log(dirname + " name not changed")

if __name__ == "__main__":
    main()
