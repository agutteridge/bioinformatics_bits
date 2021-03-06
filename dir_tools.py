import os
import sys
import re


# Produces a string of whitespace-delimited regions from a bed file
# to use with -r flag in mpileup (-l with bed file is inefficient)
def read_bed(path_to_bed):
    all_regions = list()

    with open(path_to_bed, "r") as datafile:
        for line in datafile:
            columns = line.strip().split("\t")
            region = {
                "full": "%s:%s-%s" % (columns[0], columns[1], columns[2]),
                "chrom": columns[0]
            }
            all_regions.append(region)
        datafile.close()

    return all_regions


def get_sample_info(sample_path):
    # Sample number
    sample_num = ""
    sample_name = ""
    # adding forward slash at end
    if sample_path[-1] is not "/":
        sample_path = sample_path + "/"

    data_path = sample_path + "Data/Intensities/BaseCalls/"
    # get sample number between S1-96
    for filename in os.listdir(data_path):
        # assumes no non-run files will have _R1_!
        if ".fastq" in filename and sample_num is "":
            sample_re = re.match(r"[\S*\s*]*_(S\d*)_[\S*\s*]*", filename)

            if sample_re:
                sample_num = sample_re.group(1)

    # Assumes Illumina numeric ID is at least 8 digits long
    name_re = re.match(r"[\S*\s*]*/([\S*\s*]*-\d\d\d\d\d\d\d\d+)/\Z",
                       sample_path)

    if name_re:
        sample_name = name_re.group(1)

    try:
        if sample_num is not "" and sample_name is not "":
            return {
                "num": sample_num,
                "name": sample_name
            }
        else:
            raise Exception("Sample name and/or number could not be " +
                            "retrieved.")
    except Exception as e:
        for message in e.args:
            print(message)
            print("sample number was: " + sample_num)
            print("sample name was: " + sample_name)
            sys.exit(0)


# Tests path to run directory, supplied as a command line argument.
# Extracts run name to be used as variable.
def get_run_info(run_path):
    run_name = ""

    try:
        if not os.path.isdir(run_path):
            raise OSError("Not a valid path.")

        # adding forward slash at end
        if run_path[-1] is not "/":
            run_path = run_path + "/"

        # get run name i.e. last dir
        # Assumes Illumina numeric ID is at least 8 digits long
        run_re = re.match(r"[\S*\s*]*/([\S*\s*]*-\d\d\d\d\d\d\d\d+)/\Z",
                          run_path)

        if run_re:
            run_name = run_re.group(1)
        else:
            raise AttributeError("Not a run directory.")
    except IndexError:
        print("Please give the path to the run directory as a command line" +
              " argument.")
        sys.exit(0)
    except OSError as oe:
        for e in oe.args:
            print(e)
        sys.exit(0)
    except AttributeError as ae:
        for e in ae.args:
            print(e)
        sys.exit(0)

    return {
        "name": run_name,
        "path": run_path
    }
