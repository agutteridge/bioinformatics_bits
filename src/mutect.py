import os
import sys
import subprocess
import re

import dir_tools
import varscan


def read_write(vcf_path, result_path, sample_name, sample_num):
    with open(vcf_path, 'r') as input_file:
        with open(result_path, 'a') as output_file:
            input_file.readline()  # Title
            input_file.readline()  # Column names

            for line in input_file:
                columns = line.split("\t")
                hg38_chrom = columns[0]
                hg38_coord = int(columns[1])
                ref = columns[3]
                var = columns[4]

                output_file.write(sample_name + "\t")
                output_file.write(sample_num + "\t")
                output_file.write(line.strip())

                # if "KEEP" in columns[34]:
                extra_info = varscan.liftover_oncotator(ref, var, hg38_chrom, hg38_coord)
                output_file.write("\t" + extra_info[0] + "\t")
                output_file.write(extra_info[1] + "\n")
                # else:
                #     output_file.write("\n")



# Returns True if *_call_stats.txt exists already
def check_output(data_path):
    for filename in os.listdir(data_path):
        if "_call_stats.txt" in filename:
            return True

    return False


def main():
    paths = dir_tools.get_all_paths(sys.argv)

    for dirname in os.listdir(paths["run_path"]):
        # Not running MuTect on hgDNA samples
        is_normal = False
        re_dir = re.match(r"H\d_[A|B][\S*\s*]*", dirname)
        if re_dir:
            is_normal = True

        # To check it is dir, not file
        if os.path.isdir(paths["run_path"] + dirname) and not is_normal:
            # BaseSpace (Illumina) directory structure
            data_path = (paths["run_path"] + dirname +
                         "/Data/Intensities/BaseCalls/")

            info = dir_tools.get_sample_info(paths["run_path"] + dirname)
            sample_ID = data_path + paths["run_name"] + "_" + info["sample_num"]

            if not check_output(data_path):
                return_code = subprocess.check_call([
                                    os.getcwd() + "/mutect.sh",
                                    sample_ID])

                if return_code is not 0:
                    print("Error in " + info["sample_name"])
                    # sys.exit(0)

            read_write(sample_ID + "_call_stats.txt",
                       os.getcwd() + "/mutect_results.txt",
                       info["sample_name"],
                       info["sample_num"])


if __name__ == "__main__":
    main()
