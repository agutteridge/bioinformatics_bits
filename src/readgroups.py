import os
import sys
import subprocess

import dir_tools


def main():
    paths = dir_tools.get_all_paths(sys.argv)

    for dirname in os.listdir(paths["run_path"]):
        # To check it is dir, not file
        if os.path.isdir(paths["run_path"] + dirname):
            # BaseSpace (Illumina) directory structure
            data_path = (paths["run_path"] + dirname +
                         "/Data/Intensities/BaseCalls/")

            info = dir_tools.get_sample_info(paths["run_path"] + dirname)
            sample_ID = data_path + paths["run_name"] + "_" + info["sample_num"]

            return_code = subprocess.check_call([
                                os.getcwd() + "/readgroups.sh",
                                sample_ID,
                                info["sample_name"]])

            if return_code is not 0:
                raise Exception("Error in " + info["sample_name"])
                sys.exit(0)

if __name__ == "__main__":
    main()
