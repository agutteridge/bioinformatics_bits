import os
import sys
import subprocess

import dir_tools
import config


def main():
    run = dir_tools.get_run_info(sys.argv)
    sample_list = dir_tools.load_samples()

    for dirname in os.listdir(run["path"]):
        # To check it is dir, not file
        if (os.path.isdir(run["path"] + dirname) and dirname in sample_list):
            # BaseSpace (Illumina) directory structure
            data_path = (run["path"] + dirname +
                         "/Data/Intensities/BaseCalls/")

            sample = dir_tools.get_sample_info(run["path"] + dirname)
            prefix = run["name"] + "_" + sample["num"]

            args_list = [config.bamreadcount_profile,
                         "-f",
                         config.ref_fasta,
                         "-l",
                         "IDH_hotspots.bed",
                         "-w",
                         "1",
                         data_path + prefix + "_sorted.bam"]

            p = subprocess.Popen(args_list, stdout=subprocess.PIPE)

            while 1:
                base = p.stdout.readline().decode("UTF-8")
                if base[0:3] == "chr":
                    with open(os.getcwd() + "/IDH_hotspots.txt",
                              "a") as output_file:
                        output_file.write(sample["num"] + "\t" +
                                          sample["name"] + "\t")
                        output_file.write(base.strip() + "\n")
                if not base and p.returncode is not None:
                    break
                p.poll()


if __name__ == "__main__":
    main()
