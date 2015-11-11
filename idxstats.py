import os
import sys
import subprocess

import dir_tools


def main():
    paths = dir_tools.get_all_paths(sys.argv)
    for dirname in os.listdir(paths["run_path"]):
        sample_path = paths["run_path"] + dirname + "/"
        # To check it is dir, not file
        if os.path.isdir(sample_path):
            info = dir_tools.get_sample_info(sample_path)
            data_path = sample_path + "Data/Intensities/BaseCalls/"

            if os.path.isdir(data_path):
                for filename in os.listdir(data_path):
                    if "_sorted.bam" in filename:
                        to_write = list()  # to write to .csv file
                        to_write.append(info["sample_num"])
                        to_write.append(info["sample_name"])

                        # subset!!!
                        subprocess_args = [samtools view -Shu -L cranio1.bed S34_sorted.bam > S34_sorted_subset.bam]
                        index_args = [samtools index S34_sorted_subset.bam]

                        subprocess.check_call()

                        popen_args = ["samtools",
                                      "idxstats",
                                      data_path + filename]

                        p = subprocess.Popen(popen_args, stdout=subprocess.PIPE)

                        while 1:
                            row = p.stdout.readline()
                            row_list = row.decode('UTF-8').split("\t")
                            if len(row_list) > 1:
                                if (row_list[1] is not "0"
                                        or row_list[2] is not "0"):
                                    mapped = row_list[1]
                                    unmapped = row_list[2]
                                    to_write.append(row_list[0])
                                    to_write.append(mapped)
                                    to_write.append(unmapped)
                                    percent = int(mapped) / (int(unmapped) + int(mapped))
                                    to_write.append("{0:.1f}".format(percent * 100))
                            if not row and p.returncode is not None:
                                break
                            p.poll()
                        print("done %d" % p.returncode)

                        with open(os.path.join(os.getcwd(), 'idxstats.csv'), 'a') as datafile:
                            for w in to_write:
                                datafile.write(w + "\t")

                            datafile.write("\n")
                            datafile.close()

            # Iterate through files, produce idxstats, append to file!
            # for filename in temppath:
            #     print("do something")


if __name__ == "__main__":
    main()
