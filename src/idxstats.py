import os
import sys
import subprocess

import dir_tools


def main():
    paths = dir_tools.get_all_paths(sys.argv)

    with open(os.path.join(os.getcwd(), 'noclips_idxstats.csv'), 'a') as datafile:
        # Column headers
        datafile.write("S#\tType\tchr1 (H3F3A), mapped\t" +
                       "chr1 (H3F3A), unmapped\tchr15 (IDH2), mapped\t" +
                       "chr15 (IDH2), unmapped\tchr2 (IDH1), mapped\t" +
                       "chr2 (IDH1), unmapped\tchr3 (CTNNB1), mapped\t" +
                       "chr3 (CTNNB1), unmapped\tchr6 (H3.1), mapped\t" +
                       "chr6 (H3.1), unmapped\tchr7 (BRAF), mapped\t" +
                       "chr7 (BRAF), unmapped\tother, mapped\t" +
                       "other, unmapped\t% mapped\t% mapped and on target\n")
        datafile.close()

    for dirname in os.listdir(paths["run_path"]):
        sample_path = paths["run_path"] + dirname + "/"
        # To check it is dir, not file
        if os.path.isdir(sample_path):
            info = dir_tools.get_sample_info(sample_path)
            data_path = sample_path + "Data/Intensities/BaseCalls/"

            if os.path.isdir(data_path):
                for filename in os.listdir(data_path):
                    if ("noclips.bam" in filename and
                            "noclips.bam.bai" not in filename):
                        to_write = list()  # to write to .csv file
                        to_write.append(info["sample_num"])
                        print(info["sample_num"])
                        to_write.append(info["sample_name"])

                        popen_args = ["samtools",
                                      "idxstats",
                                      data_path + filename]

                        p = subprocess.Popen(popen_args,
                                             stdout=subprocess.PIPE)

                        # For calculations
                        total = 0
                        total_mapped = 0
                        total_mapped_on_target = 0
                        other_mapped = 0
                        other_unmapped = 0

                        while 1:
                            row = p.stdout.readline()
                            rows = row.decode('UTF-8').split("\t")
                            row_list = list()
                            for r in rows:
                                row_list.append(r.strip())

                            if len(row_list) > 1:
                                if (row_list[0] == "chr1" or
                                        row_list[0] == "chr15" or
                                        row_list[0] == "chr2" or
                                        row_list[0] == "chr3" or
                                        row_list[0] == "chr6" or
                                        row_list[0] == "chr7"):
                                    to_write.append(row_list[2])
                                    to_write.append(row_list[3])
                                    total_mapped_on_target += int(row_list[2])
                                else:
                                    other_mapped += int(row_list[2])
                                    other_unmapped += int(row_list[3])
                                total_mapped += int(row_list[2])
                                total += int(row_list[2]) + int(row_list[3])

                            if not row and p.returncode is not None:
                                break
                            p.poll()
                        print("done %d" % p.returncode)

                        to_write.append(str(other_mapped))
                        to_write.append(str(other_unmapped))
                        if total is not 0:
                            to_write.append("{0:.1f}".format(
                                total_mapped / total * 100))
                            to_write.append("{0:.1f}".format(
                                total_mapped_on_target / total * 100))

                        with open(os.path.join(os.getcwd(), 'noclips_idxstats.csv'),
                                  'a') as datafile:
                            for w in to_write:
                                datafile.write(w + "\t")

                            datafile.write("\n")
                            datafile.close()


if __name__ == "__main__":
    main()
