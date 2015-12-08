import os
import sys
import subprocess
import statistics

import dir_tools
import config


def init_output(output_filename):
    with open(os.getcwd() + "/output/" + output_filename,
              "a") as output_file:
        output_file.write(
            "Sample number\tSample name\tchr1 (H3F3A)\tchr15 (IDH2)\t" +
            "chr2 (IDH1)\tchr3 (CTNNB1)\tchr6 (H3.1)\tchr7 (BRAF)\n")
        output_file.close()


def main():
    run = dir_tools.get_run_info(sys.argv[1])
    init_output(sys.argv[2])

    for dirname in os.listdir(run["path"]):
        # To check it is dir, not file
        if os.path.isdir(run["path"] + dirname):
            # BaseSpace (Illumina) directory structure
            data_path = (run["path"] + dirname +
                         "/Data/Intensities/BaseCalls/")

            sample = dir_tools.get_sample_info(run["path"] + dirname)
            prefix = run["name"] + "_" + sample["num"]

            args_list = [config.bamreadcount_path,
                         "-f",
                         config.ref_fasta,
                         "-l",
                         os.getcwd() + "/input/cranio1.bed",
                         "-w",
                         "0",
                         data_path + prefix + "_sorted.bam"]

            p = subprocess.Popen(args_list, stdout=subprocess.PIPE)
            chrom_data = dict()

            while 1:
                base = p.stdout.readline().decode("UTF-8")
                columns = base.strip().split("\t")

                if ("chr" in columns[0] and
                        columns[0] in chrom_data):
                    chrom_data[columns[0]].append(columns[3])

                elif ("chr" in columns[0] and
                        columns[0] not in chrom_data):
                    chrom_data[columns[0]] = list()
                    chrom_data[columns[0]].append(columns[3])

                if not base and p.returncode is not None:
                    break
                p.poll()

            chrom_averages = dict()

            for (key, value) in chrom_data:
                chrom_averages[key] = statistics.mean(
                    chrom_data[key])

            with open(os.getcwd() + "/output/bam-readcount.txt",
                      "a") as output_file:
                output_file.write(sample["num"] + "\t" +
                                  sample["name"] + "\t" +
                                  chrom_averages["chr1"] + "\t" +
                                  chrom_averages["chr15"] + "\t" +
                                  chrom_averages["chr2"] + "\t" +
                                  chrom_averages["chr3"] + "\t" +
                                  chrom_averages["chr6"] + "\t" +
                                  chrom_averages["chr7"] + "\n")
                output_file.write(base.strip() + "\n")
                output_file.close()

if __name__ == "__main__":
    main()
