###############################################################################
# Calls shell script with VarScan2 and Samtools mpileup
# Variants are annotated (requires connection to internet) and written to file
###############################################################################
import os
import subprocess
import argparse

import dir_tools
import config
import annotate


def annotate_results(input_data,
                     output_filename,
                     sample_num,
                     sample_name,
                     whole_path):
    final_data = list()

    for l in input_data:
        final_data.append(sample_num)  # first column is sample number
        final_data.append(sample_name)  # second column is sample name
        columns = l.strip().split("\t")
        final_data.extend(columns[:6])  # direct from VCF
        extra_info = columns[9].strip().split(":")
        final_data.extend(extra_info)

        ref = columns[3]
        var = columns[4]

        # protein change
        final_data.append(annotate.run(ref, var, columns[0], columns[1]))
        final_data.append("\n")

    write_to_file(output_filename, final_data)


def write_to_file(output_filename, data):
    with open(os.path.join(config.output_dir, output_filename),
              "a") as datafile:
        for v in data:
            if v:
                if v[-1] == "\n":
                    datafile.write(v)
                else:
                    datafile.write(v + "\t")
        datafile.close()


def run_varscan(prefix, data_path):
    varscan_args = [config.scripts_dir + "varscan.sh",
                    prefix,
                    data_path]

    p = subprocess.Popen(varscan_args, stdout=subprocess.PIPE)
    row_list = list()

    while 1:
        row = p.stdout.readline().decode("UTF-8")

        if not row and p.returncode is not None:
            break
        elif row:
            if row[0] != "#":  # VCFv4.1 comments and column header
                row_list.append(row.strip())

        p.poll()

    print("done %d" % p.returncode)

    return row_list


def execute(run, output_filename):
    for dirname in os.listdir(run["path"]):
        # To check it is dir, not file
        if os.path.isdir(run["path"] + dirname):
            # BaseSpace (Illumina) directory structure
            data_path = (run["path"] + dirname +
                         "/Data/Intensities/BaseCalls/")
            sample = dir_tools.get_sample_info(run["path"] + dirname)
            prefix = run["name"] + "_" + sample["num"]

            varscan_results = run_varscan(
                prefix,
                data_path)

            if varscan_results != "no reads.":
                annotate_results(varscan_results,
                                 output_filename,
                                 sample["num"],
                                 sample["name"],
                                 data_path + prefix)


# Copies column headers to new file
def init_output(output_filename):
    inf = config.input_dir + "varscan_template.txt"
    outf = config.output_dir + output_filename

    with open(inf, mode="r") as template, open(outf, mode="a") as output:
        for line in template:
            output.write(line)


def main(shell_args):
    run = dir_tools.get_run_info(shell_args.run_path)
    init_output(shell_args.output_filename)
    execute(run, shell_args.output_filename)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run VarScan2 using bam files for each sample.")
    parser.add_argument("--run_path", "-rp", type=str, required=True,
                        help="Path to the run directory")
    parser.add_argument("--output_filename", "-o", type=str, required=True,
                        help="Name of the output file (incl extension)")

    args = parser.parse_args()
    main(args)
