###############################################################################
# Run VarScan2
# all variants are saved to one output file
#
# See parse_those_args for required arguments.

import os
import urllib
import json
import subprocess
import argparse

from pyliftover import LiftOver

import dir_tools
import config


# Global var for caching of Oncotator results
oncotator_dict = dict()
# Global var for pyliftover results (it takes a surprisingly long time)
pyliftover_dict = dict()


def pyliftover(hg38_chrom, hg38_coord):
    hg38_key = "%s:%s" % (hg38_chrom, hg38_coord)

    if hg38_key not in pyliftover_dict:
        lo = LiftOver(config.input_dir + "hg38ToHg19.over.chain.gz")
        result = lo.convert_coordinate(hg38_chrom, int(hg38_coord))

        if result is not None:
            coords_list = result[0]

            pyliftover_dict[hg38_key] = {
                "chrom": coords_list[0],
                "coord": str(coords_list[1])
            }

    return pyliftover_dict[hg38_key]


def oncotator(ref, var, hg19_chrom, hg19_coord):
    url = ("http://www.broadinstitute.org/oncotator/mutation/%s_%s_%s_%s_%s/" %
           (hg19_chrom, hg19_coord, hg19_coord, ref, var))

    if url not in oncotator_dict:
        try:
            oncodata = json.loads(
                urllib.request
                .urlopen(url)
                .read()
                .decode("UTF-8"))

            if "protein_change" in oncodata:
                oncotator_dict[url] = oncodata["protein_change"]
            else:
                raise urllib.error.URLError("Expected JSON not received.")
        except (urllib.error.HTTPError,
                urllib.error.URLError) as e:
            print(e.reason)
            print(url)

    return oncotator_dict[url]


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

        hg19 = pyliftover(columns[0], columns[1])
        protein_change = oncotator(ref, var, hg19["chrom"], hg19["coord"])

        final_data.append(protein_change)
        final_data.append(hg19["coord"])
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


def run_varscan(prefix, data_path, region, chrom):
    varscan_args = [config.scripts_dir + "varscan.sh",
                    prefix,
                    data_path,
                    region,
                    "_" + chrom]

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


def execute(run, bedfile, output_filename):
    for dirname in os.listdir(run["path"]):
        # To check it is dir, not file
        if os.path.isdir(run["path"] + dirname):
            # BaseSpace (Illumina) directory structure
            data_path = (run["path"] + dirname +
                         "/Data/Intensities/BaseCalls/")
            sample = dir_tools.get_sample_info(run["path"] + dirname)
            prefix = run["name"] + "_" + sample["num"]
            regions = dir_tools.read_bed(bedfile)

            for r in regions:
                varscan_results = run_varscan(
                    prefix,
                    data_path,
                    r["full"],
                    r["chrom"])

                annotate_results(varscan_results,
                                 output_filename,
                                 sample["num"],
                                 sample["name"],
                                 data_path + prefix)


# Copies column headers to new file
def init_output(output_filename):
    return_code = subprocess.check_call([
        "cp",
        config.input_dir + "varscan_template.txt",
        config.output_dir + output_filename])

    if return_code is not 0:
        print("Varscan results file not initialised.")


def main(shell_args):
    run = dir_tools.get_run_info(shell_args.run_path)
    init_output(shell_args.output_filename)
    execute(run, shell_args.input_bed, shell_args.output_filename)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run VarScan2 using bam files for each sample.")
    parser.add_argument("--run_path", "-rp", type=str, required=True,
                        help="Path to the run directory")
    parser.add_argument("--output_filename", "-o", type=str, required=True,
                        help="Name of the output file.")
    parser.add_argument("--input_bed", "-ib", type=str, required=True,
                        help="Path to the input bed file")

    args = parser.parse_args()
    main(args)
