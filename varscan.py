import os
import urllib
import json
import sys
import subprocess
import re

from pyliftover import LiftOver

import dir_tools


# Create dict with sample number as key and dict of chromosome read counts as
# value
def read_idxstats():
    with open(os.path.join(os.getcwd(), 'idxstats.csv'), 'r') as datafile:
        all_data = dict()
        datafile.readline()  # discard first line

        for line in datafile:
            columns = line.split("\t")
            all_data[columns[0]] = {
                "chr1": columns[2],
                "chr15": columns[4],
                "chr2": columns[6],
                "chr3": columns[8],
                "chr6": columns[10],
                "chr7": columns[12]
            }

    datafile.close()
    return all_data


def url_request(chrom, coord, ref, var):
    url = ("http://www.broadinstitute.org/oncotator/mutation/" +
           chrom + "_" +
           coord + "_" +
           coord + "_" +
           ref + "_" +
           var + "/")

    try:
        result = json.loads(
            urllib.request
            .urlopen(url)
            .read()
            .decode('UTF-8')
        )

        if "protein_change" in result:
            return result
        else:
            raise urllib.error.URLError("Expected JSON not received.")
    except (urllib.error.HTTPError,
            urllib.error.URLError) as e:
        print(e.reason)
        print(url)


def init_varscan_csv():
    with open(os.path.join(os.getcwd(), "varscan.csv"), 'a') as datafile:
        datafile.write("S#\tSample name\tDownsample\tChr\thg38 Pos\tRef\tVar" +
                       "\tFreq\tProtein change\thg19 Pos\tdbSNP\n")
        datafile.close()


def read_from_vcf(vcf_filename):
    varscan_data = list()

    with open(os.path.join(vcf_filename), 'r') as datafile:
        datafile.readline()  # discard first line

        for l in datafile:
            columns = l.split("\t")
            varscan_data += columns[1:4]
            hg38_chrom = columns[0]
            hg38_coord = int(columns[1])
            ref = columns[2]
            var = columns[3]

            re_freq = re.match(r"[\S*\s*]*:(\d*.*\d*%):[\S*\s*]*",
                               columns[4])
            frequency = ""

            if re_freq:
                frequency = re_freq.group(1)
                varscan_data.append(frequency)
                varscan_data.extend(
                    liftover_oncotator(
                        ref,
                        var,
                        hg38_chrom,
                        hg38_coord))

    return varscan_data


def liftover_oncotator(ref, var, hg38_chrom, hg38_coord):
    results = list()
    # zip file in root folder
    lo = LiftOver("hg38ToHg19.over.chain.gz")
    result = lo.convert_coordinate(hg38_chrom, hg38_coord)

    if result is not None:
        coords_list = result[0]
        print("hg19: " + coords_list[0] + " " +
              str(coords_list[1]))
        oncodata = url_request(coords_list[0],
                               str(coords_list[1]),
                               ref,
                               var)
        results.append(oncodata["protein_change"])
        results.append(oncodata["start"])
        results.append(oncodata["dbSNP_RS"])
        results.append("\n")
    else:
        results.append("hg19 coordinate not found.")
        return results


# Returns True if *chrN.vcf exists already
def check_vcf(data_path, chrom):
    for filename in os.listdir(data_path):
        if (chrom + ".vcf") in filename:
            return True

    return False


def write_to_csv(data_path, var_info):
    with open(os.path.join(os.getcwd(), "varscan.csv"),
              'a') as datafile:
        for v in var_info:
            if len(v) > 0:
                if v[-1] == "\n":
                    datafile.write(v)
                else:
                    datafile.write(v + "\t")
            else:
                datafile.write(v + "\t")
        datafile.close()


def main():
    init_varscan_csv()
    read_data = read_idxstats()
    paths = dir_tools.get_all_paths(sys.argv)

    for dirname in os.listdir(paths["run_path"]):
        # To check it is dir, not file
        if os.path.isdir(paths["run_path"] + dirname):
            # BaseSpace (Illumina) directory structure
            data_path = (paths["run_path"] + dirname +
                         "/Data/Intensities/BaseCalls/")

            info = dir_tools.get_sample_info(paths["run_path"] + dirname)
            chrom_reads = read_data[info["sample_num"]]

            for chrom in chrom_reads:
                var_info = list()
                downsample = 8000 / int(chrom_reads[chrom])
                prefix = paths["run_name"] + "_" + info["sample_num"]
                var_info.append(info["sample_num"])
                var_info.append(info["sample_name"])
                var_info.append("{0:.2f}".format(downsample))
                var_info.append(chrom)

                if not check_vcf():
                    return_code = subprocess.check_call([
                                        os.getcwd() + "/varscan.sh",
                                        prefix,
                                        str(downsample),
                                        chrom,
                                        data_path])

                    if return_code is 0:
                        print("sample " + info["sample_num"] + " done")

                var_info.extend(
                    read_from_vcf(data_path,
                                  prefix +
                                  "_" +
                                  chrom +
                                  ".vcf"))

                write_to_csv(data_path, var_info)


if __name__ == "__main__":
    main()
