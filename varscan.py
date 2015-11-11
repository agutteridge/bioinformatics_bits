import os
import urllib
import json
import sys
import subprocess
import re

from pyliftover import LiftOver

import dir_tools


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


def init_varscan_txt():
    with open(os.path.join(os.getcwd(), "varscan.csv"), 'a') as datafile:
        datafile.write("S\#\tDownsample\tChr\thg38 Pos\tRef\tVar\tFreq\tProtein change\thg19 Pos\n")
        datafile.close()


def main():
    init_varscan_txt()
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
                downsample = 8000 / int(chrom_reads[chrom])
                prefix = paths["run_name"] + "_" + info["sample_num"]
                return_code = subprocess.check_call([
                                    os.getcwd() + "/varscan.sh",
                                    prefix,
                                    str(downsample),
                                    chrom,
                                    data_path])

                if return_code is 0:
                    print("sample " + info["sample_num"] + " done")

                var_info = list()

                with open(os.path.join(data_path, prefix + "_" + chrom + ".vcf"), 'r') as datafile:
                    datafile.readline()  # discard first line

                    for l in datafile:
                        var_info.append(info["sample_num"])
                        var_info.append("{0:.2f}".format(downsample))
                        var_info.append(chrom)
                        columns = l.split("\t")
                        var_info += columns[1:4]
                        hg38_chrom = columns[0]
                        hg38_coord = int(columns[1])
                        ref = columns[2]
                        var = columns[3]

                        re_freq = re.match(r"[\S*\s*]*:(\d*.*\d*%):[\S*\s*]*", columns[4])
                        frequency = ""

                        if re_freq:
                            frequency = re_freq.group(1)
                            var_info.append(frequency)

                        # zip file in root folder
                        lo = LiftOver("hg38ToHg19.over.chain.gz")
                        result = lo.convert_coordinate(hg38_chrom, hg38_coord)

                        if result is not None:
                            coords_list = result[0]
                            print("hg19: " + coords_list[0] + " " + str(coords_list[1]))
                            oncodata = url_request(coords_list[0],
                                                   str(coords_list[1]),
                                                   ref,
                                                   var)
                            var_info.append(oncodata["protein_change"])
                            var_info.append(oncodata["start"])
                            var_info.append("\n")
                        else:
                            print("hg19 coordinate not found.")

                with open(os.path.join(os.getcwd(), "varscan.csv"), 'a') as datafile:
                    length = len(var_info)

                    i = 0
                    for v in var_info:
                        if i is length - 1:
                            datafile.write(v)
                        else:
                            datafile.write(v + "\t")
                    datafile.close()


# returns None or a list
if __name__ == "__main__":
    main()
