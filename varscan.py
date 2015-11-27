import os
import urllib
import json
import sys
import subprocess
import re

from pyliftover import LiftOver

import dir_tools
import config


# Global var for caching of Oncotator results
oncotator_dict = dict()


def url_request(chrom, coord, ref, var):
    url = ("http://www.broadinstitute.org/oncotator/mutation/" +
           chrom + "_" +
           coord + "_" +
           coord + "_" +
           ref + "_" +
           var + "/")

    if url in oncotator_dict:
        return oncotator_dict[url]
    else:
        try:
            result = json.loads(
                urllib.request
                .urlopen(url)
                .read()
                .decode('UTF-8')
            )

            if "protein_change" in result:
                oncotator_dict[url] = {
                    "protein_change": result["protein_change"],
                    "start": result["start"]
                }
                return result
            else:
                raise urllib.error.URLError("Expected JSON not received.")
        except (urllib.error.HTTPError,
                urllib.error.URLError) as e:
            print(e.reason)
            print(url)


def init_varscan_csv():
    with open(os.path.join(os.getcwd(), "varscan_no_downsampling.csv"),
              'a') as datafile:
        datafile.write("S#\tSample name\tChr\thg38 Pos\tdbSNP\t" +
                       "Ref\tVar\tQuality\tFilter\tInfo\tFormat\tSamples\t" +
                       "Protein change\thg19 Pos\tFreq\n")
        datafile.close()


def read_from_vcf(vcf_filename,
                  sample_num,
                  sample_name,
                  chrom):
    varscan_data = list()

    with open(os.path.join(vcf_filename), 'r') as datafile:

        for l in datafile:
            if l[0] != "#" and l != "":  # VCFv4.1 comments and column header
                varscan_data.append(sample_num)
                varscan_data.append(sample_name)
                columns = l.strip().split("\t")
                varscan_data.append(columns[0])
                varscan_data.extend(columns[1:])
                hg38_chrom = columns[0]
                hg38_coord = int(columns[1])
                ref = columns[3]
                var = columns[4]

                varscan_data.extend(
                    liftover_oncotator(
                        ref,
                        var,
                        hg38_chrom,
                        hg38_coord))

                re_freq = re.match(r"[\S*\s*]*:(\d*.*\d*%):[\S*\s*]*",
                                   columns[9])
                if re_freq:
                    frequency = re_freq.group(1)
                    varscan_data.append(frequency)

                varscan_data.append("\n")

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

        if oncodata["protein_change"] is not None:
            results.append(oncodata["protein_change"])
        else:
            results.append("\t")

        results.append(oncodata["start"])
    else:
        results.append("hg19 coordinate not found.")

    return results


# Returns True if *chrN.vcf exists already
def check_vcf(data_path, chrom):
    for filename in os.listdir(data_path):
        if (chrom + "_nd.vcf") in filename:
            return True

    return False


def write_to_csv(data_path, var_info):
    with open(os.path.join(os.getcwd(), "varscan_no_downsampling.csv"),
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
    run = dir_tools.get_run_info(sys.argv)

    for dirname in os.listdir(run["path"]):
        # To check it is dir, not file
        if os.path.isdir(run["path"] + dirname):
            # BaseSpace (Illumina) directory structure
            data_path = (run["path"] + dirname +
                         "/Data/Intensities/BaseCalls/")

            sample = dir_tools.get_sample_info(run["path"] + dirname)

            prefix = run["name"] + "_" + sample["num"]

            # if not check_vcf(data_path, chrom):
            return_code = subprocess.check_call([
                                os.getcwd() + "/bash/varscan.sh",
                                prefix,
                                data_path,
                                config.input_dir])

            if return_code is not 0:
                raise Exception("Error in varscan.sh")
                sys.exit(0)

            vcf_results = read_from_vcf(data_path +
                                        prefix +
                                        "_nd.vcf",
                                        sample["num"],
                                        sample["name"])

            if vcf_results:
                write_to_csv(data_path, vcf_results)


if __name__ == "__main__":
    main()
