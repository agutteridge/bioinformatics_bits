import os
import urllib
import json
import sys
import subprocess
import re

from pyliftover import LiftOver

import dir_tools


# Global var for caching of Oncotator results
oncotator_dict = dict()


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
    with open(os.path.join(os.getcwd(), "varscan_no_downsampling.csv"), 'a') as datafile:
        datafile.write("S#\tSample name\tDownsample\tChr\thg38 Pos\tdbSNP\t" +
                       "Ref\tVar\tQuality\tFilter\tInfo\tFormat\tSamples\t" +
                       "Protein change\thg19 Pos\tFreq\n")
        datafile.close()


def read_from_vcf(vcf_filename,
                  sample_num,
                  sample_name,
                  downsample,
                  chrom):
    varscan_data = list()

    with open(os.path.join(vcf_filename), 'r') as datafile:
        datafile.readline()  # discard first line

        for l in datafile:
            if l[0] != "#" and l != "":  # VCFv4.1 comments and column header
                varscan_data.append(sample_num)
                varscan_data.append(sample_name)
                varscan_data.append(downsample)
                varscan_data.append(chrom)
                columns = l.strip().split("\t")
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

                # run coverage.sh

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
        if (chrom + ".vcf") in filename:
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


def get_region(chrom):
    if chrom == "chr6":
        return "chr6:26031956-26032018"
    elif chrom == "chr2":
        return "chr2:208248374-208248400"
    elif chrom == "chr3":
        return "chr3:41224593-41224650"
    elif chrom == "chr7":
        return "chr7:140753327-140753392"
    elif chrom == "chr1":
        return "chr1:226064415-226064499"
    elif chrom == "chr15":
        return "chr15:90088582-90088636"
    else:
        raise Exception("Chromosome not in correct format.")


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
                # downsample = 8000 / int(chrom_reads[chrom])
                # if downsample > 1:
                downsample = 1
                prefix = paths["run_name"] + "_" + info["sample_num"]

                # if not check_vcf(data_path, chrom):
                return_code = subprocess.check_call([
                                    os.getcwd() + "/varscan.sh",
                                    prefix,
                                    str(downsample),
                                    chrom,
                                    get_region(chrom),
                                    data_path])

                if return_code is not 0:
                    raise Exception("Error in varscan.sh")
                    sys.exit(0)

                vcf_results = read_from_vcf(data_path +
                                            prefix +
                                            "_" +
                                            chrom +
                                            ".vcf",
                                            info["sample_num"],
                                            info["sample_name"],
                                            "{0:.2f}".format(downsample),
                                            chrom)

                if vcf_results:
                    write_to_csv(data_path, vcf_results)
                else:
                    print(info["sample_name"] +
                          ", " +
                          chrom +
                          " downsampled to " +
                          str(downsample * 100) +
                          "%")


if __name__ == "__main__":
    main()
