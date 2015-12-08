import os
import sys
import subprocess
import re

import dir_tools


def main():
    run = dir_tools.get_run_info(sys.argv[1])
    number_alignments = 0

    for dirname in os.listdir(run["path"]):

        # To check it is dir, not file
        if os.path.isdir(run["path"] + dirname):
            # BaseSpace (Illumina) directory structure
            data_path = (run["path"] + dirname +
                         "/Data/Intensities/BaseCalls/")

            # So that alignments are not done twice
            alignment_done = False
            for filename in os.listdir(data_path):
                if ".bam.bai" in filename:
                    alignment_done = True

            if not alignment_done:
                try:
                    # Whole names of downloaded, compressed FastQ files
                    file1 = ""
                    file2 = ""
                    # Sample number
                    sample = ""

                    for filename in os.listdir(data_path):
                        # assumes no non-fastq files will have _R1_ or _R2_ !
                        if "_R1_" in filename:

                            if sample is "":
                                # get sample number between S1-96
                                sample_re = re.match(
                                    r"[\S*\s*]*_(S\d+)_[\S*\s*]*",
                                    filename)

                                if sample_re:
                                    sample = sample_re.group(1)
                                else:
                                    log("Regex error in " + filename)

                            if "fastq.gz" in filename:
                                subprocess.check_call(["gunzip",
                                                       data_path + filename])
                                # removes .gz from strings
                                file1 = filename[:-3]
                            else:
                                file1 = filename

                        elif "_R2_" in filename:
                            if "fastq.gz" in filename:
                                subprocess.check_call(["gunzip",
                                                       data_path + filename])
                                file2 = filename[:-3]
                            else:
                                file2 = filename

                    # SAM and BAM files will have prefix of run#_S#
                    run_sample = run["name"] + "_" + sample

                    if (file1 is not "" and
                            file2 is not "" and
                            run_sample is not ""):

                        return_code = subprocess.check_call([
                                        os.getcwd() + "/shell/batch_align.sh",
                                        file1,
                                        file2,
                                        run_sample,
                                        data_path])

                        if return_code is 0:
                            number_alignments = number_alignments + 1
                            log(dirname + " aligned successfully.")

                        else:
                            log("error no. %s: %s, subprocess.check_call."
                                % (str(return_code), sample))

                    else:
                        log("unnamed variable in %s\nfile1: %s\nfile2: %s\nfull sample name: %s" %
                            (sample, file1, file2, run_sample))

                except Exception as e:
                    log("Exception in directory " + dirname)
                    log(e)
                    sys.exit(0)

    log("Alignments so far: " + str(number_alignments))


def log(message):
    with open(os.path.join(os.getcwd(), 'batch_align_log.txt'),
              'a') as datafile:
        datafile.write(message + "\n")
        datafile.close()

if __name__ == "__main__":
    main()
