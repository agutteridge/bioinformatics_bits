import os
import sys
import subprocess
import re


def main():
    number_alignments = 0

    path_str = ""

    try:
        path_str = sys.argv[1]
    except IndexError:
        print("""Please give the path to the run directory as a command line
               argument.""")
        sys.exit(0)

    # omitting forward slash at beginning
    if path_str[0] is "/":
        path_str = path_str[1:]

    # adding forward slash at end
    if path_str[-1] is not "/":
        path_str = path_str + "/"

    # get run name
    run_name = ""
    run_re = re.match(r"[\S*\s*]*/([\S*\s*]*)/\Z", path_str)

    if run_re:
        run_name = run_re.group(1)

    bio = os.getcwd() + "/"
    root = bio + path_str
    print("Path set to " + root)

    for dirname in os.listdir(root):

        # To check it is dir, not file
        if os.path.isdir(root + dirname):
            # BaseSpace (Illumina) directory structure
            temppath = root + dirname + "/Data/Intensities/BaseCalls/"

            # So that alignments are not done twice
            alignment_done = False
            for filename in os.listdir(temppath):
                if ".bam.bai" in filename:
                    alignment_done = True

            if not alignment_done:
                try:
                    # Whole names of downloaded, compressed FastQ files
                    file1 = ""
                    file2 = ""
                    # Sample number
                    sample = ""

                    for filename in os.listdir(temppath):
                        # assumes no non-fastq files will have _R1_ or _R2_ !
                        if "_R1_" in filename:

                            if sample is "":
                                # get sample number between S1-96
                                sample_re = re.match(
                                    r"[\S*\s*]*_(S\d*)_[\S*\s*]*",
                                    filename)

                                if sample_re:
                                    sample = sample_re.group(1)

                            if "fastq.gz" in filename:
                                subprocess.check_call(["gunzip",
                                                       temppath + filename])
                                # removes .gz from strings
                                file1 = filename[:-3]
                            else:
                                file1 = filename

                        elif "_R2_" in filename:
                            if "fastq.gz" in filename:
                                subprocess.check_call(["gunzip",
                                                       temppath + filename])
                                file2 = filename[:-3]
                            else:
                                file2 = filename

                    # SAM and BAM files will have prefix of run#_S#
                    run_sample = run_name + "_" + sample

                    if (file1 is not "" and
                            file2 is not "" and
                            run_sample is not ""):

                        return_code = subprocess.check_call([
                                            bio + "batch_align.sh",
                                            file1,
                                            file2,
                                            run_sample,
                                            temppath,
                                            bio + "ref/hg38bwaidx"])

                        if return_code is 0:
                            number_alignments = number_alignments + 1

                            with open(os.path.join('./', 'log.txt'),
                                      'a') as datafile:
                                datafile.write(dirname +
                                               "aligned and indexed.\n")
                                datafile.write("Alignments so far: " +
                                               str(number_alignments))
                                datafile.close()

                        else:
                            problem(dirname + " with subprocess.check_call.")

                    else:
                        problem(dirname +
                                " with getting sample and file names.")

                except Exception as e:
                    problem(dirname)
                    print(e)
                    sys.exit(0)


def problem(dirname):
    with open(os.path.join('./', 'log.txt'), 'a') as datafile:
        datafile.write("agh! problem in " + dirname + "\n")
        datafile.close()

if __name__ == "__main__":
    main()
