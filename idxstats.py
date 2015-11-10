import os
import sys
import subprocess
import re


def main():
    all_stats = list()
    root = os.getcwd()

    for dirname in os.listdir(root):

        # To check it is dir, not file
        if os.path.isdir(root + dirname):
            temppath = root + dirname + "/Data/Intensities/BaseCalls/"

            # Sample number
            sample = ""

            both_gz = 0

            for filename in os.listdir(temppath):

                # assumes no non-run files will have _R1_ or _R2_ !
                if "_R1_" in filename:

                    if sample is "":
                        # get sample number between S1-96
                        sample_re = re.match(r"[\S*\s*]*_(S\d*)_[\S*\s*]*", filename)
 
                        if sample_re:
                            sample = sample_re.group(1)                        

                    if "fastq.gz" in filename:
                        subprocess.check_call(["gunzip", temppath + filename])
                        # removes .gz from strings
                        file1 = filename[:-3]
                    else:
                        file1=filename
                            
                elif "_R2_" in filename:
                    if "fastq.gz" in filename:
                        subprocess.check_call(["gunzip", temppath + filename])
                        file2 = filename[:-3]
                    else:
                        file2=filename

            # SAM and BAM files will have prefix of run_S#
            run_sample = run_name + "_" + sample

            if file1 is not "" and file2 is not "" and run_sample is not "":
                return_code = subprocess.check_call([
                    bio + "batch_align.sh",
                    file1, 
                    file2,
                    run_sample,
                    temppath,
                    bio + "ref/hg38bwaidx"])

                if return_code is 0:
                    number_alignments = number_alignments + 1

                    with open(os.path.join('./', 'log.txt'), 'a') as datafile:
                        datafile.write(dirname + "aligned and indexed.\n")
                        datafile.write("Alignments so far: " + str(number_alignments))
                        datafile.close()

                else:
                    problem(dirname)

            else:
                problem(dirname)

def problem(dirname):
    with open(os.path.join('./', 'log.txt'), 'a') as datafile:
        datafile.write("agh! problem in " + dirname + "\n")
        datafile.close()

if __name__ == "__main__":
    main()