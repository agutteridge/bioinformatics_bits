import os
import sys
import subprocess

import dir_tools


def main():
    paths = dir_tools.get_all_paths(sys.argv)

    with open(os.getcwd() + '/varscan_anomalous.csv', 'r') as input_file:
        input_file.readline()  # ignore header
        all_columns = list()

        for varscan_line in input_file:
            varscan_columns = varscan_line.strip().split("\t")
            region = (varscan_columns[3] + ":" + varscan_columns[4] +
                      "-" + varscan_columns[4])

            if os.path.isdir(paths["run_path"] + varscan_columns[1]):
                print("processing " + varscan_columns[1] + ", " + region)
                popen_args = [os.getcwd() + '/coverage.sh',
                              paths["run_name"],
                              varscan_columns[1],
                              varscan_columns[0],
                              region]

                p = subprocess.Popen(popen_args, stdout=subprocess.PIPE)

                while 1:
                    coverage = p.stdout.readline().decode('UTF-8')
                    varscan_columns.append(coverage.strip())
                    if not coverage and p.returncode is not None:
                        break
                    p.poll()

                all_columns.append('\t'.join(varscan_columns))

        with open(os.getcwd() + '/varscan_coverage.csv', 'w')as output_file:
            for c in all_columns:
                output_file.write(c + "\n")
            output_file.close()

if __name__ == "__main__":
    main()
