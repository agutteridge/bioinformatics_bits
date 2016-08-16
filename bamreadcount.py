import os
import subprocess
import argparse

import dir_tools
import config
from inputfiles import columns
import varscan

MODE = ''
OUTPUT_FILENAME = ''


def init_output():
    with open(config.output_dir + OUTPUT_FILENAME, 'w') as o:
        if MODE == 'READS':
            o.write(columns.reads)
            o.close()
        elif MODE == 'VARIANTS':
            o.write(columns.variants)
            o.close()
        else:
            raise Exception('invalid mode')


def run_reads(data_path, sample, prefix):
    args_list = [config.bamreadcount_path,
                 '-b',
                 '15',
                 '-f',
                 config.ref_fasta,
                 '-l',
                 config.input_dir + 'hotspots.bed',
                 '-w',
                 '0',
                 data_path + prefix + '.bam']

    p = subprocess.Popen(args_list, stdout=subprocess.PIPE)
    chrom_data = dict()

    while 1:
        base = p.stdout.readline().decode('UTF-8')
        base_data = base.strip().split('\t')

        if 'chr' in base_data[0]:
            chrom_data[base_data[0]] = int(base_data[3])

        if not base and p.returncode is not None:
            break
        p.poll()

    with open(config.output_dir + OUTPUT_FILENAME, 'a') as o:
        o.write(sample['num'] + '\t' +
                sample['name'] + '\t')

        chrom_list = ['chr1', 'chr2', 'chr3', 'chr6', 'chr7', 'chr15']

        for c in chrom_list:
            if c in chrom_data:
                o.write(str(chrom_data[c]) + '\t')
            else:
                o.write('0\t')

        o.write('\n')
        o.close()


def per_base(base_data, sample, var=None, ofile=OUTPUT_FILENAME):
    ref_results = list()
    ref = base_data[2]

    # find reference base
    for base in base_data[5:9]:
        quals = base.split(':')
        if quals[0] == ref:
            ref_results.extend(quals[2:])

    results = list()
    depth_all = int(base_data[3])

    # potential var bases (ACGT only)
    for base in base_data[5:9]:
        quals = base.split(':')
        if quals[0] != ref:
            depth_var = int(quals[1])
            freq = (depth_var / depth_all) * 100

            if (not var and freq > 0.5 and depth_var > 2) or var == quals[0]:
                print(ref + ' > ' + quals[0])
                results.extend(base_data[0:3])
                results.append(quals[0])
                results.append(varscan.external(
                    base_data[2],
                    quals[0],
                    base_data[0],
                    base_data[1]))
                results.append(str(depth_all))
                results.append(str(depth_var))
                results.append('{0:.2f}'.format(freq))
                results.extend(quals[2:])
                results.extend(ref_results)

                with open(config.output_dir + ofile, 'a') as o:
                    o.write(sample['num'] + '\t' +
                            sample['name'] + '\t')

                    for r in results:
                        o.write(r + '\t')

                    o.write('\n')
                    o.close()

                results = list()


def run_variants(data_path, sample, prefix):
    args_list = [config.bamreadcount_path,
                 '-f',
                 config.ref_fasta,
                 '-l',
                 config.input_dir + 'cranio1.bed',
                 '-w',
                 '0',
                 data_path + prefix + '.bam']

    p = subprocess.Popen(args_list, stdout=subprocess.PIPE)

    while 1:
        base = p.stdout.readline().decode('UTF-8')
        base_data = base.strip().split('\t')

        if 'chr' in base_data[0]:
            per_base(base_data, sample)

        if not base and p.returncode is not None:
            break
        p.poll()


def main(shell_args):
    global MODE, OUTPUT_FILENAME
    MODE = str.upper(shell_args.mode)
    OUTPUT_FILENAME = shell_args.output_filename
    run = dir_tools.get_run_info(shell_args.run_path)
    init_output()

    for dirname in os.listdir(run['path']):
        # To check it is dir, not file
        if os.path.isdir(run['path'] + dirname):
            # BaseSpace (Illumina) directory structure
            data_path = (run['path'] + dirname +
                         '/Data/Intensities/BaseCalls/')

            sample = dir_tools.get_sample_info(run['path'] + dirname)
            prefix = run['name'] + '_' + sample['num']

            if MODE == 'READS':
                run_reads(data_path, sample, prefix)
            elif MODE == 'VARIANTS':
                run_variants(data_path, sample, prefix)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Run bamreadcount.')
    parser.add_argument('--run_path', '-rp', type=str, required=True,
                        help='Path to the run directory')
    parser.add_argument('--output_filename', '-o', type=str, required=True,
                        help='Name of the output file.')
    parser.add_argument('--mode', '-m', type=str, required=True,
                        help='READS or VARIANTS')

    args = parser.parse_args()
    main(args)
