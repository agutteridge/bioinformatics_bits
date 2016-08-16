import os
import socket
import urllib
import json
from datetime import datetime

from pyliftover import LiftOver

import config

# Global var for caching of Oncotator results
oncotator_dict = dict()
# Global var for pyliftover results (it takes a surprisingly long time)
pyliftover_dict = dict()
# Output file is in csv format
# with the nomenclature e.g. varscan_160812-165400.csv
fname = 'varscan_' + datetime.today().strftime('%y%m%d-%H%M%S') + '.csv'


def write_to_file(data):
    with open(os.path.join(config.output_dir, fname), 'a') as df:
        for v in data:
            if v:
                if v[-1] == '\n':
                    df.write(v)
                else:
                    df.write(v + '\t')
        df.close()


# Returns true if a connection to the internet is present
def is_connected():
    try:
        host = socket.gethostbyname('www.google.com')
        socket.create_connection((host, 80), 2)
        return True
    except:
        return False


def pyliftover(hg38_chrom, hg38_coord):
    hg38_key = '%s:%s' % (hg38_chrom, hg38_coord)

    if hg38_key not in pyliftover_dict:
        lo = LiftOver(config.input_dir + 'hg38ToHg19.over.chain.gz')
        result = lo.convert_coordinate(hg38_chrom, int(hg38_coord))

        if result is not None:
            coords_list = result[0]

            pyliftover_dict[hg38_key] = {
                'chrom': coords_list[0],
                'coord': str(coords_list[1])
            }

    return pyliftover_dict[hg38_key]


def oncotator(ref, var, hg19_chrom, hg19_coord):
    url = ('http://www.broadinstitute.org/oncotator/mutation/%s_%s_%s_%s_%s/' %
           (hg19_chrom, hg19_coord, hg19_coord, ref, var))

    if url not in oncotator_dict:
        try:
            oncodata = json.loads(
                urllib.request
                .urlopen(url)
                .read()
                .decode('UTF-8'))

            if 'protein_change' in oncodata:
                if oncodata['protein_change'] is '':
                    oncotator_dict[url] = 'n/a'
                else:
                    oncotator_dict[url] = (oncodata['gene'] +
                                           '-' +
                                           oncodata['protein_change'][2:])
                    print(oncotator_dict[url])
            else:
                raise urllib.error.URLError('Expected JSON not received.')
        except (urllib.error.HTTPError,
                urllib.error.URLError) as e:
            print(e.reason)
            print(url)

    return oncotator_dict[url]


def annotate(ref, var, chrom, coord):
    if is_connected():
        hg19 = pyliftover(chrom, coord)
        protein_change = oncotator(ref, var, hg19['chrom'], hg19['coord'])
        return protein_change
    else:
        return 'offline mode'


def read_vcf(path):
    sample_name = 'not sure what your nomenclature is, so regex here'
    sample_num = 'same here'
    # let me know if you need help with regex!
    # or do str.split('_') if this works/is easier

    with open(path, mode='r') as v:
        # each line of vcf
        for l in v:
            if l[0] != "#":  # VCFv4.1 comments and column header
                final_data = list()

                final_data.append(sample_num)  # first column is sample number
                final_data.append(sample_name)  # second column is sample name
                columns = l.strip().split('\t')
                final_data.extend(columns[:6])  # direct from VCF

                ref = columns[3]
                var = columns[4]
                # protein change
                final_data.append(annotate(ref, var, columns[0], columns[1]))

                extra_info = columns[9].strip().split(':')
                final_data.extend(extra_info)
                final_data.append('\n')

        write_to_file(final_data)


# Copies column headers to new file
def init_output():
    inf = config.input_dir + 'varscan_template.txt'
    outf = config.output_dir + fname

    with open(inf, mode='r') as template, open(outf, mode='a') as output:
        for line in template:
            output.write(line)


def main():
    init_output()
    for fi in os.listdir(os.getcwd()):
        if os.path.isfile(os.getcwd() + '/' + fi) and '.vcf' in fi:
            read_vcf(fi)


if __name__ == '__main__':
    main()
