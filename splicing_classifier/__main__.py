import logging
import argparse
import pathlib as pl

# TODO: Implement pathlib to handle file creation and pathing

logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(description='Toolkit designed to help with \
                                              quick splicing classification \
                                              and analysis of long-read \
                                              sequencing data. More \
                                              information in the GitHub \
                                              readme.')

parser.add_argument('--input-file', '-i', type=str, required=True,
                    dest='input_file_path',
                    help='Input long-read bed/bam file')
parser.add_argument('--output-dir', '-o', type=str, required=True,
                    dest='output_dir_path',
                    help='')
parser.add_argument('--genome_annotation', '-g', type=str, required=True,
                    dest='genome_anno_path',
                    help='GFF3 file containing information on the genome to\
                          which the long-read input is aligned')
parser.add_argument('--mode', '-m', type=str, required=True,
                    dest='action_mode',
                    help='What mode of should the splicing classifier engage')

args = parser.parse_args()

print(args.input_file_path)
