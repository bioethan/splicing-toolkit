import pathlib as pl
import pandas as pd

# Functions to process GFF3 files and other inputs.
# Currently, this module only contains the process_gff3 and
# process_gff3_utr functions.This function reads in a GFF3 file
# and returns the transcript, exon and intron information in
# dataframe format. Future iterations could extend this module
# to include functions to process other input file types.


def process_gff(path_gff3=None):
    """
    Function to process a GFF3 file and return the transcript,
    exon and intron information in dataframe format. This includes
    information on the orginal transcript ID associated with the exons
    and introns to make splicing classification straightforward.

    Args:
      path_gff3: The path to the GFF3 file to process.

    Returns:
        transcript_df: A pandas dataframe containing transcript coordinate
        information.
        exon_df: A pandas dataframe containing exon coordinate information.
        intron_df: A pandas dataframe containing intron coordinate information.
    """
    # Read in the GFF3 file
    gff3_df = pd.read_csv(path_gff3, sep='\t', comment='#',
                          header=None,
                          names=['#chrom', 'source', 'type',
                                 'start', 'end', 'score',
                                 'strand', 'phase', 'attributes'])

    # Now process mRNA transcript info
    transcript_df = gff3_df.loc[gff3_df['type'] == 'mRNA'][['#chrom', 'start',
                                                            'end', 'strand',
                                                            'attributes',
                                                            'score']]
    transcript_df['transcript_id'] = transcript_df['attributes'].str.extract(
        r'ID=(.*?);')

    # Now onto processing the exons and intron info
    exon_df = gff3_df.loc[gff3_df['type'] == 'CDS'][['#chrom', 'start',
                                                     'end', 'strand',
                                                     'attributes', 'score']]

    intron_df = gff3_df.loc[gff3_df['type'] == 'intron'][['#chrom', 'start',
                                                          'end', 'strand',
                                                          'attributes',
                                                          'score']]

    # Extract the parent transcript ID and exon/intron ID
    exon_df['parent_transcript_id'] = exon_df['attributes'].str.extract(
        r'Parent=(.*)$')
    exon_df['exon_id'] = exon_df['attributes'].str.extract(r'ID=(.*?);')

    intron_df['parent_transcript_id'] = intron_df['attributes'].str.extract(
        r'Parent=(.*)$')
    intron_df['intron_id'] = intron_df['attributes'].str.extract(r'ID=(.*?);')

    # Format final returned dataframes
    transcript_df = transcript_df[['#chrom', 'start',
                                   'end', 'transcript_id',
                                   'score', 'strand']]

    exon_df = exon_df[['#chrom', 'start',
                       'end', 'exon_id',
                       'score', 'strand',
                       'parent_transcript_id']].set_index(
                           'parent_transcript_id')

    intron_df = intron_df[['#chrom', 'start',
                           'end', 'intron_id',
                           'score', 'strand',
                           'parent_transcript_id']].set_index(
                               'parent_transcript_id')

    return transcript_df, exon_df, intron_df


def process_gff_utrs(path_gff3=None):
    """
    Function to process a GFF3 file and return the 3' and 5' UTR
    regions in dataframe format. This is an optional function and is
    designed to be used in conjunction with the process_gff3 function,
    but especially if the UTR information is being interrogated.

    Args:
      path_gff3: The path to the GFF3 file to process.

    Returns:
        five_prime_utr_df: A pandas dataframe containing 5' UTR coordinate
        information.
        three_prime_utr_df: A pandas dataframe containing 3' UTR coordinate
        information.
    """
    # Read in the GFF3 file
    gff3_df = pd.read_csv(path_gff3, sep='\t', comment='#',
                          header=None,
                          names=['#chrom', 'source', 'type',
                                 'start', 'end', 'score',
                                 'strand', 'phase', 'attributes'])

    # Now onto processing the exons and intron info
    five_prime_utr_df = gff3_df.loc[gff3_df['type'] == 'five_prime_UTR'][[
                                                                  '#chrom',
                                                                  'start',
                                                                  'end',
                                                                  'strand',
                                                                  'attributes',
                                                                  'score']]

    three_prime_utr_df = gff3_df.loc[gff3_df['type'] == 'three_prime_UTR'][[
                                                                '#chrom',
                                                                'start',
                                                                'end',
                                                                'strand',
                                                                'attributes',
                                                                'score']]

    # Extract the parent transcript ID and UTR IDs
    five_prime_utr_df['parent_transcript_id'] = five_prime_utr_df[
        'attributes'].str.extract(r'Parent=(.*)$')
    five_prime_utr_df['five_prime_utr_id'] = five_prime_utr_df[
        'attributes'].str.extract(r'ID=(.*?);')

    three_prime_utr_df['parent_transcript_id'] = three_prime_utr_df[
        'attributes'].str.extract(r'Parent=(.*)$')
    three_prime_utr_df['three_prime_utr_id'] = three_prime_utr_df[
        'attributes'].str.extract(r'ID=(.*?);')

    # Format final returned UTR dataframes
    five_prime_utr_df = five_prime_utr_df[['#chrom', 'start',
                                           'end', 'five_prime_utr_id',
                                           'score', 'strand',
                                           'parent_transcript_id']].set_index('parent_transcript_id')

    three_prime_utr_df = three_prime_utr_df[['#chrom', 'start',
                                             'end', 'three_prime_utr_id',
                                             'score', 'strand',
                                             'parent_transcript_id']].set_index('parent_transcript_id')

    return five_prime_utr_df, three_prime_utr_df


# TODO: Write code to parse out BED12 rows into the
# consistuent exons and introns for splicing classification
def parse_long_read_introns_exons(long_read_bed12):
    """
    Function to parse a BED12 row containing a single long-read into its
    constituent exons and introns.

    Args:
        long_read_bed12: A pandas dataframe containing a single row of a BED12
        file.

    Returns:
        long_read_exons: A pandas dataframe containing information on the exons
        present in the current long-read data.
        long_read_introns: A pandas dataframe containing information on the
        introns present in the current long-read data.
    """

    columns = long_read_bed12.iloc[0].values.flatten().tolist()
    chrom = columns[0]
    start = int(columns[1])
    name = columns[3]
    strand = columns[5]
    block_count = int(columns[9])
    block_sizes = list(map(int, columns[10].split(',')))
    block_starts = list(map(int, columns[11].split(',')))

    # Generating exon start and end coordinates
    exon_starts = [start + block_starts[i] for i in range(block_count)]
    exon_ends = [exon_starts[i] + block_sizes[i] for i in range(block_count)]

    # Numbering for exons based on the strand (back to front or front to back)
    exon_numbers = [i + 1 if strand == '+' else block_count - i
                    for i in range(block_count)]
    exon_names = [f"{name}_exon_{exon_numbers[i]}" for i in range(block_count)]

    # Final exon list for the current long-read
    exons = [[chrom,
              exon_starts[i],
              exon_ends[i],
              exon_names[i],
              ".",
              strand]
             for i in range(block_count)]

    # Generating intron start and end coordinates
    intron_starts = [exon_ends[i] for i in range(block_count - 1)]
    intron_ends = [start + block_starts[i + 1] for i in range(block_count - 1)]

    # Numbering for introns between exons based on the strand
    intron_numbers = [i + 1 if strand == '+' else block_count - i - 1 for i in
                      range(block_count - 1)]
    intron_names = [f"{name}_intron_{intron_numbers[i]}" for i in
                    range(block_count - 1)]

    # Final intron list for the current long-read
    introns = [[chrom,
                intron_starts[i],
                intron_ends[i],
                intron_names[i],
                ".",
                strand] for i in range(block_count - 1)]

    exon_df = pd.DataFrame(data=(exons),
                           columns=['chrom',
                                    'start',
                                    'end',
                                    'name',
                                    'score',
                                    'strand'])
    intron_df = pd.DataFrame(data=(introns),
                             columns=['chrom',
                                      'start',
                                      'end',
                                      'name',
                                      'score',
                                      'strand'])

    return exon_df, intron_df


def main():
    return 0


if __name__ == "__main__":
    main()
