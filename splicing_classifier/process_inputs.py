import pathlib as pl
import pandas as pd

# TODO: Implement path handling for the GFF3 file location
# TODO: Implement validation for GFF3 file (if you have two distinct exons
# in a gene, you necessarily have an intron). If the GFF3 file is not
# formatted correctly, raise an error/correct the error by adding an intron to
# the ref_intron_df, flag to the user in some log file, and continue. This is
# being done in the validate_and_correct_exon_and_intron_dfs function.


def process_gff(path_gff3=None):
    """
    Function to process a GFF3 file and return the transcript,
    exon and intron information in dataframe format. This includes
    information on the orginal transcript ID associated with the exons
    and introns to make splicing classification straightforward.

    Args:
      path_gff3: The path to the GFF3 file to process.

    Returns:
        transcript_df (pd.DataFrame): A dataframe containing transcript
        coordinate information.
        exon_df (pd.DataFrame): A dataframe containing exon coordinate
        information.
        intron_df (pd.DataFrame): A dataframe containing intron coordinate
        information.
    """
    # Read in the GFF3 file
    gff3_df = pd.read_csv(path_gff3, sep='\t', comment='#',
                          header=None,
                          names=['chrom',
                                 'source',
                                 'type',
                                 'start',
                                 'end',
                                 'score',
                                 'strand',
                                 'phase',
                                 'attributes'])

    # Now process mRNA transcript info
    transcript_df = gff3_df.loc[gff3_df['type'] == 'mRNA'][['chrom',
                                                            'start',
                                                            'end',
                                                            'strand',
                                                            'attributes',
                                                            'score']]
    transcript_df['transcript_id'] = transcript_df['attributes'].str.extract(
        r'ID=(.*?);')

    # Now onto processing the exons and intron info
    exon_df = gff3_df.loc[gff3_df['type'] == 'CDS'][['chrom',
                                                     'start',
                                                     'end',
                                                     'strand',
                                                     'attributes',
                                                     'score']]

    intron_df = gff3_df.loc[gff3_df['type'] == 'intron'][['chrom',
                                                          'start',
                                                          'end',
                                                          'strand',
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
    transcript_df = transcript_df[['chrom',
                                   'start',
                                   'end',
                                   'transcript_id',
                                   'score',
                                   'strand']]

    exon_df = exon_df[['chrom',
                       'start',
                       'end',
                       'exon_id',
                       'score',
                       'strand',
                       'parent_transcript_id']].set_index(
                           'parent_transcript_id')

    intron_df = intron_df[['chrom',
                           'start',
                           'end',
                           'intron_id',
                           'score',
                           'strand',
                           'parent_transcript_id']].set_index(
                               'parent_transcript_id')

    return transcript_df, exon_df, intron_df


def process_gff_utrs(path_gff3=None):
    """
    Function to process a GFF3 file and return the 3' and 5' UTR
    regions in dataframe format. This is an optional function and is
    designed to be used in conjunction with the process_gff3 function,
    especially if the UTR information is being interrogated.

    Args:
      path_gff3: The path to the GFF3 file to process.

    Returns:
        five_prime_utr_df (pd.DataFrame): A dataframe containing 5' UTR
        coordinate information.
        three_prime_utr_df (pd.DataFrame): A dataframe containing 3' UTR
        coordinate information.
    """
    # Read in the GFF3 file
    gff3_df = pd.read_csv(path_gff3, sep='\t', comment='#',
                          header=None,
                          names=['chrom',
                                 'source',
                                 'type',
                                 'start',
                                 'end',
                                 'score',
                                 'strand',
                                 'phase',
                                 'attributes'])

    # Now onto processing the exons and intron info
    five_prime_utr_df = gff3_df.loc[gff3_df['type'] == 'five_prime_UTR'][[
                                                                  'chrom',
                                                                  'start',
                                                                  'end',
                                                                  'strand',
                                                                  'attributes',
                                                                  'score']]

    three_prime_utr_df = gff3_df.loc[gff3_df['type'] == 'three_prime_UTR'][[
                                                                'chrom',
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
    five_prime_utr_df = five_prime_utr_df[['chrom',
                                           'start',
                                           'end',
                                           'five_prime_utr_id',
                                           'score',
                                           'strand',
                                           'parent_transcript_id']
                                          ].set_index('parent_transcript_id')

    three_prime_utr_df = three_prime_utr_df[['chrom',
                                             'start',
                                             'end',
                                             'three_prime_utr_id',
                                             'score',
                                             'strand',
                                             'parent_transcript_id']
                                            ].set_index('parent_transcript_id')

    return five_prime_utr_df, three_prime_utr_df


def parse_long_read_introns_exons(long_read_bed12):
    """
    Function to parse a BED12 row containing a single long-read into its
    constituent exons and introns.

    Args:
        long_read_bed12 : A pandas dataframe containing a single row of a BED12
        file.

    Returns:
        long_read_exons (pd.DataFrame): A dataframe containing information on
        the exons present in the current long-read data.
        long_read_introns (pd.DataFrame): A dataframe containing information
        on the introns present in the current long-read data.
    """

    # May need to alter this depending on what I pass in as the long read row
    chrom = long_read_bed12['lr_chrom']
    start = long_read_bed12['lr_start']
    name = long_read_bed12['lr_name']
    strand = long_read_bed12['lr_strand']
    block_count = int(long_read_bed12['lr_blocks'])
    block_sizes = list(map(int, long_read_bed12['lr_block_lengths'].split(',')))
    block_starts = list(map(int, long_read_bed12['lr_block_starts'].split(',')))

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

def validate_and_correct_exon_and_intron_dfs(exon_df, intron_df):
    '''
    Function to validate the exon and intron dataframes generated from the
    GFF3 file and correct them if necessary. This function will check if
    the assumptions about the exon-intron relationship are correct in the GFF3.
    For every two exons, there should be an intron. If this is not the case,
    the function will add an intron entry to the intron_df and flag this to the
    user in a log file.

    Args:
        exon_df (pd.DataFrame): A dataframe containing exon coordinate
        information from a GFF3.
        intron_df (pd.DataFrame): A dataframe containing intron coordinate
        information from a GFF3.

    Returns:
        exon_df (pd.DataFrame): A dataframe containing exon coordinate
        information with any necessary corrections.
        intron_df (pd.DataFrame): A dataframe containing intron coordinate
        information with any necessary corrections.
    '''

    return 0


def main():
    return 0


if __name__ == "__main__":
    main()
