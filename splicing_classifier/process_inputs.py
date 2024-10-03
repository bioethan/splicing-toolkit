import pandas as pd


class GFFProcessor:
    """
    Class to process a GFF3 file and store the transcript,
    exon and intron information in dataframe format. This includes
    information on the orginal transcript ID associated with the exons
    and introns to make splicing classification straightforward.

    Args:
      path_gff3: The path to the GFF3 file to process.
      exon_type: The type of feature to be considered as an exon.
      intron_type: The type of feature to be considered as an intron.
      transcript_type: The type of feature to be considered as a transcript.

    Attributes:
        path_gff3 (path): The path to the GFF3 file to process.
        exon_type (str): The type of feature to be considered as an exon.
        intron_type (str): The type of feature to be considered as an intron.
        transcript_type (str): The type of feature to be considered as a
        transcript.
        transcript_df (pd.DataFrame): A dataframe containing transcript
        coordinate information.
        exon_df (pd.DataFrame): A dataframe containing exon coordinate
        information.
        intron_df (pd.DataFrame): A dataframe containing intron coordinate
        information.
    """
    def __init__(self, path_gff3=None, exon_type='CDS',
                 intron_type='intron', transcript_type='mRNA'):

        # Recording the path to the GFF file and the feature keys
        self.path_gff3 = path_gff3
        self.exon_type = exon_type
        self.intron_type = intron_type
        self.transcript_type = transcript_type

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
        transcript_df = gff3_df.loc[gff3_df['type'] == transcript_type][[
            'chrom', 'start', 'end', 'strand', 'attributes', 'score']]
        transcript_df['transcript_id'] = transcript_df['attributes'].str.extract(
        r'ID=(.*?);')

        # Now onto processing the exons and intron info
        exon_df = gff3_df.loc[gff3_df['type'] == exon_type][[
            'chrom', 'start', 'end', 'strand', 'attributes', 'score']]

        intron_df = gff3_df.loc[gff3_df['type'] == intron_type][[
            'chrom', 'start', 'end', 'strand', 'attributes', 'score']]

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

        # Now storing the dataframes in the class
        self.transcript_df = transcript_df.copy()
        self.exon_df = exon_df.copy()
        self.intron_df = intron_df.copy()
        self.gff3_df = gff3_df.copy()
        self.five_prime_utr_df = None
        self.three_prime_utr_df = None

    def validate_exon_and_intron_dfs(self, correct_mode=False):
        """
        Function to validate the exon and intron dataframes generated from the
        GFF3 file and correct them if necessary. This function will check if
        the assumptions about the exon-intron relationship are correct in the
        GFF3. For every two exons, there should be an intron. If this is not
        the case, the function will add an intron entry to the intron_df and
        flag this to the user in a log file.

        Args:
            correct_mode (bool): A flag to indicate if the function should
            automatically correct the exon and intron dataframes when
            necessary.
        """
        # Gathering exon and intron counts per transcript ID
        exon_df_count = self.exon_df.groupby('parent_transcript_id').size()
        intron_df_count = self.intron_df.groupby('parent_transcript_id').size()

        temp_exon_df = self.exon_df.copy(deep=True)
        temp_intron_df = self.intron_df.copy(deep=True)

        # Now checking if each transcript has n exons and n-1 introns
        for transcript_id in exon_df_count.index:
            exon_count = exon_df_count[transcript_id]
            intron_count = intron_df_count.get(transcript_id, 0)
            # Default to 0 if transcript_id not in intron_df_count

            if exon_count != intron_count + 1:
                print(f"Discrepancy found for transcript {transcript_id}: {exon_count} exons, {intron_count} introns")
                if correct_mode:
                    # Implement correction logic here
                    print("Correcting exon and intron dataframes")
                    print("Need to implement correction logic here")
                    continue
                else:
                    print(f"""Correct mode not enabled. Removing transcript {transcript_id} from exon and intron dataframes.""")
                    temp_exon_df = temp_exon_df.drop(transcript_id)
                    temp_intron_df = temp_intron_df.drop(transcript_id)

        self.exon_df = temp_exon_df
        self.intron_df = temp_intron_df
    
    def process_utrs(self, five_prime_utr_type='five_prime_UTR', three_prime_utr_type='three_prime_UTR'):
        """
        Function to process the 5' and 3' UTR information from the GFF3 file stored in the class. Saves
        two dataframes containing the 5' and 3' UTR information respectively in the class.

        Args:
        five_prime_utr_type (str): The type of feature to be considered as a 5' UTR.
        three_prime_utr_type (str): The type of feature to be considered as a 3' UTR.

        Returns:
            five_prime_utr_df (pd.DataFrame): A dataframe containing 5' UTR
            coordinate information.
            three_prime_utr_df (pd.DataFrame): A dataframe containing 3' UTR
            coordinate information.
        """
        # Read in the 
        five_prime_utr_df = self.gff3_df.loc[self.gff3_df['type'] == five_prime_utr_type][[
                                                                    'chrom',
                                                                    'start',
                                                                    'end',
                                                                    'strand',
                                                                    'attributes',
                                                                    'score']]

        three_prime_utr_df = self.gff3_df.loc[self.gff3_df['type'] == three_prime_utr_type][[
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

        self.five_prime_utr_df = five_prime_utr_df.copy()
        self.three_prime_utr_df = three_prime_utr_df.copy()


def parse_long_read_introns_exons(long_read_bed12):
    """
    Function to parse a BED12 row containing a single long-read into its
    constituent exons and introns.

    Args:
        long_read_bed12 : A pandas dataframe containing a single row of a BED12 file.

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


if __name__ == "__main__":
    path_to_sample_gff3 = 'data/S_pombe_all_chr.gff3'
    processor = GFFProcessor(path_to_sample_gff3)
    processor.validate_exon_and_intron_dfs(correct_mode=False)
