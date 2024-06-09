import pathlib as pl
import pandas as pd


def process_gff3(path_gff3=None):
    """
    Function to process a GFF3 file and return the transcript, exon and intron information
    Args:
      path_gff3: The path to the GFF3 file to process.

    Returns:
        transcript_df: A pandas dataframe containing transcript coordinate information.
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
                                                            'attributes', 'score']]
    transcript_df['transcript_id'] = transcript_df['attributes'].str.extract(r'ID=(.*?);')

    # Now onto processing the exons and intron info
    exon_df = gff3_df.loc[gff3_df['type'] == 'CDS'][['#chrom', 'start',
                                                     'end', 'strand',
                                                     'attributes', 'score']]

    intron_df = gff3_df.loc[gff3_df['type'] == 'intron'][['#chrom', 'start',
                                                          'end', 'strand',
                                                          'attributes', 'score']]

    # Extract the parent transcript ID and exon/intron ID
    exon_df['parent_transcript_id'] = exon_df['attributes'].str.extract(r'Parent=(.*)$')
    exon_df['exon_id'] = exon_df['attributes'].str.extract(r'ID=(.*?);')

    intron_df['parent_transcript_id'] = intron_df['attributes'].str.extract(r'Parent=(.*)$')
    intron_df['intron_id'] = intron_df['attributes'].str.extract(r'ID=(.*?);')

    # Format final returned dataframes
    transcript_df = transcript_df[['#chrom', 'start',
                                   'end', 'transcript_id',
                                   'score', 'strand']]

    exon_df = exon_df[['#chrom', 'start',
                       'end', 'exon_id',
                       'score', 'strand',
                       'parent_transcript_id']].set_index('parent_transcript_id')

    intron_df = intron_df[['#chrom', 'start',
                           'end', 'intron_id',
                           'score', 'strand',
                           'parent_transcript_id']].set_index('parent_transcript_id')

    return transcript_df, exon_df, intron_df


def main():
    print('hello')
    x, y, z = process_gff3(path_gff3='tests/test_documents/S_pombe_all_chr.gff3')
    print(x)
    print(y)
    print(z)


if __name__ == "__main__":
    main()
