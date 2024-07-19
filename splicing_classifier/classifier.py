import pandas as pd
import numpy as np
import pybedtools as pybed
import multiprocessing as mp
from process_inputs import parse_long_read_introns_exons

# TODO: Allow for more specification in multiprocessing.
# TODO: Add more flags for different types of splicing events
# TODO: Generate more test scenarios for the classifier
# TODO: Correct flake8 error messages (somehow...)
# TODO: Implement main functionality (or remove, tbd)
# TODO: Change implementation of pybedtools to clean up files more quickly
#       Need to find implementation where they use with statement to
#       automatically clean up files (and before python exits)


def classify_reads(lr_bed_row_df, ref_exon_df, ref_intron_df):
    """
    Classify long reads based on their overlap with reference transcripts.

    Args:
    lr_bed_row_df (pd.DataFrame): A DataFrame containing long read data
    that has been intersected against the transcript_df taken from
    processing the GFF.
    ref_exon_df (pd.DataFrame): A DataFrame containing reference exon data
    generated from the GFF.
    ref_intron_df (pd.DataFrame): A DataFrame containing reference intron data
    generated from the GFF.

    Returns:
    A list containing the long read name, the most likely transcript ID,
    the splicing status, the number of introns spliced, and any flags.
    """
    # String containing all relevant flags thrown during processing
    flag_string = []
    if lr_bed_row_df.shape[0] > 1:
        # Adding information about the multiple gene overlaps and
        # lengths of gene/overlaps
        flag_string.append(
            f'multiple_gene_overlap:{",".join(lr_bed_row_df.transcript_name)}')
        flag_string.append(
            f'gene_lengths:{",".join((lr_bed_row_df.transcript_end - lr_bed_row_df.transcript_start).apply(str))}')
        flag_string.append(
            f'gene_overlap_bases:{",".join(lr_bed_row_df.overlap_bases.apply(str))}')

    # Getting the correct row from the lr_bed12 and the transcript
    # that has the greatest overlap with the long_read
    lr_bed_row = lr_bed_row_df.loc[lr_bed_row_df.overlap_bases.idxmax()]
    overlapping_transcript_name = lr_bed_row.transcript_name

    # Getting the exon and intron info for that transcript
    transcript_exons = ref_exon_df.loc[[overlapping_transcript_name]]
    transcript_introns = ref_intron_df.loc[[overlapping_transcript_name]]

    introns_exist = (transcript_exons.shape[0] > 1)

    # Entering the intron classification workflow
    if introns_exist:

        # Flagging single_intron_genes
        if transcript_introns.shape[0] == 1:
            flag_string.append('single_intron_gene')

        long_read_exons, long_read_introns = parse_long_read_introns_exons(lr_bed_row)

        # Checking assumptions about the 5-prime and 3-prime positions
        # of LR. Should contain all exons and introns between those two
        # points in the read
        lr_start = lr_bed_row['lr_start']
        lr_end = lr_bed_row['lr_end']

        # Building logic to grab expected introns and exons for a given LR
        expected_introns_df = transcript_introns.loc[
                                (transcript_introns['start'] >= lr_start) &
                                (transcript_introns['end'] <= lr_end)]

        expected_exons_df = transcript_exons.loc[
                            (transcript_exons['end'] >= lr_start) &
                            (transcript_exons['start'] <= lr_end)]

        # Checking to see if 5-prime end of LR is not in exon 1
        # TODO, check this logic with the new df naming scheme
        if np.min(expected_exons_df['exon_id'].str.split(':').apply(lambda x: int(x[-1]))) > 1:
            flag_string.append('5_prime_starts_after_first_exon')

        # If a LR doesn't cover a single intron completely
        if expected_introns_df.empty:
            splicing_status = 'single_exon_unspliced'
            flag_string.append('read_covers_no_introns')
            num_introns_spliced = 0
            return [lr_bed_row.lr_name,
                    overlapping_transcript_name,
                    splicing_status,
                    num_introns_spliced,
                    ';'.join(flag_string)]

        # If there are no LR introns, but there should be, it is unspliced
        if long_read_introns.empty:
            splicing_status = 'unspliced'
            num_introns_spliced = 0
            return [lr_bed_row.lr_name,
                    overlapping_transcript_name,
                    splicing_status,
                    num_introns_spliced,
                    ';'.join(flag_string)]

        # Bed files for intron and exon
        expected_intron_bed = pybed.BedTool.from_dataframe(expected_introns_df)
        expected_exon_bed = pybed.BedTool.from_dataframe(expected_exons_df)
        long_read_exons_bed = pybed.BedTool.from_dataframe(long_read_exons)
        long_read_introns_bed = pybed.BedTool.from_dataframe(long_read_introns)

        # We want to define a tolerance (about how many bases are required
        # for an overlap to officially occur?)
        # Current threshold is greater than 5 bases
        lr_exons_gene_introns = long_read_exons_bed.intersect(
                                expected_intron_bed, wo=True)
        lr_introns_gene_exons = long_read_introns_bed.intersect(
                                expected_exon_bed, wo=True)
        lr_introns_gene_introns = long_read_introns_bed.intersect(
                                  expected_intron_bed, wo=True)

        # Now time to define splicing status for various types of events
        lr_exons_gene_introns_df = lr_exons_gene_introns.to_dataframe(
                                    header=None,
                                    names=['lr_chrom',
                                           'lr_start',
                                           'lr_end',
                                           'lr_name',
                                           'lr_score',
                                           'lr_strand',
                                           'chrom_intron',
                                           'start_intron',
                                           'end_intron',
                                           'name_intron',
                                           'score_intron',
                                           'strand_intron',
                                           'overlap_bases'])

        lr_introns_gene_exons_df = lr_introns_gene_exons.to_dataframe(
                                   header=None,
                                   names=['lr_chrom',
                                          'lr_start',
                                          'lr_end',
                                          'lr_name',
                                          'lr_score',
                                          'lr_strand',
                                          'exon_chrom',
                                          'exon_start',
                                          'exon_end',
                                          'exon_name',
                                          'exon_score',
                                          'exon_strand',
                                          'overlap_bases'])

        lr_introns_gene_introns_df = lr_introns_gene_introns.to_dataframe(
                                     header=None,
                                     names=['lr_chrom',
                                            'lr_start',
                                            'lr_end',
                                            'lr_name',
                                            'lr_score',
                                            'lr_strand',
                                            'intron_chrom',
                                            'intron_start',
                                            'intron_end',
                                            'intron_name',
                                            'intron_score',
                                            'intron_strand',
                                            'overlap_bases'])

        # Now filtering for more than 7 overlapping bases to allow for
        # some tolerance
        if not lr_exons_gene_introns_df.empty:
            lr_exons_gene_introns_df = lr_exons_gene_introns_df.loc[
                lr_exons_gene_introns_df['overlap_bases'] > 7]

        if not lr_introns_gene_exons_df.empty:
            lr_introns_gene_exons_df = lr_introns_gene_exons_df.loc[
                lr_introns_gene_exons_df['overlap_bases'] > 7]

        if not lr_introns_gene_introns_df.empty:
            lr_introns_gene_introns_df = lr_introns_gene_introns_df.loc[
                lr_introns_gene_introns_df['overlap_bases'] > 7]

        num_introns_spliced = long_read_introns.shape[0]

        # Splicing status logic!

        # If there is a LR intron over a gene exon
        # (exon skipping or noncanonical splice sites)
        # OR if there is a LR exon over a gene intron but the number of
        # introns over introns matches the expected intron count
        # (internal intron splice site)
        if lr_introns_gene_exons_df.shape[0] > 0\
            or (lr_exons_gene_introns_df.shape[0] > 0 and
                (lr_introns_gene_introns_df.shape[0] == expected_introns_df.shape[0])):
            splicing_status = 'alternatively_spliced'
            flag_string.append('noncanonical_splice_sites_or_exon_skipping')
        elif lr_exons_gene_introns_df.shape[0] > 0:
            splicing_status = 'partially_spliced'
        elif lr_exons_gene_introns_df.empty and lr_introns_gene_exons_df.empty:
            splicing_status = 'fully_spliced'
        else:
            splicing_status = 'check_gene'

    # Workflow for single exon genes
    else:
        flag_string.append('single_exon_gene')
        # Checking for alternative splicing
        # (splicing in a single exon gene)
        if lr_bed_row['lr_blocks'][0] > 1:
            # Append the number of splicing events to the flag string
            flag_string.append(
                f'count_splicing_events:{lr_bed_row["lr_blocks"]}')
            # Also append the num introns spliced
            num_introns_spliced = lr_bed_row['lr_blocks']
            splicing_status = 'single_exon_alternatively_spliced'
        else:
            splicing_status = 'single_exon_unspliced'
            num_introns_spliced = 0

    return [lr_bed_row.lr_name,
            overlapping_transcript_name,
            splicing_status,
            num_introns_spliced,
            ';'.join(flag_string)]


def process_long_reads(lr_data_df, ref_transcript_df,
                       ref_exon_df, ref_intron_df):
    """
    Process the long read sequencing data from lr_data_df and assign splicing
    status to each read based on the reference transcript data in
    ref_transcript_df, ref_exon_df, and ref_intron_df.

    Args:
        lr_data_df (pd.DataFrame): A DataFrame containing long read sequencing
        data.
        ref_transcript_df (pd.DataFrame): A DataFrame containing reference
        transcript data.
        ref_exon_df (pd.DataFrame): A DataFrame containing reference exon data.
        ref_intron_df (pd.DataFrame): A DataFrame containing reference intron
        data.

    Returns:
        A DataFrame containing the results of the splicing analysis.
    """
    # First grabbing the gene and data beds
    lr_data_bed = pybed.BedTool.from_dataframe(lr_data_df)
    ref_transcript_bed = pybed.BedTool.from_dataframe(ref_transcript_df)

    data_gene_intersect = lr_data_bed.intersect(ref_transcript_bed,
                                                s=True,
                                                wo=True)

    # Defining the specific dfs from the gene intersect
    overlap_df = data_gene_intersect.to_dataframe(header=None,
                                                  names=['lr_chrom',
                                                         'lr_start',
                                                         'lr_end',
                                                         'lr_name',
                                                         'lr_score',
                                                         'lr_strand',
                                                         'lr_thick_start',
                                                         'lr_thick_end',
                                                         'lr_rgb',
                                                         'lr_blocks',
                                                         'lr_block_lengths',
                                                         'lr_block_starts',
                                                         'transcript_chrom',
                                                         'transcript_start',
                                                         'transcript_end',
                                                         'transcript_name',
                                                         'transcript_score',
                                                         'transcript_strand',
                                                         'overlap_bases']
                                                  ).drop_duplicates()

    # Want to check that this method of chunking does not
    # generate a lot of IO overhead
    lr_groups = list(overlap_df.groupby('lr_name'))

    # Create a multiprocessing Pool and apply the worker function to each group
    with mp.Pool() as pool:
        results = pool.starmap(classify_reads, [(lr_group,
                                                ref_exon_df,
                                                ref_intron_df)
                                                for _, lr_group in lr_groups])

    # Concatenate the results back into a single DataFrame
    classified_lr_data = pd.DataFrame(results,
                                      columns=['lr_name',
                                               'most_likely_transcript_id',
                                               'splicing_status',
                                               'num_introns_spliced',
                                               'flags'])

    # Return the classified long read data
    return classified_lr_data


def main():
    return 0


if __name__ == '__main__':
    main()
