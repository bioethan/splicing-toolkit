import pandas as pd
import numpy as np
import pybedtools as pybed
import multiprocessing as mp
from process_inputs import parse_long_read_introns_exons


def classify_reads(lr_bed12_row, ref_exon_df, ref_intron_df):
    flag_string = []

    # First check, determine if there are more than gene for a given LR
    if lr_bed12_row.shape[0] > 1:

        # Flag string multigene overlap
        flag_string.append(
            f'multiple_gene_overlap:{",".join(lr_bed12_row["transcript_name"])}')
        flag_string.append(
            f'gene_lengths:{",".join(int(lr_bed12_row["transcript_end"]) - int(lr_bed12_row["transcript_start"]))}')
        flag_string.append(
            f'gene_overlap_bases:{",".join(lr_bed12_row["overlap_bases"].apply(str))}')

        # Final update to name_df for finding gene name
        lr_bed12_row = lr_bed12_row.iloc[[lr_bed12_row['overlap_bases'].idxmax()]].reset_index(drop=True)

    # Grabbing subsequent gene name for splicing analysis
    transcript_name = lr_bed12_row['transcript_name'][0]

    # TODO: If less than 25 bases of overlap with the gene of interest,
    # trash and move on

    # Grab intron and exon information for the transcript associated 
    # with the long-read of interest
    associated_exon_df = ref_exon_df.loc[ref_exon_df.parent_transcript_id == transcript_name]
    associated_intron_df = ref_intron_df.loc[ref_intron_df.parent_transcript_id == transcript_name]

    introns_exist = (ref_exon_df.loc[ref_exon_df.parent_transcript_id == transcript_name].shape[0] > 1)

    # Workflow for multi-exon genes
    if introns_exist:

        # Flagging single intron genes
        if associated_intron_df.shape[0] == 1:
            flag_string.append('single_intron_gene')

        long_read_exons, long_read_introns = parse_long_read_introns_exons(lr_bed12_row)

        # Checking assumptions about 5-prime and 3-prime positions of LR
        # Should contain all exons and introns between those two points in
        # the read
        lr_start = lr_bed12_row['lr_start'][0]
        lr_end = lr_bed12_row['lr_end'][0]

        # Building logic to grab expected introns and exons for a given LR
        expected_introns_df = associated_intron_df.loc[
                             (associated_intron_df['start'] >= lr_start) &
                             (associated_intron_df['end'] <= lr_end)]
        expected_exons_df = associated_exon_df.loc[
                            (associated_exon_df['end'] >= lr_start) &
                            (associated_exon_df['start'] <= lr_end)]

        # Checking to see if 5-prime end of LR is not in exon 1
        if np.min(expected_exons_df['exon_id'].str.split(':').apply(lambda x: int(x[2]))) > 1:
            flag_string.append('5_prime_starts_after_first_exon')

        # If a LR doesn't cover a single intron completely
        if expected_introns_df.empty:
            splicing_status = 'single_exon_unspliced'
            flag_string.append('read_covers_no_introns')
            num_introns_spliced = 0
            return [lr_bed12_row.iloc[0].lr_name,
                    transcript_name,
                    splicing_status,
                    num_introns_spliced,
                    ';'.join(flag_string)]

        # If there are no LR introns, but there should be, it is unspliced
        if long_read_introns.empty:
            splicing_status = 'unspliced'
            num_introns_spliced = 0
            return [lr_bed12_row.iloc[0].lr_name,
                    transcript_name,
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
                                    names=['chromLR',
                                           'startLR',
                                           'endLR',
                                           'nameLR',
                                           'scoreLR',
                                           'strandLR',
                                           'chromIntron',
                                           'startIntron',
                                           'endIntron',
                                           'nameIntron',
                                           'scoreIntron',
                                           'strandIntron',
                                           'overlapBases'])

        lr_introns_gene_exons_df = lr_introns_gene_exons.to_dataframe(
                                   header=None,
                                   names=['chromLR',
                                          'startLR',
                                          'endLR',
                                          'nameLR',
                                          'scoreLR',
                                          'strandLR',
                                          'chromExon',
                                          'startExon',
                                          'endExon',
                                          'nameExon',
                                          'scoreExon',
                                          'strandExon',
                                          'overlapBases'])

        lr_introns_gene_introns_df = lr_introns_gene_introns.to_dataframe(
                                     header=None,
                                     names=['chromLR',
                                            'startLR',
                                            'endLR',
                                            'nameLR',
                                            'scoreLR',
                                            'strandLR',
                                            'chromIntron',
                                            'startIntron',
                                            'endIntron',
                                            'nameIntron',
                                            'scoreIntron',
                                            'strandIntron',
                                            'overlapBases'])

        # Now filtering for more than 7 overlapping bases to allow for
        # some tolerance
        if not lr_exons_gene_introns_df.empty:
            lr_exons_gene_introns_df = lr_exons_gene_introns_df.loc[
                lr_exons_gene_introns_df['overlapBases'] > 7]

        if not lr_introns_gene_exons_df.empty:
            lr_introns_gene_exons_df = lr_introns_gene_exons_df.loc[
                lr_introns_gene_exons_df['overlapBases'] > 7]

        if not lr_introns_gene_introns_df.empty:
            lr_introns_gene_introns_df = lr_introns_gene_introns_df.loc[
                lr_introns_gene_introns_df['overlapBases'] > 7]

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
        if lr_bed12_row['lr_blocks'][0] > 1:

            # Append the number of splicing events to the flag string
            flag_string.append(
                f'count_splicing_events:{lr_bed12_row["lr_blocks"][0]}')

            # Also append the num introns spliced
            num_introns_spliced = lr_bed12_row['lr_blocks'][0]
            splicing_status = 'single_exon_alternatively_spliced'
        else:
            splicing_status = 'single_exon_unspliced'
            num_introns_spliced = 0

    return [lr_bed12_row.iloc[0].lr_name,
            transcript_name,
            splicing_status,
            num_introns_spliced,
            ';'.join(flag_string)]


def worker(chunk, ref_exon_df, ref_intron_df):
    """
    Worker function for multiprocessing. Applies classify_reads to each row in
    a chunk of a DataFrame.

    Args:
        chunk (pd.DataFrame): A chunk of a DataFrame containing long read data.
        ref_exon_df (pd.DataFrame): A DataFrame containing reference exon data.
        ref_intron_df (pd.DataFrame): A DataFrame containing reference intron
        data.

    Returns:
        A DataFrame containing the results of classify_reads applied to each
        row in the chunk.
    """
    return chunk.apply(classify_reads,
                       args=(ref_exon_df, ref_intron_df),
                       axis=1)


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
        results = pool.starmap(worker, [(lr_group,
                                         ref_exon_df,
                                         ref_intron_df)
                                        for _, lr_group in lr_groups])

    # Concatenate the results back into a single DataFrame
    classified_lr_data = pd.concat(results)

    # Return the classified long read data
    return classified_lr_data


def main():
    return 0


if __name__ == '__main__':
    main()
