import pathlib as pl
import pandas as pd


def gtf_to_bed12(path_gtf=None, out_dir="~/tmp"):
    """
    Converts the genes present in a GTF file to a BED12 for easy processing of
    transcriptomic information.

    Inputs:
    path_gtf: A path to the GTF file for generation of the bed12
        Input: str

    path_out: Output directory for the resulting gene, exon, and intron bed
    files.
        Input: str
    """

    # Making the output directory
    out_path = pl.Path(out_dir)
    out_path.mkdir(mode=0o755, exist_ok=True)

    # Grabbing the contents of the GTF in a pd df
    gtf_contents = pd.read_csv(path_gtf, comment='#', sep='\t', header=0,
                               names=['chrom', 'source', 'feature', 'start',
                                      'end', 'score', 'strand', 'frame',
                                      'attributes'])

    # Gathering list of genes from the attributes column

    # List of genes
    coding_genes = gtf_contents[gtf_contents['feature'] == 'gene']
    coding_genes = coding_genes[coding_genes['attributes']
                                .str.contains('protein_coding')]
    coding_genes['gene_id'] = coding_genes['attributes']\
        .str.extract(r'gene_id "(.*?)";', expand=False)

    # List of exons
    coding_exons = gtf_contents[gtf_contents['feature'] == 'exon']
    coding_exons = coding_exons[coding_exons['attributes']
                                .str.contains('protein_coding')]
    coding_exons['gene_id'] = coding_exons['attributes']\
        .str.extract(r'gene_id "(.*?)";', expand=False)
    coding_exons['exon_num'] = coding_exons['attributes']\
        .str.extract(r'exon_number "(\d+)"')
    coding_exons['exon_name'] = coding_exons['gene_id'] \
        + '_exon_' \
        + coding_exons['exon_num']

    coding_exons['exon_num'] = coding_exons['exon_num'].apply(int)
    coding_exons['length'] = coding_exons['end'] - coding_exons['start']
    coding_exons['length'] = coding_exons['length'].apply(str)
    coding_exons['start'] = coding_exons['start'].apply(str)
    coding_exons['end'] = coding_exons['end'].apply(str)

    # Generatin coding exon information
    coding_exons_info = coding_exons.groupby(by=['chrom',
                                                 'gene_id',
                                                 'strand']).agg({'exon_num': 'max',
                                                                 'start': ','.join,
                                                                 'end': ','.join,
                                                                 'length': ','.join}).reset_index()

    coding_genes_info = coding_genes.join(coding_exons_info.set_index('gene_id'), on='gene_id', rsuffix= '_exons')
    coding_genes_info['block_starts'] = coding_genes_info['start_exons'].str.split(',').apply(lambda x: [int(l) for l in x])
    coding_genes_info['block_starts'] = coding_genes_info.apply(lambda row: [x - int(row['start']) for x in row['block_starts']], axis=1)
    coding_genes_info['block_starts'] = coding_genes_info.apply(lambda row: [list(reversed(row['block_starts'])) if row['strand'] == '-' else row['block_starts']], axis=1)
    coding_genes_info['length'] = coding_genes_info.apply(lambda row: row['length'].split(','), axis=1)
    coding_genes_info['length'] = coding_genes_info.apply(lambda row: [list(reversed(row['length'])) if row['strand'] == '-' else row['length']], axis=1)
    coding_genes_info['block_starts'] = coding_genes_info.apply(lambda x: ','.join([str(i) for i in x['block_starts'][0]]), axis=1)
    coding_genes_info['length'] = coding_genes_info.apply(lambda x: ','.join([str(i) for i in x['length'][0]]), axis=1)

    # Exon bed file
    exons_bed12 = pd.DataFrame({
        'chrom': coding_exons['chrom'],
        'start': coding_exons['start'],
        'end': coding_exons['end'],
        'name': coding_exons['exon_name'],
        'score': coding_exons['score'],
        'strand': coding_exons['strand'],
    })

    genes_bed12 = pd.DataFrame({
        'chrom': coding_genes_info['chrom'],
        'start': coding_genes_info['start'],
        'end': coding_genes_info['end'],
        'name': coding_genes_info['gene_id'],
        'score': coding_genes_info['score'],
        'strand': coding_genes_info['strand'],
        'thick_start': coding_genes_info['start'],
        'thick_end': coding_genes_info['end'],
        'item_rgb': '0,0,255',
        'block_count': coding_genes_info['exon_num'],
        'block_sizes': coding_genes_info['length'],
        'block_starts': coding_genes_info['block_starts'],
    })

    # Generating intron coords with a special function
    intron_df = coding_exons_info[coding_exons_info['exon_num'] > 1]\
        .reset_index()
    intron_bed12 = intron_coordinate_generator(intron_df['exon_num'],
                                               intron_df['gene_id'],
                                               intron_df['chrom'],
                                               intron_df['start'],
                                               intron_df['end'],
                                               intron_df['strand'])

    # Opening exon and general gene files in
    # the outdir(out_path / 'GTF_genes.bed')
    exons_bed12.to_csv((out_path / 'GTF_exons.bed'), sep='\t', index=False,
                       header=False)
    genes_bed12.to_csv((out_path / 'GTF_genes.bed'), sep='\t', index=False,
                       header=False)
    intron_bed12.to_csv((out_path / 'GTF_introns.bed'), sep='\t', index=False,
                        header=False)

    return 0


# returns a df containing in bed12 format intron coords
def intron_coordinate_generator(num_exons,
                                gene_id,
                                chroms,
                                start_coords,
                                end_coords,
                                strand):
    """
    Tool for generating the intron dataframe, separated out in part because of
    the increased complexity required in calculating the coordinates of the
    intron bed file from the GTF

    Inputs:
    All of the following inputs must be the same size and each row must
    correspond to the same feature row in the GTF

    num_exons: A list/series of exon counts
    gene_id A list of gene_ids of
    chroms: A list of chrom identifiers
    start_coords: A list of exon start coords
    end_coords: A list of exon end coordinates
    strand: A list of exon strandedness values

    """
    # Starting info
    num_genes = len(num_exons)
    intron_starts = []
    intron_ends = []
    intron_strands = []
    intron_names = []
    intron_chroms = []
    intron_lengths = []

    # Looping through all genes
    for gene_num in range(num_genes):
        exon_starts = start_coords[gene_num].split(',')
        exon_ends = end_coords[gene_num].split(',')
        for exon_num in range(num_exons[gene_num] - 1):
            intron_chrom = chroms[gene_num]
            intron_name = gene_id[gene_num] + f'_intron{exon_num + 1}'

            # If positively stranded
            if strand[gene_num] == '+':
                intron_name = gene_id[gene_num] + f'_intron{exon_num + 1}'
                intron_start = int(exon_ends[exon_num])
                intron_end = int(exon_starts[exon_num + 1])
            # If negatively stranded
            else:

                intron_end = int(exon_starts[exon_num])
                intron_start = int(exon_ends[exon_num + 1])

            intron_length = intron_end - intron_start
            intron_starts.append(intron_start)
            intron_ends.append(intron_end)
            intron_chroms.append(intron_chrom)
            intron_names.append(intron_name)
            intron_strands.append(strand[gene_num])
            intron_lengths.append(intron_length)

    intron_df = pd.DataFrame({
        'chrom': intron_chroms,
        'start': intron_starts,
        'end': intron_ends,
        'name': intron_names,
        'score': '.',
        'strand': intron_strands,
    })

    return intron_df


# I want to process a gff3 file and return three pandas dataframes
# The gene_df is the most basic, and essentially looks like a bed file of the
# gene information
# The exon_df and intron_dfs are more complex, but effectively contain the
# coordinates of the exons and introns of the genes in the gff3 file,
# but with the index corresponding to the gene_id
def process_gff3(path_gff3=None):
    # Read in the GFF3 file
    gff3_df = pd.read_csv(path_gff3, sep='\t', comment='#',
                          header=None,
                          names=['#chrom', 'source', 'type',
                                 'start', 'end', 'score',
                                 'strand', 'phase', 'attributes'])

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

    exon_df['parent_transcript_id'] = exon_df['attributes'].str.extract(r'Parent=(.*)$')
    exon_df['exon_id'] = exon_df['attributes'].str.extract(r'ID=(.*?);')

    intron_df['parent_transcript_id'] = intron_df['attributes'].str.extract(r'Parent=(.*)$')
    intron_df['intron_id'] = intron_df['attributes'].str.extract(r'ID=(.*?);')

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
    x, y, z = process_gff3(path_gff3='tests/test_documents/Schizosaccharomyces_pombe_all_chromosomes-2.gff3')
    print(x)
    print(y)
    print(z)


if __name__ == "__main__":
    main()
