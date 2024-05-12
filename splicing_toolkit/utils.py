import pandas as pd
import numpy as np
import pathlib as pl
from tqdm import tqdm
import pybedtools as pybed


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
    coding_genes = coding_genes[coding_genes['attributes'].str.contains('protein_coding')]
    coding_genes['gene_id'] = coding_genes['attributes'].str.extract(r'gene_id "(.*?)";', expand=False)

    # List of exons
    coding_exons = gtf_contents[gtf_contents['feature'] == 'exon']
    coding_exons = coding_exons[coding_exons['attributes'].str.contains('protein_coding')]
    coding_exons['gene_id'] = coding_exons['attributes'].str.extract(r'gene_id "(.*?)";', expand=False)
    coding_exons['exon_num'] = coding_exons['attributes'].str.extract(r'exon_number "(\d+)"')
    coding_exons['exon_name'] = coding_exons['gene_id'] + '_exon_' + coding_exons['exon_num']
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
    coding_genes_info['block_starts'] = coding_genes_info['start_exons'].str\
        .split(',').apply(lambda x: [int(l) for l in x])
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


def parse_bed12_introns_exons(bed12_row):
    columns = bed12_row.iloc[0].values.flatten().tolist()
    chrom = columns[0]
    start = int(columns[1])
    name = columns[3]
    strand = columns[5]
    block_count = int(columns[9])
    block_sizes = list(map(int, columns[10].split(',')))
    block_starts = list(map(int, columns[11].split(',')))

    exons = []
    introns = []

    for i in range(block_count):
        exon_start = start + block_starts[i]
        exon_end = exon_start + block_sizes[i]

        # Numbering for exons based on the strand
        if strand == '+':
            exon_number = i + 1
        else:
            exon_number = block_count - i

        exons.append([chrom, exon_start, exon_end,
                      f"{name}_exon_{exon_number}", ".", strand])

        if i < block_count - 1:
            intron_start = exon_end
            intron_end = start + block_starts[i + 1]

            # Numbering for introns between exons based on the strand
            if strand == '+':
                intron_number = i + 1
            else:
                intron_number = (block_count - i - 1)

            introns.append([chrom, intron_start, intron_end,
                            f"{name}_intron_{intron_number}", ".", strand])

    exon_df = pd.DataFrame(data=(exons), columns=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    intron_df = pd.DataFrame(data=(introns), columns=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    return exon_df, intron_df


def splicing_classifier(LRS_url, output_save):
    # First grabbing the gene and data beds
    data_file = pybed.BedTool(LRS_url)
    genes_file = pybed.BedTool('/Users/ethanbrown/Documents/splicing_toolkit/splicing_toolkit/tmp/GTF_genes.bed')
    intron_file = pd.read_csv('/Users/ethanbrown/Documents/splicing_toolkit/splicing_toolkit/tmp/GTF_introns.bed', sep='\t', names=['chromIntron', 'startIntron', 'endIntron', 'nameIntron', 'scoreIntron', 'strandIntron'])
    exon_file = pd.read_csv('/Users/ethanbrown/Documents/splicing_toolkit/splicing_toolkit/tmp/GTF_exons.bed', sep='\t', names=['chromExon', 'startExon', 'endExon', 'nameExon', 'scoreExon', 'strandExon'])
    data_gene_intersect = data_file.intersect(genes_file, s=True, wo=True)

    # Output file var
    output_list = []

    # Defining the specific dfs from the gene intersect
    df = data_gene_intersect.to_dataframe(header=None, names=["chrom",
                                                              "start",
                                                              "end",
                                                              "name",
                                                              "score",
                                                              "strand",
                                                              'thickStart',
                                                              'thickEnd',
                                                              'rgb',
                                                              'blocks',
                                                              'blockLengths',
                                                              'blockStarts',
                                                              'geneChrom',
                                                              'geneStart',
                                                              'geneEnd',
                                                              'geneName',
                                                              'geneScore',
                                                              'geneStrand',
                                                              'geneThickStart',
                                                              'geneThickEnd',
                                                              'geneRgb',
                                                              'geneBlocks',
                                                              'geneBlockLengths',
                                                              'geneBlockStarts',
                                                              'overlapBases'])

    # overlap_df drops most of the extra info about the genes to preserve
    # geneName and basic geneinfo
    overlap_df = df[['chrom', 'start', 'end', 'name', 'score', 'strand',
                     'thickStart', 'thickEnd', 'rgb', 'blocks', 'blockLengths',
                     'blockStarts', 'geneName', 'geneStart', 'geneEnd',
                     'geneBlocks', 'geneBlockStarts', 'geneBlockLengths',
                     'overlapBases']]

    # data_df contains information about the reads that map to specific genes,
    # effectively a bed12 of those overlapping LRs
    data_df = df[['chrom', 'start', 'end', 'name', 'score', 'strand',
                  'thickStart', 'thickEnd', 'rgb', 'blocks', 'blockLengths',
                  'blockStarts']].drop_duplicates()

    # Enter the splicing analysis workflow
    print(LRS_url)
    for i in tqdm(data_df['name']):

        # Starting to create flags for the final output
        flag_string = []
        name_df = overlap_df.loc[overlap_df['name'] == i].reset_index(drop=True)
        long_read_bed_df = data_df.loc[data_df['name'] == i].reset_index(drop=True)

        # First check, determine if there are more than gene for a given LR
        if name_df.shape[0] > 1:

            # TODO: Decide what metric determines a LR belongs to a gene,
            # straight overlap or whichever gene has a 5' start site closest
            # to the LR 5' start site?
            # Current method is gene with largest amount of overlap

            # Flag string multigene overlap
            flag_string.append(f'multiple_gene_overlap:{",".join(name_df["geneName"])}')
            flag_string.append(f'gene_lengths:{",".join((name_df["geneEnd"] - name_df["geneStart"]).apply(str))}')
            flag_string.append(f'gene_overlap_bases:{",".join(name_df["overlapBases"].apply(str))}')

            # Final update to name_df for finding gene name
            name_df = name_df.iloc[[name_df['overlapBases'].idxmax()]].reset_index(drop=True)

        # Grabbing subsequent gene name for splicing analysis
        gene_name = name_df['geneName'][0]

        # If less than 25 bases of overlap with the gene of interest, trash and move on
        if name_df['overlapBases'][0] < 25:
            continue

        # Check to see if gene has more than one block (i.e. exon)
        introns_exist = (name_df['geneBlocks'][0] > 1)

        # Workflow for multi-exon genes
        if introns_exist:

            # Generating intron and exon dfs for splicing classification
            intron_df = intron_file.loc[intron_file['nameIntron'].str.contains(gene_name)]
            exon_df = exon_file.loc[exon_file['nameExon'].str.contains(gene_name)]

            # Flagging single intron genes
            if intron_df.shape[0] == 1:
                flag_string.append('single_intron_gene')

            long_read_exons, long_read_introns = parse_bed12_introns_exons(long_read_bed_df)

            # Checking assumptions about 5-prime and 3-prime positions of LR
            # Should contain all exons and introns between those two points in
            # the read
            lr_start = long_read_bed_df['start'][0]
            lr_end = long_read_bed_df['end'][0]

            # Building logic to grab expected introns and exons for a given LR
            expected_introns_df = intron_df.loc[(intron_df['startIntron'] >= lr_start) & (intron_df['endIntron'] <= lr_end)]
            expected_exons_df = exon_df.loc[(exon_df['endExon'] >= lr_start) & (exon_df['startExon'] <= lr_end)]

            # Checking to see if 5-prime end of LR is not in exon 1
            if np.min(expected_exons_df['nameExon'].str.split('_').apply(lambda x: int(x[2]))) > 1:
                flag_string.append('5_prime_starts_after_first_exon')

            # If a LR doesn't cover a single intron completely
            if expected_introns_df.empty:
                splicing_status = 'single_exon_unspliced'
                flag_string.append('read_covers_no_introns')
                num_introns_spliced = 0
                # TODO add an output then break here
                # OUTPUT i, gene_name, splicing_status, none, flag_string
                output_list.append([i, gene_name, splicing_status, num_introns_spliced, ';'.join(flag_string)])
                continue

            # If there are no LR introns, but there should be, it is unspliced
            if long_read_introns.empty:
                splicing_status = 'unspliced'
                num_introns_spliced = 0
                # TODO add an output then continue here
                # OUTPUT i, gene_name, splicing_status, none, flag_string
                output_list.append([i, gene_name, splicing_status, num_introns_spliced, ';'.join(flag_string)])
                continue

            # Bed files for intron and exon
            expected_intron_bed = pybed.BedTool.from_dataframe(intron_df)
            expected_exon_bed = pybed.BedTool.from_dataframe(exon_df)
            long_read_exons_bed = pybed.BedTool.from_dataframe(long_read_exons)
            long_read_introns_bed = pybed.BedTool.from_dataframe(long_read_introns)

            # We want to define a tolerance (about how many bases are required
            # for an overlap to officially occur?)
            # Current threshold is greater than 5 bases
            lr_exons_gene_introns = long_read_exons_bed.intersect(expected_intron_bed, wo=True)
            lr_introns_gene_exons = long_read_introns_bed.intersect(expected_exon_bed, wo=True)
            lr_introns_gene_introns = long_read_introns_bed.intersect(expected_intron_bed, wo=True)

            # Now time to define splicing status for various types of events
            #lr_exons_gene_exons_df = lr_exons_gene_exons.to_dataframe(header=None, names=["chromLR", "startLR", "endLR", "nameLR", "scoreLR", "strandLR", 'chromExon', 'startExon', 'endExon', 'nameExon', 'scoreExon', 'strandExon', 'overlapBases'])
            lr_exons_gene_introns_df = lr_exons_gene_introns.to_dataframe(header=None, names=["chromLR", "startLR", "endLR", "nameLR", "scoreLR", "strandLR", 'chromIntron', 'startIntron', 'endIntron', 'nameIntron', 'scoreIntron', 'strandIntron', 'overlapBases'])
            lr_introns_gene_exons_df = lr_introns_gene_exons.to_dataframe(header=None, names=["chromLR", "startLR", "endLR", "nameLR", "scoreLR", "strandLR", 'chromExon', 'startExon', 'endExon', 'nameExon', 'scoreExon', 'strandExon', 'overlapBases'])
            lr_introns_gene_introns_df = lr_introns_gene_introns.to_dataframe(header=None, names=["chromLR", "startLR", "endLR", "nameLR", "scoreLR", "strandLR", 'chromIntron', 'startIntron', 'endIntron', 'nameIntron', 'scoreIntron', 'strandIntron', 'overlapBases'])

            # Now filtering for more than 7 overlapping bases to allow for some tolerance
            #lr_exons_gene_exons_df = lr_exons_gene_exons_df.loc[lr_exons_gene_exons_df['overlapBases'] > 7]
            if not lr_exons_gene_introns_df.empty:
                lr_exons_gene_introns_df = lr_exons_gene_introns_df.loc[lr_exons_gene_introns_df['overlapBases'] > 7]

            if not lr_introns_gene_exons_df.empty:
                lr_introns_gene_exons_df = lr_introns_gene_exons_df.loc[lr_introns_gene_exons_df['overlapBases'] > 7]

            if not lr_introns_gene_introns_df.empty:
                lr_introns_gene_introns_df = lr_introns_gene_introns_df.loc[lr_introns_gene_introns_df['overlapBases'] > 7]

            num_introns_spliced = long_read_introns.shape[0]

            # Splicing status logic!

            # If there is a LR intron over a gene exon
            # (exon skipping or noncanonical splice sites)
            # OR if there is a LR exon over a gene intron but the number of
            # introns over introns matches the expected intron count
            # (internal intron splice site)
            if lr_introns_gene_exons_df.shape[0] > 0 \
               or (lr_exons_gene_introns_df.shape[0] > 0
               and (lr_introns_gene_introns_df.shape[0] == expected_introns_df.shape[0])):
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
            # Checking for alternative splicing (splicing in a single exon gene)
            if long_read_bed_df['blocks'][0] > 1:
                flag_string.append(f'count_splicing_events:{long_read_bed_df["blocks"][0]}')
                num_introns_spliced = long_read_bed_df['blocks'][0]
                splicing_status = 'single_exon_alternatively_spliced'
            else:
                splicing_status = 'single_exon_unspliced'
                num_introns_spliced = 0

        output_list.append([i, gene_name, splicing_status, num_introns_spliced, ';'.join(flag_string)])

    final_df = pd.DataFrame(output_list, columns=['LR_name', 'gene_name', 'splicing_status', 'num_introns', 'flagged_reasons'])
    final_df.to_csv(output_save, index=False, sep='\t')
