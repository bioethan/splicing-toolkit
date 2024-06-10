from splicing_classifier.input_module.process_gff import process_gff, \
    process_gff_utrs, parse_long_read_introns_exons


def test_process_gff():
    test_url = 'test_documents/S_pombe_all_chr.gff3'
    transcript_df, exon_df, intron_df = process_gff(test_url)
    return 0


def test_process_gff_utrs():
    test_url = 'test_documents/S_pombe_all_chr.gff3'
    five_prime_utr, three_prime_utr = process_gff_utrs(test_url)
    return 0


def test_parse_long_read_introns_exons():
    test_url = 'test_documents/S_pombe_all_chr.gff3'
    lr_exons, lr_introns = parse_long_read_introns_exons(test_url)
    return 0
