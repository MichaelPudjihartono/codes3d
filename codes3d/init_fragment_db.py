import pybedtools
import pandas as pd
import sys
from tqdm import tqdm
import argparse
import os
import time
from sqlalchemy import create_engine
import codes3d


def find_variant_fragments(
        query_fp,
        fragment_fp,
        output_fp,
        db_url,
        table):
    """ Create table of variant fragment IDs.
    """
    fragment_bed = pybedtools.BedTool(fragment_fp)
    desc = 'Finding variant fragments...'
    bar_format = '{desc}: {n_fmt} {unit}'
    t = tqdm(total=0, unit='variants', desc=desc, disable=False,
             bar_format=bar_format)
    chunksize = 200000
    cols = ['chr', 'frag_id', 'id', 'variant_id']
    pybedtools.set_tempdir(os.path.dirname(output_fp))
    db = create_engine(db_url, echo=False)
    idx = 1
    for query_df in pd.read_csv(
            query_fp, sep='\t', compression='gzip', chunksize=chunksize,
            usecols=['chr', 'variant_pos', 'variant_id']):
        query_df = query_df.rename(columns={'variant_pos': 'end'})
        query_df['id'] = range(idx, idx+len(query_df))
        query_df['start'] = query_df['end'] - 1
        query_bed = pybedtools.BedTool.from_dataframe(
            query_df[['chr', 'start', 'end', 'id', 'variant_id']])
        df = fragment_bed.intersect(query_bed, wb=True)
        df = df.to_dataframe(
            names=['frag_chr', 'frag_start', 'frag_end', 'frag_id',
                   'chrom', 'start', 'end', 'id', 'variant_id'])
        cols = ['frag_id', 'chrom', 'id']
        mode = 'w' if t.total == 0 else 'a'
        header = 'True' if t.total == 0 else None
        df[cols].to_csv(output_fp, sep='\t', header=header, mode=mode)
        if table:
            if_exists = 'replace' if t.total == 0 else 'append'
            df[cols].to_sql(table, con=db, if_exists=if_exists, index=False)
        idx += len(query_df)
        t.total += len(query_df)
        t.update(len(query_df))
    t.close()
    pybedtools.cleanup(remove_all=True)
    if not table:
        return
    create_index(table, db)


def find_gene_fragments(
        query_fp,
        fragment_fp,
        output_fp,
        db_url,
        table):
    """ Create table of gene fragment IDs.
    """
    fragment_bed = pybedtools.BedTool(fragment_fp)
    desc = 'Finding gene fragments...'
    bar_format = '{desc}: {n_fmt} {unit}'
    t = tqdm(total=0, unit='genes', desc=desc, disable=False,
             bar_format=bar_format)
    chunksize = 2000
    pybedtools.set_tempdir(os.path.dirname(output_fp))
    db = create_engine(db_url, echo=False)
    for query_df in pd.read_csv(
            query_fp, sep='\t', compression='infer', chunksize=chunksize
    ):
        #query_df.columns = ['id', 'chr', 'start', 'end', 'gene', 'gencode_id']
        query_bed = pybedtools.BedTool.from_dataframe(
            query_df[['chrom', 'start', 'end', 'id', 'gencode_id']])
        #    print(query_bed)
        df = fragment_bed.intersect(query_bed, wb=True)
        df = df.to_dataframe(
            names=['frag_chr', 'frag_start', 'frag_end', 'frag_id',
                   'chrom', 'start', 'end', 'id', 'gencode_id'])
        cols = ['frag_id', 'chrom', 'id']
        mode = 'w' if t.total == 0 else 'a'
        header = 'True' if t.total == 0 else None
        df[cols].to_csv(output_fp, sep='\t', header=header,
                        mode=mode, index=False)
        if table:
            if_exists = 'replace' if t.total == 0 else 'append'
            df[cols].to_sql(table, con=db, if_exists=if_exists, index=False)
        t.total = len(query_df)
        t.update(len(query_df))
    t.close()
    pybedtools.cleanup(remove_all=True)
    if not table:
        return
    create_index(table, db)


def create_index(table, db):
    print('Creating index for {}...'.format(table))
    db.execute(
        '''CREATE INDEX idx_{}_frag_chrom ON {}(frag_id, chrom) 
        '''.format(table, table))
    db.execute(
        '''CREATE INDEX idx_{}_id ON {}(id) 
        '''.format(table, table))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Create s PostgreSQL table of variant or gene fragment IDs.")
    parser.add_argument(
        '-f', '--fragments', required=True,
        help='The filepath to the fragment BED file with chrom, start, end, fragment ID.')
    parser.add_argument(
        '-v', '--variants',
        help='''The filepath to the variant lookup file. Must contain chr,
        variant_pos, and variant_id columns''')
    parser.add_argument(
        '-g', '--genes',
        help='''The filepath to a sorted gene reference file containing 
        id, chrom, start, end, name, and gencode_id columns.
        Note that the id indicates the sorting (chrom, start) 
        order of the file starting at 1.''')
    parser.add_argument(
        '-t', '--table',
        help='Name of table e.g. variant_lookup_mboi.')
    parser.add_argument(
        '-u', '--db-url', required=True,
        help='URL of database e.g postgresql://user:password@hostname/database')
    parser.add_argument(
        '-c', '--config',
        default=os.path.join(os.path.dirname(__file__),
                             '../docs/codes3d.conf'),
        help='''The configuration file to be used for this  
        run (default: codes3d.conf)''')

    parser.add_argument(
        '-o', '--output', required=True,
        help='The filepath to output file.')
    args = parser.parse_args()
    c = codes3d.CODES3D(args.config)
    start_time = time.time()
    if not args.variants and not args.genes:
        sys.exit('Variants or gene bed file missing.')
    if os.path.exists(args.output):
        os.remove(args.output)
    if args.variants:
        find_variant_fragments(
            args.variants, args.fragments, args.output,
            args.db_url, args.table)
    elif args.genes:
        find_gene_fragments(
            args.genes, args.fragments, args.output,
            args.db_url, args.table)
    print('Time elapsed: {:.2f} mins'.format(
        (time.time() - start_time)/60))
