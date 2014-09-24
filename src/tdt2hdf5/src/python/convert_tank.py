import argparse

from tdt2hdf5.convert import import_tank,write_tank

if __name__ == '__main__':

    ap = argparse.ArgumentParser(description='Convert a TDT tank exported from MATLAB to HDF5')
    ap.add_argument('--tank_dir', type=str, required=True, help='The directory of your exported TDT tank.')
    ap.add_argument('--stim_file', type=str, required=True, help='The directory of your exported TDT tank.')
    ap.add_argument('--output_dir', type=str, required=True, help='The directory to write the HDF5 files to.')
    ap.add_argument('--block_mapping_file', type=str, default=None, help='A file that maps your block names to block names of the format <Site>_<Protocol>_L<depth>R<depth>')
    ap.add_argument('--check_zscore', type=int, default=1, help='Whether to filter out cells based on zscore')
    ap.add_argument('--sort_codes', type=str, default='*', help='A comma-separated list of spike sort codes. * means all sort codes.')

    argvals = ap.parse_args()

    tank = import_tank(argvals.tank_dir, argvals.stim_file, block_mapping_file=argvals.block_mapping_file)
    write_tank(tank, argvals.output_dir, sort_codes=argvals.sort_codes, check_zscores=argvals.check_zscore, pval_thresh=0.05)






