#!/usr/bin/env python3

from metaprogramming_r import SalmonResultSheetWriter

import click

@click.command()
@click.option('--directory', '-d', help='Where should the files be saved?')
@click.option('--file_location', '-f', help="Where are the files?")
@click.option('--all_replicates', '-r', is_flag=True, help="include all replicates? or just 1-3?")

def create_files(directory, file_location, all_replicates):
    file_writer = SalmonResultSheetWriter(directory, file_location, all_replicates)
    file_writer.write_result_creation_sheet()

if __name__ == '__main__':
    create_files()


# python3 metap_script.py -d "python/src" -t "ZM" -r