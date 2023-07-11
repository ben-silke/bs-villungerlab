from metaprogramming_r import RFileWriter


import click

@click.command()
@click.option('--treatment', '-t', help='treatment (e.g. ZM)')
@click.option('--directory', '-d', help='Where should the files be saved?')
@click.option('--file_location', '-f', default="data_directory = file.path(getwd(), glue('data/organised/{treatment}/output_salmon'))", help="Where are the files?")
@click.option('--all_replicates', '-r', is_flag=True, help="include all replicates? or just 1-3?")

def create_files(treatment, directory, file_location, all_replicates):
    file_writer = RFileWriter(treatment, directory, file_location, all_replicates)
    file_writer.write_r_file()
    file_writer.write_markdown_file()

if __name__ == '__main__':
    create_files()


# python3 metap_script.py -d "python/src" -t "ZM" -r