from metaprogramming_r import SalmonRFileWriter


import click

@click.command()
@click.option('--treatment', '-t', help='treatment (e.g. ZM)')
@click.option('--directory', '-d', help='Where should the files be saved?')
@click.option('--file_location', '-f', default="data_directory = file.path('/nobackup/lab_villunger/bsilke', glue('organised/{treatment}/output_salmon'))", help="Where are the files?")
@click.option('--all_replicates', '-r', is_flag=True, help="include all replicates? or just 1-3?")
@click.option('--all_treatments', '-a', is_flag=True, help='run for all replicates?')
@click.option('--include_data_create', '-mdf', is_flag=True, help='run for all replicates?')

def create_files(treatment, directory, file_location, all_replicates, all_treatments, include_data_create=False):
    if all_treatments:
        treatments = [
        'ZM',
        'Noc',
        'DHCB',
        'Nutl',
        'Etop' 
    ]
    for treatment in treatments:
        file_writer = SalmonRFileWriter(treatment, directory, file_location, all_replicates)
        file_writer.write_r_file()
        file_writer.write_markdown_file(include_data_create)
    
    else:
        file_writer = SalmonRFileWriter(treatment, directory, file_location, all_replicates)
        file_writer.write_r_file()
        file_writer.write_markdown_file(include_data_create)

if __name__ == '__main__':
    create_files()


# python3 metap_script.py -d "python/src" -t "ZM" -r