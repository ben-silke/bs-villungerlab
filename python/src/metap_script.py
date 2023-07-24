from metaprogramming_r import SalmonRFileWriter, StarRFileWriter

import click

@click.command()
@click.option('--treatment', '-t', help='treatment (e.g. ZM)')
@click.option('--directory', '-d', help='Where should the files be saved?')
@click.option('--file_location', '-f', help="Where are the files?")
@click.option('--all_replicates', '-r', is_flag=True, help="include all replicates? or just 1-3?")
@click.option('--all_treatments', '-a', is_flag=True, help='run for all replicates?')
@click.option('--markdown', '-md', is_flag=True, help='run for all replicates?')
@click.option('--star', '-s', is_flag=True, help="run for star?")

def create_files(treatment, directory, file_location, all_replicates, all_treatments, markdown=False, star=False):
    if star:
        file_writer = StarRFileWriter(treatment, directory, file_location, all_replicates)
    else:
        file_writer = SalmonRFileWriter(treatment, directory, file_location, all_replicates)

    if all_treatments:
        treatments = [
        'ZM',
        'Noc',
        'DHCB',
        'Nutl',
        'Etop' 
    ]
    for treatment in treatments:
        file_writer.treatment = treatment
        file_writer.write_r_file()
        if markdown:
            file_writer.write_markdown_file()
    
    else:
        file_writer.write_r_file()
        if markdown:
            file_writer.write_markdown_file()

if __name__ == '__main__':
    create_files()


# Run to create star scripts:
# python3 metap_script.py -d "python/src" -a -r -s -md -f "data/"