
import click
import pandas as pd
import os
@click.command()
@click.option('--file', '-f', help='what file are we changing??')
@click.option('--header', '-h', help="What column are we extracting?")

def read_file(file, header):
    print(os.getcwd())
    os.chdir('../../')
    df = pd.read_excel(file)
    print(df.head())
    print(df.iloc[0])
    headers = df.iloc[0]
    



if __name__ == '__main__':
    read_file()


# python3 metap_script.py -d "python/src" -t "ZM" -r

# python3 save_text_file.py -f 'lukas_proteomics/Omics/Proteomics/PD_exported/Reprocessed_3-pat1-vs-ctrl2&4_pat2-vs-ctrl1&3/CC_AVi_172-M1158-F1-F5-P13031-1-OTITOT-SumInd_TotalPeptideNorm-16plex-pat2_vs_control13-L2FC0.4_upreg_AR-20230418-V3.xlsx' -h 'symbol'