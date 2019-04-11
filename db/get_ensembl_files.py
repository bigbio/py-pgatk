"""
Get species related files using the ENSEMBL REST API 
"""

import os
from json import loads

import click
from requests import get
from shutil import move
from urllib import request, error


def get_species(server='http://rest.ensembl.org', endpoint='/info/species'):
    """
    Get species info from the server
    """
    species_info = loads(get(server + endpoint, headers={"Content-Type": 'application/json'}).text)
    return species_info['species']


def download_file(file_url: str, file_name: str) -> str:
    """
    Download file_url and move it to file_name, 
    do nothing if file_name already exists
    """
    if os.path.isfile(file_name):
        return file_name
    try:
        downloaded_file = request.urlretrieve(file_url)[0]
    except error.URLError:
        print("Incorrect URL or file not found: ", file_url)
        return None

    # move the pep file to the desired name
    if os.path.isfile(downloaded_file):
        if os.stat(downloaded_file).st_size > 1000:
            move(downloaded_file, file_name)
            return file_name
        else:
            print("Corrupt File (size<1kb): ", file_url)
            os.remove(downloaded_file)
            return None
    else:
        print("Failed to download the file: ", file_url)
        return None


def get_pep_files(species_info: dict, output_dir: str,
                  ftp_server="ftp://ftp.ensembl.org/pub") -> list:
    """
    Generate Pep file name from the species info and download the Pep file
    """
    files = []
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    for species in species_info:
        try:
            file_name = '{}.{}.pep.all.fa.gz'.format(species['name'][0].upper()
                                                     + species['name'][1:],
                                                     species['assembly'])
            file_url = '{}/release-{}/fasta/{}/pep/{}'.format(ftp_server,
                                                              species['release'],
                                                              species['name'],
                                                              file_name)
        except KeyError:
            print("No valid info is available species: ", species)
            continue

        files.append(download_file(file_url, output_dir + '/' + file_name))

    return files


def get_gtf_files(species_info: dict, output_dir: str,
                  ftp_server="ftp://ftp.ensembl.org/pub") -> list:
    """
    Generate GTF file name from the species info and download the GTF file
    """
    files = []
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    for species in species_info:
        try:
            file_name = '{}.{}.{}.gtf.gz'.format(species['name'][0].upper()
                                                 + species['name'][1:],
                                                 species['assembly'],
                                                 species['release'], )
            file_url = '{}/release-{}/gtf/{}/{}'.format(ftp_server,
                                                        species['release'],
                                                        species['name'],
                                                        file_name)
        except KeyError:
            print("No valid info is available species: ", species)
            continue

        files.append(download_file(file_url, output_dir + '/' + file_name))

    return files


@click.command()
@click.option('--output-dir', '-o', help='Output directory for the peptide databases', default="./")
@click.pass_context
def main(output_dir):
    species_info = get_species()
    pep_files = get_pep_files(species_info, output_dir='{}pep-release{}'.format(output_dir, species_info[0]['release']))
    gtf_files = get_gtf_files(species_info, output_dir='{}gtf-release{}'.format(output_dir, species_info[0]['release']))


if __name__ == '__main__':
    main()
