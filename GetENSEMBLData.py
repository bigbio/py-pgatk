"""
Use the ENSEMBL REST API to get species and pep fasta file
"""
import os
import json
import requests
import shutil
from urllib import request,error

def get_species(server='http://rest.ensembl.org', 
                endpoint='/info/species', 
                content_type = 'application/json'):
    #get all info for the species
    species_info = json.loads(requests.get(server+endpoint, headers={ "Content-Type" : content_type}).text)
    return species_info['species']
    
def get_pep_files(species_info, ftp_server = "ftp://ftp.ensembl.org/pub/"):
    """
    Generate the file name from the species info and download the file
    """
    pep_files = []
    for species in species_info:
        gz_file_name = '{}.{}.pep.all.fa.gz'.format(species['name'][0].upper()+species['name'][1:], species['assembly'])
        pep_file_url = '{}/release-{}/fasta/{}/pep/{}'.format(ftp_server, species['release'], species['name'], gz_file_name)
        print(pep_file_url)
        if os.path.isfile(gz_file_name):
            pep_files.append(gz_file_name)
            continue
        #download the pep file
        try:
            downloaded_file = request.urlretrieve(pep_file_url)[0]
        except error.URLError:
            print("Failed to download: ", pep_file_url)
            continue
        #move the pep file to a desired name
        if os.path.isfile(downloaded_file):
            if os.stat(downloaded_file).st_size>1000:
                shutil.move(downloaded_file, gz_file_name)
                pep_files.append(gz_file_name)
            else:
                print("Failed to get the pep file: ", pep_file_url)
                os.remove(downloaded_file)
    return pep_files
    
if __name__ == '__main__':
    species_info = get_species()
    get_pep_files(species_info)
