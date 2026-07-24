#Author: Adam Zelezny
#Date: 2026-07-06

import os
import urllib.request

def download_geo_processed_data(gse_accession: str, file_name: str, download_dir: str='./data') -> str:
    if not os.path.exists(download_dir):
        os.makedirs(download_dir)
    prefix = gse_accession[:6] + 'nnn'
    base_url = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse_accession}/suppl/'
    download_url = base_url + file_name
    local_file_path = os.path.join(download_dir, file_name)
    print(f'Constructed URL: {download_url}')
    if not os.path.exists(local_file_path):
        print(f'Downloading {file_name}... This is a large file (approx 2.3GB).')
        urllib.request.urlretrieve(download_url, local_file_path)
        print('Download complete.')
    else:
        print(f'File {file_name} already exists locally. Skipping download.')
    return local_file_path
gse_id = 'GSE231800'
target_file_name = 'GSE231800_S4F_50.rds.gz'
try:
    file_path = download_geo_processed_data(gse_id, target_file_name)
    print(f'\nData successfully downloaded to: {file_path}')
except Exception as e:
    print(f'\nError during download: {e}')
