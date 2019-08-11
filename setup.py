#!/usr/bin/env python

# ============================================================================ #
# setup.py: download the database
#
# Author: Alessio Milanese (milanese@embl.de)
#
# ============================================================================ #

motus_version = "2.5.0"
link_db = "https://zenodo.org/record/3364101/files/db_mOTU_v2.5.0.tar.gz"
md5_db = "f533a7b55fc133589f08f50648548b42"
DOI_db = "10.5281/zenodo.3364101"

import os
import sys
import tempfile
import shutil
import subprocess
import hashlib
import time

#function that detect the python version
def python_version():
    if(sys.version_info >= (3,0,0)):
        return(3)
    else:
        return(2)

# load correct library
type_download = ""
if python_version() == 2:
    import urllib2
    type_download = "python2"
else:
    import urllib.request
    type_download = "python3"

# function to print progress bar for python 3
def reporthook(count, block_size, total_size):
    global start_time
    if count == 0:
        start_time = time.time()
        return
    duration = time.time() - start_time
    progress_size = int(count * block_size)
    speed = int(progress_size / (1024 * duration))
    percent = int(count * block_size * 100 / total_size)
    sys.stdout.write("\r %d%%, %d MB, %d KB/s, %d seconds passed" %
                    (percent, progress_size / (1024 * 1024), speed, duration))
    sys.stdout.flush()

def save_f(url, filename):
    if "--no-download-progress" in sys.argv:
        urllib.request.urlretrieve(url, filename)
    else:
        urllib.request.urlretrieve(url, filename, reporthook)


# function to check md5
def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()



# position of the script -------------------------------------------------------
path_mOTUs = os.path.realpath(__file__)
path_array = path_mOTUs.split("/")
relative_path = "/".join(path_array[0:-1])
relative_path = relative_path + "/"

# ------------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------------
def main(argv=None):
    sys.stderr.write(" ------------------------------------------------------------------------------\n")
    sys.stderr.write("|                           SETUP CLASSIFY-GENOMES                             |\n")
    sys.stderr.write(" ------------------------------------------------------------------------------\n")
    # download the files -------------------------------------------------------
    path_versions = relative_path + "db_mOTU/versions"
    if os.path.isfile(path_versions) and "--force-redownload" not in sys.argv:
        sys.stderr.write("Database already downloaded. Not doing anything.\n"
                         "Use --force-redownload to download again.\n")
        sys.exit(0)

    sys.stderr.write("Download the compressed motus database (~1Gb)\n")
    db_name = relative_path+"db_mOTU.tar.gz"

    if type_download == "python2":
        u = urllib2.urlopen(link_db)
        f = open(db_name, 'wb')
        meta = u.info()
        file_size = int(meta.getheaders("Content-Length")[0])

        file_size_dl = 0
        block_sz = 100000
        while True:
            buffer = u.read(block_sz)
            if not buffer:
                break

            file_size_dl += len(buffer)
            f.write(buffer)
            status = r"%10d  [%3.2f%%]" % (file_size_dl, file_size_dl * 100. / file_size)
            status = status + chr(8)*(len(status)+1)
            if "--no-download-progress" not in sys.argv:
                sys.stderr.write(status)

        f.close()
        sys.stderr.write("\n")

    if type_download == "python3":
        save_f(link_db, db_name)
        sys.stderr.write("\n")

    # check md5 ----------------------------------------------------------------
    sys.stderr.write("\nCheck md5: ")
    current_md5 = md5(db_name)

    if current_md5 == md5_db:
        sys.stderr.write("MD5 verified\n")
    else:
        sys.stderr.write("MD5 verification failed!\n")
        os.remove(db_name)
        sys.exit(1)


    # extract files ------------------------------------------------------------
    sys.stderr.write("Extract files from the archive...")
    extract_cmd = "tar -zxvf "+db_name+" -C "+relative_path
    try:
        FNULL = open(os.devnull, 'w')
        process = subprocess.Popen(extract_cmd.split(),stderr=FNULL,stdout=FNULL)
        output, error = process.communicate()
    except:
        sys.stderr.write("Error: failed to extract files\n")
        sys.exit(1)
    if process.returncode:
        sys.stderr.write("Error: failed to extract files\n")
        sys.exit(1)
    else:
        sys.stderr.write("done\n")

    # --- remove db file
    sys.stderr.write("Remove zipped file...")
    os.remove(db_name)
    sys.stderr.write("done\n")

    # remove files that are not used
    os.remove(os.path.abspath(os.path.join(relative_path, 'db_mOTU/db_mOTU_bam_header_CEN')))
    os.remove(os.path.abspath(os.path.join(relative_path, 'db_mOTU/db_mOTU_bam_header_NR')))
    os.remove(os.path.abspath(os.path.join(relative_path, 'db_mOTU/db_mOTU_taxonomy_CAMI.tsv')))
    os.remove(os.path.abspath(os.path.join(relative_path, 'db_mOTU/db_mOTU_MAP_MGCs_to_mOTUs_in-line.tsv')))
    os.remove(os.path.abspath(os.path.join(relative_path, 'db_mOTU/db_mOTU_DB_CEN.fasta.amb')))
    os.remove(os.path.abspath(os.path.join(relative_path, 'db_mOTU/db_mOTU_DB_CEN.fasta.ann')))
    os.remove(os.path.abspath(os.path.join(relative_path, 'db_mOTU/db_mOTU_DB_CEN.fasta.annotations')))
    os.remove(os.path.abspath(os.path.join(relative_path, 'db_mOTU/db_mOTU_DB_CEN.fasta.bwt')))
    os.remove(os.path.abspath(os.path.join(relative_path, 'db_mOTU/db_mOTU_DB_CEN.fasta.pac')))
    os.remove(os.path.abspath(os.path.join(relative_path, 'db_mOTU/db_mOTU_DB_CEN.fasta.sa')))
    os.remove(os.path.abspath(os.path.join(relative_path, 'db_mOTU/db_mOTU_genes_length_NR')))
    os.remove(os.path.abspath(os.path.join(relative_path, 'db_mOTU/db_mOTU_DB_NR.fasta.amb')))
    os.remove(os.path.abspath(os.path.join(relative_path, 'db_mOTU/db_mOTU_DB_NR.fasta.ann')))
    os.remove(os.path.abspath(os.path.join(relative_path, 'db_mOTU/db_mOTU_DB_NR.fasta.bwt')))
    os.remove(os.path.abspath(os.path.join(relative_path, 'db_mOTU/db_mOTU_DB_NR.fasta.pac')))
    os.remove(os.path.abspath(os.path.join(relative_path, 'db_mOTU/db_mOTU_DB_NR.fasta.sa')))
    os.remove(os.path.abspath(os.path.join(relative_path, 'db_mOTU/db_mOTU_padding_coordinates_NR.tsv')))
    os.remove(os.path.abspath(os.path.join(relative_path, 'db_mOTU/db_mOTU_padding_coordinates_CEN.tsv')))
    os.remove(os.path.abspath(os.path.join(relative_path, 'db_mOTU/db_mOTU_taxonomy_ref-mOTUs_short_names.tsv')))
    shutil.rmtree(os.path.abspath(os.path.join(relative_path, 'db_mOTU/db_mOTU_test')))

    return 0        # success

#-------------------------------- run main -------------------------------------
if __name__ == '__main__':
    status = main()
    sys.exit(status)
