import socket
import time
import sys
import traceback
import urllib.error
import urllib.parse
import urllib.request
import requests
import json
from Bio import Entrez

from mavetools.models.utils import parseFasta

def is_connected():

    """
    Tests if there is an active internet connection.
    """

    try:
        # connect to the host -- tells us if the host is actually
        # reachable
        socket.create_connection(("1.1.1.1", 53))
        return True
    except OSError:
        pass
    return False


def connection_sleep_cycle():

    """
    A waiting routine for short internet outages.
    """

    while not is_connected():
        print('No connection, sleeping a bit and then try again')
        time.sleep(30)


def retrieve_transcript_sequences(transcript_ids, recursed = False):
    rest_url = 'https://rest.ensembl.org/sequence/id'

    headers = { "Content-Type" : "application/json", "Accept" : "application/json"}
    transcript_seq_dict = {}
    chunksize = 45
    transcript_ids = list(transcript_ids)
    chunk = transcript_ids[:chunksize]
    i = 0
    try_number = 0
    empty_seqs = {}
    while len(chunk) > 0:

        data = json.dumps({'ids':chunk})

        try:
            r = requests.post(rest_url, headers = headers, data = data, params = {'type' : 'protein'})
        except:
            try_number += 1
            if try_number == 4:
                [e, f, g] = sys.exc_info()
                g = traceback.format_exc()
                print(f'POST request failed: {e}\n{f}\n{g}')
                break
            time.sleep(try_number*2)
            continue

        if not r.ok:
            try:
                r.raise_for_status()
            except requests.exceptions.HTTPError as errh:
                print(f'Transcript Sequence Retrieval failed: {data}\nError: {errh.response.status_code}\n{errh.response.text}')
            if recursed:
                return {}, {transcript_ids[0]:None}
            for transcript_id in chunk:
                seq_sub_dict, sub_empty_seqs = retrieve_transcript_sequences([transcript_id], recursed=True)
                transcript_seq_dict.update(seq_sub_dict)
                empty_seqs.update(sub_empty_seqs)
            i += 1
            chunk = transcript_ids[(chunksize*i):(chunksize*(i+1))]
            continue
            #r.raise_for_status()
            #return None

        decoded = r.json()
        try_number = 0

        for entry in decoded:
            transcript_id = entry['query']
            nt_seq = entry['seq']
            transcript_seq_dict[transcript_id] = nt_seq

        i += 1
        chunk = transcript_ids[(chunksize*i):(chunksize*(i+1))]

    return transcript_seq_dict, empty_seqs

def getUniprotSequence(uniprot_ac, tries=0):

    """
    Fetches sequence data from Uniprot database.

    Parameters
    ----------

    uniprot_ac
        Uniprot accesion identifier of the protein sequence to fetch.

    tries
        A count parameter for the trial and error loop for short time internet disconnections.

    Returns
    -------

    wildtype_sequence
        Amino acid sequence of the querried protein.
    """

    if uniprot_ac is None:
        print(f'Uniprot Ac is None')
        return None

    if not uniprot_ac[0:3] == 'UPI':
        url = 'https://www.uniprot.org/uniprot/%s.fasta' % uniprot_ac
    else:
        url = 'https://www.uniprot.org/uniparc/%s.fasta' % uniprot_ac
    connection_sleep_cycle()
    try:
        request = urllib.request.Request(url)
        response = urllib.request.urlopen(request, timeout=(tries + 1) * 10)
        page = response.read(9000000).decode('utf-8')
    except:
        if tries < 3:
            return getSequence(uniprot_ac, tries=tries + 1)
        else:
            return None

    lines = page.split("\n")

    wildtype_sequences = []

    for line in lines:
        if line == '':
            continue
        if line[0] == '>':
            continue
        wildtype_sequences.append(line)

    wildtype_sequence = ("".join(wildtype_sequences)).replace(" ", "").replace("\n", "")

    if wildtype_sequence == '':
        return None

    return wildtype_sequence

def get_refseq_sequences(refseqs, seq_type='protein'):

    """
    Fetches sequence data from NCBI refseq database.
    Uses the biopython library.

    Parameters
    ----------

    refseqs
        comma-separated string of refseq identifiers

    seq_type
        type of sequence to be retrieved, either 'nucleotide' or 'protein'

    Returns
    -------

    seq_map
        dictionary of {refseq_identifier:sequence}
    """

    Entrez.email = ''

    ret_type = 'fasta_cds_aa'
    if seq_type == 'protein' or seq_type == 'nucleotide':
        ret_type = 'fasta'

    try:
        net_handle = Entrez.efetch(
        db=seq_type, id=refseqs, rettype=ret_type, retmode="text"
        )
    except:
        [e, f, g] = sys.exc_info()
        g = traceback.format_exc()
        print(f'Sequence retrieval failed for {refseqs=}\n{e}\n{f}\n{g}')
        return {}
    page = net_handle.read()
    net_handle.close()
    right_split = '_prot'
    left_split = '|'
    if seq_type == 'protein':
        right_split = ' '
        left_split = None    
    elif seq_type == 'nucleotide':
        left_split = None
        right_split = '.'

    seq_map = parseFasta(page=page, left_split=left_split, right_split=right_split)

    return seq_map
