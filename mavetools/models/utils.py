import attr
import collections
import os
from typing import BinaryIO

from .licence import Licence


ONE_TO_THREE = {
    'C': 'CYS',
    'D': 'ASP',
    'S': 'SER',
    'V': 'VAL',
    'Q': 'GLN',
    'K': 'LYS',
    'P': 'PRO',
    'T': 'THR',
    'F': 'PHE',
    'A': 'ALA',
    'H': 'HIS',
    'G': 'GLY',
    'I': 'ILE',
    'L': 'LEU',
    'R': 'ARG',
    'W': 'TRP',
    'N': 'ASN',
    'Y': 'TYR',
    'M': 'MET',
    'E': 'GLU',
    'X': 'UNK'
}


def attrs_filter(attr, value):
    """
    ???
    Parameters
    ----------
    attr
    value

    Returns
    -------

    """
    return value is not None


def attrs_serializer(inst, field, value):
    """
    ???
    Parameters
    ----------
    inst
    field
    value

    Returns
    -------

    """
    if isinstance(value, str):
        if os.path.isfile(value):
            ext = os.path.splitext(value)[1]
            return (f"{field.name}{ext}", open(value, "rb"), "application/octet-stream")
        return value
    if value is not None:
        return value


def get_variant_type(hgvs_pro):

    """
    Routine to deduce the type of variation from hgvs_pro identifier.

    Parameters
    ----------

    hgvs_pro
        hgvs_pro identifier

    Returns
    -------

    variant_type
        As string. If variant types are used more in a future implementation, then one should consider to create corresponding classes and datastructures.
    """

    if hgvs_pro == '_wt' or hgvs_pro[-1] == '=' or hgvs_pro == 'p.(=)' or hgvs_pro == '_sy':
        return 'synonymous'

    if hgvs_pro[-1] == '?' or hgvs_pro[-1] == '*':
        return 'unknown'

    if hgvs_pro[-3:] == 'Ter' or hgvs_pro[:5] == 'p.Ter':
        return 'nonsense'

    if hgvs_pro[-1] == ']':
        return 'multi'

    if hgvs_pro.count('fs*') > 0 or hgvs_pro[-2:] == 'fs':
        return 'frameshift'

    if hgvs_pro.count('delins') > 0:
        return 'indel'

    if hgvs_pro.count('del') > 0:
        return 'deletion'

    if hgvs_pro.count('ins') > 0:
        return 'insertion'

    return 'sav'


def median(l):

    """
    Calculates the median of a list.

    Parameters
    ----------

    l
        A list.

    med
        Median value of the given list.
    """

    n = len(l)
    l = sorted(l)
    if n == 1:
        return l[0]
    if n == 0:
        return None
    if n % 2 == 0:
        med = (l[(n // 2) - 1] + l[n // 2]) / 2.0
    else:
        med = l[(n - 1) // 2]
    return med


def score_scale_function(reported_score, median_synonymous, bottom_perc_median):
    """
    Scaling function for variant effect scores.
    Based on 'Analysis of Large-Scale Mutagenesis Data To Assess the Impact of Single Amino Acid Substitutions' doi: 10.1534/genetics.117.300064
    Slightly improved to work on a broader spectrum of score distributions.

    Parameters
    ----------

    reported_score
        Effect value reported in the scoreset table.

    median_synonymous
        Typical neutral effect value of the corresponding experiment.

    bottom_perc_median
        Typical strong effect value of the corresponding experiment.

    Returns
    -------

    scaled_score
        Scaled effect score.
    """

    #Determine the direction of the effect scores.
    if bottom_perc_median > median_synonymous:
        sign = 1
    else:
        sign = -1

    #Addapted scaling function
    scaled_score = sign * ((reported_score - median_synonymous) / (-1*abs(bottom_perc_median - median_synonymous))) + 1

    return scaled_score


def parseFasta(path=None, new_file=None, lines=None, page=None, left_split=None, right_split=' '):

    """
    Parses a fasta file.

    Parameters
    ----------

    path
        Path to the fasta file. Can be omitted, if a page or lines is given.

    new_file
        Path to where the fasta file should be written. Can be used, when not giving an actual input path.

    lines
        A list of line strings representing the fasta file. Can be omitted, if path or page is given.

    page
        A string representing the fasta file. Can be omitted, if path or lines is given.

    left_split
        A character, which can be used to split the Identifier string given after the '>' symbol in the fasta file format.

    right_split
        A character, which can be used to split the Identifier string given after the '>' symbol in the fasta file format.

    Returns
    -------

    seq_map
        A dictionary of sequence identifiers mapping to sequences.
    """

    if lines is None and page is None:
        f = open(path, 'r')
        lines = f.read().split('\n')
        f.close()
    elif lines is None:
        lines = page.split('\n')

    seq_map = {}
    n = 0

    if new_file is not None:
        new_lines = []

    for line in lines:
        if len(line) == 0:
            continue
        if line[0] == '>':
            entry_id = line[1:]
            if left_split is not None:
                entry_id = entry_id.split(left_split, 1)[1]
            if right_split is not None:
                entry_id = entry_id.split(right_split, 1)[0]
            seq_map[entry_id] = ''
            n += 1
            if new_file is not None:
                new_lines.append(line)
        else:
            seq_map[entry_id] += line
            if new_file is not None:
                new_lines.append(line)

    if new_file is not None:
        f = open(new_file, 'w')
        f.write('\n'.join(new_lines))
        f.close()

    return seq_map


def disect_hgvs_pro(hgvs_pro):

    """
    Routine to split a hgvs_pro into its parts.
    Works only for SAV-like hgvs_pro identifiers!!!

    Parameters
    ----------

    hgvs_pro
        hgvs_pro identifier.

    Returns
    -------

    left_part
        Wildtype amino acid of the SAV.

    right_part
        Mutant amino acid of the SAV.

    pos
        Position of the SAV in the corresponding protein sequence.
    """

    left_part = hgvs_pro[:5]
    right_part = hgvs_pro[-3:]
    pos = int(hgvs_pro[5:-3])
    return left_part, right_part, pos

def apply_offset_to_hgvs_pro(hgvs_pro, offset):

    """
    Routine to add an offset to a hgvs_pro identifier.

    Parameters
    ----------

    hgvs_pro
        hgvs_pro identifier.

    offset
        An integer for that to position in the identifier has to be shifted.

    Returns
    -------

    new_hgvs_pro
        Shifted hgvs_pro identifier.
    """

    if offset == 0:
        return hgvs_pro

    left_part, right_part, old_pos = disect_hgvs_pro(hgvs_pro)

    new_pos = old_pos + offset
    new_hgvs_pro = f'{left_part}{new_pos}{right_part}'
    return new_hgvs_pro


def offset_loop(seq, offset, hgvs_pros, urn, verbosity = 0):

    """
    Identifies problematic offsets and calculates the rate of mismatched amino acids.

    Parameters
    ----------

    seq
        Amino acid sequence of the protein.

    offset
        Offset given in the scoreset metadata.

    hgvs_pros
        List of hgvs_pro identifiers of the corresponding scoreset.

    urn
        MaveDB urn identifier of the corresponding scoreset.

    verbosity
        Verbosity level of the function.

    Returns
    -------

    mismatch_found
        Boolean, True if at least one mismatched amino acid got detected.

    hit_rate
        Rate of correctly matched amino acids.
    """

    mismatch_found = False
    correct_count = 0
    incorrect_count = 0

    for hgvs_pro in hgvs_pros:

        left_part, right_part, old_pos = disect_hgvs_pro(hgvs_pro)
        wt_aa_three_letter = left_part[2:].upper()

        new_pos = old_pos + offset
        try:
            seq_wt_aa_three_letter = ONE_TO_THREE[seq[new_pos-1]]
        except:
            incorrect_count += 1
            if incorrect_count < 2 and verbosity >= 1:
                print(f'Offset failure: seq_pos too large, {urn}: {old_pos} + {offset} = {new_pos} > {len(seq)}')
            mismatch_found = True
            continue

        if wt_aa_three_letter != seq_wt_aa_three_letter:
            incorrect_count += 1
            if incorrect_count < 2 and verbosity >= 1:
                print(f'Offset failure: AA not matched, {urn}: Pos and Offset: {old_pos} + {offset} = {new_pos}; Seq AA: {seq_wt_aa_three_letter}, HGVS AA: {wt_aa_three_letter}')
            mismatch_found = True
            continue
        correct_count += 1
    if (correct_count + incorrect_count) > 0:
        hit_rate = correct_count / (correct_count + incorrect_count)
    else:
        hit_rate = 0
    return mismatch_found, hit_rate

def check_offset(seq, offset, hgvs_pros, urn, verbosity = 0):

    """
    Sanity check of the offset value given in the scoreset metadata.
    When the sanity check fails, tries to find the correct offset.

    Parameters
    ----------

    seq
        Amino acid sequence of the target protein of the corresponding scoreset.

    offset
        Offset value given in the scoreset metadata.

    hgvs_pros
        List of all hgvs_pro identifiers given in the scoreset.

    urn
        MaveDB urn identifier of the scoreset.

    verbosity
        Verbosity level of the function.

    Returns
    -------

    best_offset
        The best offset found, can be None, if none got at least 90% correctly matched amino acids.
    """

    hit_rate_thresh = 0.9
    best_offset = None

    mismatch_found, hit_rate = offset_loop(seq, offset, hgvs_pros, urn, verbosity = verbosity)

    if not mismatch_found:
        return offset

    if hit_rate > hit_rate_thresh:
        hit_rate_thresh = hit_rate
        best_offset = offset

    original_offset = offset
    offset = 0
    mismatch_found, hit_rate = offset_loop(seq, offset, hgvs_pros, urn, verbosity = verbosity)

    if not mismatch_found:
        return offset

    if hit_rate > hit_rate_thresh:
        hit_rate_thresh = hit_rate
        best_offset = offset

    offset = original_offset - 1
    mismatch_found, hit_rate = offset_loop(seq, offset, hgvs_pros, urn, verbosity = verbosity)

    if not mismatch_found:
        return offset

    if hit_rate > hit_rate_thresh:
        hit_rate_thresh = hit_rate
        best_offset = offset

    offset = original_offset + 1
    mismatch_found, hit_rate = offset_loop(seq, offset, hgvs_pros, urn, verbosity = verbosity)

    if not mismatch_found:
        return offset

    if hit_rate > hit_rate_thresh:
        hit_rate_thresh = hit_rate
        best_offset = offset

    offset = -1
    mismatch_found, hit_rate = offset_loop(seq, offset, hgvs_pros, urn, verbosity = verbosity)

    if not mismatch_found:
        return offset

    if hit_rate > hit_rate_thresh:
        hit_rate_thresh = hit_rate
        best_offset = offset

    return best_offset

def prepare_for_encoding(nested_dict):
    """
    Prepares data for encoding by converting the data in the provided nested_dict into
    a json_dict and file_dict
    Parameters
    ----------
    nested_dict (dictionary): data to be converted

    Returns
    -------
    json_dict
    file_dict
    """
    json_dict = {}
    file_dict = {}
    for k, v in nested_dict.items():
        if isinstance(v, tuple):
            file_dict[k] = v
        elif isinstance(v, collections.MutableMapping):
            j, f = prepare_for_encoding(v)
            # Keep the original keys for nested dicts, but flatten the files dict
            json_dict[k] = j
            file_dict.update(f)
        else:
            json_dict[k] = v

    return json_dict, file_dict
