import attr
import collections
import os
import math
from typing import BinaryIO

from .licence import Licence
from structman.lib.globalAlignment import (
    init_bp_aligner_class,
    call_biopython_alignment,
)

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

ONE_TO_THREE_LC = {
    'C': 'Cys',
    'D': 'Asp',
    'S': 'Ser',
    'V': 'Val',
    'Q': 'Gln',
    'K': 'Lys',
    'P': 'Pro',
    'T': 'Thr',
    'F': 'Phe',
    'A': 'Ala',
    'H': 'His',
    'G': 'Gly',
    'I': 'Ile',
    'L': 'Leu',
    'R': 'Arg',
    'W': 'Trp',
    'N': 'Asn',
    'Y': 'Tyr',
    'M': 'Met',
    'E': 'Glu',
    'X': 'Unk'
}

THREE_TO_ONE = {
    "00C": "C", "01W": "X", "02K": "A", "03Y": "C", "07O": "C",
    "08P": "C", "0A0": "D", "0A1": "Y", "0A2": "K", "0A8": "C",
    "0AA": "V", "0AB": "V", "0AC": "G", "0AD": "G", "0AF": "W",
    "0AG": "L", "0AH": "S", "0AK": "D", "0AM": "A", "0AP": "C",
    "0AU": "U", "0AV": "A", "0AZ": "P", "0BN": "F", "0C ": "C",
    "0CS": "A", "0DC": "C", "0DG": "G", "0DT": "T", "0FL": "A",
    "0G ": "G", "0NC": "A", "0SP": "A", "0U ": "U", "0YG": "YG",
    "10C": "C", "125": "U", "126": "U", "127": "U", "128": "N",
    "12A": "A", "143": "C", "175": "ASG", "193": "X", "1AP": "A",
    "1MA": "A", "1MG": "G", "1PA": "F", "1PI": "A", "1PR": "N",
    "1SC": "C", "1TQ": "W", "1TY": "Y", "1X6": "S", "200": "F",
    "23F": "F", "23S": "X", "26B": "T", "2AD": "X", "2AG": "A",
    "2AO": "X", "2AR": "A", "2AS": "X", "2AT": "T", "2AU": "U",
    "2BD": "I", "2BT": "T", "2BU": "A", "2CO": "C", "2DA": "A",
    "2DF": "N", "2DM": "N", "2DO": "X", "2DT": "T", "2EG": "G",
    "2FE": "N", "2FI": "N", "2FM": "M", "2GT": "T", "2HF": "H",
    "2LU": "L", "2MA": "A", "2MG": "G", "2ML": "L", "2MR": "R",
    "2MT": "P", "2MU": "U", "2NT": "T", "2OM": "U", "2OT": "T",
    "2PI": "X", "2PR": "G", "2SA": "N", "2SI": "X", "2ST": "T",
    "2TL": "T", "2TY": "Y", "2VA": "V", "2XA": "C", "32S": "X",
    "32T": "X", "3AH": "H", "3AR": "X", "3CF": "F", "3DA": "A",
    "3DR": "N", "3GA": "A", "3MD": "D", "3ME": "U", "3NF": "Y",
    "3QN": "K", "3TY": "X", "3XH": "G", "4AC": "N", "4BF": "Y",
    "4CF": "F", "4CY": "M", "4DP": "W", "4F3": "GYG", "4FB": "P",
    "4FW": "W", "4HT": "W", "4IN": "W", "4MF": "N", "4MM": "X",
    "4OC": "C", "4PC": "C", "4PD": "C", "4PE": "C", "4PH": "F",
    "4SC": "C", "4SU": "U", "4TA": "N", "4U7": "A", "56A": "H",
    "5AA": "A", "5AB": "A", "5AT": "T", "5BU": "U", "5CG": "G",
    "5CM": "C", "5CS": "C", "5FA": "A", "5FC": "C", "5FU": "U",
    "5HP": "E", "5HT": "T", "5HU": "U", "5IC": "C", "5IT": "T",
    "5IU": "U", "5MC": "C", "5MD": "N", "5MU": "U", "5NC": "C",
    "5PC": "C", "5PY": "T", "5SE": "U", "5ZA": "TWG", "64T": "T",
    "6CL": "K", "6CT": "T", "6CW": "W", "6HA": "A", "6HC": "C",
    "6HG": "G", "6HN": "K", "6HT": "T", "6IA": "A", "6MA": "A",
    "6MC": "A", "6MI": "N", "6MT": "A", "6MZ": "N", "6OG": "G",
    "70U": "U", "7DA": "A", "7GU": "G", "7JA": "I", "7MG": "G",
    "8AN": "A", "8FG": "G", "8MG": "G", "8OG": "G", "9NE": "E",
    "9NF": "F", "9NR": "R", "9NV": "V", "A  ": "A", "A1P": "N",
    "A23": "A", "A2L": "A", "A2M": "A", "A34": "A", "A35": "A",
    "A38": "A", "A39": "A", "A3A": "A", "A3P": "A", "A40": "A",
    "A43": "A", "A44": "A", "A47": "A", "A5L": "A", "A5M": "C",
    "A5N": "N", "A5O": "A", "A66": "X", "AA3": "A", "AA4": "A",
    "AAR": "R", "AB7": "X", "ABA": "A", "ABR": "A", "ABS": "A",
    "ABT": "N", "ACB": "D", "ACL": "R", "AD2": "A", "ADD": "X",
    "ADX": "N", "AEA": "X", "AEI": "D", "AET": "A", "AFA": "N",
    "AFF": "N", "AFG": "G", "AGM": "R", "AGT": "C", "AHB": "N",
    "AHH": "X", "AHO": "A", "AHP": "A", "AHS": "X", "AHT": "X",
    "AIB": "A", "AKL": "D", "AKZ": "D", "ALA": "A", "ALC": "A",
    "ALM": "A", "ALN": "A", "ALO": "T", "ALQ": "X", "ALS": "A",
    "ALT": "A", "ALV": "A", "ALY": "K", "AN8": "A", "AP7": "A",
    "APE": "X", "APH": "A", "API": "K", "APK": "K", "APM": "X",
    "APP": "X", "AR2": "R", "AR4": "E", "AR7": "R", "ARG": "R",
    "ARM": "R", "ARO": "R", "ARV": "X", "AS ": "A", "AS2": "D",
    "AS9": "X", "ASA": "D", "ASB": "D", "ASI": "D", "ASK": "D",
    "ASL": "D", "ASM": "X", "ASN": "N", "ASP": "D", "ASQ": "D",
    "ASU": "N", "ASX": "B", "ATD": "T", "ATL": "T", "ATM": "T",
    "AVC": "A", "AVN": "X", "AYA": "A", "AYG": "AYG", "AZK": "K",
    "AZS": "S", "AZY": "Y", "B1F": "F", "B1P": "N", "B2A": "A",
    "B2F": "F", "B2I": "I", "B2V": "V", "B3A": "A", "B3D": "D",
    "B3E": "E", "B3K": "K", "B3L": "X", "B3M": "X", "B3Q": "X",
    "B3S": "S", "B3T": "X", "B3U": "H", "B3X": "N", "B3Y": "Y",
    "BB6": "C", "BB7": "C", "BB8": "F", "BB9": "C", "BBC": "C",
    "BCS": "C", "BE2": "X", "BFD": "D", "BG1": "S", "BGM": "G",
    "BH2": "D", "BHD": "D", "BIF": "F", "BIL": "X", "BIU": "I",
    "BJH": "X", "BLE": "L", "BLY": "K", "BMP": "N", "BMT": "T",
    "BNN": "F", "BNO": "X", "BOE": "T", "BOR": "R", "BPE": "C",
    "BRU": "U", "BSE": "S", "BT5": "N", "BTA": "L", "BTC": "C",
    "BTR": "W", "BUC": "C", "BUG": "V", "BVP": "U", "BZG": "N",
    "C  ": "C", "C12": "TYG", "C1X": "K", "C25": "C", "C2L": "C",
    "C2S": "C", "C31": "C", "C32": "C", "C34": "C", "C36": "C",
    "C37": "C", "C38": "C", "C3Y": "C", "C42": "C", "C43": "C",
    "C45": "C", "C46": "C", "C49": "C", "C4R": "C", "C4S": "C",
    "C5C": "C", "C66": "X", "C6C": "C", "C99": "TFG", "CAF": "C",
    "CAL": "X", "CAR": "C", "CAS": "C", "CAV": "X", "CAY": "C",
    "CB2": "C", "CBR": "C", "CBV": "C", "CCC": "C", "CCL": "K",
    "CCS": "C", "CCY": "CYG", "CDE": "X", "CDV": "X", "CDW": "C",
    "CEA": "C", "CFL": "C", "CFY": "FCYG", "CG1": "G", "CGA": "E",
    "CGU": "E", "CH ": "C", "CH6": "MYG", "CH7": "KYG", "CHF": "X",
    "CHG": "X", "CHP": "G", "CHS": "X", "CIR": "R", "CJO": "GYG",
    "CLE": "L", "CLG": "K", "CLH": "K", "CLV": "AFG", "CM0": "N",
    "CME": "C", "CMH": "C", "CML": "C", "CMR": "C", "CMT": "C",
    "CNU": "U", "CP1": "C", "CPC": "X", "CPI": "X", "CQR": "GYG",
    "CR0": "TLG", "CR2": "GYG", "CR5": "G", "CR7": "KYG", "CR8": "HYG",
    "CRF": "TWG", "CRG": "THG", "CRK": "MYG", "CRO": "GYG", "CRQ": "QYG",
    "CRU": "EYG", "CRW": "ASG", "CRX": "ASG", "CS0": "C", "CS1": "C",
    "CS3": "C", "CS4": "C", "CS8": "N", "CSA": "C", "CSB": "C",
    "CSD": "C", "CSE": "C", "CSF": "C", "CSH": "SHG", "CSI": "G",
    "CSJ": "C", "CSL": "C", "CSO": "C", "CSP": "C", "CSR": "C",
    "CSS": "C", "CSU": "C", "CSW": "C", "CSX": "C", "CSY": "SYG",
    "CSZ": "C", "CTE": "W", "CTG": "T", "CTH": "T", "CUC": "X",
    "CWR": "S", "CXM": "M", "CY0": "C", "CY1": "C", "CY3": "C",
    "CY4": "C", "CYA": "C", "CYD": "C", "CYF": "C", "CYG": "C",
    "CYJ": "X", "CYM": "C", "CYQ": "C", "CYR": "C", "CYS": "C",
    "CZ2": "C", "CZO": "GYG", "CZZ": "C", "D11": "T", "D1P": "N",
    "D3 ": "N", "D33": "N", "D3P": "G", "D3T": "T", "D4M": "T",
    "D4P": "X", "DA ": "A", "DA2": "X", "DAB": "A", "DAH": "F",
    "DAL": "A", "DAR": "R", "DAS": "D", "DBB": "T", "DBM": "N",
    "DBS": "S", "DBU": "T", "DBY": "Y", "DBZ": "A", "DC ": "C",
    "DC2": "C", "DCG": "G", "DCI": "X", "DCL": "X", "DCT": "C",
    "DCY": "C", "DDE": "H", "DDG": "G", "DDN": "U", "DDX": "N",
    "DFC": "C", "DFG": "G", "DFI": "X", "DFO": "X", "DFT": "N",
    "DG ": "G", "DGH": "G", "DGI": "G", "DGL": "E", "DGN": "Q",
    "DHA": "S", "DHI": "H", "DHL": "X", "DHN": "V", "DHP": "X",
    "DHU": "U", "DHV": "V", "DI ": "I", "DIL": "I", "DIR": "R",
    "DIV": "V", "DLE": "L", "DLS": "K", "DLY": "K", "DM0": "K",
    "DMH": "N", "DMK": "D", "DMT": "X", "DN ": "N", "DNE": "L",
    "DNG": "L", "DNL": "K", "DNM": "L", "DNP": "A", "DNR": "C",
    "DNS": "K", "DOA": "X", "DOC": "C", "DOH": "D", "DON": "L",
    "DPB": "T", "DPH": "F", "DPL": "P", "DPP": "A", "DPQ": "Y",
    "DPR": "P", "DPY": "N", "DRM": "U", "DRP": "N", "DRT": "T",
    "DRZ": "N", "DSE": "S", "DSG": "N", "DSN": "S", "DSP": "D",
    "DT ": "T", "DTH": "T", "DTR": "W", "DTY": "Y", "DU ": "U",
    "DVA": "V", "DXD": "N", "DXN": "N", "DYG": "DYG", "DYS": "C",
    "DZM": "A", "E  ": "A", "E1X": "A", "ECC": "Q", "EDA": "A",
    "EFC": "C", "EHP": "F", "EIT": "T", "ENP": "N", "ESB": "Y",
    "ESC": "M", "EXB": "X", "EXY": "L", "EY5": "N", "EYS": "X",
    "F2F": "F", "FA2": "A", "FA5": "N", "FAG": "N", "FAI": "N",
    "FB5": "A", "FB6": "A", "FCL": "F", "FFD": "N", "FGA": "E",
    "FGL": "G", "FGP": "S", "FHL": "X", "FHO": "K", "FHU": "U",
    "FLA": "A", "FLE": "L", "FLT": "Y", "FME": "M", "FMG": "G",
    "FMU": "N", "FOE": "C", "FOX": "G", "FP9": "P", "FPA": "F",
    "FRD": "X", "FT6": "W", "FTR": "W", "FTY": "Y", "FVA": "V",
    "FZN": "K", "G  ": "G", "G25": "G", "G2L": "G", "G2S": "G",
    "G31": "G", "G32": "G", "G33": "G", "G36": "G", "G38": "G",
    "G42": "G", "G46": "G", "G47": "G", "G48": "G", "G49": "G",
    "G4P": "N", "G7M": "G", "GAO": "G", "GAU": "E", "GCK": "C",
    "GCM": "X", "GDP": "G", "GDR": "G", "GFL": "G", "GGL": "E",
    "GH3": "G", "GHG": "Q", "GHP": "G", "GL3": "G", "GLH": "Q",
    "GLJ": "E", "GLK": "E", "GLM": "X", "GLN": "Q", "GLQ": "E",
    "GLU": "E", "GLX": "Z", "GLY": "G", "GLZ": "G", "GMA": "E",
    "GMS": "G", "GMU": "U", "GN7": "G", "GND": "X", "GNE": "N",
    "GOM": "G", "GPL": "K", "GS ": "G", "GSC": "G", "GSR": "G",
    "GSS": "G", "GSU": "E", "GT9": "C", "GTP": "G", "GVL": "X",
    "GYC": "CYG", "GYS": "SYG", "H2U": "U", "H5M": "P", "HAC": "A",
    "HAR": "R", "HBN": "H", "HCS": "X", "HDP": "U", "HEU": "U",
    "HFA": "X", "HGL": "X", "HHI": "H", "HHK": "AK", "HIA": "H",
    "HIC": "H", "HIP": "H", "HIQ": "H", "HIS": "H", "HL2": "L",
    "HLU": "L", "HMR": "R", "HOL": "N", "HPC": "F", "HPE": "F",
    "HPH": "F", "HPQ": "F", "HQA": "A", "HRG": "R", "HRP": "W",
    "HS8": "H", "HS9": "H", "HSE": "S", "HSL": "S", "HSO": "H",
    "HTI": "C", "HTN": "N", "HTR": "W", "HV5": "A", "HVA": "V",
    "HY3": "P", "HYP": "P", "HZP": "P", "I  ": "I", "I2M": "I",
    "I58": "K", "I5C": "C", "IAM": "A", "IAR": "R", "IAS": "D",
    "IC ": "C", "IEL": "K", "IEY": "HYG", "IG ": "G", "IGL": "G",
    "IGU": "G", "IIC": "SHG", "IIL": "I", "ILE": "I", "ILG": "E",
    "ILX": "I", "IMC": "C", "IML": "I", "IOY": "F", "IPG": "G",
    "IPN": "N", "IRN": "N", "IT1": "K", "IU ": "U", "IYR": "Y",
    "IYT": "T", "IZO": "M", "JJJ": "C", "JJK": "C", "JJL": "C",
    "JW5": "N", "K1R": "C", "KAG": "G", "KCX": "K", "KGC": "K",
    "KNB": "A", "KOR": "M", "KPI": "K", "KST": "K", "KYQ": "K",
    "L2A": "X", "LA2": "K", "LAA": "D", "LAL": "A", "LBY": "K",
    "LC ": "C", "LCA": "A", "LCC": "N", "LCG": "G", "LCH": "N",
    "LCK": "K", "LCX": "K", "LDH": "K", "LED": "L", "LEF": "L",
    "LEH": "L", "LEI": "V", "LEM": "L", "LEN": "L", "LET": "X",
    "LEU": "L", "LEX": "L", "LG ": "G", "LGP": "G", "LHC": "X",
    "LHU": "U", "LKC": "N", "LLP": "K", "LLY": "K", "LME": "E",
    "LMF": "K", "LMQ": "Q", "LMS": "N", "LP6": "K", "LPD": "P",
    "LPG": "G", "LPL": "X", "LPS": "S", "LSO": "X", "LTA": "X",
    "LTR": "W", "LVG": "G", "LVN": "V", "LYF": "K", "LYK": "K",
    "LYM": "K", "LYN": "K", "LYR": "K", "LYS": "K", "LYX": "K",
    "LYZ": "K", "M0H": "C", "M1G": "G", "M2G": "G", "M2L": "K",
    "M2S": "M", "M30": "G", "M3L": "K", "M5M": "C", "MA ": "A",
    "MA6": "A", "MA7": "A", "MAA": "A", "MAD": "A", "MAI": "R",
    "MBQ": "Y", "MBZ": "N", "MC1": "S", "MCG": "X", "MCL": "K",
    "MCS": "C", "MCY": "C", "MD3": "C", "MD6": "G", "MDH": "X",
    "MDO": "ASG", "MDR": "N", "MEA": "F", "MED": "M", "MEG": "E",
    "MEN": "N", "MEP": "U", "MEQ": "Q", "MET": "M", "MEU": "G",
    "MF3": "X", "MFC": "GYG", "MG1": "G", "MGG": "R", "MGN": "Q",
    "MGQ": "A", "MGV": "G", "MGY": "G", "MHL": "L", "MHO": "M",
    "MHS": "H", "MIA": "A", "MIS": "S", "MK8": "L", "ML3": "K",
    "MLE": "L", "MLL": "L", "MLY": "K", "MLZ": "K", "MME": "M",
    "MMO": "R", "MMT": "T", "MND": "N", "MNL": "L", "MNU": "U",
    "MNV": "V", "MOD": "X", "MP8": "P", "MPH": "X", "MPJ": "X",
    "MPQ": "G", "MRG": "G", "MSA": "G", "MSE": "M", "MSL": "M",
    "MSO": "M", "MSP": "X", "MT2": "M", "MTR": "T", "MTU": "A",
    "MTY": "Y", "MVA": "V", "N  ": "N", "N10": "S", "N2C": "X",
    "N5I": "N", "N5M": "C", "N6G": "G", "N7P": "P", "NA8": "A",
    "NAL": "A", "NAM": "A", "NB8": "N", "NBQ": "Y", "NC1": "S",
    "NCB": "A", "NCX": "N", "NCY": "X", "NDF": "F", "NDN": "U",
    "NEM": "H", "NEP": "H", "NF2": "N", "NFA": "F", "NHL": "E",
    "NIT": "X", "NIY": "Y", "NLE": "L", "NLN": "L", "NLO": "L",
    "NLP": "L", "NLQ": "Q", "NMC": "G", "NMM": "R", "NMS": "T",
    "NMT": "T", "NNH": "R", "NP3": "N", "NPH": "C", "NPI": "A",
    "NRP": "LYG", "NRQ": "MYG", "NSK": "X", "NTY": "Y", "NVA": "V",
    "NYC": "TWG", "NYG": "NYG", "NYM": "N", "NYS": "C", "NZH": "H",
    "O12": "X", "O2C": "N", "O2G": "G", "OAD": "N", "OAS": "S",
    "OBF": "X", "OBS": "X", "OCS": "C", "OCY": "C", "ODP": "N",
    "OHI": "H", "OHS": "D", "OIC": "X", "OIP": "I", "OLE": "X",
    "OLT": "T", "OLZ": "S", "OMC": "C", "OMG": "G", "OMT": "M",
    "OMU": "U", "ONE": "U", "ONH": "A", "ONL": "X", "OPR": "R",
    "ORN": "A", "ORQ": "R", "OSE": "S", "OTB": "X", "OTH": "T",
    "OTY": "Y", "OXX": "D", "P  ": "G", "P1L": "C", "P1P": "N",
    "P2T": "T", "P2U": "U", "P2Y": "P", "P5P": "A", "PAQ": "Y",
    "PAS": "D", "PAT": "W", "PAU": "A", "PBB": "C", "PBF": "F",
    "PBT": "N", "PCA": "E", "PCC": "P", "PCE": "X", "PCS": "F",
    "PDL": "X", "PDU": "U", "PEC": "C", "PF5": "F", "PFF": "F",
    "PFX": "X", "PG1": "S", "PG7": "G", "PG9": "G", "PGL": "X",
    "PGN": "G", "PGP": "G", "PGY": "G", "PHA": "F", "PHD": "D",
    "PHE": "F", "PHI": "F", "PHL": "F", "PHM": "F", "PIA": "AYG",
    "PIV": "X", "PLE": "L", "PM3": "F", "PMT": "C", "POM": "P",
    "PPN": "F", "PPU": "A", "PPW": "G", "PQ1": "N", "PR3": "C",
    "PR5": "A", "PR9": "P", "PRN": "A", "PRO": "P", "PRS": "P",
    "PSA": "F", "PSH": "H", "PST": "T", "PSU": "U", "PSW": "C",
    "PTA": "X", "PTH": "Y", "PTM": "Y", "PTR": "Y", "PU ": "A",
    "PUY": "N", "PVH": "H", "PVL": "X", "PYA": "A", "PYO": "U",
    "PYX": "C", "PYY": "N", "QLG": "QLG", "QMM": "Q", "QPA": "C",
    "QPH": "F", "QUO": "G", "R  ": "A", "R1A": "C", "R4K": "W",
    "RC7": "HYG", "RE0": "W", "RE3": "W", "RIA": "A", "RMP": "A",
    "RON": "X", "RT ": "T", "RTP": "N", "S1H": "S", "S2C": "C",
    "S2D": "A", "S2M": "T", "S2P": "A", "S4A": "A", "S4C": "C",
    "S4G": "G", "S4U": "U", "S6G": "G", "SAC": "S", "SAH": "C",
    "SAR": "G", "SBL": "S", "SC ": "C", "SCH": "C", "SCS": "C",
    "SCY": "C", "SD2": "X", "SDG": "G", "SDP": "S", "SEB": "S",
    "SEC": "A", "SEG": "A", "SEL": "S", "SEM": "S", "SEN": "S",
    "SEP": "S", "SER": "S", "SET": "S", "SGB": "S", "SHC": "C",
    "SHP": "G", "SHR": "K", "SIB": "C", "SIC": "DC", "SLA": "P",
    "SLR": "P", "SLZ": "K", "SMC": "C", "SME": "M", "SMF": "F",
    "SMP": "A", "SMT": "T", "SNC": "C", "SNN": "N", "SOC": "C",
    "SOS": "N", "SOY": "S", "SPT": "T", "SRA": "A", "SSU": "U",
    "STY": "Y", "SUB": "X", "SUI": "DG", "SUN": "S", "SUR": "U",
    "SVA": "S", "SVV": "S", "SVW": "S", "SVX": "S", "SVY": "S",
    "SVZ": "X", "SWG": "SWG", "SYS": "C", "T  ": "T", "T11": "F",
    "T23": "T", "T2S": "T", "T2T": "N", "T31": "U", "T32": "T",
    "T36": "T", "T37": "T", "T38": "T", "T39": "T", "T3P": "T",
    "T41": "T", "T48": "T", "T49": "T", "T4S": "T", "T5O": "U",
    "T5S": "T", "T66": "X", "T6A": "A", "TA3": "T", "TA4": "X",
    "TAF": "T", "TAL": "N", "TAV": "D", "TBG": "V", "TBM": "T",
    "TC1": "C", "TCP": "T", "TCQ": "Y", "TCR": "W", "TCY": "A",
    "TDD": "L", "TDY": "T", "TFE": "T", "TFO": "A", "TFQ": "F",
    "TFT": "T", "TGP": "G", "TH6": "T", "THC": "T", "THO": "X",
    "THR": "T", "THX": "N", "THZ": "R", "TIH": "A", "TLB": "N",
    "TLC": "T", "TLN": "U", "TMB": "T", "TMD": "T", "TNB": "C",
    "TNR": "S", "TOX": "W", "TP1": "T", "TPC": "C", "TPG": "G",
    "TPH": "X", "TPL": "W", "TPO": "T", "TPQ": "Y", "TQI": "W",
    "TQQ": "W", "TRF": "W", "TRG": "K", "TRN": "W", "TRO": "W",
    "TRP": "W", "TRQ": "W", "TRW": "W", "TRX": "W", "TS ": "N",
    "TST": "X", "TT ": "N", "TTD": "T", "TTI": "U", "TTM": "T",
    "TTQ": "W", "TTS": "Y", "TY1": "Y", "TY2": "Y", "TY3": "Y",
    "TY5": "Y", "TYB": "Y", "TYI": "Y", "TYJ": "Y", "TYN": "Y",
    "TYO": "Y", "TYQ": "Y", "TYR": "Y", "TYS": "Y", "TYT": "Y",
    "TYU": "N", "TYW": "Y", "TYX": "X", "TYY": "Y", "TZB": "X",
    "TZO": "X", "U  ": "U", "U25": "U", "U2L": "U", "U2N": "U",
    "U2P": "U", "U31": "U", "U33": "U", "U34": "U", "U36": "U",
    "U37": "U", "U8U": "U", "UAR": "U", "UCL": "U", "UD5": "U",
    "UDP": "N", "UFP": "N", "UFR": "U", "UFT": "U", "UMA": "A",
    "UMP": "U", "UMS": "U", "UN1": "X", "UN2": "X", "UNK": "X",
    "UR3": "U", "URD": "U", "US1": "U", "US2": "U", "US3": "T",
    "US5": "U", "USM": "U", "VAD": "V", "VAF": "V", "VAL": "V",
    "VB1": "K", "VDL": "X", "VLL": "X", "VLM": "X", "VMS": "X",
    "VOL": "X", "WCR": "GYG", "X  ": "G", "X2W": "E", "X4A": "N",
    "X9Q": "AFG", "XAD": "A", "XAE": "N", "XAL": "A", "XAR": "N",
    "XCL": "C", "XCN": "C", "XCP": "X", "XCR": "C", "XCS": "N",
    "XCT": "C", "XCY": "C", "XGA": "N", "XGL": "G", "XGR": "G",
    "XGU": "G", "XPR": "P", "XSN": "N", "XTH": "T", "XTL": "T",
    "XTR": "T", "XTS": "G", "XTY": "N", "XUA": "A", "XUG": "G",
    "XX1": "K", "XXY": "THG", "XYG": "DYG", "Y  ": "A", "YCM": "C",
    "YG ": "G", "YOF": "Y", "YRR": "N", "YYG": "G", "Z  ": "C",
    "Z01": "A", "ZAD": "A", "ZAL": "A", "ZBC": "C", "ZBU": "U",
    "ZCL": "F", "ZCY": "C", "ZDU": "U", "ZFB": "X", "ZGU": "G",
    "ZHP": "N", "ZTH": "T", "ZU0": "T", "ZZJ": "A", "AME": "M",
    "KFP": "L", "7O5": "A", "AC5": "L", "ZXW": "X", "LGY": "L",
    "PRJ": "P", "G5G": "L", "CNT": "T", "69P": "X", "STA": "X",
    "DNW": "A", "8MC": "X", "BBS": "X", "L3O": "L", "823": "N",
    "ADS": "X", "DPN": "F", "2PP": "X", "PRR": "X", "0UO": "W",
    "HOA": "X", "PLF": "X", "OCE": "X", "E6F": "X", "T8L": "T",
    "JMH": "X", "34H": "X", "6NA": "X", "A5R": "X", "0MG": "X",
    "VLT": "X", "8DT": "X", "81S": "X", "0EA": "Y", "0E5": "T",
    "FBE": "X", "6J9": "X", "RNG": "X", "3EG": "X", "0QE": "X",
    "IL0": "X", "SUJ": "X", "B8N": "X", "DKA": "X", "9MN": "X",
    "ACE": "X", "BAL": "X", "7XC": "F", "OSL": "X", "BTK": "L",
    "IVA": "X", "3WU": "X", "3WT": "X", "3FB": "X", "HT7": "W",
    "V9M": "X"
}

def aac_to_hgvs_pro(aac):
    aa1 = aac[0]
    aa2 = aac[-1]
    pos = aac[1:-1]

    aa1_three_letter = ONE_TO_THREE_LC[aa1]
    aa2_three_letter = ONE_TO_THREE_LC[aa2]

    hgvs_pro = f'p.{aa1_three_letter}{pos}{aa2_three_letter}'

    return hgvs_pro

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
    
    if hgvs_pro.count('dup') > 0:
        return 'duplication'

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


def score_scale_function(reported_score, median_synonymous, bottom_perc_median, top_perc_median):
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
    if (sign * (reported_score - median_synonymous)) > 0:
        stretch = -1*abs(bottom_perc_median - median_synonymous)
    else:
        stretch = -1*abs(top_perc_median - median_synonymous)

    if stretch == 0:
        scaled_score = 1
    else:
        scaled_score = sign * ((reported_score - median_synonymous) / stretch) + 1

    return scaled_score

def second_scale(score, min_value, max_value, shift):
    if shift is None:
        return score
    if shift < 0.:
        return score
    zero_one_score = (score - min_value)

    zero_one_scaled_score = zero_one_score ** (1 / shift)

    scaled_score = zero_one_scaled_score + min_value

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
        #print(line)
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


def disect_hgvs(hgvs):

    """
    Routine to split a hgvs into its parts.

    Parameters
    ----------

    hgvs
        hgvs identifier.

    Returns
    -------

    precursor
        hgvs precursor symbol (for example "p." for protein hgvs_pro)

    parts
        list with disected hgvs mutation identifiers
    """

    precursor = hgvs[:2]
    hgvs_body = hgvs[2:]

    if len(hgvs_body) < 4:
        muts = [('invalid', hgvs_body)]

    elif hgvs_body[0] == '[':
        muts = []

        for mut in hgvs_body[1:-1].split(';'):
            muts.append(disect_hgvs_mut(mut))

    else:
        muts = [disect_hgvs_mut(hgvs_body)]

    return precursor, muts


def disect_hgvs_pro_sav(hgvs_pro):

    """
    Routine to split a SAV hgvs_pro into its parts.

    Parameters
    ----------

    hgvs_pro
        hgvs_pro identifier.

    Returns
    -------

    precursor
        hgvs precursor symbol (for example "p." for protein hgvs_pro)

    left_part
        Wildtype amino acid of the SAV.

    right_part
        Mutant amino acid of the SAV.

    pos
        Position of the SAV in the corresponding protein sequence.
    """

    precursor = hgvs_pro[:2]
    hgvs_body = hgvs_pro[2:]

    left_part, right_part, pos = disect_hgvs_single_pos(hgvs_body)
    return precursor, left_part, right_part, pos


def disect_hgvs_single_pos(hgvs_single_pos):

    if hgvs_single_pos.count(':') > 0:
        hgvs_single_pos = hgvs_single_pos.rsplit(':',1)[1]
    if hgvs_single_pos[:2] == 'p.':
        hgvs_single_pos = hgvs_single_pos[2:]

    left_part = hgvs_single_pos[:3]
    if hgvs_single_pos[-1] == '=':
        right_part = hgvs_single_pos[-1]
        pos = int(hgvs_single_pos[3:-1])
    else:
        right_part = hgvs_single_pos[-3:]
        pos = int(hgvs_single_pos[3:-3])
    return left_part, right_part, pos

def disect_hgvs_mut(hgvs_mut):

    if len(hgvs_mut) < 4 or hgvs_mut[-1] == '?':
        return 'invalid', hgvs_mut

    if hgvs_mut.count('_') > 0:

        if hgvs_mut.count('delins') > 0:
            separator = 'delins'
        elif hgvs_mut.count('ins') > 0:
            separator = 'ins'
        elif hgvs_mut.count('del') > 0:
            separator = 'del'
        else:
            print(hgvs_mut)
            sys.exit(1)

        left_tuple, right_body = hgvs_mut.split('_')
        left_part = left_tuple[:3]
        left_pos = int(left_tuple[3:])

        if right_body.count('delins') > 0:
            separator = 'delins'
        elif right_body.count('del') > 0:
            separator = 'del'
        elif right_body.count('ins') > 0:
            separator = 'ins'
        else:
            print(hgvs_mut)
            sys.exit(1)
        right_tuple,tail = right_body.split(separator)
        right_part = right_tuple[:3]
        right_pos = int(right_tuple[3:])

        return 'indel', left_part, left_pos, right_part, right_pos, f'{separator}{tail}'
    
    elif hgvs_mut.count('del') > 0 or hgvs_mut.count('ins') > 0:
        """TODO
        if hgvs_mut.count('delins') > 0:
            separator = 'delins'
        elif hgvs_mut.count('ins') > 0:
            separator = 'ins'
        elif hgvs_mut.count('del') > 0:
            separator = 'del'
        """
        return 'invalid', hgvs_mut

    elif hgvs_mut.count('fs') > 0:
        if hgvs_mut.count('fs*') > 0:
            mut_body, fs_number = hgvs_mut.split('fs*')
            left_part, right_part, pos = disect_hgvs_single_pos(mut_body)
            return 'frameshift', left_part, right_part, pos, f'fs*{fs_number}'
        else:
            left_part = hgvs_mut[:3]
            pos = int(hgvs_mut[3:-2])
            right_part = 'fs'
            return 'frameshift', left_part, right_part, pos, ''

    else:
        left_part, right_part, pos = disect_hgvs_single_pos(hgvs_mut)
        return 'sav', left_part, right_part, pos


def apply_offset_to_hgvs_mut_parts(parts, offset):
    if parts[0] == 'invalid':
        return parts[1]
    elif parts[0] == 'sav':
        left_part, right_part, pos = parts[1:]
        new_pos = pos + offset
        return f'{left_part}{new_pos}{right_part}'
    elif parts[0] == 'frameshift':
        left_part, right_part, pos, tail = parts[1:]
        new_pos = pos + offset
        return f'{left_part}{new_pos}{right_part}{tail}'
    else:
        left_part, left_pos, right_part, right_pos, tail = parts[1:]
        new_left_pos = left_pos + offset
        new_right_pos = right_pos + offset

        return f'{left_part}{new_left_pos}_{right_part}{new_right_pos}{tail}'


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

    precursor, muts = disect_hgvs(hgvs_pro)
    if len(muts) == 1:
        new_mut = apply_offset_to_hgvs_mut_parts(muts[0], offset)
        new_hgvs_pro = f'{precursor}{new_mut}'
    else:
        new_muts = []
        for mut in muts:
            new_muts.append(apply_offset_to_hgvs_mut_parts(mut, offset))
        internal_string = ';'.join(new_muts)
        new_hgvs_pro = f'{precursor}[{internal_string}]'
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

        precursor, left_part, right_part, old_pos = disect_hgvs_pro_sav(hgvs_pro)
        wt_aa_three_letter = left_part.upper()

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


def extract_seq(hgvs_pros):
    chars = [None]
    for hgvs_pro in hgvs_pros:
        precursor, left_part, right_part, pos = disect_hgvs_pro_sav(hgvs_pro)
        wt_aa_three_letter = left_part.upper()
        wt_aa = THREE_TO_ONE[wt_aa_three_letter]

        while len(chars) <= pos:
            chars.append('G')
        
        chars[pos] = wt_aa

    seq = ''.join(chars[1:])
    return seq


def get_offset_from_aligned_seq(aligned_extracted_sequence):
    offset = 0
    for char in aligned_extracted_sequence:
        if char != '-':
            break
        offset += 1

    return offset



def repair_seq(seq, offset, hgvs_pros):
    pos_dict = {}
    for hgvs_pro in hgvs_pros:
        precursor, left_part, right_part, old_pos = disect_hgvs_pro_sav(hgvs_pro)
        new_pos = old_pos + offset
        wt_aa_three_letter = left_part.upper()
        wt_aa = THREE_TO_ONE[wt_aa_three_letter]
        if new_pos not in pos_dict:
            pos_dict[new_pos] = wt_aa
        elif wt_aa != pos_dict[new_pos]:
            return False, f'Mismatching WT AA {wt_aa=} {new_pos=} {pos_dict[new_pos]=}'
    repair_seq = ''
    for pos, wt_aa in enumerate(seq):
        seq_pos = pos+1
        if seq_pos in pos_dict:
            repair_seq += pos_dict[seq_pos]
        else:
            repair_seq += wt_aa
    return True, repair_seq

def check_offset(seq, offset, hgvs_pros, urn, fallback_seq, verbosity = 0):

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

    if verbosity > 0:
        print(f'Call of check_offset {urn=}:\n{seq=}\n{fallback_seq=}')

    hit_rate_thresh = 0.5
    hit_rate = None
    best_offset = None

    mismatch_found, hit_rate = offset_loop(seq, offset, hgvs_pros, urn, verbosity = verbosity)

    if not mismatch_found:
        return offset, False, hit_rate, seq

    if hit_rate > hit_rate_thresh:
        hit_rate_thresh = hit_rate
        best_offset = offset

    original_offset = offset

    extracted_seq = extract_seq(hgvs_pros)

    if verbosity > 0:
        print(f'In check_offset {urn=}:\n{extracted_seq=}')

    (aligned_sequence, aligned_extracted_sequence) = (
        call_biopython_alignment(seq, extracted_seq)
    )

    offset = get_offset_from_aligned_seq(aligned_extracted_sequence)
    mismatch_found, hit_rate = offset_loop(seq, offset, hgvs_pros, urn, verbosity = verbosity)

    if not mismatch_found:
        return offset, False, hit_rate, seq

    if hit_rate > hit_rate_thresh:
        hit_rate_thresh = hit_rate
        best_offset = offset
    
    offset = 0
    mismatch_found, hit_rate = offset_loop(seq, offset, hgvs_pros, urn, verbosity = verbosity)

    if not mismatch_found:
        return offset, False, hit_rate, seq

    if hit_rate > hit_rate_thresh:
        hit_rate_thresh = hit_rate
        best_offset = offset

    offset = original_offset - 1
    mismatch_found, hit_rate = offset_loop(seq, offset, hgvs_pros, urn, verbosity = verbosity)

    if not mismatch_found:
        return offset, False, hit_rate, seq

    if hit_rate > hit_rate_thresh:
        hit_rate_thresh = hit_rate
        best_offset = offset

    offset = original_offset + 1
    mismatch_found, hit_rate = offset_loop(seq, offset, hgvs_pros, urn, verbosity = verbosity)

    if not mismatch_found:
        return offset, False, hit_rate, seq

    if hit_rate > hit_rate_thresh:
        hit_rate_thresh = hit_rate
        best_offset = offset

    offset = -1
    mismatch_found, hit_rate = offset_loop(seq, offset, hgvs_pros, urn, verbosity = verbosity)

    if not mismatch_found:
        return offset, False, hit_rate, seq

    if hit_rate > hit_rate_thresh:
        hit_rate_thresh = hit_rate
        best_offset = offset

    if best_offset is None:
        mismatch_found, hit_rate = offset_loop(fallback_seq, 0, hgvs_pros, urn, verbosity = verbosity)
        if hit_rate > hit_rate_thresh:
            best_offset = 0
            return best_offset, True, hit_rate, fallback_seq

    elif hit_rate_thresh < 1.0:
        successful, repaired_seq = repair_seq(seq, best_offset, hgvs_pros)
        if not successful:
            print(f'{repaired_seq=}')
            return None, False, hit_rate_thresh, None
        else:
            return best_offset, False, hit_rate_thresh, repaired_seq
    return best_offset, False, hit_rate_thresh, seq

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
