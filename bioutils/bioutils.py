#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

Description:
    ...

Usage:
    bioutils.py

Options:
    ...

"""

import os
import sys
import shutil
import subprocess
import gzip
import numpy

from docopt import docopt
    
# Standard 20 amino acids.
AAS = \
    [("A","ALA"),("C","CYS"),("D","ASP"),("E","GLU"),
     ("F","PHE"),("G","GLY"),("H","HIS"),("I","ILE"),
     ("K","LYS"),("L","LEU"),("M","MET"),("N","ASN"),
     ("P","PRO"),("Q","GLN"),("R","ARG"),("S","SER"),
     ("T","THR"),("V","VAL"),("W","TRP"),("Y","TYR")]

# One-to-three character mapping.
AA_123 = dict(AAS)

# Three-to-one character mapping.
AA_321 = dict(((y,x) for x,y in AAS))

# Maximal accessible solvent area (ASA) of amino acids, used for computing the
# relative solvent accessibility (RSA), as published by
#
#   Sander & Rost, (1994), Proteins, 20:216-226
#
# These lines were modified from the `DSSP.py` file in the BioPython library.
MAX_ASA = {}
MAX_ASA["ALA"] = 106.0
MAX_ASA["CYS"] = 135.0
MAX_ASA["ASP"] = 163.0
MAX_ASA["GLU"] = 194.0
MAX_ASA["PHE"] = 197.0
MAX_ASA["GLY"] = 84.0
MAX_ASA["HIS"] = 184.0
MAX_ASA["ILE"] = 169.0
MAX_ASA["LYS"] = 205.0
MAX_ASA["LEU"] = 164.0
MAX_ASA["MET"] = 188.0
MAX_ASA["ASN"] = 157.0
MAX_ASA["PRO"] = 136.0
MAX_ASA["GLN"] = 198.0
MAX_ASA["ARG"] = 248.0
MAX_ASA["SER"] = 130.0
MAX_ASA["THR"] = 142.0
MAX_ASA["VAL"] = 142.0
MAX_ASA["TRP"] = 227.0
MAX_ASA["TYR"] = 222.0


# Use the mean value in MAX_ASA for non-standard amino acids.
DEFAULT_MAX_ASA = numpy.mean(MAX_ASA.values())


def get_gunzipped_fn(fn, tmp_dir="/tmp", verbose=False):
    """Generated temporary path to a gunzipped copy of a file. The original file
    path is returned if the file is not gzipped. If the file is gunzipped to a
    temporary file, it is up to the user to remove it.

    Parameters
    ----------
    fn : str
        Path to the gzip file (may not be gzipped).
    tmp_dir : str
        Path to the temporary directory (optional).

    Returns
    -------
    str
        Path to the temporary gunzipped file.
    bool
        True if the original file was a gzip file, False otherwise.

    """

    if not fn.endswith(".gz"):
        return fn, False

    # temporary gunzipped file location.
    gunzipped_fn = \
        os.path.join(tmp_dir, os.path.basename(os.path.splitext(fn)[0]))

    # pythonically gunzip the file into a temporary file.
    with open(gunzipped_fn, "wb") as tmp:
        shutil.copyfileobj(gzip.open(fn), tmp)

    assert(os.path.exists(gunzipped_fn))

    if verbose:
        sys.stderr.write("Created {}".format(gunzipped_fn) + os.linesep)

    return gunzipped_fn, True


def get_aa1(aa):
    if len(aa) == 1:
        return aa if aa in AA_123 else "X"

    assert(len(aa) == 3)

    return AA_321.get(aa, "X")


def get_aa3(aa):
    if len(aa) == 3:
        return aa if aa in AA_321 else "UNK"

    assert(len(aa) == 1)

    return AA_123.get(aa, "UNK")


def get_max_asa(aa):
    """Get the max. ASA value for a given one- or three-letter amino acid. A
    default ASA value is returned for non-standard amino acids.

    """

    if len(aa) == 1:
        max_asa = MAX_ASA.get(get_aa3(aa), DEFAULT_MAX_ASA)
    elif len(aa) == 3:
        max_asa = MAX_ASA.get(aa, DEFAULT_MAX_ASA)
    else:
        raise NotImplementedError

    return max_asa


def get_dssp_fields(pdb_path, chain_id):
    """Run DSSP and parse the output to obtain the various secondary structure
    and amino sequence properties for a single protein chain.  Gzipped PDB files
    are transparently gunzipped.

    Parameters
    ----------
    pdb_path : str
        Path to the PDB file.

    chain_id : str
        Chain identifier.

    Returns
    -------
    pdb_id : str
        The PDB ID of this protein chain. This may also be some other accession
        code (e.g. SCOP domain identifiers).

    aa_seq : str
        Single-letter amino acid sequence of this chain.

    ss_seq : str
        Secondary structure sequence, as determined by DSSP.

    bp1s : list of int
        Ordinal indices of bridge-partners (BP1 in DSSP).

    bp2s : list of int
        Ordinal indices of bridge-partners (BP2 in DSSP).

    asas : list of float
        Absolute solvent accessibilities, directly from DSSP.

    rsas : list of float
        Relative solvent accessibilities.

    """

    pdb_fn, delete_pdb_fn = get_gunzipped_fn(pdb_path)

    pdb_id = None
    dssp_cmd = "dssp {}".format(pdb_fn)

    with subprocess.Popen(dssp_cmd, shell=True,
            stdout=subprocess.PIPE).stdout as f:
        aa_seq = []
        ss_seq = []
        bp1s = []
        bp2s = []
        asas = []

        start = False
        for line in ( ln.rstrip() for ln in f ):
            if line.startswith("  #  RESIDUE"):
                start = True
                continue
            elif line.startswith("HEADER"):
                pdb_id = line[62:66].lower()
                continue

            if not start:
                continue

            # Chain break.
            if line[13] == "!":
                continue

            if line[11] != chain_id:
                continue

            # Remember, lower-case amino acids are cysteines.
            aa = get_aa1("C" if line[13].islower() else line[13])
            ss = "-" if line[16] == " " else line[16]

            aa_seq.append(aa)
            ss_seq.append(ss)

            bp1s.append(int(line[25:29].strip()))
            bp2s.append(int(line[29:33].strip()))
            asas.append(float(line[35:38].strip()))

    # Convert some lists to numpy arrays.
    bp1s = numpy.array(bp1s)
    bp2s = numpy.array(bp2s)
    asas = numpy.array(asas)

    # Compute the RSA values.
    rsas = asas / numpy.array(map(get_max_asa, aa_seq))

    if os.path.exists(pdb_fn) and delete_pdb_fn:
        os.remove(pdb_fn)
        assert(not os.path.exists(pdb_fn))

    assert(pdb_id is not None)

    return pdb_id, "".join(aa_seq), "".join(ss_seq), bp1s, bp2s, asas, rsas


if __name__ == '__main__':
    opts = docopt(__doc__)

    print get_dssp_fields("1UBQ.pdb.gz", "A")
