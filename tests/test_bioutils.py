#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

Description:
    Unit tests for the bioutils library.

"""

import os
import sys
import unittest
import filecmp
import bioutils

from distutils import spawn

PWD = os.path.dirname(os.path.realpath(__file__))
PDB_DIR = os.path.join(PWD, "pdb_files")


class MainTestSuite(unittest.TestCase):
    """The main test suite for bioutils.

    """

    def test_dssp_installed(self):
        """Check if DSSP is installed and exists in the system PATH.

        """

        # Check if DSSP is installed
        self.assertTrue(spawn.find_executable("dssp"),
                "DSSP needs to be installed/linked in PATH"
                " and renamed/aliased to 'dssp'.")

        return

    def test_get_gunzipped_fn(self):
        """Check if gunzipped PDB files are successfully gunzipped and deleted.

        """

        from bioutils import get_gunzipped_fn

        pdb_fn_01 = os.path.join(PDB_DIR, "1ubq.pdb")
        pdb_fn_02 = os.path.join(PDB_DIR, "1ubq.pdb.gz")

        tmp_fn_01, delete_pdb_fn_01 = get_gunzipped_fn(pdb_fn_01)
        tmp_fn_02, delete_pdb_fn_02 = get_gunzipped_fn(pdb_fn_02)

        # Check if the temporary file deletion flag was correctly set.
        self.assertFalse(delete_pdb_fn_01)
        self.assertTrue(delete_pdb_fn_02)

        # The .pdb file shouldn't have changed.
        self.assertEqual(pdb_fn_01, tmp_fn_01)
        self.assertNotEqual(pdb_fn_02, tmp_fn_02)

        # Check if the gunzipping was successful.
        self.assertTrue(os.path.exists(tmp_fn_01))
        self.assertTrue(os.path.exists(tmp_fn_02))

        # Check if gunzipping didn't corrupt the files.
        self.assertTrue(filecmp.cmp(pdb_fn_01, tmp_fn_01))
        self.assertTrue(filecmp.cmp(pdb_fn_01, tmp_fn_02))

        if os.path.exists(tmp_fn_02):
            os.remove(tmp_fn_02)

        self.assertFalse(os.path.exists(tmp_fn_02))

        # Check if original pdb files weren't wrongly deleted.
        self.assertTrue(os.path.exists(pdb_fn_01))
        self.assertTrue(os.path.exists(pdb_fn_02))

        return
