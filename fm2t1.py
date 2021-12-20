#!/usr/bin/env python3
#
# fm2t1.py: Register field map to T1-weighted brain
# jfreg2: Joint fMRI Registration (Version 2)
#
# Copyright (c) 2021, Jeffrey M. Engelmann
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors
#   may be used to endorse or promote products derived from this software
#   without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# jfreg2 uses the FMRI Software Library (FSL), which is developed by the
# Analysis Group at the Wellcome Center for Integrative Neuroimaging at the
# University of Oxford. FSL is released for non-commercial use under an open-
# source license. https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Licence
#

from jfreg2 import __version__
from jfreg2 import cmd
from jfreg2 import dset_exists
from jfreg2 import strip_ext
from jfreg2 import delete
from jfreg2 import check_flirt
import sys
import os
import argparse


def fm2t1(argv):
    '''
Register field map to T1-weighted brain

Output datasets and files:
PREFIX_fm_mag_brain_to_t1_brain      Fieldmap (magnitude) registered to T1
PREFIX_fm_rads_brain_to_t1_brain     Fieldmap (phase) registered to T1
PREFIX_fm_mag_brain_to_t1_brain.mat  Transformation matrix (fieldmap to T1)'''

    check_flirt()

    parser = argparse.ArgumentParser(
            prog='fm2t1.py',
            add_help=False,
            allow_abbrev=False,
            description=fm2t1.__doc__,
            usage='%(prog)s [options] <switches>',
            formatter_class=argparse.RawTextHelpFormatter)

    g1 = parser.add_argument_group('required switches')

    g1.add_argument('--t1-brain',
            required=True,
            metavar='T1_BRAIN',
            help='T1-weighted dataset in subject space (brain-extracted)')

    g1.add_argument('--fm-rads',
            required=True,
            metavar='FM_RADS',
            help='Field map in units of radians/second (not brain-extracted)')

    g1.add_argument('--fm-mag-brain',
            required=True,
            metavar='FM_MAG_BRAIN',
            help='Magnitude-reconstructed field map (brain-extracted)')

    g1.add_argument('--unwarp-dir',
        required=True,
        choices=['x', 'y', 'z', 'x-', 'y-', 'z-'],
        metavar='DIR',
        help='Unwarp direction (%(choices)s)')

    g1.add_argument('--prefix',
            required=True,
            metavar='PREFIX',
            help='Output prefix')

    g2 = parser.add_argument_group('options')

    g2.add_argument('--search',
        type=int,
        default=90,
        choices=[0, 90, 180],
        metavar='DEG',
        help='Search angle in degrees (%(choices)s; default: %(default)d)')

    g2.add_argument('--overwrite',
            action='store_true',
            help='Replace output files if they already exist')

    g2.add_argument('--keep-all',
            action='store_true',
            help='Keep intermediate files')

    g2.add_argument('--version',
            action='version',
            version=__version__,
            help='Show version number and exit')

    g2.add_argument('--help',
            action='help',
            help='Show this help message and exit')

    opts = parser.parse_args(argv)

    print('jfreg2 begins....')

    # Look for input datasets
    t1_brain = strip_ext(opts.t1_brain)
    if not dset_exists(t1_brain):
        raise IOError('Could not find --t1-brain: %s' % t1_brain)
    print('--t1-brain: %s' % t1_brain)

    fm_rads = strip_ext(opts.fm_rads)
    if not dset_exists(fm_rads):
        raise IOError('Could not find --fm-rads: %s' % fm_rads)
    print('--fm-rads: %s' % fm_rads)

    fm_mag_brain = strip_ext(opts.fm_mag_brain)
    if not dset_exists(fm_mag_brain):
        raise IOError('Could not find --fm-mag-brain: %s' % fm_mag_brain)
    print('--fm-mag-brain: %s' % fm_mag_brain)

    fm_mag_brain_to_t1_brain_dset = opts.prefix + '_fm_mag_brain_to_t1_brain'
    fm_mag_brain_to_t1_brain_mat = fm_mag_brain_to_t1_brain_dset + '.mat'
    fm_mag_brain_mask = opts.prefix + '_fm_mag_brain_mask'
    fm_mag_brain_mask_inv = fm_mag_brain_mask + '_inv'
    fm_mag_brain_mask_idx = fm_mag_brain_mask + '_idx'
    fm_mag_brain_mask50 = fm_mag_brain_mask + '50'
    fm_mag_brain_mask_ero = fm_mag_brain_mask + '_ero'
    fm_mag_brain_masked = fm_mag_brain_mask + 'ed'
    fm_rads_brain = opts.prefix + '_fm_rads_brain'
    fm_rads_brain_tmp_fmapfilt = fm_rads_brain + '_tmp_fmapfilt'
    fm_rads_brain_mask = fm_rads_brain + '_mask'
    fm_rads_brain_unmasked = fm_rads_brain + '_unmasked'
    fm_rads_brain_to_t1_brain = fm_rads_brain + '_to_t1_brain'
    fm_rads_brain_to_t1_brain_pad0 = fm_rads_brain_to_t1_brain + '_pad0'
    fm_rads_brain_to_t1_brain_inner_mask = (fm_rads_brain_to_t1_brain +
            '_inner_mask')

    # Look for output datasets and overwrite if requested
    out = [
            fm_mag_brain_to_t1_brain_dset,
            fm_mag_brain_to_t1_brain_mat,
            fm_rads_brain_to_t1_brain
    ]

    for o in out:
        print('Output: %s' % o)
        if dset_exists(o):
            if opts.overwrite:
                cmd('imrm', o, echo=False)
            else:
                raise IOError('Output dataset already exists: %s' % o)
        elif os.path.isfile(o):
            if opts.overwrite:
                os.remove(o)
            else:
                raise IOError('Output dataset already exists: %s' % o)

    # Temporary datasets and files
    tmp = [
            fm_mag_brain_mask,
            fm_mag_brain_mask_inv,
            fm_mag_brain_mask_idx,
            fm_mag_brain_mask50,
            fm_mag_brain_mask_ero,
            fm_mag_brain_masked,
            fm_rads_brain,
            fm_rads_brain_tmp_fmapfilt,
            fm_rads_brain_mask,
            fm_rads_brain_unmasked,
            fm_rads_brain_to_t1_brain_pad0,
            fm_rads_brain_to_t1_brain_inner_mask
    ]

    for t in tmp:
        delete(t)

    try:
        # Preprocess field map
        # These steps come from mainfeatreg's preprocessFieldmaps subroutine

        print('Preprocessing fieldmap....')

        cmd('fslmaths',
                fm_mag_brain,
                '-bin',
                fm_mag_brain_mask,
                '-odt', 'short')

        cmd('fslmaths',
                fm_rads,
                '-abs',
                '-bin',
                '-mas', fm_mag_brain_mask,
                '-mul', str(-1),
                '-add', str(1),
                '-bin',
                fm_mag_brain_mask_inv)

        cmd('cluster',
                '-i', fm_mag_brain_mask_inv,
                '-t', str(0.5),
                '--no_table',
                '-o', fm_mag_brain_mask_idx)

        capture = cmd('fslstats',
                fm_mag_brain_mask_idx,
                '-R',
                capture_output=True)
        maxidx = float(capture.stdout.decode().split()[1])

        cmd('fslmaths',
                fm_mag_brain_mask_idx,
                '-thr', str(maxidx),
                '-bin',
                '-mul', str(-1),
                '-add', str(1),
                '-bin',
                '-mas', fm_mag_brain_mask,
                fm_mag_brain_mask)

        capture = cmd('fslstats',
                fm_rads,
                '-k', fm_mag_brain_mask,
                '-P', str(50),
                capture_output=True)
        mean = float(capture.stdout.decode())

        cmd('fslmaths',
                fm_rads,
                '-sub', str(mean),
                '-mas', fm_mag_brain_mask,
                fm_rads_brain)

        capture = cmd('fslstats',
                fm_mag_brain,
                '-P', str(98),
                capture_output=True)
        thresh50 = float(capture.stdout.decode())/2.0

        cmd('fslmaths',
                fm_mag_brain,
                '-thr', str(thresh50),
                fm_mag_brain_mask50)

        cmd('fslmaths',
                fm_mag_brain_mask,
                '-ero',
                fm_mag_brain_mask_ero)

        cmd('fslmaths',
                fm_mag_brain_mask_ero,
                '-add', fm_mag_brain_mask50,
                '-thr', str(0.5),
                '-bin',
                fm_mag_brain_mask)

        cmd('fslmaths',
                fm_rads_brain,
                '-mas', fm_mag_brain_mask,
                fm_rads_brain)

        cmd('fslmaths',
                fm_mag_brain,
                '-mas', fm_mag_brain_mask,
                fm_mag_brain_masked)

        cmd('fslmaths',
                fm_mag_brain_mask,
                '-ero',
                fm_mag_brain_mask_ero)

        cmd('fugue',
                '--loadfmap=%s' % fm_rads_brain,
                '--savefmap=%s' % fm_rads_brain_tmp_fmapfilt,
                '--mask=%s' % fm_mag_brain_mask,
                '--despike',
                '--despikethreshold=2.1')

        cmd('fslmaths',
                fm_rads_brain,
                '-sub', fm_rads_brain_tmp_fmapfilt,
                '-mas', fm_mag_brain_mask_ero,
                '-add', fm_rads_brain_tmp_fmapfilt,
                fm_rads_brain)

        capture = cmd('fslstats',
                fm_rads_brain,
                '-k', fm_mag_brain_mask,
                '-P', str(50),
                capture_output=True)
        median = float(capture.stdout.decode())

        cmd('fslmaths',
                fm_rads_brain,
                '-sub', str(median),
                '-mas', fm_mag_brain_mask,
                fm_rads_brain)

        # Register field map to T1-weighted image
        # These steps are based on FSL's epi_reg
        #
        # Brain-extracted field map to brain-extracted T1 is used as the final
        # fieldmap to T1 registration, rather than an initial registration
        # followed by registration of the field map head to T1 head. This is
        # because of inconsistent performance of the fieldmap head to T1 head
        # registration.
        # See https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=FSL;360f54df.1910

        print('Registering fieldmap to T1-weighted dataset....')

        cmd('flirt',
                '-in', fm_mag_brain_masked,
                '-ref', t1_brain,
                '-out', fm_mag_brain_to_t1_brain_dset,
                '-omat', fm_mag_brain_to_t1_brain_mat,
                '-cost', 'corratio',
                '-dof', str(6),
                '-searchrx', str(-opts.search), str(opts.search),
                '-searchry', str(-opts.search), str(opts.search),
                '-searchrz', str(-opts.search), str(opts.search))

        cmd('fslmaths',
                fm_mag_brain_masked,
                '-abs',
                '-bin',
                fm_rads_brain_mask)

        cmd('fslmaths',
                fm_rads_brain,
                '-abs',
                '-bin',
                '-mul', fm_rads_brain_mask,
                fm_rads_brain_mask)

        cmd('fugue',
                '--loadfmap=%s' % fm_rads_brain,
                '--mask=%s' % fm_rads_brain_mask,
                '--unmaskfmap',
                '--savefmap=%s' % fm_rads_brain_unmasked,
                '--unwarpdir=%s' % opts.unwarp_dir)

        cmd('applywarp',
                '-i', fm_rads_brain_unmasked,
                '-r', t1_brain,
                '--premat=%s' % fm_mag_brain_to_t1_brain_mat,
                '-o', fm_rads_brain_to_t1_brain_pad0)

        cmd('fslmaths',
                fm_rads_brain_to_t1_brain_pad0,
                '-abs',
                '-bin',
                fm_rads_brain_to_t1_brain_inner_mask)

        cmd('fugue',
                '--loadfmap=%s' % fm_rads_brain_to_t1_brain_pad0,
                '--mask=%s' % fm_rads_brain_to_t1_brain_inner_mask,
                '--unmaskfmap',
                '--savefmap=%s' % fm_rads_brain_to_t1_brain,
                '--unwarpdir=%s' % opts.unwarp_dir)

    finally:
        if not opts.keep_all:
            for t in tmp:
                delete(t)

    print('fm2t1 ends....')


if __name__ == '__main__':
    fm2t1(sys.argv[1:])
