#!/usr/bin/env python3
#
# epi2std.py: Warp EPI to standard (MNI152) space
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
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.
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
from jfreg2 import check_prefix
from jfreg2 import dset_exists
from jfreg2 import strip_ext
from jfreg2 import delete
from jfreg2 import check_flirt
import sys
import os
import argparse


def epi2std(argv):
    '''Warp EPI dataset to standard (MNI152) space'''

    check_flirt()

    # Look for the standard template (MNI152_T1_2mm_brain)
    standard_brain = os.path.join(os.environ['FSLDIR'], 'data', 'standard',
            'MNI152_T1_2mm_brain')
    if not dset_exists(standard_brain):
        raise IOError('Could not find standard template: %s' % standard_brain)

    parser = argparse.ArgumentParser(
            add_help=False,
            allow_abbrev=False,
            description=epi2std.__doc__)

    g1 = parser.add_argument_group('required switches')

    g1.add_argument('--epi',
            required=True,
            metavar='EPI',
            help='EPI dataset in subject space')

    g1.add_argument('--epi-base-to-t1-mat',
            required=True,
            metavar='MAT1',
            help='''Transformation matrix for EPI base volume to T1-weighted
            dataset in subject space (from epi2t1.py)''')

    g1.add_argument('--t1-to-mni152-mat',
            required=True,
            metavar='MAT2',
            help='''Transformation matrix for T1-weighted dataset to
            MNI152_T1_2mm (from flirt linear registration)''')

    g1.add_argument('--prefix',
            required=True,
            metavar='PREFIX',
            help='Output prefix')

    g2 = parser.add_argument_group('distortion correction')

    g2.add_argument('--epi-base-to-t1-shift',
            metavar='SHIFT',
            help='''Shift map for EPI dataset to T1-weighted dataset in subject
            space (from epi2t1.py)''')

    g2.add_argument('--unwarp-dir',
            choices=['x', 'y', 'z', 'x-', 'y-', 'z-'],
            metavar='DIR',
            help='Unwarp direction (%(choices)s)')

    g3 = parser.add_argument_group('options')

    g3.add_argument('--t1-to-mni152-warp',
            required=False,
            metavar='WARP',
            help='''Warp field for T1-weighted dataset to MNI152_T1_2mm (from
            fnirt nonlinear registration)''')

    g3.add_argument('--mc-cat',
            required=False,
            metavar='MC',
            help='Concatenated motion correction matrices (from epi2t1.py)')

    g3.add_argument('--final-interp',
            default='trilinear',
            choices=['nn', 'trilinear', 'sinc', 'spline'],
            metavar='INTERP',
            help='''Interpolation method for functional to standard warp
            (%(choices)s; default: %(default)s)''')

    g3.add_argument('--output-resolution',
            type=float,
            default=2.0,
            metavar='MM',
            help='''Output resolution for functional dataset (mm;
            default: %(default)1.1f)''')

    g3.add_argument('--overwrite',
            action='store_true',
            help='Replace output files if they already exist')

    g3.add_argument('--keep-all',
            action='store_true',
            help='Keep intermediate files')

    g3.add_argument('--version',
            action='version',
            version=__version__,
            help='Show version number and exit')

    g3.add_argument('--help',
            action='help',
            help='Show this help message and exit')

    opts = parser.parse_args(argv)

    print('epi2std begins....')

    check_prefix(opts.prefix)

    # Look for input datasets and files and set output dataset names
    epi = strip_ext(opts.epi)
    if not dset_exists(epi):
        raise IOError('Could not find --epi: %s' % epi)
    print('--epi: %s' % epi)

    epi_base_to_t1_mat = opts.epi_base_to_t1_mat
    if not os.path.isfile(epi_base_to_t1_mat):
        raise IOError('Could not find --epi-base-to-t1-mat: %s' %
                epi_base_to_t1_mat)
    print('--epi-base-to-t1-mat: %s' % epi_base_to_t1_mat)

    t1_to_mni152_mat = opts.t1_to_mni152_mat
    if not os.path.isfile(t1_to_mni152_mat):
        raise IOError('Could not find --t1-to-mni152-mat: %s' %
                t1_to_mni152_mat)
    print('--t1-to-mni152-mat: %s' % t1_to_mni152_mat)

    epi_to_mni152_dset = opts.prefix
    epi_base_to_mni152_mat = epi_to_mni152_dset + '.mat'
    epi_base_to_mni152_warp = epi_to_mni152_dset + '_warp'
    standard_brain_ores = opts.prefix + '_tmp_ores'

    out = [epi_to_mni152_dset, epi_base_to_mni152_mat, epi_base_to_mni152_warp]
    tmp = [standard_brain_ores]

    # Nonlinear T1 to MNI152?
    t1_to_mni152_warp = opts.t1_to_mni152_warp
    if t1_to_mni152_warp:
        t1_to_mni152_warp = strip_ext(t1_to_mni152_warp)
        if not dset_exists(t1_to_mni152_warp):
            raise IOError('Could not find --t1-to-mni152-warp: %s' %
                    t1_to_mni152_warp)
        print('--t1-to-mni152-warp: %s' % t1_to_mni152_warp)

    # Motion correction premat?
    mc_cat = opts.mc_cat
    if mc_cat:
        if not os.path.isfile(mc_cat):
            raise IOError('Could not find --mc-cat: %s' % mc_cat)
        print('--mc-cat: %s' % mc_cat)

    # Distortion correction with shift map?
    epi_base_to_t1_shift = None
    if(opts.epi_base_to_t1_shift or opts.unwarp_dir):
        shift = True
    else:
        shift = False

    if shift:
        if not opts.epi_base_to_t1_shift:
            parser.error(
                    'Distortion correction requires --epi-base-to-t1-shift')
        if not opts.unwarp_dir:
            parser.error('Distortion correction requires --unwarp-dir')

        epi_base_to_t1_shift = strip_ext(opts.epi_base_to_t1_shift)
        if not dset_exists(epi_base_to_t1_shift):
            raise IOError('Could not find --epi-base-to-t1-shift: %s' %
                    epi_base_to_t1_shift)
        print('--epi-base-to-t1-shift: %s' % epi_base_to_t1_shift)

    # Look for output datasets and overwrite if requested
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

    # Delete temporary datasets and files
    for t in tmp:
        delete(t)

    try:
        print('Generating standard space template at output resolution....')

        cmd('flirt',
                '-in', standard_brain,
                '-ref', standard_brain,
                '-out', standard_brain_ores,
                '-applyisoxfm', str(opts.output_resolution),
                '-interp', 'trilinear')

        print('Generating EPI to MNI152 warp field....')

        # Concatenate EPI base to T1 and T1 to MNI152 transformation matrices
        cmd('convert_xfm',
                '-omat', epi_base_to_mni152_mat,
                '-concat', t1_to_mni152_mat,
                epi_base_to_t1_mat)

        # Generate EPI to standard warp field
        convertwarp = [
                'convertwarp',
                '--ref=%s' % standard_brain_ores,
                '--out=%s' % epi_base_to_mni152_warp,
                '--relout'
        ]

        if t1_to_mni152_warp:
            convertwarp += [
                    '--premat=%s' % epi_base_to_t1_mat,
                    '--warp1=%s' % t1_to_mni152_warp
            ]
        else:
            convertwarp += ['--premat=%s' % epi_base_to_mni152_mat]

        if shift:
            convertwarp += [
                    '--shiftmap=%s' % epi_base_to_t1_shift,
                    '--shiftdir=%s' % opts.unwarp_dir
            ]

        cmd(*convertwarp)

        print('Warping EPI dataset to standard space....')

        applywarp = [
                'applywarp',
                '-i', epi,
                '-r', standard_brain_ores,
                '-o', epi_to_mni152_dset,
                '-w', epi_base_to_mni152_warp,
                '--rel',
                '--paddingsize=%s' % str(1),
                '--interp=%s' % opts.final_interp
        ]

        if mc_cat:
            applywarp += ['--premat=%s' % mc_cat]

        cmd(*applywarp)

    finally:
        if not opts.keep_all:
            for t in tmp:
                delete(t)

    print('epi2std ends....')


if __name__ == '__main__':
    epi2std(sys.argv[1:])
