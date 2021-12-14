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

import sys
import os
import re
import shutil
import argparse
import subprocess
import glob

__version__ = '2.0.0 (pre-release)'


def cmd(*cmd, echo=True, capture_output=False):
    """Run a shell command (echoing it to stdout by default)"""
    if(echo):
        print(' '.join(cmd))
    return subprocess.run(cmd, check=True, capture_output=capture_output)


def dset_exists(path):
    """Test whether a dataset exists"""
    # If there is a file extension, check that the file exists
    if(path.endswith('.nii') or path.endswith('.nii.gz')):
        return os.path.isfile(path)
    elif(path.endswith('.hdr')):
        if not os.path.isfile(path):
            return False
        if os.path.isfile(path.replace('.hdr', '.img')):
            return True
        if os.path.isfile(path.replace('.hdr', '.img.gz')):
            return True
        return False
    elif(path.endswith('.hdr.gz')):
        if not os.path.isfile(path):
            return False
        if os.path.isfile(path.replace('.hdr.gz'), '.img.gz'):
            return True
        return False
    else:
        # If there is no file extension, try adding valid extensions
        if dset_exists(path + '.nii'):
            return True
        if dset_exists(path + '.nii.gz'):
            return True
        if dset_exists(path + '.hdr'):
            return True
        if dset_exists(path + '.hdr.gz'):
            return True
        return False


def strip_ext(path):
    """Strip NIfTI extensions"""
    path = path.removesuffix('.hdr')
    path = path.removesuffix('.hdr.gz')
    path = path.removesuffix('.nii')
    path = path.removesuffix('.nii.gz')
    return path


def delete(x):
    if dset_exists(x):
        cmd('imrm', x, echo=False)
    elif os.path.isdir(x):
        shutil.rmtree(x)
    elif os.path.isfile(x):
        os.remove(x)


def check_flirt():
    """Check FLIRT version"""
    FLIRT_MIN = (6, 0)
    capture = cmd('flirt', '-version', echo=False, capture_output=True)
    flirt_verstr = capture.stdout.decode().splitlines()[0]
    pattern = r'^FLIRT\s+version\s+(\d+)\.(\d+).*$'
    match = re.search(pattern, flirt_verstr)
    assert match is not None
    flirt_ver = tuple(map(int, match.groups()))
    if flirt_ver < FLIRT_MIN:
        raise RuntimeError('FLIRT %d.%d or newer is required' % FLIRT_MIN)

def placeholder(argv):
    # Look for the standard template (MNI152_T1_2mm_brain)
    standard_brain = os.path.abspath(os.path.join(os.environ['FSLDIR'], 'data',
            'standard', 'MNI152_T1_2mm_brain'))
    if not dset_exists(standard_brain):
        raise IOError('Could not find standard template: %s' % standard_brain)

    # Get FLIRT schedule for BBR
    bbr_schedule = os.path.abspath(os.path.join(os.environ['FSLDIR'],
        'etc', 'flirtsch', 'bbr.sch'))
    if not os.path.isfile(bbr_schedule):
        raise IOError('Could not find FLIRT BBR schedule: %s' % bbr_schedule)

    # Generate dict that maps unwarp direction to flirt phase encode direction
    pe_dir = {
        'x': 1,
        'y': 2,
        'z': 3,
        'x-': -1,
        'y-': -2,
        'z-': -3
    }

    parser = argparse.ArgumentParser(
        prog='jfreg2',
        add_help=False,
        allow_abbrev=False,
        description='Joint fMRI Registration',
        usage='%(prog)s [options] <switches> FUNC [FUNC ...]')

    g1 = parser.add_argument_group('required switches')

    g1.add_argument('--outdir',
        required=True,
        metavar='OUTDIR',
        help='Output directory (will be created)')

    g1.add_argument('--t1-brain',
        required=True,
        metavar='T1_BRAIN',
        help='T1-weighted dataset in subject space (brain-extracted)')

    g1.add_argument('--t1-head',
        required=True,
        metavar='T1_HEAD',
        help='T1-weighted dataset in subject space (not brain-extracted)')

    g1.add_argument('--t1-brain-wmseg',
        required=True,
        metavar='WMSEG',
        help='White matter map from FAST segmentation of T1-weighted dataset')

    g1.add_argument('--t1-to-standard-mat',
        required=True,
        metavar='MAT',
        help='''Transformation matrix for T1-weighted dataset to MNI152_T1_2mm
        (from flirt linear registration)''')

    g1.add_argument('--fieldmap-rads',
        required=True,
        metavar='FM_RADS',
        help='Field map in units of radians/second (phase-reconstructed)')

    g1.add_argument('--fieldmap-mag-brain',
        required=True,
        metavar='FM_MAG_BRAIN',
        help='Magnitude-reconstructed field map (brain-extracted)')

    g1.add_argument('--fieldmap-mag-head',
        required=True,
        metavar='FM_MAG_HEAD',
        help='Magnitude-reconstructed field map (not brain-extracted)')

    g1.add_argument('--echo-time',
        type=float,
        required=True,
        metavar='TE',
        help='Echo time (in seconds) of EPI used for functional MRI')

    g1.add_argument('--echo-spacing',
        type=float,
        required=True,
        metavar='ES',
        help='Effective echo spacing (in seconds) of the EPI used for fMRI')

    g1.add_argument('--unwarp-dir',
        required=True,
        choices=['x', 'y', 'z', 'x-', 'y-', 'z-'],
        metavar='UNWARP_DIR',
        help='Unwarp direction (%(choices)s)')

    g2 = parser.add_argument_group('options')

    g2.add_argument('--t1-to-standard-warp',
        required=False,
        metavar='WARP',
        help='''Warp field for T1-weighted dataset to MNI152_T1_2mm (from fnirt
        nonlinear registration)''')

    g2.add_argument('--base-volume',
        type=int,
        default=0,
        metavar='VOL',
        help='Functional base volume number (0-indexed); default: %(default)d')

    g2.add_argument('--output-resolution',
        type=float,
        default=2.0,
        metavar='MM',
        help='''Output resolution for functional dataset (mm;
        default: %(default)1.1f)''')

    g2.add_argument('--search',
        type=int,
        default=90,
        choices=[0, 90, 180],
        metavar='DEG',
        help='Search angle in degrees (%(choices)s; default: %(default)d)')

    g2.add_argument('--final-interp',
        default='trilinear',
        choices=['nn', 'trilinear', 'sinc', 'spline'],
        metavar='INTERP',
        help='''Interpolation method for functional to standard warp
        (%(choices)s; default: %(default)s)''')

    g2.add_argument('--overwrite',
        action='store_true',
        help='Replace output directory if it already exists')

    g2.add_argument('--keep-all',
        action='store_true',
        help='Keep intermediate files')

    g2.add_argument('--version',
        action='version',
        version=version,
        help='Show version number and exit')

    g2.add_argument('--help',
        action='help',
        help='Show this help message and exit')

    parser.add_argument('func',
        nargs='+',
        metavar='FUNC',
        help='Functional MRI dataset(s)')

    opts = parser.parse_args(argv)

    print('jfreg2 begins....')

    # Look for T1-to-standard transformation matrix
    t1_to_standard_mat = os.path.abspath(opts.t1_to_standard_mat)
    if not os.path.isfile(t1_to_standard_mat):
        raise IOError('Could not find --t1-to-standard-mat file: %s' %
            t1_to_standard_mat)
    print('--t1-to-standard-mat: %s' % t1_to_standard_mat)

    # Look for T1-to-standard warp field
    t1_to_standard_warp = opts.t1_to_standard_warp
    if t1_to_standard_warp:
        t1_to_standard_warp = os.path.abspath(t1_to_standard_warp)
        if not dset_exists(t1_to_standard_warp):
            raise IOError('Could not find --t1-to-standard-warp file: %s' %
                t1_to_standard_warp)
    print('--t1-to-standard-warp: %s' % t1_to_standard_warp)

    # Look for input datasets
    t1_brain = os.path.abspath(opts.t1_brain)
    if not dset_exists(t1_brain):
        raise IOError('Could not find --t1-brain dataset: %s' % t1_brain)
    print('--t1-brain: %s' % t1_brain)

    t1_head = os.path.abspath(opts.t1_head)
    if not dset_exists(t1_head):
        raise IOError('Could not find --t1-head dataset: %s' % t1_head)
    print('--t1-head: %s' % t1_head)

    wmseg = os.path.abspath(opts.t1_brain_wmseg)
    if not dset_exists(wmseg):
        raise IOError('Could not find --t1-brain-wmseg dataset: %s' % wmseg)
    print('--t1-brain-wmseg: %s' % wmseg)

    fm_rads = os.path.abspath(opts.fieldmap_rads)
    if not dset_exists(fm_rads):
        raise IOError('Could not find --fieldmap-rads dataset: %s' % fm_rads)
    print('--fieldmap-rads: %s' % fm_rads)

    fm_mag_brain = os.path.abspath(opts.fieldmap_mag_brain)
    if not dset_exists(fm_mag_brain):
        raise IOError('Could not find --fieldmap-mag-brain dataset: %s' %
            fm_mag_brain)
    print('--fieldmap-mag-brain: %s' % fm_mag_brain)

    fm_mag_head = os.path.abspath(opts.fieldmap_mag_head)
    if not dset_exists(fm_mag_head):
        raise IOError('Could not find --fieldmap-mag-head dataset: %s' %
            fm_mag_head)
    print('--fieldmap-mag-head: %s' % fm_mag_head)

    for f in opts.func:
        if not dset_exists(f):
            raise IOError('Could not find functional dataset: %s' % f)
        print('functional dataset: %s' % f)

    # Warn if --fieldmap-mag-brain appears not to be brain-extracted
    print('Checking fieldmap....')
    capture = cmd('fslstats', fm_mag_brain, '-V', capture_output=True)
    nonzero_voxels = float(capture.stdout.decode().split()[0])
    print('Non-zero voxels: %d' % int(nonzero_voxels))
    capture = cmd('fslstats', fm_mag_brain, '-v', capture_output=True)
    total_voxels = float(capture.stdout.decode().split()[0])
    print('Total voxels: %d' % int(total_voxels))
    prop_nonzero = nonzero_voxels / total_voxels;
    print('Proportion non-zero: %f' % prop_nonzero)
    if(prop_nonzero > 0.9):
        print('WARNING: High proportion of non-zero voxels in brain-extracted '
            'field map')

    # Create output directory
    outdir = os.path.abspath(opts.outdir)
    if os.path.isdir(outdir):
        if opts.overwrite:
            shutil.rmtree(outdir)
        else:
            raise IOError('Output directory already exists: %s' % outdir)
    os.mkdir(outdir)

    # Temporary datasets and files
    tmp = []

    # Create a dataset in standard space at the desired output resolution
    standard_brain_ores = os.path.join(outdir,
        os.path.basename(standard_brain)) + '_ores'

    tmp += [ standard_brain_ores ]

    cmd('flirt',
        '-in', standard_brain,
        '-ref', standard_brain,
        '-out', standard_brain_ores,
        '-applyisoxfm', str(opts.output_resolution),
        '-interp', 'trilinear')

    # ###############################
    # LOOP ON FUNCTIONAL MRI DATASETS
    # ###############################

    for f in opts.func:
        print('Processing functional dataset: %s....' % f)

        f_base = os.path.join(outdir, strip_ext(os.path.basename(f)) + '_base')

        # Extract base volume and convert to floating-point
        cmd('fslroi',
            f,
            f_base,
            str(opts.base_volume),
            str(1))

        cmd('fslmaths',
            f_base,
            f_base,
            '-odt',
            'float')

        # ##########################################
        # FUNCTIONAL TO STRUCTURAL (T1) REGISTRATION
        # ##########################################

        print('Registering functional base volume to T1-weighted dataset....')

        f_base_to_t1_init_mat = f_base + '_to_t1_init.mat'
        f_base_to_t1_dset = f_base + '_to_t1'
        f_base_to_t1_warp = f_base_to_t1_dset + '_warp'
        f_base_to_t1_mat = f_base_to_t1_dset + '.mat'
        f_base_to_t1_inv_mat = f_base_to_t1_dset + '_inv.mat'
        f_base_fm_rads_to_base_dset = f_base + '_fm_rads_to_base'
        f_base_fm_rads_to_base_mat = f_base_fm_rads_to_base_dset + '.mat'
        f_base_fm_rads_to_base_mask = f_base_fm_rads_to_base_dset + '_mask'
        f_base_fm_rads_to_base_shift = f_base_fm_rads_to_base_dset + '_shift'

        tmp += [
            f_base_to_t1_init_mat,
            f_base_to_t1_inv_mat,
            f_base_fm_rads_to_base_dset,
            f_base_fm_rads_to_base_mat,
            f_base_fm_rads_to_base_mask
        ]

        # Initial functional base to T1 (brain) alignment using 6 DOF
        cmd('flirt',
            '-in', f_base,
            '-ref', t1_brain,
            '-omat', f_base_to_t1_init_mat,
            '-cost', 'corratio',
            '-dof', str(6),
            '-searchrx', str(-opts.search), str(opts.search),
            '-searchry', str(-opts.search), str(opts.search),
            '-searchrz', str(-opts.search), str(opts.search))

        # Register functional base volume to T1 (head) using BBR and field map
        cmd('flirt',
            '-in', f_base,
            '-ref', t1_head,
            '-omat', f_base_to_t1_mat,
            '-cost', 'bbr',
            '-dof', str(6),
            '-wmseg', wmseg,
            '-init', f_base_to_t1_init_mat,
            '-nosearch',
            '-schedule', bbr_schedule,
            '-echospacing', str(opts.echo_spacing),
            '-pedir', str(pe_dir[opts.unwarp_dir]),
            '-fieldmap', fm_rads_brain_to_t1_head)

        # Generate warp fields for use with other registrations
        cmd('convert_xfm',
            '-omat', f_base_to_t1_inv_mat,
            '-inverse', f_base_to_t1_mat)

        cmd('convert_xfm',
            '-omat', f_base_fm_rads_to_base_mat,
            '-concat',
            f_base_to_t1_inv_mat,
            fm_to_t1_head_mat)

        cmd('applywarp',
            '-i', fm_rads_brain_unmasked,
            '-r', f_base,
            '--premat=%s' % f_base_fm_rads_to_base_mat,
            '-o', f_base_fm_rads_to_base_dset)

        cmd('fslmaths',
            f_base_fm_rads_to_base_dset,
            '-abs',
            '-bin',
            f_base_fm_rads_to_base_mask)

        cmd('fugue',
            '--loadfmap=%s' % f_base_fm_rads_to_base_dset,
            '--mask=%s' % f_base_fm_rads_to_base_mask,
            '--saveshift=%s' % f_base_fm_rads_to_base_shift,
            '--unmaskshift',
            '--dwell=%s' % str(opts.echo_spacing),
            '--unwarpdir=%s' % opts.unwarp_dir)

        cmd('convertwarp',
            '-r', t1_head,
            '-s', f_base_fm_rads_to_base_shift,
            '--postmat=%s' % f_base_to_t1_mat,
            '-o', f_base_to_t1_warp,
            '--shiftdir=%s' % opts.unwarp_dir,
            '--relout')

        # Apply the warp field: This does the final base to T1 registration
        cmd('applywarp',
            '-i', f_base,
            '-r', t1_head,
            '-o', f_base_to_t1_dset,
            '-w', f_base_to_t1_warp,
            '--interp=spline',
            '--rel')

        # #################
        # MOTION CORRECTION
        # #################

        print('Motion correcting the functional dataset....')

        f_mc = os.path.join(outdir, strip_ext(os.path.basename(f)) + '_mc')
        f_mc_mat = f_mc + '.mat'
        f_mc_cat = f_mc + '.cat'

        tmp += [ f_mc, f_mc_mat ]

        cmd('mcflirt',
            '-in', f,
            '-out', f_mc,
            '-reffile', f_base,
            '-mats',
            '-plots')

        # Concatenate the transformation matrices into one big file
        with open(f_mc_cat, 'wb') as catfile:
            for i in sorted(glob.glob('%s/*' % f_mc_mat)):
                with open(i, 'rb') as matfile:
                    shutil.copyfileobj(matfile, catfile)

        # #########################################
        # WARP FUNCTIONAL DATASET TO STANDARD SPACE
        # #########################################

        print('Warping functional dataset to standard space....')

        f_to_standard_dset = os.path.join(outdir,
            strip_ext(os.path.basename(f)) + '_to_standard')
        f_base_to_standard_mat = f_base + '_to_standard.mat'
        f_base_to_standard_warp = f_base + '_to_standard_warp'

        # Concatenate base to T1 and T1 to standard transformation matrices
        cmd('convert_xfm',
            '-omat', f_base_to_standard_mat,
            '-concat', t1_to_standard_mat,
            f_base_to_t1_mat)

        # Generate functional to standard warp field
        if t1_to_standard_warp:
            cmd('convertwarp',
                '--ref=%s' % standard_brain_ores,
                '--shiftmap=%s' % f_base_fm_rads_to_base_shift,
                '--premat=%s' % f_base_to_t1_mat,
                '--warp1=%s' % t1_to_standard_warp,
                '--out=%s' % f_base_to_standard_warp,
                '--relout',
                '--shiftdir=%s' % opts.unwarp_dir)
        else:
            cmd('convertwarp',
                '--ref=%s' % standard_brain_ores,
                '--shiftmap=%s' % f_base_fm_rads_to_base_shift,
                '--premat=%s' % f_base_to_standard_mat,
                '--out=%s' % f_base_to_standard_warp,
                '--relout',
                '--shiftdir=%s' % opts.unwarp_dir)

        # Apply the warp. Use motion correction matrices as premat.
        cmd('applywarp',
            '-i', f,
            '--premat=%s' % f_mc_cat,
            '-r', standard_brain_ores,
            '-o', f_to_standard_dset,
            '-w', f_base_to_standard_warp,
            '--rel',
            '--paddingsize=%s' % str(1),
            '--interp=%s' % opts.final_interp)

    if not opts.keep_all:
        for t in tmp:
            if dset_exists(t):
                cmd('imrm', t, echo=False)
            elif os.path.isdir(t):
                shutil.rmtree(t)
            elif os.path.isfile(t):
                os.remove(t)
            else:
                raise RuntimeError('Could not remove %s' % t)

    print('jfreg2 ends....')
