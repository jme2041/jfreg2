#!/usr/bin/env python3
#
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

def main(argv):
    """Main jfreg2 routine"""
    version = '2.0.0 (pre-release)'

    # Check FLIRT version
    FLIRT_MIN = (6, 0)
    capture = cmd('flirt', '-version', echo=False, capture_output=True)
    flirt_verstr = capture.stdout.decode().splitlines()[0]
    pattern = r'^FLIRT\s+version\s+(\d+)\.(\d+).*$'
    match = re.search(pattern, flirt_verstr)
    assert match is not None
    flirt_ver = tuple(map(int, match.groups()))
    if flirt_ver < FLIRT_MIN:
        raise RuntimeError('FLIRT %d.%d or newer is required' % FLIRT_MIN)

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

    g2 = parser.add_argument_group('options')

    g2.add_argument('--standard-brain',
        default='MNI152_T1_2mm_brain.nii.gz',
        metavar='STANDARD_BRAIN',
        help='''Standard space template (T1-weighted and brain-extracted);
        default: %(default)s''')

    g2.add_argument('--search',
        type=int,
        default=90,
        choices=[0, 90, 180],
        metavar='DEG',
        help='Search angle in degrees (%(choices)s; default: %(default)d)')

    g2.add_argument('--sigloss-thresh',
        type=float,
        default=10.0,
        metavar='SIGLOSS',
        help='Signal loss threshold (percent); default: %(default).2f')

    g2.add_argument('--overwrite',
        action='store_true',
        help='Replace output directory if it already exists')

    g2.add_argument('--version',
        action='version',
        version=version,
        help='Show version number and exit')

    g2.add_argument('--help',
        action='help',
        help='Show this help message and exit')

    opts = parser.parse_args(argv)

    # Look for the standard template (try FSLDIR as fallback)
    standard_brain = opts.standard_brain
    if not dset_exists(standard_brain):
        st2 = os.path.abspath(os.path.join(os.environ['FSLDIR'], 'data',
            'standard', standard_brain))
    if not dset_exists(st2):
        raise IOError('Could not find --standard-brain dataset: %s' %
            standard_brain)
    standard_brain = st2

    print('jfreg2 begins....')

    # Look for input datasets
    t1_brain = os.path.abspath(opts.t1_brain)
    if not dset_exists(t1_brain):
        raise IOError('Could not find --t1-brain dataset: %s' % t1_brain)
    print('--t1-brain: %s' % t1_brain)

    t1_head = os.path.abspath(opts.t1_head)
    if not dset_exists(t1_head):
        raise IOError('Could not find --t1-head dataset: %s' % t1_head)
    print('--t1-head: %s' % t1_head)

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

    # ####################
    # PREPROCESS FIELD MAP
    # ####################

    # These steps come from mainfeatreg's preprocessFieldmaps subroutine

    print('Preprocessing fieldmap....')

    fm_mag_brain_mask = os.path.join(outdir,
        strip_ext(os.path.basename(fm_mag_brain)) + '_mask')
    fm_mag_brain_mask_inv = fm_mag_brain_mask + '_inv'
    fm_mag_brain_mask_idx = fm_mag_brain_mask + '_idx'
    fm_mag_brain_mask50 = fm_mag_brain_mask + '50'
    fm_mag_brain_mask_ero = fm_mag_brain_mask + '_ero'
    fm_mag_brain_masked = fm_mag_brain_mask + 'ed'
    fm_mag_brain_masked_siglossed = fm_mag_brain_masked + '_siglossed'
    fm_rads_brain = os.path.join(outdir,
        strip_ext(os.path.basename(fm_rads)) + '_brain')
    fm_rads_brain_tmp_fmapfilt = fm_rads_brain + '_tmp_fmapfilt'
    fm_rads_brain_sigloss = fm_rads_brain + '_sigloss'


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

    # Get signal loss estimate
    cmd('sigloss',
        '-i', fm_rads_brain,
        '--te=%f' % opts.echo_time,
        '-m', fm_mag_brain_mask,
        '-s', fm_rads_brain_sigloss)

    sigloss_thresh = 1 - opts.sigloss_thresh/100.0
    cmd('fslmaths',
        fm_rads_brain_sigloss,
        '-mul', fm_mag_brain_masked,
        fm_mag_brain_masked_siglossed,
        '-odt', 'float')

    cmd('fslmaths',
        fm_rads_brain_sigloss,
        '-thr', str(sigloss_thresh),
        fm_rads_brain_sigloss,
        '-odt', 'float')

    cmd('imrm', fm_rads_brain_tmp_fmapfilt)
    cmd('imrm', fm_mag_brain_mask_ero)
    cmd('imrm', fm_mag_brain_mask50)
    cmd('imrm', fm_mag_brain_mask_idx)
    cmd('imrm', fm_mag_brain_mask_inv)

    # ################################################
    # SPATIAL NORMALIZATION OF THE T1-WEIGHTED DATASET
    # ################################################

    print('Spatially normalizing T1-weighted dataset....')

    t1_brain_to_standard_dset = os.path.join(outdir,
        strip_ext(os.path.basename(t1_brain)) + '_to_standard')
    t1_brain_to_standard_mat = t1_brain_to_standard_dset + '.mat'
    t1_head_to_standard_dset = os.path.join(outdir,
        strip_ext(os.path.basename(t1_head)) + '_to_standard')

    # Register the brain-extracted T1-weighted dataset to the standard
    # template. Use 12-parameter affine transformation.
    cmd('flirt',
        '-in', t1_brain,
        '-ref', standard_brain,
        '-out', t1_brain_to_standard_dset,
        '-omat', t1_brain_to_standard_mat,
        '-cost', 'corratio',
        '-dof', str(12),
        '-searchrx', str(-opts.search), str(opts.search),
        '-searchry', str(-opts.search), str(opts.search),
        '-searchrz', str(-opts.search), str(opts.search),
        '-interp', 'trilinear')

    assert(dset_exists(t1_brain_to_standard_dset))
    assert(os.path.isfile(t1_brain_to_standard_mat))

    # Apply the resulting transformation to the non-brain-extracted dataset
    cmd('flirt',
        '-in', t1_head,
        '-ref', standard_brain,
        '-out', t1_head_to_standard_dset,
        '-init', t1_brain_to_standard_mat,
        '-applyxfm',
        '-interp', 'trilinear')

    assert(dset_exists(t1_head_to_standard_dset))

if __name__ == '__main__':
    main(sys.argv[1:])
