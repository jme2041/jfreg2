#!/usr/bin/env python3
#
# epi2t1.py: Register EPI to T1-weighted brain (optionally with field map)
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
import shutil
import argparse
import glob
import csv


def epi2t1(argv):
    '''Motion-correct and register EPI dataset to T1-weighted dataset'''

    check_flirt()

    parser = argparse.ArgumentParser(
            add_help=False,
            allow_abbrev=False,
            description=epi2t1.__doc__)

    g1 = parser.add_argument_group('required switches')

    g1.add_argument('--t1-brain',
            required=True,
            metavar='T1_BRAIN',
            help='T1-weighted dataset in subject space (brain-extracted)')

    g1.add_argument('--t1-head',
            required=True,
            metavar='T1_HEAD',
            help='T1-weighted dataset in subject space (not brain-extracted)')

    g1.add_argument('--epi',
            required=True,
            metavar='EPI',
            help='EPI dataset in subject space (not brain-extracted)')

    g1.add_argument('--prefix',
            required=True,
            metavar='PREFIX',
            help='Output prefix')

    g2 = parser.add_argument_group('boundary-based registration (BBR)')

    g2.add_argument('--t1-brain-wmseg',
            metavar='WMSEG',
            help='White matter mask from FAST segmentation of T1 brain')

    g3 = parser.add_argument_group('distortion correction (requires BBR)')

    g3.add_argument('--fm-rads-unmasked',
            metavar='FM_RADS_UNMASKED',
            help='Unmasked field map in units of radians/s (from fm2t1.py)')

    g3.add_argument('--fm-rads-brain-t1',
            metavar='FM_RADS_BRAIN_T1',
            help='Field map (radians/s) aligned to T1 brain (from fm2t1.py)')

    g3.add_argument('--fm-to-t1-mat',
            metavar='FM_TO_T1_MAT',
            help='Field map to T1 brain transformation matrix (from fm2t1.py)')

    g3.add_argument('--echo-spacing',
            type=float,
            metavar='ES',
            help='Effective echo spacing (seconds) of the EPI used for fMRI')

    g3.add_argument('--unwarp-dir',
            choices=['x', 'y', 'z', 'x-', 'y-', 'z-'],
            metavar='DIR',
            help='Unwarp direction (%(choices)s)')

    g4 = parser.add_argument_group('options')

    g4.add_argument('--base-volume',
            type=int,
            default=0,
            metavar='VOL',
            help='Functional base volume (0-indexed); default: %(default)d')

    g4.add_argument('--base-dset',
            metavar='BASE',
            help='EPI base volume dataset (overrides --base-volume)')

    g4.add_argument('--search',
            type=int,
            default=90,
            choices=[0, 90, 180],
            metavar='DEG',
            help='Search angle in degrees (%(choices)s; default: %(default)d)')

    g4.add_argument('--overwrite',
            action='store_true',
            help='Replace output files if they already exist')

    g4.add_argument('--keep-all',
            action='store_true',
            help='Keep intermediate files')

    g4.add_argument('--version',
            action='version',
            version=__version__,
            help='Show version number and exit')

    g4.add_argument('--help',
            action='help',
            help='Show this help message and exit')

    opts = parser.parse_args(argv)

    print('epi2t1 begins....')

    check_prefix(opts.prefix)

    # Look for input datasets and set output dataset names
    t1_brain = strip_ext(opts.t1_brain)
    if not dset_exists(t1_brain):
        raise IOError('Could not find --t1-brain: %s' % t1_brain)
    print('--t1-brain: %s' % t1_brain)

    t1_head = strip_ext(opts.t1_head)
    if not dset_exists(t1_head):
        raise IOError('Could not find --t1-head: %s' % t1_head)
    print('--t1-head: %s' % t1_head)

    epi = strip_ext(opts.epi)
    if not dset_exists(epi):
        raise IOError('Could not find --epi: %s' % epi)
    print('--epi: %s' % epi)

    out = []
    create_base = True
    base_needs_float = False
    if opts.base_dset:
        pre_base = strip_ext(opts.base_dset)
        if not dset_exists(pre_base):
            raise IOError('Could not find --base-dset: %s' % pre_base)
        create_base = False
        print('--base-dset: %s' % pre_base)
        capture = cmd('fslval', pre_base, 'datatype', capture_output=True)
        datatype = int(capture.stdout.decode())
        if datatype == 16:
            # If floating point, use the provided --base-dset as epi_base
            epi_base = pre_base
        else:
            # If not floating point, need to create epi_base
            epi_base = pre_base + '_float'
            base_needs_float = True
            out += [epi_base]
    else:
        # If --base-dset is not set, need to create epi_base
        epi_base = opts.prefix + '_base'
        out += [epi_base]

    epi_base_to_t1_dset = epi_base + '_to_t1'
    epi_base_to_t1_warp = epi_base_to_t1_dset + '_warp'
    epi_base_to_t1_mat = epi_base_to_t1_dset + '.mat'
    epi_mc = opts.prefix + '_mc'
    mc_cat = epi_mc + '.cat'
    mc_par = epi_mc + '.par'
    mc_mat = epi_mc + '.mat'
    mc_abs = epi_mc + '_abs.rms'
    mc_abs_mean = epi_mc + '_abs_mean.rms'
    mc_rel = epi_mc + '_rel.rms'
    mc_rel_mean = epi_mc + '_rel_mean.rms'
    epi_mc_mask = epi_mc + '_mask'
    epi_mc_res_mse = epi_mc + '_res_mse'
    epi_mc_res_mse0 = epi_mc_res_mse + '0'
    epi_mc_res_mse1 = epi_mc_res_mse + '1'
    epi_mc_res_mse_diff = epi_mc_res_mse + '_diff'
    epi_mc_outliers = epi_mc + '_outliers'
    epi_mc_one = epi_mc + '_one'
    epi_mc_zero = epi_mc + '_zero'
    epi_mc_stp = epi_mc + '_stp'
    epi_mc_singleev = epi_mc + '_singleev.txt'
    epi_mc_censor = epi_mc + '_censor.txt'
    epi_mc_confounds = epi_mc + '_confounds.txt'

    out += [
            epi_base_to_t1_dset,
            epi_base_to_t1_warp,
            epi_base_to_t1_mat,
            mc_cat,
            mc_par,
            mc_abs,
            mc_abs_mean,
            mc_rel,
            mc_rel_mean,
            epi_mc_censor,
            epi_mc_confounds
    ]

    tmp = [
            epi_mc,
            mc_mat,
            epi_mc_mask,
            epi_mc_res_mse,
            epi_mc_res_mse0,
            epi_mc_res_mse1,
            epi_mc_res_mse_diff,
            epi_mc_outliers,
            epi_mc_one,
            epi_mc_zero,
            epi_mc_stp,
            epi_mc_singleev
    ]

    # BBR?
    epi_base_to_t1_init_mat = None
    schedule = None
    bbr = True if opts.t1_brain_wmseg else False
    if bbr:
        wmseg = strip_ext(opts.t1_brain_wmseg)
        if not dset_exists(wmseg):
            raise IOError('Could not find --t1-brain-wmseg: %s' % wmseg)
        print('--t1-brain-wmseg: %s' % wmseg)

        epi_base_to_t1_init_mat = epi_base + '_to_t1_init.mat'
        tmp += [epi_base_to_t1_init_mat]

        # Get FLIRT schedule for BBR
        schedule = os.path.abspath(os.path.join(os.environ['FSLDIR'],
            'etc', 'flirtsch', 'bbr.sch'))
        if not os.path.isfile(schedule):
            raise IOError('Could not find FLIRT BBR schedule: %s' % schedule)

    # Distortion correction with field map?
    epi_base_to_t1_shift = None
    epi_base_to_t1_inv_mat = None
    fm_rads_brain_to_base_dset = None
    fm_rads_brain_to_base_mat = None
    fm_rads_brain_to_base_mask = None
    pe_dir = None
    if (opts.fm_rads_unmasked or
            opts.fm_rads_brain_t1 or
            opts.fm_to_t1_mat or
            opts.echo_spacing or
            opts.unwarp_dir):
        fm = True
    else:
        fm = False

    if fm:
        if not opts.fm_rads_unmasked:
            parser.error('Distortion correction requires --fm-rads-unmasked')
        if not opts.fm_rads_brain_t1:
            parser.error('Distortion correction requires --fm-rads-brain-t1')
        elif not opts.fm_to_t1_mat:
            parser.error('Distortion correction requires --fm-to-t1-mat')
        elif not opts.echo_spacing:
            parser.error('Distortion correction requires --echo-spacing')
        elif not opts.unwarp_dir:
            parser.error('Distortion correction requires --unwarp-dir')
        elif not bbr:
            parser.error('Distortion correction requires --t1-brain-wmseg')

        fm_rads_unmasked = strip_ext(opts.fm_rads_unmasked)
        if not dset_exists(fm_rads_unmasked):
            raise IOError('Could not find --fm-rads-unmasked: %s' %
                    fm_rads_unmasked)

        fm_rads_brain_t1 = strip_ext(opts.fm_rads_brain_t1)
        if not dset_exists(fm_rads_brain_t1):
            raise IOError('Could not find --fm-rads-brain-t1: %s' %
                    fm_rads_brain_t1)
        print('--fm-rads-brain-t1: %s' % fm_rads_brain_t1)

        fm_to_t1_mat = opts.fm_to_t1_mat
        if not os.path.isfile(fm_to_t1_mat):
            raise IOError('Could not find --fm-to-t1-mat: %s' % fm_to_t1_mat)
        print('--fm-to-t1-mat: %s' % fm_to_t1_mat)

        epi_base_to_t1_shift = epi_base_to_t1_dset + '_shift'
        epi_base_to_t1_inv_mat = epi_base_to_t1_dset + '_inv.mat'
        fm_rads_brain_to_base_dset = opts.prefix + '_fm_rads_brain_to_base'
        fm_rads_brain_to_base_mat = fm_rads_brain_to_base_dset + '.mat'
        fm_rads_brain_to_base_mask = fm_rads_brain_to_base_dset + '_mask'
        out += [epi_base_to_t1_shift]
        tmp += [
                epi_base_to_t1_inv_mat,
                fm_rads_brain_to_base_dset,
                fm_rads_brain_to_base_mat,
                fm_rads_brain_to_base_mask
        ]

        # Map unwarp direction to flirt phase encode direction
        pe_dir = {'x': 1, 'y': 2, 'z': 3, 'x-': -1, 'y-': -2, 'z-': -3}

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

    if opts.keep_all:
        for t in tmp:
            print('Temporary: %s' % t)

    # Delete temporary datasets and files
    for t in tmp:
        delete(t)

    try:
        if create_base:
            print('Extracting EPI base volume....')

            cmd('fslroi',
                    epi,
                    epi_base,
                    str(opts.base_volume),
                    str(1))

            cmd('fslmaths',
                    epi_base,
                    epi_base,
                    '-odt',
                    'float')

        else:
            print('Checking EPI base volume....')

            capture = cmd('fslnvols', pre_base, capture_output=True)
            bvols = int(capture.stdout.decode())

            if bvols != 1:
                raise RuntimeError('EPI base dataset must have one volume')

            if base_needs_float:
                cmd('fslmaths', pre_base, epi_base, '-odt', 'float')

        # Initial 6 DOF registration with correlation ratio cost function
        print('Registering EPI base volume to T1-weighted dataset....')

        cmd('flirt',
                '-in', epi_base,
                '-ref', t1_brain,
                '-omat', epi_base_to_t1_init_mat
                        if epi_base_to_t1_init_mat else epi_base_to_t1_mat,
                '-cost', 'corratio',
                '-dof', str(6),
                '-searchrx', str(-opts.search), str(opts.search),
                '-searchry', str(-opts.search), str(opts.search),
                '-searchrz', str(-opts.search), str(opts.search))

        if bbr:
            flirt = [
                    'flirt',
                    '-in', epi_base,
                    '-ref', t1_head,
                    '-omat', epi_base_to_t1_mat,
                    '-cost', 'bbr',
                    '-dof', str(6),
                    '-wmseg', wmseg,
                    '-init', epi_base_to_t1_init_mat,
                    '-nosearch',
                    '-schedule', schedule
            ]

            if fm:
                flirt += [
                        '-echospacing', str(opts.echo_spacing),
                        '-pedir', str(pe_dir[opts.unwarp_dir]),
                        '-fieldmap', fm_rads_brain_t1
                ]

            cmd(*flirt)

        # Generate warp fields for use with other registrations
        # Without distortion correction, a warp field isn't strictly necessary,
        # but we generate it anyway to allow using applywarp to go to standard
        # space, which is easier when incorporating motion parameters.

        print('Generating EPI base to T1 warp....')

        convertwarp = [
                'convertwarp',
                '-r', t1_head,
                '--postmat=%s' % epi_base_to_t1_mat,
                '-o', epi_base_to_t1_warp,
                '--relout'
        ]

        if fm:
            cmd('convert_xfm',
                    '-omat', epi_base_to_t1_inv_mat,
                    '-inverse', epi_base_to_t1_mat)

            cmd('convert_xfm',
                    '-omat', fm_rads_brain_to_base_mat,
                    '-concat',
                    epi_base_to_t1_inv_mat,
                    fm_to_t1_mat)

            cmd('applywarp',
                    '-i', fm_rads_unmasked,
                    '-r', epi_base,
                    '--premat=%s' % fm_rads_brain_to_base_mat,
                    '-o', fm_rads_brain_to_base_dset)

            cmd('fslmaths',
                    fm_rads_brain_to_base_dset,
                    '-abs',
                    '-bin',
                    fm_rads_brain_to_base_mask)

            cmd('fugue',
                    '--loadfmap=%s' % fm_rads_brain_to_base_dset,
                    '--mask=%s' % fm_rads_brain_to_base_mask,
                    '--saveshift=%s' % epi_base_to_t1_shift,
                    '--unmaskshift',
                    '--dwell=%s' % str(opts.echo_spacing),
                    '--unwarpdir=%s' % opts.unwarp_dir)

            convertwarp += [
                    '-s', epi_base_to_t1_shift,
                    '--shiftdir=%s' % opts.unwarp_dir
            ]

        cmd(*convertwarp)

        # Apply the warp field: This does the final EPI base to T1 registration
        print('Warping EPI base volume to T1-weighted dataset....')

        cmd('applywarp',
                '-i', epi_base,
                '-r', t1_head,
                '-o', epi_base_to_t1_dset,
                '-w', epi_base_to_t1_warp,
                '--interp=spline',
                '--rel')

        # Motion correction
        print('Motion correcting EPI dataset....')

        mcflirt = [
                'mcflirt',
                '-in', epi,
                '-out', epi_mc,
                '-mats',
                '-plots',
                '-rmsrel',
                '-rmsabs'
        ]

        if create_base:
            mcflirt += ['-refvol', str(opts.base_volume)]
        else:
            mcflirt += ['-reffile', epi_base]

        cmd(*mcflirt)

        # Concatenate the transformation matrices into one big file
        with open(mc_cat, 'wb') as catfile:
            for i in sorted(glob.glob('%s/*' % mc_mat)):
                with open(i, 'rb') as matfile:
                    shutil.copyfileobj(matfile, catfile)

        # The following detects motion outliers and produces a censor list that
        # can be passed to AFNI's 3dDeconvolve via the CENSORTR option and a
        # design matrix that can be used to add motion confound regressors to
        # an AFNI or FSL GLM. The procedure is that used in fsl_motion_outliers
        # with the rmsrel criterion for determining outliers, with the
        # exception that the base volume is selected via opts.base_volume
        # rather than automatically using the middle volume.

        print('Detecting motion outliers....')

        capture = cmd('fslnvols', epi_mc, capture_output=True)
        tmax = int(capture.stdout.decode())
        tmax1 = tmax-1

        capture = cmd('fslstats',
                epi_mc,
                '-P', str(2),
                '-P', str(98),
                capture_output=True)
        thr2, thr98 = capture.stdout.decode().split()
        thr2 = float(thr2)
        thr98 = float(thr98)
        robthr = thr2 + 0.1 * (thr98 - thr2)

        cmd('fslmaths',
                epi_mc,
                '-Tmean',
                '-thr', str(robthr),
                '-bin',
                epi_mc_mask)

        capture = cmd('fslstats',
                epi_mc,
                '-k', epi_mc_mask,
                '-P', str(50),
                capture_output=True)
        brainmed = float(capture.stdout.decode())

        capture = cmd('fslstats', epi_mc_mask, '-m', capture_output=True)
        maskmean = float(capture.stdout.decode())

        # Compute residual mean square error (--rmsrel in fsl_motion_outliers)
        cmd('fslmaths',
                epi_mc,
                '-sub', epi_base,
                '-mas', epi_mc_mask,
                '-div', str(brainmed),
                '-sqr',
                '-Xmean',
                '-Ymean',
                '-Zmean',
                '-div', str(maskmean),
                '-sqrt',
                epi_mc_res_mse,
                '-odt', 'float')

        # Compute difference
        cmd('fslroi',
                epi_mc_res_mse,
                epi_mc_res_mse0,
                str(0), str(1),
                str(0), str(1),
                str(0), str(1),
                str(0), str(tmax1))

        cmd('fslroi',
                epi_mc_res_mse,
                epi_mc_res_mse1,
                str(0), str(1),
                str(0), str(1),
                str(0), str(1),
                str(1), str(tmax1))

        cmd('fslmaths',
                epi_mc_res_mse1,
                '-sub', epi_mc_res_mse0,
                '-abs',
                epi_mc_res_mse_diff)

        # Compute outlier threshold
        capture = cmd('fslstats',
                epi_mc_res_mse_diff,
                '-p', str(25),
                '-p', str(75),
                capture_output=True)
        p25, p75 = capture.stdout.decode().split()
        p25 = float(p25)
        p75 = float(p75)
        threshv = p75 + 1.5 * (p75 - p25)

        # Find outliers
        cmd('fslmaths',
                epi_mc_res_mse_diff,
                '-thr', str(threshv),
                '-bin',
                epi_mc_outliers)

        cmd('fslroi',
                epi_mc_outliers,
                epi_mc_one,
                str(0), str(1),
                str(0), str(1),
                str(0), str(1),
                str(0), str(1))

        cmd('fslmaths',
                epi_mc_one,
                '-mul', str(0),
                epi_mc_zero)

        cmd('fslmerge',
                '-t',
                epi_mc_outliers,
                epi_mc_zero,
                epi_mc_outliers)

        capture = cmd('fslstats',
                epi_mc_outliers,
                '-V',
                capture_output=True)
        nmax = int(capture.stdout.decode().split()[0])

        # Get volume indices of outliers and create FSL design matrix
        outliers = []
        dm = []
        for i in range(tmax):
            cmd('fslmaths',
                    epi_mc_outliers,
                    '-roi',
                    str(0), str(1),
                    str(0), str(1),
                    str(0), str(1),
                    str(i), str(1),
                    epi_mc_stp)

            capture = cmd('fslstats', epi_mc_stp, '-V', capture_output=True)
            val = int(capture.stdout.decode().split()[0])
            if val > 0:
                outliers.append(i)
                cmd('fslmeants',
                        '-i', epi_mc_stp,
                        '-o', epi_mc_singleev)
                sev = None
                with open(epi_mc_singleev) as f:
                    sev = [line.strip() for line in f]
                dm.append(sev)

        print('Found %d outliers over %f' % (nmax, threshv))
        print(*outliers, sep=' ')

        # Write design matrix to file
        dm = zip(*dm)
        with open(epi_mc_confounds, 'w') as f:
            writer = csv.writer(f, delimiter=' ')
            writer.writerows(dm)

        # Write (0-indexed) list of censored TRs to text file
        with open(epi_mc_censor, 'w') as f:
            f.write(' '.join(str(i) for i in outliers))

    finally:
        if not opts.keep_all:
            for t in tmp:
                delete(t)

    print('epi2t1 ends....')


if __name__ == '__main__':
    epi2t1(sys.argv[1:])
