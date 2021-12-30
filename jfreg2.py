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

import sys
import os
import re
import shutil
import argparse
import subprocess


__version__ = '2.0.0 (pre-release)'


def cmd(*cmd, echo=True, capture_output=False):
    '''Run a shell command (echoing it to stdout by default)'''
    if(echo):
        print(' '.join(cmd))
    return subprocess.run(cmd, check=True, capture_output=capture_output)


def dset_exists(path):
    '''Test whether a dataset exists'''
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
    '''Strip NIfTI extensions'''
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
    '''Check FLIRT version'''
    FLIRT_MIN = (6, 0)
    capture = cmd('flirt', '-version', echo=False, capture_output=True)
    flirt_verstr = capture.stdout.decode().splitlines()[0]
    pattern = r'^FLIRT\s+version\s+(\d+)\.(\d+).*$'
    match = re.search(pattern, flirt_verstr)
    assert match is not None
    flirt_ver = tuple(map(int, match.groups()))
    if flirt_ver < FLIRT_MIN:
        raise RuntimeError('FLIRT %d.%d or newer is required' % FLIRT_MIN)
