#!/usr/bin/env python2

################################################################################
## Common functions for bad
## Copyright (C) 2014  Rachel Domagalski: domagalski@berkeley.edu
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## ## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

import os as _os
import sys as _sys
import numpy as _np
import os.path as _op
import aipy.cal as _ac
import numpy.random as _npr
import omnical.calib as _oc

def get_aa(calfile):
    """
    Simple wrapper around aipy's get_aa function.
    """
    if calfile[-3:] == '.py':
        calfile = calfile[:-3]
    try:
        return _ac.get_aa(calfile, 0.1, 0.1, 1)
    except SyntaxError:
        _sys.path.insert(0, _op.dirname(_op.abspath(calfile)))
        return _ac.get_aa(_op.basename(calfile), 0.1, 0.1, 1)

def gen_corr(gains, reds, visdata=None, noise=None, bad_ants=[]):
    """
    Generate a set of random visibilities and use them with the antenna
    gains to return a set of fake measured correlations. Alternatively,
    visibility data can be provided if needed.
    """
    # Generate random visibility data.
    nreds = len(reds)
    gen_vis = False
    if visdata is None:
        gen_vis = True
        visdata = _npr.randn(nreds) + 1j * _npr.randn(nreds)
        visdata = visdata.astype(_np.complex64)

    # Generate random gaussian noise for each antenna
    gen_noise = False
    if noise is None:
        noise = _np.zeros(len(gains)) # Default to no noise.
    elif not hasattr(noise, '__len__'):
        gen_noise = True
        nant = len(gains)
        noise = (_npr.randn(nant) + 1j * _npr.randn(nant)) * noise / _np.sqrt(2)

    # Generate correlator data
    corrdata = {}
    for i, redlist in enumerate(reds):
        for a1, a2 in redlist:
            corrdata[(a1, a2)] = _np.conj(gains[a1])*gains[a2]*visdata[i]
            corrdata[(a1, a2)] += _np.conj(noise[a1])*noise[a2]
            if a1 in bad_ants:
                corrdata[(a1, a2)] += _npr.randn() + 1j * _npr.randn()
            if a2 in bad_ants:
                corrdata[(a1, a2)] += _npr.randn() + 1j * _npr.randn()
            corrdata[(a1,a2)] = _np.array([[corrdata[(a1,a2)]]], _np.complex64)

    # Return the data.
    output = [corrdata]
    if gen_vis:
        output.insert(0, visdata)
    if gen_noise:
        output.append(noise)
    if len(output) == 1:
        return corrdata
    else:
        return output

def mkgrid(gridding):
    """
    Ncols,Nrows,dx,dy
    """
    # Get the information to generate the grid from
    # XXX make it hexagonal
    ncols, nrows, delta_x, delta_y = gridding.split(',')
    ncols = int(ncols)
    nrows = int(nrows)
    delta_x = float(delta_x)
    delta_y = float(delta_y)
    nants = nrows * ncols

    # simple rectangular grid
    posfunc = lambda i: [delta_x*(i%nrows), delta_y*(i/nrows), 0.]
    antpos = _np.array([posfunc(i) for i in range(nants)])
    return antpos

def print_progress(step,
                   total,
                   prog_str='Percent complete:',
                   quiet=False,
                   progfile='progress'):
    """
    Print the progress of some iteration through data. The step is the
    current i for i in range(total). This function can also display
    progress quietly by writing it to a file.

    Input:

    - ``step``: The iteration, starting at 0, of progress.
    - ``total``: Total number of iterations completed.
    - ``prog_str``: Message to print with the progress number.
    - ``quiet``: Setting this to true saves the progress to a file.
    - ``progfile``: The file to save progress to when in quiet mode.
    """
    progress = round(100 * float(step+1) / total, 2)
    progress = '\r' + prog_str + ' ' + str(progress) + '%\t\t'
    if quiet:
        if step + 1 == total:
            _os.system('rm -f ' + progfile)
        else:
            with open(progfile, 'w') as f:
                f.write(progress[1:-2] + '\n')
    else:
        print progress,
        if step == total - 1:
            print
        else:
            _sys.stdout.flush()
