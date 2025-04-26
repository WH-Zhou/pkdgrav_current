#!/usr/bin/env python3
'''
Automate movie generation from pkdgrav output.

Invoke with '-h' or '--help' for a usage message.
Assumes ffmpeg, povray, rastoppm, and ssdraw are in the search path.
Written by Derek C. Richardson 7/20/18.
'''

import argparse
import glob
import os
import subprocess
import sys


def get_params(filename):
    '''Return lines of text read from parameter file.'''
    try:
        return open(filename, 'r').readlines()
    except IOError:
        return None


def get_param(params, target, pos):
    '''Return parameter value from given position on target line.'''
    matches = [x.split()[pos] for x in params if x.startswith(target)]
    if len(matches) != 1:
        print(f'Unable to find unique match for "{target}" in parameters.')
        exit(1)
    return matches[0]


def ffmpeg_clean(filetype):
    '''Remove softlinks used to create the movie.'''
    for file in glob.glob('ffmpeg*.' + filetype):
        os.remove(file)


def ffmpeg_links(frames, filetype):
    '''Create consecutively numbered softlinks to frames for ffmpeg.'''
    ffmpeg_clean(filetype)
    count = 0
    # ffmpeg requires numbers in filenames be consecutive integers...
    for frame in frames:
        os.symlink(f'{frame}.{filetype}', f'ffmpeg{count:012d}.{filetype}')
        count += 1


def status_start(args, message):
    '''Begin active status message.'''
    if args.quiet:
        return
    # make room for backspaces if this is an addressable terminal...
    if sys.stdout.isatty():
        print(message + '    ', end='')  # 4 spaces for '100%'
        status_incr.count = 0
    else:
        print(message)


def status_incr(args, max):
    '''Show current progress and increment counter.'''
    if args.quiet or not sys.stdout.isatty():
        return
    percent = 100 * status_incr.count // max
    print(f'\b\b\b\b{percent:3d}%', end='', flush=True)
    status_incr.count += 1


def status_end(args):
    '''Finish active status message, showing 100% completion.'''
    if args.quiet or not sys.stdout.isatty():
        return
    print('\b\b\b\b100%')


def frame_skip(args, frame, filetype):
    '''Return True if current frame already drawn and should be skipped.'''
    if filetype == 'pov':
        # check for gzip'd version...
        if frame_skip(args, frame, filetype + '.gz'):
            return True
    if os.path.isfile(frame + '.' + filetype):
        if args.keep:
            return True
        if args.force:
            return False
        print(f'\n{frame}.{filetype} exists.'
              ' Use --force or --keep to override.')
        exit(1)
    else:
        return False


def post_ssdraw(frame, do_povray, width, height, ifiletype):
    '''Run POV-Ray/rastoppm to convert ssdraw output frame to ffmpeg input.'''
    if do_povray:
        subprocess.run(f'gzip -d {frame}.pov.gz', shell=True)  # gunzip first
        result = subprocess.run(
                f'povray -i{frame}.pov -w{width} -h{height}'
                ' +A -J -D >& /dev/null', shell=True
                )
        if result.returncode != 0:
            print(f'\nError running POV-Ray (frame {frame}).')
            exit(1)
    else:
        assert width is None and height is None
        result = subprocess.run(f'rastoppm {frame}.ras > {frame}.ppm',
                                shell=True)
        if result.returncode != 0:
            print(f'\nError running rastoppm (frame {frame}).')
            exit(1)
    os.remove(frame + '.' + ifiletype)


def main():

    # parse arguments...

    parser = argparse.ArgumentParser(
            description='Make a movie from pkdgrav files.'
            )
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-f', '--force', action='store_true',
                       help='force overwriting render frames')
    group.add_argument('-k', '--keep', action='store_true',
                       help='preserve render frames')
    parser.add_argument('-b', '--batch', action='store_true',
                        help='process frames with ssdraw all at once '
                        'instead of one at a time (needed for certain '
                        'camera tracking options)')
    parser.add_argument('-c', '--color24', action='store_true',
                        help='use 24-bit color (POV-Ray only)')
    parser.add_argument('-o', '--orient', action='store_true',
                        help='use particle orientations (POV-Ray only)')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='disable informational messages')
    parser.add_argument('-n', '--nth', default=1, type=int,
                        help='render every n-th frame (default: %(default)d)')
    parser.add_argument('-r', '--rate', default='25',
                        help='movie frame rate in frames per second '
                        '(default: %(default)s)')
    args = parser.parse_args()

    # read parameter files...

    pkdparams = get_params('ss.par')
    if not pkdparams:
        if not args.quiet:
            print('Unable to read from "ss.par"...assuming "ss" for basename.')
        basename = 'ss'
    else:
        basename = get_param(pkdparams, 'achOutName', 2)

    ssdparams = get_params('ssdraw.par')
    if not ssdparams:
        print('Unable to read from "ssdraw.par"...aborting.')
        exit(1)

    # full vs. reduced output determination below done solely on basis
    # of filenames, not by looking at magic numbers...

    outtemp = glob.glob(basename + '.[0-9]*[0-9]')
    outputs = glob.glob(basename + '.[0-9]*[0-9].r')
    compare = ' '.join(outputs)

    using_full_output = False
    for output in outtemp:
        if output not in compare:
            outputs.append(output)
            using_full_output = True

    if not args.quiet:
        if using_full_output:
            if len(compare) > 0:
                print('Using a mix of full and reduced outputs.')
            else:
                print('Using full outputs.')
        else:
            if len(outtemp) > 0:
                print('Using reduced outputs in place of full outputs.')
            elif len(outputs) > 0:
                print('Using reduced outputs.')
            else:
                print('No outputs found!')

    outputs.sort()
    outputs = outputs[::args.nth]

    do_povray = (get_param(ssdparams, 'Particle shape', 2) == '2')

    # set file types and get frame dimensions if needed...

    if do_povray:
        if not os.path.isfile('povray.inc'):
            print('Missing povray.inc file (get it from pkdgrav etc folder).')
            exit(1)
        ifiletype = 'pov'
        ofiletype = 'png'
        width = int(get_param(ssdparams, 'Frame size', 2))
        # note use of eval() is a possible security risk...see
        # https://stackoverflow.com/questions/2371436/evaluating-a-mathematical-expression-in-a-string for alternatives # noqa
        height = int(width /
                     eval(f"eval({get_param(ssdparams, 'Aspect ratio', 2)})"))
    else:
        ifiletype = 'ras'
        ofiletype = 'ppm'
        width = None  # not used
        height = None  # ditto

    # add optional ssdraw arguments...

    ssdraw_args = ''
    if args.color24:
        if not do_povray:
            print('24-bit color requires POV-Ray frames...aborting.')
            exit(1)
        ssdraw_args += ' -c'
    if args.orient:
        if not do_povray:
            print('Orientations option requires POV-Ray frames...aborting.')
            exit(1)
        ssdraw_args += ' -o'

    # process the frames, using batch mode if requested...

    if args.batch:
        print('Running ssdraw in batch mode...', end='', flush=True)
        # send all the needed frames to ssdraw at once...
        frames = [f for f in outputs if not frame_skip(args, f, ofiletype)]
        result = subprocess.run(f'ssdraw {ssdraw_args} {" ".join(frames)}',
                                capture_output=True, shell=True)
        if result.returncode != 0:
            raise Exception('Error running ssdraw. Output:\n',
                            result.stdout.decode('utf-8'))
        print('done!')
        # now process ssdraw output one frame at a time...
        status_start(args, 'Post-processing frames...')
        for frame in frames:
            status_incr(args, len(frames))
            post_ssdraw(frame, do_povray, width, height, ifiletype)
    else:
        status_start(args, 'Processing frames...')
        # do everything one frame at a time in this case...
        for frame in outputs:
            status_incr(args, len(outputs))
            if frame_skip(args, frame, ofiletype):
                continue
            result = subprocess.run('ssdraw' + ssdraw_args + ' ' + frame,
                                    capture_output=True, shell=True)
            if result.returncode != 0:
                raise Exception('Error running ssdraw. Output:\n',
                                result.stdout.decode('utf-8'))
            post_ssdraw(frame, do_povray, width, height, ifiletype)

    status_end(args)

    # make the movie!...

    if not args.quiet:
        print('Making movie...')

    ffmpeg_links(outputs, ofiletype)
    result = subprocess.run(f'ffmpeg -framerate {args.rate} '
                            f'-i ffmpeg%012d.{ofiletype} '
                            '-loglevel warning -pix_fmt yuv420p -y movie.mp4',
                            shell=True)
    if result.returncode != 0:
        print('WARNING: non-zero return code from ffmpeg.')

    ffmpeg_clean(ofiletype)

    # clean up if needed...

    if not args.keep:
        for frame in outputs:
            os.remove(frame + '.' + ofiletype)


if __name__ == '__main__':
    # subprocess.run(..., capture_output=True, ...) needs Python 3.7...
    if sys.version_info.major == 3 and sys.version_info.minor < 7:
        raise Exception('Python version 3.7 or later required.')
    # intercept keyboard interrupt to make friendlier message...
    try:
        main()
    except KeyboardInterrupt:
        print('Interrupted.')
