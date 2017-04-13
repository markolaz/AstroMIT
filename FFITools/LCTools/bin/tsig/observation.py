#
# Copyright (C) 2017 - Massachusetts Institute of Technology (MIT)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""
The observation is a combination of the catalog and the exposure to
be fed into a camera to produce the desired output location and effect.
"""

import glob
import os
import math
import logging
logger = logging.getLogger(__name__)

from astropy.io import fits
import numpy as np

from .util import to_bool, to_float
from .util.configurable import ConfigurableObject
from .spacecraft.transformer import Transformer
from .effects import *


# This array defines the order in which effects will be applied
EFFECTS = [
    FlatField,
    Smear,
    ShotNoise,
    CosmicRays,
    ReadoutNoise,
    DarkCurrent,
    Saturation,
    Undershoot,
    LineRinging,
    TransferEfficiency,
    BlackSky,
    BadPixels,
]


class Observation(ConfigurableObject):
    """
    The observation executes a series of exposures to generate a set of
    images.  In its simplest form, it does a set of exposures for each camera
    for a single pointing of the spacecraft with no spacecraft motion.

    1 exposure set at 2m cadence requires 60 2s exposures (2*60/2)
    1 exposure set at 30m cadence requires 900 2s exposures (30*60/2)
    """
    def __init__(self, num_exposures=0, cadence=2, camera_fov_radius=0.0,
                 apply_cosmic_mitigation=True, apply_lightcurves=True,
                 output_directory=None, save_combined_ccds=False,
                 retain_2s_images=False, apply_ccd_markup=False,
                 **kw):
        super(Observation, self).__init__()
        self.spacecraft = None
        self.mission = None
        self.catalog = None

        self.single_read = 2 # length of a single exposure, in seconds
        self.num_exposures = int(num_exposures)
        # FIXME: cadence should be a list, not a scaler and should contain either 2min ior 30min
        self.cadence = int(cadence) # how often to stack exposures, in minutes
        self.cadence_s = self.cadence * 60

        if to_bool(apply_cosmic_mitigation):
            cmcfg = kw.get('CosmicMitigation', {})
            self.stacker = CosmicMitigation(**cmcfg)
        else:
            self.stacker = ExposureStacker()

        self.lcg = None
        if to_bool(apply_lightcurves):
            lccfg = kw.get('LightCurves', {})
            self.lcg = LightCurveGenerator(**lccfg)

        # where to put everything
        self.output_directory = output_directory
        # optionally combine images from all 4 CCDs into a single image
        self.save_combined_ccds = to_bool(save_combined_ccds)
        # optionally save every two-second exposure
        self.retain_2s_images = to_bool(retain_2s_images)
        # optional override to camera field of view radius
        self.radius = float(camera_fov_radius)
        # place marks on the CCD image to verify position and reference frame
        self.apply_ccd_markup = to_bool(apply_ccd_markup)
        # see which effects should be applied
        # FIXME: Instan objects for effects here instead of later.
        self.effects = []
        ecfg = kw.get('Effects', {})
        for e in EFFECTS:
            name = e.__name__
            if name in ecfg and not to_bool(ecfg[name].pop('enable', True)):
                logger.debug("effect %s: disabled" % name)
            else:
                self.effects.append(e)
                logger.debug("effect %s: enabled" % name)
        self.effects_cfg = ecfg

    def set_spacecraft(self, spacecraft):
        self.spacecraft = spacecraft

    def set_catalog(self, catalog):
        self.catalog = catalog

    def set_mission(self, mission):
        self.mission = mission

    def observe(self):
        """
        Make an observation.  Loop through time to take multiple exposures,
        moving the spacecraft through a series of locations at each point in
        time.  Save an image for each CCD at each time step.  Periodically
        stack those exposures into the cadence (2m or 30m) images.
        """

        epoch = self.mission.epoch
        epoch_s = 0 # epoch in seconds
        end_t = epoch_s + self.num_exposures * self.cadence_s # end time

        # use the specified output directory...
        outdir = self.output_directory
        if outdir is None:
            # ... or default to something sane if nothing specified
            outdir = self.get_dirname(self.mission.pointing_ra,
                                      self.mission.pointing_dec, epoch)
        # 2s exposures go in a subdirectory within the output directory
        outdir_e = "%s/exposures" % outdir
        try:
            os.mkdir(outdir)
            os.mkdir(outdir_e)
        except os.error as e:
            raise ValueError("observation aborted: %s" % e)
        logger.info("output directory is %s" % outdir)
        logger.info("output directory for 2s exposures is %s" % outdir_e)

        # Point the spacecraft for the duration of the observation.
        self.spacecraft.set_pointing(
            self.mission.pointing_ra,
            self.mission.pointing_dec,
            self.mission.pointing_roll)

        # Provide some feedback about what will happen
        logger.info("num exposures: %s" % self.num_exposures)
        logger.info("single read: %ss" % self.single_read)
        logger.info("cadence: %ss" % self.cadence_s)
        logger.info("epoch: %s (start=%s end=%s)" % (epoch, epoch_s, end_t))
        logger.info("effects: %s" % ''.join([e.__name__ for e in self.effects]))
        logger.info("save combined: %s" % self.save_combined_ccds)
        logger.info("retain 2s images: %s" % self.retain_2s_images)
        logger.info("radius override: %s" % self.radius)
        logger.info("apply CCD markup: %s" % self.apply_ccd_markup)

        # Transformer knows the relationships between frames of reference
        tr = Transformer()

        # count of cadence images
        self.c_count = 0
        # count of two-second exposures
        self.e_count = 0
        first_exposure = self.e_count

        t = epoch_s
        while t <= end_t:
            logger.debug("time=%s (end_t=%s)" % (t, end_t))

            # Move spacecraft to the correct position and velocity
            (position, velocity) = self.mission.get_spacecraft_state(t)
            self.spacecraft.set_position(*position)
            self.spacecraft.set_velocity(*velocity)

            # do a 2-second exposure for each CCD in each camera
            for (cam_i, camera) in enumerate(self.spacecraft.camera):
                cam_id = cam_i + 1
                logger.debug("%s" % camera.label)
                radius = self.radius or camera.fov_radius
                # one subpixel (super-resolution) buffer for each CCD
                super_buffer = dict()
                for (ccd_i, ccd) in enumerate(camera.ccd):
                    ccd_id = ccd_i + 1
                    super_buffer[ccd_id] = SuperBuffer(ccd)
                # query for the celestial objects
                # FIXME: make the catalog return a queryresults object
                # FIXME: do the query+lightcurves only once for each camera
                ra_cam, dec_cam, roll_cam = camera.get_pointing()
                self.catalog.query(ra_cam, dec_cam, radius)
                if self.lcg:
                    self.lcg.apply_lightcurves(self.catalog)
                ra, dec, tmag, teff = self.catalog.stars_snapshot(epoch=epoch)
                # FIXME: use vectors in celestial_to_pixel instead of looping
                for idx in range(len(tmag)):
                    # map the celestial objects onto each subpixel buffer
                    x, y, ccd_n, fa_x, fa_y = tr.celestial_to_pixel(
                        ra[idx], dec[idx],
                        self.mission.pointing_ra,
                        self.mission.pointing_dec,
                        self.mission.pointing_roll,
                        camera.get_geometry(),
                        camera.get_ccd_geometries(),
                        subpixels_per_pixel=SuperBuffer.RESOLUTION)
                    self.apply_electrons(super_buffer[ccd_n[0]],
                                         x[0], y[0], fa_x[0], fa_y[0],
                                         tmag[idx], teff[idx])
                for (ccd_i, ccd) in enumerate(camera.ccd):
                    ccd_id = ccd_i + 1
                    logger.debug("%s %s" % (camera.label, ccd.label))
                    jitter = None # FIXME: get jitter
                    pixels = self.downsample(super_buffer[ccd_id], jitter)
                    hdulist = self.make_hdulist(
                        pixels, self.spacecraft, camera, ccd)
                    hdulist = self.apply_effects(
                        hdulist, self.effects, self.effects_cfg)
                    # FIXME: This should be moved to an effect
                    if self.apply_ccd_markup:
                        self.apply_marks(hdulist[0].data.T, ccd_id)
                    # save each exposure to disk so it can be stacked later
                    self.save_image(hdulist, outdir_e, cam_id, ccd_id,
                                    self.e_count, self.single_read)

            if (t - epoch_s) % self.cadence_s == 0 and t != epoch_s:
                # it is time to stack the exposures into a cadence image
                logger.debug("stack images at t=%s" % t)
                for (cam_i, camera) in enumerate(self.spacecraft.camera):
                    cam_id = cam_i + 1
                    composite = None
                    if self.save_combined_ccds:
                        composite = dict()
                    for (ccd_i, ccd) in enumerate(camera.ccd):
                        ccd_id = ccd_i + 1
                        hdulist = self.make_hdulist(
                            None, self.spacecraft, camera, ccd)
                        hdulist = self.stacker.stack_exposures(
                            hdulist, outdir_e, cam_id, ccd_id,
                            first_exposure, self.e_count)
                        # FIXME: don't do apply ccd markup here
                        if self.apply_ccd_markup:
                            self.apply_marks(hdulist[0].data.T, ccd_id)
                        self.save_image(hdulist, outdir, cam_id, ccd_id,
                                        self.e_count, self.cadence_s)
                        if composite is not None:
                            composite[ccd_id] = hdulist[0].data.T
                    if composite is not None:
                        hdulist = self.create_composite(composite)
                        self.save_image(hdulist, outdir, cam_id, None,
                                        self.e_count, self.cadence_s)
                    # FIXME: get postage stamps from full image

                # clear the exposures to make way for the next image
                first_exposure = self.e_count + 1
                if not self.retain_2s_images:
                    self.delete_2s_exposures(outdir_e)
                self.c_count += 1

            t += self.single_read
            self.e_count += 1

    @staticmethod
    def apply_electrons(super_buffer, x, y, fa_x, fa_y, tmag, teff):
        """
        Put the electrons into each pixel.  Note that the pixels of the
        super-resolution buffer count electrons, not ADUs.

        super_buffer - subpixel buffer with same reference frame as the CCD
                       but higher resolution
        """
        # Provisional conversion of tess magnitude to number of electrons as
        # per Al Levine's notes from formula for zodiacal light.
        num_electrons = 2.4 * 10.0 ** (-0.4 * (tmag - 22.0))
        # FIXME: apply psf, 0.5 is a fudge factor
        i = int(x + 0.5 * super_buffer.num_buffer_pixels) # FIXME
        j = y + super_buffer.num_buffer_pixels # FIXME
        try:
            super_buffer.pixels[i, j] += num_electrons
        except IndexError:
            # ignore anything that is not within the buffer
            pass

    @staticmethod
    def downsample(super_buffer, jitter):
        # FIXME: apply jitter
        return super_buffer.get_ccd_pixels()

# FIXME: Make the observation another object that has get_headers
# and turn this make_hdulist into a static method OR move it to
# utils.
    def make_hdulist(self, image, *objs):
        """Produce an hdulist and add the headers from each object"""
        head = fits.Header()
        head['EXPTIME'] = (self.single_read, '[s] single exposure time')
        head['CADTIME'] = (self.cadence_s, '[s] stacked exposure time')
        head['NREADS'] = (
            np.round(self.cadence_s / self.single_read).astype(np.int),
            'number of stacked exposures')
        head['COUNTER'] = (
            self.e_count, 'number of exposures since start')
        head['CCOUNTER'] = (
            self.c_count, 'number of stacked exposures since start')

#        head['READTIME'] = (
#            readout_time, '[s] time to transfer to frame store')
#        head['BJD0'] = (bjd0, '[day] base time subtracted from all BJD')
#        head['BJD'] = (bjd - bjd0, '[day] mid-exposure time - BJD0')
#        head['BJD_TDB'] = (bjd, '[day] BJD_TDB')
#        head['ANTISUN'] = (bjd_antisun - bjd0, '[day] time of antisun - BJD0')
#        head['EPOCH'] = (epoch, '[year] epoch of mid-exposure time')

# from the header ICD draft by chelsea and andras
# JD, JDTAI, JDTB - julian day of midexpo in some uniform time
# EXPTIME - gross exposure time
# BJD DX, BJD DY, BJD DZ - coordinates of spacecraft wrt solar system
#   barycenter.  if units are speed of light times SI day, then BJD-JD
#   correction can be retrieved by computing the scalar product of this and the
#   normal vector of the star.  J2000 coordinates are preferred.
# TIMESLICA, TIMELICB, TIMESLICC, TIMESLICD - in kepler it was TIMSLICE
# TIERRELA - uncertainty and/or precision about onboard timekeeping

# TESS QA, TESS QB, TESS QC, TESS QD - spacecraft quaternion
# RA NOM, DEC NOM, ROLL NOM, RAC(CRVAL1), DECC (CRVAL2) - nominal approximate
#   J2000 positions of the optical axes (including field rotation) and the ccd
# centers
# RA (CRVAL1) - ra for ccd center
# DECC (CRVAL2) - dec for ccd center

# pointing jitter reported by the ADCS
# position of moon and earth - moon vector wrt spacecraft?
#   earth, moon, and sun vectors?

        for obj in objs:
            if hasattr(obj, "get_headers"):
                for (name, value) in obj.get_headers():
                    head[name] = value
        if image is None:
            width = head['CCDCOLS']
            height = head['CCDROWS']
            image = np.zeros(shape=(width, height), dtype=np.int32)
        hdu = fits.PrimaryHDU(image.T, header=head)
        hdulist = fits.HDUList([hdu])
        return hdulist

    @staticmethod
    def create_composite(image_dict):
        """
        Combine 4 CCD images into a single image.  Flip CCD 1 and 2.
        """
        # assume that all 4 CCD images have the same dimensions
        # FIXME: insert the non-science pixels
        (width, height) = image_dict[1].shape
        image = np.zeros(shape=(2 * width, 2 * height), dtype=np.int32)
        image[width:, height:] = image_dict[1][::-1, ::-1]
        image[0:width, height:] = image_dict[2][::-1, ::-1]
        image[0:width, 0:height] = image_dict[3]
        image[width:, 0:height] = image_dict[4]
        head = fits.Header()
        # FIXME: insert the headers
        # FIXME TOO: Change CCD to CCDn for each CCD HEADER
        # ERROR if the header name is more than 7 letters
        hdu = fits.PrimaryHDU(image.T, header=head)
        hdulist = fits.HDUList([hdu])
        return hdulist

    @staticmethod
    def apply_effects(hdulist, effects, effects_cfg):
        """Take each effect in order and apply it to the FITS object"""
        for effect in effects:
            try:
                ename = effect.__name__
                cfg = effects_cfg.get(ename, {})
                logger.debug("apply %s" % ename)
                hdulist = effect(**cfg).apply(hdulist)
            except (NotImplementedError, ImportError) as err:
                logger.error(str(err))
        return hdulist

    @staticmethod
    def apply_marks(pixels, ccd_id):
        """Apply a grid in corner of CCD to mark origin and ID"""
        # draw the origin block
        pixels[10:200, 10:200] = 1000
        # draw smaller blocks in origin to indicate the CCD number
        if ccd_id >= 1:
            pixels[120:180, 120:180] = 0
        if ccd_id >= 2:
            pixels[30:90, 120:180] = 0
        if ccd_id >= 3:
            pixels[30:90, 30:90] = 0
        if ccd_id >= 4:
            pixels[120:180, 30:90] = 0
        # draw a line in the direction of increasing x
        pixels[250:500, 20:30] = 500
        # draw a box around the CCD
        pixels[0:, 1:2] = 300
        pixels[1:2, 0:] = 300
        pixels[-2:-1, 0:] = 300
        pixels[0:, -2:-1] = 300

    @staticmethod
    def get_dirname(ra, dec, epoch):
        """
        Default file path is the CWD plus:

        $CWD/%(ra)2dh%(dec)2dd%(epoch)4de/ 

        e.g., /tmp/01h45d2019e/
        """
        return "%02dh%02dd%04de" % (ra, dec, epoch)

    @staticmethod
    def get_filename(outdir, cam_id, ccd_id, frame, cadence, imgtype='full'):
        """
        The filenames for each time should be:
        
        cam%(camera)d/
          ccd%(ccd)d/
            tsig-%(frame)7d-cam%1d-ccd%1d-%(cadence)-%(full).fits

        cadence: in seconds
                 maps to: llc (30mins), slc (2mins), mlc (2sec)
        full: imagette or full
              ffi (full image), tps (postage stamp)

        e.g., cam1/ccd1/tsig-0000001-cam1-ccd1-llc-ffi.fits
        """
        if imgtype == 'imagette':
            imgtype_label = 'tps'
        elif imgtype == 'full':
            imgtype_label = 'ffi'
        else:
            raise TypeError("Unknown image type '%s'" % imgtype)

        # FIXME: if imagette, include identifier in filename

        parts = [outdir]
        if ccd_id is not None:
            parts.append("tsig-%07d-cam%1d-ccd%1d-%d-%s.fits" % (
                frame, cam_id, ccd_id, cadence, imgtype_label))
        else:
            parts.append("tsig-%07d-cam%1d-%s-%s.fits" % (
                frame, cam_id, cadence, imgtype_label))
        return '/'.join(parts)

    @staticmethod
    def save_image(hdulist, outdir, cam_id, ccd_id, frame, cadence):
        filename = Observation.get_filename(
            outdir, cam_id, ccd_id, frame, cadence)
        logger.debug("save %s" % filename)
        hdulist.writeto(filename)

    @staticmethod
    def delete_2s_exposures(outdir):
        pattern = "%s/*-2-ffi.fits"
        for filename in glob.glob(pattern):
            try:
                logger.debug("delete 2s image %s" % filename)
                os.remove(filename)
            except IOError as e:
                logger.error("delete failed: %s" % e)


class MissionProfile(ConfigurableObject):
    """
    Contains all the spacecraft positioning for each of the times required
    for the observation.
    """
    def __init__(self, filename=None,
                 pointing=None, position=None, velocity=None, epoch=2018):
        """
        If a filename is specified, get the position and velocity for each
        point in time from the file.  Otherwise, use a fixed position and
        velocity for the entire mission.

        The pointing is constant for the duration of the mission.
        """
        if pointing:
            self.pointing_ra = float(pointing[0])
            self.pointing_dec = float(pointing[1])
            self.pointing_roll = float(pointing[2])
        elif not filename:
            self.pointing_ra = 30.0
            self.pointing_dec = 45.0
            self.pointing_roll = 0.0

        self.position = position or (0, 0, 0)
        self.velocity = velocity or (0, 0, 0)
        self.epoch    = epoch or 2018

        self.timed_mission = None

        if filename:
            if not os.path.isfile(filename):
                raise IOError("Cannot read mission file: %s" % str(filename))
            # XXX Load timed_mission list from filename
            # XXX load epoch and pointing here
            #self.timed_mission = ...

    def get_spacecraft_state(self, t):
        # FIXME: interpolate position and velocity as needed
        return self.position, self.velocity


class ExposureStacker(object):
    def __init__(self):
        logger.debug("generic stacker")

    def stack_exposures(self, hdulist,
                        imgdir, cam_id, ccd_id, start_frame, end_frame):
        """
        Put data from multiple two-second exposures into a single image.
        Two-second exposures are in files on disk, so we have to read each
        file and add it to the final image.
        """
        for i in range(start_frame, end_frame + 1):
            filename = Observation.get_filename(imgdir, cam_id, ccd_id, i, 2)
            logger.debug("read exposure %s" % filename)
            try:
                new_hdulist = fits.open(filename)
                assert(hdulist[0].shape == new_hdulist[0].shape)
                hdulist[0].data += new_hdulist[0].data
            except IOError as e:
                logger.error("read failed: %s" % e)
        return hdulist


class CosmicMitigation(ConfigurableObject):
    """
    Cosmic mitigation is an algorithm to reduce the effect of cosmic rays.
    The algorithm is applied while stacking two-second images.
    """
    def __init__(self, block_size=10):
        super(CosmicMitigation, self).__init__()
        """
        count - number of 2s images to consider as a block
        """
        logger.debug("cosmic mitigation stacker: block_size=%s" % block_size)
        self.block_size = int(block_size)

    def stack_exposures(self, hdulist,
                        imgdir, cam_id, ccd_id, start_frame, end_frame):
        # FIXME: this is not a precise implementation of the algorithm, since
        # the block_size might not be a multiple of the number of frames and
        # since it considers only one highest and one lowest for each pixel
        minval = np.zeros_like(hdulist[0].data) # FIXME: set to high value
        maxval = np.zeros_like(hdulist[0].data)
        idx = 0
        for i in range(start_frame, end_frame + 1):
            filename = Observation.get_filename(imgdir, cam_id, ccd_id, i, 2)
            logger.debug("read exposure %s" % filename)
            try:
                new_hdulist = fits.open(filename)
                assert(hdulist[0].shape == new_hdulist[0].shape)
                minval = np.minimum(minval, new_hdulist[0].data)
                maxval = np.maximum(maxval, new_hdulist[0].data)
                idx += 1
                hdulist[0].data += new_hdulist[0].data
                if idx % self.block_size == 0:
                    hdulist[0].data -= minval
                    hdulist[0].data -= maxval
                    minval = np.zeros_like(hdulist[0].data)
                    maxval = np.zeros_like(hdulist[0].data)
            except IOError as e:
                logger.error("read failed: %s" % e)
        hdulist[0].data -= minval
        hdulist[0].data -= maxval
        return hdulist


class LightCurveGenerator(ConfigurableObject):
    """
    The light curve generator controls how light curves are assigned to objects.
    """
    def __init__(self, types=None, fraction=0.0, min_brightness=None):
        super(LightCurveGenerator, self).__init__()
        """
        types - which types of lightcurves should be considered
        fraction - how many of the objects should get a light curve [0.0, 1.0]
        min_brightness - apply lightcurves only to this magnitude or brighter
        """
        logger.debug("light curve generator: fraction=%s" % fraction)
        if types is None:
            types = ['ConstantLightCurve']
        self.types = types
        self.fraction = float(fraction)
        self.min_brightness = to_float(min_brightness)

    def apply_lightcurves(self, catalog):
        catalog.add_lightcurves(faintest_star=self.min_brightness,
                                stars_with_random_curves=self.fraction)


class SuperBuffer(object):
    """
    Pixel buffer with sub-pixel resolution relative to a CCD.
    """
    RESOLUTION = 1.0 # Maybe 11, 12, or 15 or 22

    def __init__(self, ccd, num_buffer_pixels=10):
        # num_buffer_pixels - number of pixels beyond the CCD pixels
        self.num_buffer_pixels = num_buffer_pixels
        self.width = int((ccd.cols + 2 * num_buffer_pixels) * self.RESOLUTION)
        self.height = int((ccd.rows + 2 * num_buffer_pixels) * self.RESOLUTION)
        self.pixels = np.zeros(shape=(self.width, self.height), dtype=np.int32)

    def get_ccd_pixels(self):
        # transform from superbuffer to CCD pixels
        # FIXME: downsample based on resolution
        return self.pixels[self.num_buffer_pixels: -self.num_buffer_pixels,
                           self.num_buffer_pixels: -self.num_buffer_pixels]
