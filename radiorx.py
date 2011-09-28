#! /usr/bin/env python
# 
# Copyright (C) 2009  The University of Texas at Austin.
# 
# This file is part of Hydra: A wireless multihop testbed.
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

from gnuradio import gr, gru
from gnuradio import uhd
from gnuradio import eng_notation
from optparse import OptionParser

from hydra import rf as rflib
from hydra.PyHydra import Callable

import sys

# FIXME: Need to determine actual max speed for USB 2.0
MAX_USB_RATE = 128e6       #128MB/sec

class default_radiorx_setup:
    d_options = {'verbose':0, 'fake_rf':False,
                 'which_board':0, 'nchannels':1, 'subdev_spec':None,
                 'sample_rate':1.0e6, 'freq':2400.0e6,
                 'rx_gain': 40.0, 'master_serialno':None}

    def get_options():
        dopt = default_radiorx_setup()
        for s in default_radiorx_setup.d_options:
            if s not in dir(dopt):
                setattr(dopt, s, default_radiorx_setup.d_options[s] )
        return dopt
    
    def check_options(options):
        for s in default_radiorx_setup.d_options:
            if s not in dir(options):
                setattr(options, s, default_radiorx_setup.d_options[s] )

    get_options = Callable (get_options)
    check_options = Callable (check_options)

class RadioRx(gr.hier_block2):
    
    def __init__(self, options):
        # check for default options before continuing with constructor
        default_radiorx_setup.check_options(options)
        self.verbose = options.verbose
        self.fake_rf = options.fake_rf

        self.rt = None
        nusrp = int((options.nchannels+1)/2.0 )
        nchan_per_usrp = 1
        if (options.nchannels>1): nchan_per_usrp = 2;
        do_rx = True
        do_tx = False
        
        # instantiate blocks
        self.src = self.make_src (options)

        # set options
        self.set_nchannels (options.nchannels)
        sys.stderr.write("[radiorx] set n_channels = %s\n"%(options.nchannels) )
        self.set_sample_rate (options.sample_rate)
        sys.stderr.write("[radiorx] set sample_rate = %.1f MHz\n"%(options.sample_rate/1e6) )
        self.set_subdev (options.subdev_spec)
        sys.stderr.write("[radiorx] set subdev = %s\n"%(options.subdev_spec) )
        self.set_freq (options.freq)
        sys.stderr.write("[radiorx] set freq = %.1f MHz\n"%(options.freq/1e6) )
        self.set_rx_gain (options.rx_gain)
        sys.stderr.write("[radiorx] set rx_gain = %s dB\n"%(options.rx_gain) )

        # connect blocks and call hier_block2 constructor
        gr.hier_block2.__init__(self, "RadioRx",
                                gr.io_signature(0, 0, 0), 
                                gr.io_signature(1, 1, gr.sizeof_gr_complex))
        self.connect(self.src, self)

    def make_src (self, options):
        if options.fake_rf:
            return rflib.fake_rx()
        else:
            return usrp.source_c(options.which_board)

    def shutdown(self):
        sys.stderr.write("radiorx.shutdown called ...\n")
        if self.fake_rf: self.src.shutdown()

    def set_nchannels (self, n):
        """ call set_sample_rate and set_subdev after setting nchannels """
        self.nchannels = n
        return self.src.set_nchannels(n)

    def set_sample_rate (self, s):
        """ assumes nchannels has been set """
        sizeof_sample = 4*self.nchannels  # 4 bytes/complex sample = 2 shorts
        max_rate = MAX_USB_RATE / (2 *sizeof_sample) # half-duplex constraint
        self.sample_rate = min (s, max_rate)
        
        # get decimation rates (lower bound on sample rate)
        MAX_DECIM = 256
        ulist = [self.src]
        for u in ulist:
            rdecim = min (int(u.adc_rate()/self.sample_rate), MAX_DECIM)
            sys.stderr.write("[radiorx]: setting decimation rate to %d on %s\n"%(rdecim, str(u) ) )
            
            # configure src - this "should" correspond to radiotx sample rate
            self.sample_rate = u.adc_rate()/rdecim   # these "should" be
            sys.stderr.write("[radiorx]: sample rate = %1.1f MHz\n"%(self.sample_rate*1e-6) )
            r = u.set_decim_rate(rdecim)
        return r

    def set_subdev (self, spec=None):
        """ assumes nchannels has been set """
        """ call set_freq after subdev is set """
        self.subdev = ()
        if self.fake_rf: return     # no subdev for fake rf
        """ get subdev spec """
        if spec is None:
            srx1 = usrp.pick_rx_subdevice(self.src)
        else:
            srx1 = spec

        """ configure USRP mux and set rx/tx subdev """
        ulist = [self.src]
        if self.nchannels == 1:
            assert len(ulist) == 1, "[radiorx]: Invalid number of USRP boards!"
            src = ulist[0]
            src.set_mux (usrp.determine_rx_mux_value(src, srx1) )
            self.subdev += (usrp.selected_subdev(src, srx1), )
            for s in self.subdev: exec("if not hasattr(s, '_u'): s._u = src")
        elif self.nchannels == 2:
            assert len(ulist) == 1, "[radiorx]: Invalid number of USRP boards!"
            srx2 = ((srx1[0]+1)%2, srx1[1])
            src = ulist[0]
            src.set_mux(gru.hexint(0x23012301) )
            #src.set_mux(gru.hexint(0x32103210) )
            self.subdev += (usrp.selected_subdev(src, srx1), )
            self.subdev += (usrp.selected_subdev(src, srx2), )
            for s in self.subdev: exec("if not hasattr(s, '_u'): s._u = src")
        elif self.nchannels == 4:
            assert len(ulist) == 2, "[radiorx]: Invalid number of USRP boards!"
            srx2 = ((srx1[0]+1)%2, srx1[1])
            for src in ulist:
                src.set_mux(gru.hexint(0x23012301) )
                self.subdev += (usrp.selected_subdev(src, srx1), )
                self.subdev += (usrp.selected_subdev(src, srx2), )
                for s in self.subdev: exec("if not hasattr(s, '_u'): s._u = src")
        else:
            return self.error("Unable to set subdev; invalid value for " \
                              + "nchannels=%d"%(self.nchannels) )
        self.set_auto_tr(True)  # enable Automatic Tx/Rx Switching

    def set_freq (self, f):
        """ assumes subdev has been set """
        if abs(f) < 1e6: f = f*1e6
        self.freq = f
        if self.fake_rf: return     # no subdev for fake rf
        for s in self.subdev:
            r = usrp.tune(s._u, s.which(), s, self.freq)
            if r and (self.verbose > 0):
                sys.stderr.write("setting frequency of %s:\n"%(str(s) ) )
                sys.stderr.write("   baseband frequency  = %s\n"%(eng_notation.num_to_str(r.baseband_freq)) )
                sys.stderr.write("   DUC/DDC offset      = %s\n"%(eng_notation.num_to_str(r.dxc_freq)) )
                sys.stderr.write("   residual frequency  = %s\n"%(eng_notation.num_to_str(r.residual_freq)) )
                sys.stderr.write("   inverted            = %s\n"%(r.inverted) )
            elif not r:
                self.error("Unable to set frequency of " \
                           + "%s to %g MHz"%(str(s),self.freq/1.0e6) )

    def set_rx_gain (self, g):
        self.rx_gain = default_radiorx_setup.d_options['rx_gain']
        if self.fake_rf: return     # no subdev for fake rf
        for s in self.subdev:
            gain_range = s.gain_range()     # [min_gain, max_gain]
            self.rx_gain = max(min(g, gain_range[1]), gain_range[0] )
            s.set_gain(self.rx_gain)
            if self.verbose>0:
                sys.stderr.write("[radiorx]: setting rx gain = %.1f dB on %s\n"%(self.rx_gain, s) )
                sys.stderr.write("           gain range = (%.1f dB, %.1f dB)\n"%(gain_range[0], gain_range[1]) )

    def set_auto_tr(self, enable):
        if self.fake_rf: return     # no subdev for fake rf
        for s in self.subdev:
            s.set_auto_tr(enable)

    def error (self, msg, level=0):
        if self.verbose >= level: sys.stderr.write("RFRX ERROR: "+str(msg)+"\n")
    
    """ Add parser options to an OptionParser """
    def add_parser_options (parser):
        if not parser.has_option("-v"):
            parser.add_option ("-v", "--verbose", type="int", \
                    default=default_radiorx_setup.d_options['verbose'], \
                    help="set verbose level of output [default=%default]")
        if not parser.has_option("--fake-rf"):
            parser.add_option("", "--fake-rf", action="store_true", \
                    default=default_radiorx_setup.d_options['fake_rf'],
                    help="enable \"fake\" RF for emulation [default=%default]")
        if not parser.has_option("-w"):
            parser.add_option ("-w", "--which-board", type="int", \
                    default=default_radiorx_setup.d_options['which_board'], \
                    help="select which USRP board to use [default=%default]")
        if not parser.has_option("--master-serialno"):
            parser.add_option ("", "--master-serialno", type="str", \
                    default=default_radiorx_setup.d_options['master_serialno'], \
                    help="specify serial no. for master USRP board [default=%default]")
        if not parser.has_option("-n"):
            parser.add_option ("-n", "--nchannels", type="int", \
                    default=default_radiorx_setup.d_options['nchannels'], \
                    help="set number of channels (or antennas) on USRP board [default=%default]")
        if not parser.has_option("-S"):
            parser.add_option ("-S", "--subdev-spec", type="subdev", \
                    default=default_radiorx_setup.d_options['subdev_spec'], \
                    help="select USRP Tx/RX side A or B")
        if not parser.has_option("-s"):
            parser.add_option ("-s", "--sample-rate", type="eng_float", \
                    default=default_radiorx_setup.d_options['sample_rate'], \
                    help="set usrp sample rate [default=%default]")
        if not parser.has_option("-f"):
            parser.add_option ("-f", "--freq", type="eng_float", \
                    default=default_radiorx_setup.d_options['freq'], \
                    help="set carrier frequency [default=%default]")
        if not parser.has_option("-G"):
            parser.add_option ("-G", "--rx-gain", type="eng_float", \
                    default=default_radiorx_setup.d_options['rx_gain'], \
                    help="set usrp receive gain in dB [default=%default]")

    add_parser_options = Callable (add_parser_options)
