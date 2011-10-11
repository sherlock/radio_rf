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

MAX_USB_RATE = 32e6       #32MB/sec

class default_radiotx_setup:
    d_options = {'verbose':0, 'debug':1, 'fake_rf':False, 'disable_tx':False,
                 'which_board':0, 'nchannels':1, 'subdev_spec':"A:0",
                 'sample_rate':1.0e6, 'freq':2400.0e6,
                 'tx_gain': 10.0, 'master_serialno':None, 'address':"type=usrp1"}

    def get_options():
        dopt = default_radiotx_setup()
        for s in default_radiotx_setup.d_options:
            if s not in dir(dopt):
                setattr(dopt, s, default_radiotx_setup.d_options[s] )
        return dopt
    
    def check_options(options):
        for s in default_radiotx_setup.d_options:
            if s not in dir(options):
                setattr(options, s, default_radiotx_setup.d_options[s] )

    get_options = Callable (get_options)
    check_options = Callable (check_options)

class RadioTx(gr.hier_block2):
    def __init__(self, options):
        # check for default options before continuing with constructor
        default_radiotx_setup.check_options(options)
        self.verbose = options.verbose
        self.debug = options.debug
        self.fake_rf = options.fake_rf or options.disable_tx

        self.rt = None
        nusrp = int((options.nchannels+1)/2.0 )
        nchan_per_usrp = 1
        if (options.nchannels>1): nchan_per_usrp = 2;
        do_rx = False
        do_tx = True

        # instantiate blocks
        self.pad = rflib.pad()
        #self.pad.set_npre(4)    # increase pre-zero padding
        #self.pad.set_npost(4)   # increase post-zero padding
        self.sink = self.make_tx (options)
        
        # set debug option
        if hasattr(self.pad, 'set_debug'): self.pad.set_debug(options.debug)

        # set other options
        #self.set_nchannels (options.nchannels)
        self.set_sample_rate (options.sample_rate)
        self.set_subdev (options.subdev_spec)
        self.set_freq (options.freq)
        self.set_tx_gain (options.tx_gain)

        # connect blocks and call hier_block2 constructor
        gr.hier_block2.__init__(self, "RadioTx",
                                gr.io_signature(0, 0, 0), 
                                gr.io_signature(0, 0, 0)) 
        self.connect(self.pad, self.sink)

    def make_tx (self, options):
        self.nchannels = options.nchannels
        if self.fake_rf:
            return rflib.fake_tx()
        else:
            #return usrp.sink_c(options.which_board)
            return uhd.usrp_sink(device_addr=options.address,\
                                   io_type=uhd.io_type.COMPLEX_FLOAT32,\
                                   num_channels=options.nchannels)

    def shutdown(self):
        sys.stderr.write("[radiotx] shutdown called ...\n")
        if self.pad: self.pad.shutdown()

    '''
    def set_nchannels (self, n):
        """ call set_sample_rate and set_subdev after setting nchannels """
        self.nchannels = n
        return self.sink.set_nchannels(n)
    '''

    def set_sample_rate (self, s):
        """ assumes nchannels has been set """
        sizeof_sample = 4 * self.nchannels  # 4 bytes/complex sample = 2 shorts
        max_rate = MAX_USB_RATE / (2 *sizeof_sample) # half-duplex constraint
        self.sample_rate = min (s, max_rate)
        self.sink.set_samp_rate(self.sample_rate)
        
        '''
        # get interpolation rate (lower bound on sample rate)
        MAX_INTERP = 512
        ulist = [self.sink]
        for u in ulist:
          rinterp = min (int(u.dac_rate()/self.sample_rate), MAX_INTERP)

          # configure sink - this "should" correspond to radiorx sample rate
          self.sample_rate = u.dac_rate()/rinterp
          r = u.set_interp_rate(rinterp)
          sys.stderr.write("[radiotx]: setting interpolation rate to %d on %s\n"%(rinterp, str(u) ) )
          sys.stderr.write("[radiotx]: sample rate = %1.1f MHz\n"%(self.sample_rate*1e-6) )
        return r
        '''
    def set_subdev (self, spec=None):
        """ assumes nchannels has been set """
        """ call set_freq after subdev is set """
        if self.fake_rf: return
        if self.nchannels ==1:
            if spec is None:
                spec = "A:0"   #default is DBoard A
        elif  self.nchannels ==2:
            spec = "A:0 B:0"
        self.sink.set_subdev_spec(spec)

    '''
    def set_subdev (self, spec=None):
        """ assumes nchannels has been set """
        """ call set_freq after subdev is set """
        self.subdev = ()
        if self.fake_rf: return     # no subdev for fake rf
        """ get subdev spec """
        if spec is None:
            stx1 = usrp.pick_tx_subdevice(self.sink)
        else:
            stx1 = spec
        
        """ configure USRP mux and set tx subdev """
        ulist = [self.sink]
        if self.nchannels == 1:
            assert len(ulist) == 1, "[radiotx]: Invalid number of USRP boards!"
            sink = ulist[0]
            sink.set_mux (usrp.determine_tx_mux_value(sink, stx1) )
            self.subdev += (usrp.selected_subdev(sink, stx1), )
            for s in self.subdev: exec("if not hasattr(s, '_u'): s._u = sink")
        elif self.nchannels == 2:
            assert len(ulist) == 1, "[radiotx]: Invalid number of USRP boards!"
            stx2 = ((stx1[0]+1)%2, stx1[1])
            sink = ulist[0]
            sink.set_mux(gru.hexint(0xBA98) )
            self.subdev += (usrp.selected_subdev(sink, stx1), )
            self.subdev += (usrp.selected_subdev(sink, stx2), )
            for s in self.subdev: exec("if not hasattr(s, '_u'): s._u = sink")
        elif self.nchannels == 4:
            assert len(ulist) == 2, "[radiorx]: Invalid number of USRP boards!"
            stx2 = ((stx1[0]+1)%2, stx1[1])
            for sink in ulist:
                sink.set_mux(gru.hexint(0xBA98) )
                self.subdev += (usrp.selected_subdev(sink, stx1), )
                self.subdev += (usrp.selected_subdev(sink, stx2), )
                for s in self.subdev: exec("if not hasattr(s, '_u'): s._u = sink")
        else:
            return self.error("Unable to set subdev; invalid value for " \
                    + "nchannels=%d"%(self.nchannels) )
        #self.set_auto_tr(True)  # enable Automatic  
        '''

    def set_freq (self, f):
        """ assumes subdev has been set """
        if abs(f) < 1e6: f = f*1e6
        self.freq = f
        if self.fake_rf: return     # no subdev for fake rf
        for i in range(self.nchannels):
            r = self.sink.set_center_freq(self.freq, i)
        
        '''
        for s in self.subdev:
            r = usrp.tune(s._u, s.which(), s, self.freq)
            if r and (self.verbose > 4):
                sys.stderr.write("setting frequency of %s:\n"%(str(s) ) )
                sys.stderr.write("   baseband frequency  = %s\n"%(eng_notation.num_to_str(r.baseband_freq) ) )
                sys.stderr.write("   DUC/DDC offset      = %s\n"%(eng_notation.num_to_str(r.dxc_freq) ) )
                sys.stderr.write("   residual frequency  = %s\n"%(eng_notation.num_to_str(r.residual_freq) ) )
                sys.stderr.write("   inverted            = %s\n"%(r.inverted) )
            elif not r:
                self.error("Unable to set frequency of " \
                         + "%s to %g MHz"%(str(s),self.freq/1.0e6) 
        '''
        
    

    def set_tx_gain (self, g):
        self.tx_gain = g
        if self.fake_rf: return     # no subdev for fake rf
        
        for i in range(self.nchannels):
            gain_range = self.sink.get_gain_range(i).to_pp_string()
            gain_range = gain_range[1:(gain_range.find(')'))]
            ulist = gain_range.split(',')   #[min_gain, max_gain]
            self.tx_gain = max(min(g, ulist[1]), ulist[0] )
            
            self.sink.set_gain(self.tx_gain, i)
        
        
        
        # no tx gain in usrp
        #pass

    '''
    def set_auto_tr(self, enable):
        if self.fake_rf: return     # no subdev for fake rf
        for s in self.subdev:
            s.set_auto_tr(enable)
            #s.set_enable(enable)
    '''

    def error (self, msg, level=0):
        if self.verbose >= level: sys.stderr.write("RFTX ERROR: "+str(msg)+"\n")
    
    """ Add parser options to an OptionParser """
    def add_parser_options (parser):
        if not parser.has_option("-v"):
            parser.add_option ("-v", "--verbose", type="int", \
                    default=default_radiotx_setup.d_options['verbose'], \
                    help="set verbose level of output [default=%default]")
        if not parser.has_option("--debug"):
            parser.add_option ("", "--debug", type="int", \
                    default=default_radiotx_setup.d_options['debug'], \
                    help="set debug level of system [default=%default]")
        if not parser.has_option("--fake-rf"):
            parser.add_option("", "--fake-rf", action="store_true", \
                    default=default_radiotx_setup.d_options['fake_rf'],
                    help="enable \"fake\" RF for emulation [default=%default]")
        if not parser.has_option("--disable-tx"):
            parser.add_option("", "--disable-tx", action="store_true", \
                    default=default_radiotx_setup.d_options['disable_tx'],
                    help="replace transmitter with a fake rf [default=%default]")
        if not parser.has_option("-w"):
            parser.add_option ("-w", "--which-board", type="int", \
                    default=default_radiotx_setup.d_options['which_board'], \
                    help="select which USRP board to use [default=%default]")
        if not parser.has_option("--master-serialno"):
            parser.add_option ("", "--master-serialno", type="str", \
                    default=default_radiotx_setup.d_options['master_serialno'], \
                    help="specify serial no. for master USRP board [default=%default]")
        if not parser.has_option("-n"):
            parser.add_option ("-n", "--nchannels", type="int", \
                    default=default_radiotx_setup.d_options['nchannels'], \
                    help="set number of channels (or antennas) on USRP board [default=%default]")
        if not parser.has_option("-S"):
            parser.add_option ("-S", "--subdev-spec", type="string", \
                    default=default_radiotx_setup.d_options['subdev_spec'], \
                    help="select USRP Tx/RX side A or B")
        if not parser.has_option("-s"):
            parser.add_option ("-s", "--sample-rate", type="eng_float", \
                    default=default_radiotx_setup.d_options['sample_rate'], \
                    help="set usrp sample rate [default=%default]")
        if not parser.has_option("-f"):
            parser.add_option ("-f", "--freq", type="eng_float", \
                    default=default_radiotx_setup.d_options['freq'], \
                    help="set carrier frequency [default=%default]")
        if not parser.has_option("-g"):
            parser.add_option ("-g", "--tx-gain", type="eng_float", \
                    default=default_radiotx_setup.d_options['tx_gain'], \
                    help="set software transmit gain in dB [default=%default]")
        if not parser.has_option("-a"):
            parser.add_option("-a","--address",type="string", \
                    default=default_radiotx_setup.d_options['address'],\
                    help="Address of UHD device, [default=%default]")

    add_parser_options = Callable (add_parser_options)
