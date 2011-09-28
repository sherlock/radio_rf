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

from gnuradio.eng_option import eng_option
from optparse import OptionParser

from gnuradio import gr
from hydra import PyHydra
from hydra.PyHydra import RadioTx, RadioRx, FuncThread
import numpy
import struct
import time

SEND_MIMO_TRAINING = 1
if SEND_MIMO_TRAINING:
    from hydra import mimo, phy

class Signal(object):
    def __init__(self, data=None):
        self.setData(data)
    def setData(self,data):
        if data is None: self.data = numpy.array([])
        else: self.data = numpy.array(data) * complex(1,0)
        self.sdata = ""
        for x in self.data:
            self.sdata += struct.pack('<ff', x.real, x.imag)
    def string(self, k=1):
        return self.sdata*k
    def vector(self):
        return self.data

def send_stoprf(rf, d=0.0):
    if (d>0.0): time.sleep(d)
    rf.shutdown()

def send_msg(m, dst, d=0.0):
    # countdown
    print "starting transmission in "
    for k in range(3,0,-1):
        print "%d ..."%(k)
        time.sleep(0.25)
    if SEND_MIMO_TRAINING: R = 100
    else: R = 1
    # start
    for k in range(R):
        if (d>0.0): time.sleep(d)
        print "sending msg ..."
        dst.insert_tail(m)

def make_exponential(A, T, f, fs=1.0e6):
    t = numpy.arange(int(T*fs) )
    x = A*numpy.cos(2*numpy.pi*(f/fs)*t) \
        + complex(0,1)*A*numpy.sin(2*numpy.pi*(f/fs)*t)
    return x

def TestTx():
    """ Parser Options """
    parser = OptionParser (option_class=eng_option)
    RadioTx.add_parser_options (parser)
    (options, args) = parser.parse_args ()
    if len(args)!=0:
        parser.print_help()
        sys.exit(1)

    print "TestTx\n=========="

    tb = gr.top_block()
    tx = RadioTx(options)
    tb.connect(tx)

    #T, w, C = 10.0, 100.0e3, 100
    T, w, C = 10.0, 40.0e3, 100
    A = 4000
    if SEND_MIMO_TRAINING:
        T = 10.0
        s = "hi"*4
        txvector = mimo.txvector()
        txvector.set_length(len(s) )
        phy = mimo.tx(options.nchannels)
        phy.configure(txvector)
        phy.set_gain(10000)
        phy.set_upsampling_factor(2)
        phy.transmit(s)
        m = phy.outputQ().delete_head()
        print "created message ... %d samples"%(m.length()/gr.sizeof_gr_complex/options.nchannels)
    else:
        P = min(T, C/w)
        k = max(1, int(T/P) )
        s = Signal()
        s.setData(make_exponential(A, P, w, options.sample_rate) )
        m = gr.message_from_string(s.string(k), 0, 1, len(s.vector())*k )
        print "created message ... %d samples"%(int(k*P*options.sample_rate) )

    interval = 0.2
    g = FuncThread(send_msg, "send msg", m, tx.pad.inputQ(), interval)
    g.setDaemon(1)
    g.start()
    
    stoptime = T+5.0
    f = FuncThread(send_stoprf, "stop rf", tx, stoptime)
    f.setDaemon(1)
    f.start()
    
    print "running ..."
    tb.run()
    print "stoppping ..."

def stop_tb(tb, d=0.0):
    if (d>0.0): time.sleep(d)
    tb.stop()

from pylab import *

def TestRx():
    """ Parser Options """
    parser = OptionParser (option_class=eng_option)
    RadioRx.add_parser_options (parser)
    (options, args) = parser.parse_args ()
    if len(args)!=0:
        parser.print_help()
        sys.exit(1)

    print "TestRx\n=========="

    tb = gr.top_block()
    rx = gr.vector_sink_c()
    rf = RadioRx(options)
    tb.connect(rf,rx)
    
    stoptime = 2.0
    if options.fake_rf:
        f = FuncThread(send_stoprf, "stop rf", rf, stoptime)
    else:
        f = FuncThread(stop_tb, "stop tb", tb, stoptime)
    f.setDaemon(1)
    f.start()
    
    print "running ..."
    tb.run()
    print "stoppping ..."

    x = numpy.array(rx.data(), numpy.complex64)
    #start, limit = 200, len(x)
    start, limit = 200, 1000
    L = min(len(x), start+limit)
    xdata = x[start:L]
    plot(range(len(xdata)), xdata.real ) 
    hold(True)
    plot(range(len(xdata)), xdata.imag, 'r')
    #plot(numpy.abs(x[start:L] ) )
    grid()
    print "plotting ..."
    show()
    print "quitting ..."
    fname = "/tmp/uuu"
    FILE = open(fname, "w")
    format = '<'+"ff"*len(x)
    args = tuple(numpy.array((x.real,x.imag)).transpose().flatten() )
    s = struct.pack(format, *args )
    FILE.write(s)
    FILE.close()
    

def TestDC():
    """ Parser Options """
    parser = OptionParser (option_class=eng_option)
    RadioRx.add_parser_options (parser)
    (options, args) = parser.parse_args ()
    if len(args)!=0:
        parser.print_help()
        sys.exit(1)

    print "TestRx\n=========="

    tb = gr.top_block()
    rx = gr.vector_sink_c()
    rf = RadioRx(options)
    tb.connect(rf,rx)
    
    stoptime = 20.0
    if options.fake_rf:
        f = FuncThread(send_stoprf, "stop rf", rf, stoptime)
    else:
        f = FuncThread(stop_tb, "stop tb", tb, stoptime)
    f.setDaemon(1)
    f.start()
    
    print "running ..."
    tb.run()
    print "stoppping ..."

    n = numpy.array(rx.data() )
    n = n[len(n)/2:]
    print "noise mean =", n.mean()
    print "noise var  =", n.var()
    print "noise std  =", n.std()

def main():
    #TestTx()
    #TestRx()
    TestDC()

if __name__=='__main__':
    main()
