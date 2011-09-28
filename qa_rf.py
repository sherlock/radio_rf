#!/usr/bin/env python
#
# Copyright (C) 2009  The University of Texas at Austin.
# Copyright 2004 Free Software Foundation, Inc.
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
# 

from gnuradio import gr, gr_unittest
import rf
import phy
import struct
import Numeric

class qa_rf (gr_unittest.TestCase):
    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    def make_waveform(self, s):
        BPSK_1 = complex(1,0)
        BPSK_0 = complex(-1,0)
        waveform = []
        for k in range(len(s) ):
            mask = 1<<7
            for n in range(8):
                x = (ord(s[k]) & mask)> 0
                mask = mask >> 1
                if (x): waveform.append(BPSK_1)
                else:   waveform.append(BPSK_0)
        return Numeric.array(waveform, Numeric.Complex32)

    def tobytes(self, c):
        s = ""
        for x in c:
            s += struct.pack('<ff', x.real, x.imag)
        return s

    def no_test_001_pad(self):
        s = "abcd"
        x = self.make_waveform(s)
        m = self.tobytes(x)
        msg = gr.message_from_string(m)
        Ntx = 2
        msg.set_arg1(Ntx)
        msg.set_arg2(len(x)/Ntx)

        pad = rf.pad()
        dst = gr.vector_sink_c()

        pad.inputQ().insert_tail(msg)
        pad.shutdown()  # schedules shutdown

        self.tb.connect (pad, dst)
        pad.set_npost(3)
        self.tb.run()
        result_data = dst.data ()
        
        block_size = pad.block_size()
        # prepad:(1+npre)*Ntx, postpad:(1+npost)*Ntx
        pre_pad = int(1+pad.npre())*Ntx*block_size
        post_pad = int(1+pad.npost())*Ntx*block_size
        pad_size = pre_pad + post_pad - len(x)%block_size
        self.assertComplexTuplesAlmostEqual (tuple(x), \
                   result_data[pre_pad:pre_pad+len(x)], len(x)+1 )
        self.assertEqual(len(result_data), len(x)+pad_size )
        self.assertEqual(phy.phyglobal.ctrlQ().count(), 1)
        
    def no_test_002_fake_tx(self):
        s = "abcd"
        x = tuple(self.make_waveform(s).tolist() )
        src = gr.vector_source_c(x)

        dst = rf.fake_tx()
        dst.set_nchannels(1)
        dst.set_interp_rate(128)
        
        self.tb.connect (src, dst)
        self.tb.run()

        expected_data = self.tobytes(x)
        self.assertEqual(dst.outputQ().count(), 1)
        msg = dst.outputQ().delete_head()
        
        type, arg1, arg2 = 0, 1.0, 1.0e6
        self.assertAlmostEqual(msg.type(), type)
        self.assertAlmostEqual(msg.arg1(), arg1)
        self.assertAlmostEqual(msg.arg1(), arg1)
        self.assertEqual (msg.to_string(), expected_data)

    def no_test_003_fake_rx(self):
        s = "abcd"
        x = tuple(self.make_waveform(s).tolist() )
        src = gr.vector_source_c(x)
        tx = rf.fake_tx()

        ntb = gr.top_block()
        ntb.connect(src, tx)
        ntb.run()

        rx = rf.fake_rx()
        dst = gr.vector_sink_c()
        msg = tx.outputQ().delete_head()

        rx.inputQ().insert_tail(msg)
        rx.shutdown()  # schedules shutdown

        self.tb.connect(rx, dst)
        self.tb.run()
        result_data = dst.data ()
        
        self.assertComplexTuplesAlmostEqual (x, result_data, len(x)+1 )

if __name__ == '__main__':
    gr_unittest.main ()
