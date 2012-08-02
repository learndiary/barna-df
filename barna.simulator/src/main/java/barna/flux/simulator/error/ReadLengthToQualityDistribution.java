/*
 * Copyright (c) 2010, Micha Sammeth
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * The names of its contributors may be not used to endorse or promote
 *       products derived from this software without specific prior written
 *       permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL MICHA SAMMETH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package barna.flux.simulator.error;

/**
 * Distribution of the quality per read position
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class ReadLengthToQualityDistribution extends Distribution {


    //BufferedWriter w;

    public ReadLengthToQualityDistribution(int size) {
        super(size);

//        try {
//            w = new BufferedWriter(new FileWriter("/home/thasso/Desktop/readQuals_35.dat"));
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
    }

    public void addRead(Read read) {
        reads++;
        int[] q = read.getQualities();
        for (int i = 0; i < size; i++) {
            values[i] += q[i];
//            try {
//                w.write(Integer.toString(q[i]));
//                if(i<size-1){
//                    w.write("\t");
//                }
//            } catch (IOException e) {
//                e.printStackTrace();
//            }
        }

//        try {
//            w.write(barna.commons.system.OSChecker.NEW_LINE);
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
    }
}
