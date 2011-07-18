/*
 * This file is part of the Flux Library.
 *
 * The code of the Flux Library may be freely distributed and modified under the terms of the
 * European Union Public Licence (EUPL) published on the web site <http://www.osor.eu/eupl/european-union-public-licence-eupl-v.1.1>.
 * Copyright for the code is held jointly by the individual authors, who should be listed
 * in @author doc comments. According to Article 5 and Article 11 of the EUPL, publications that
 * include results produced by the Flux Library are liable to reference the Work,
 * see the Flux Library homepage <http://flux.sammeth.net> for more information.
 */

package fbi.genome.errormodel;

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
//            w.write("\n");
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
    }
}
