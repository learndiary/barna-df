package fbi.genome.errormodel;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Distribution of the quality per read position
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class ReadLengthToQualityDistribution extends Distribution{


    BufferedWriter w;

    public ReadLengthToQualityDistribution(int size) {
        super(size);

        try {
            w = new BufferedWriter(new FileWriter("/home/thasso/Desktop/readQuals_35.dat"));
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void addRead(Read read){
        reads++;
        int[] q = read.getQualities();
        for (int i = 0; i < size; i++) {
            values[i] += q[i];
            try {
                w.write(Integer.toString(q[i]));
                if(i<size-1){
                    w.write("\t");
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        try {
            w.write("\n");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
