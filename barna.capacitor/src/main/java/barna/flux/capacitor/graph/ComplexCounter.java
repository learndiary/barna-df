package barna.flux.capacitor.graph;

import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 11/21/12
 * Time: 2:46 PM
 */
public class ComplexCounter {

    static final byte COUNTER_PAIRED= 0;
    static final byte COUNTER_SIMPLE= 1;
    static final byte COUNTER_STRANDED= 2;
    static final byte COUNTER_COMBINED= 3;
    static final byte COUNTER_INTRON= 4;

    static final int DEFAULT_NR_BINS_INTRONS= 10;

    int nrBinsIntron= DEFAULT_NR_BINS_INTRONS;

    Hashtable<String,int[]> map= new Hashtable<String, int[]>();
    int nrCounter= -1;
    HashMap<Byte,Byte> mapType= null;

    public ComplexCounter(byte[] counterIDs, int nrBinsIntron) {
        this.nrBinsIntron= nrBinsIntron;
        nrCounter= counterIDs.length;
        for (int i = 0; i < counterIDs.length; i++) {
            if (counterIDs[i]== COUNTER_INTRON)
                nrCounter+= nrBinsIntron- 1;
        }

        mapType= new HashMap<Byte, Byte>(counterIDs.length);
        byte pos= 0;
        for (int i = 0; i < counterIDs.length; i++) {
            mapType.put(counterIDs[i], pos);
            pos+= (counterIDs[i]== COUNTER_INTRON? nrBinsIntron: 1);
        }
    }

    public ComplexCounter(byte[] counterIDs) {
        this(counterIDs, DEFAULT_NR_BINS_INTRONS);
    }


    public void increment(String ID, byte type) {
        increment(ID, type, 1, 0);
    }

    public void increment(String ID, byte type, int value) {
        increment(ID, type, value, 0);
    }

    private int[] getBins(String ID) {
        // get counter array
        int[] ctr= null;
        if (map.contains(ID)) {
            ctr= map.get(ID);
        } else {
            ctr= new int[nrCounter];
            for (int i = 0; i < ctr.length; i++)
                ctr[i]= 0;

            map.put(ID, ctr);
        }

        return ctr;
    }

    public void increment(String ID, byte type, int value, int offset) {

        int[] ctr= getBins(ID);

        // get index to increment
        ctr[mapType.get(type)+ offset]+= value;
    }


    public int get(String ID, byte type, int offset) {
        int[] ctr= getBins(ID);

        return ctr[mapType.get(type)+ offset];
    }

    public Iterator<String> getIDs() {
        return map.keySet().iterator();
    }

}
