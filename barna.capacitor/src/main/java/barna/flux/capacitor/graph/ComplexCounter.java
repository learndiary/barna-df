package barna.flux.capacitor.graph;

import java.util.EnumSet;
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

    static enum CounterType {SIMPLE, PAIRED, STRANDED, COMBINED, INTRON};

    static final int DEFAULT_NR_BINS_INTRONS = 10;

    int nrBinsIntron;

    Hashtable<String,int[]> map= new Hashtable<String, int[]>();
    int nrCounter= -1;
    HashMap<CounterType,Byte> mapType= null;

    public ComplexCounter(EnumSet<CounterType> counterIDs, int nrBinsIntron) {
        this.nrBinsIntron = nrBinsIntron;
        nrCounter = counterIDs.size();
        if (counterIDs.contains(CounterType.INTRON))
            nrCounter += nrBinsIntron - 1;

        mapType= new HashMap<CounterType, Byte>(counterIDs.size());
        byte pos= 0;
        for (CounterType counterID : counterIDs) {
            mapType.put(counterID, pos);
            pos += (counterID.equals(CounterType.INTRON) ? nrBinsIntron : 1);
        }
    }

    public ComplexCounter(EnumSet<CounterType> counterIDs) {
        this(counterIDs, DEFAULT_NR_BINS_INTRONS);
    }


    public void increment(String ID, CounterType type) {
        increment(ID, type, 1, 0);
    }

    public void increment(String ID, CounterType type, int value) {
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

    public void increment(String ID, CounterType type, int value, int offset) {

        int[] ctr = getBins(ID);

        // get index to increment
        ctr[mapType.get(type) + offset] += value;
    }


    public int get(String ID, CounterType type, int offset) {
        int[] ctr= getBins(ID);

        return ctr[mapType.get(type) + offset];
    }

    public Iterator<String> getIDs() {
        return map.keySet().iterator();
    }

}
