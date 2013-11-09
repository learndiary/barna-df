package barna.io.sam;

import net.sf.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;

/**
 * A set of <code>SAMRecord</code> instances that can be dynamically extended. Filter criteria can be applied to
 * determine the records that are kept.
 * User: micha
 */
public class FilteredSAMRecordSet {

    /**
     * Comparator to compare two <code>SAMRecord</code> instances by their number of mismatches
     */
    public static class NMComparator implements Comparator<SAMRecord> {

        /**
         * Compares two mappings by the number of mismatches they show with the reference sequence. Note that the
         * <code>NM</code> tag has to be present in both of the <code>SAMRecord</code> instances passed to the method,
         * otherwise a <code>RuntimeException</code> is triggered.
         * @param r1 a mapping
         * @param r2 another mapping
         * @return a negative, <code>0</code> or a positive value corresponding to the first argument having less,
         * equally many or more mismatches than the second argument.
         */
        @Override
        public int compare(SAMRecord r1, SAMRecord r2) {

            Integer nm1= r1.getIntegerAttribute(SAMConstants.SAM_OPTION_NM),
                    nm2= r2.getIntegerAttribute(SAMConstants.SAM_OPTION_NM);

            if (nm1== null|| nm2== null)
                throw new IllegalArgumentException("Mapping missing tag "+ SAMConstants.SAM_OPTION_NM+ ":\n"
                        +(nm1== null?r1.toString():"")+ (nm2== null?r2.toString():""));

            return nm1.compareTo(nm2);
        }
    }

    /**
     * Default number of strata that are kept, starting with the minimum.
     */
    public static final byte DEFAULT_MAX_STRATA= 2;

    /**
     * Singleton instance of the NMComparator.
     */
    protected static Comparator defaultNMComparator= null;

    /**
     * Number of strata that are kept, starting with the minimum.
     */
    protected byte maxStrata = DEFAULT_MAX_STRATA;

    /**
     * The mapping instances currently stored in the set.
     */
    ArrayList<SAMRecord> set= null;

    /**
     * Strata of the mappings currently stored in the set.
     */
    int[] currStrata= null;

    int currStrataLength= -1;

    /**
     * Generates a singleton of the <code>NMComparator</code>.
     * @return the <code>NMComparator</code> singleton
     */
    public static Comparator getDefaultNMComparator() {
        if (defaultNMComparator== null) {
            defaultNMComparator= new NMComparator();
        }
        return defaultNMComparator;
    }

    /**
     * Creates a default set, filtering <code>DEFAULT_MAX_STRATA</code> strata while new <code>SAMRecord</code>
     * instances are added.
     */
    public FilteredSAMRecordSet() {
        this(6, DEFAULT_MAX_STRATA);
    }

    /**
     * Constructor to specify the number of strata and initial size for the set.
     * @param initialSize the initial size of the internal vector used to sore the mappings of this set. Per default,
     *                    the datastructure will extend by 50% each time the capacity needs to be increased.
     * @param maxStrata the number of strata (starting with the lowest one) stored in the set
     */
    public FilteredSAMRecordSet(int initialSize, byte maxStrata) {
        set= new ArrayList<SAMRecord>(initialSize);
        this.maxStrata= maxStrata;
        currStrata= new int[maxStrata];
        Arrays.fill(currStrata, Integer.MAX_VALUE);
        currStrataLength= 0;
    }

    /**
     * Checks whether a given mapping is added to the set, and in case adapts the distribution of strata currently
     * stored in <code>this</code> set.
     * @param mapping a mapping
     * @return the change in the strata that is stored: 0= no change (mapping discarded or describes a stratum that is
     * already stored), 1= mapping added new stratum, (-1)= a previously stored mapping was lost due to the newly
     * inserted stratum described by the mapping.
     */
    public int add(SAMRecord mapping) {

        int nm= mapping.getIntegerAttribute(SAMConstants.SAM_OPTION_NM);
        if (currStrataLength== currStrata.length && nm> currStrata[currStrata.length- 1])
            return 0;   // discard mapping

        // otherwise add mapping
        int p= Collections.binarySearch(set, mapping, getDefaultNMComparator());
        if (p< 0)
            p= -(p+ 1);
        set.add(p, mapping);

        // check for strata overflows
        int q= Arrays.binarySearch(currStrata, nm);
        if (q> 0)
            return 0;   // no change to strata
        q= -(q+ 1);
        assert(q< currStrataLength);    // condition above
        if (currStrataLength== currStrata.length) {
            // loose the highest stratum
            for (int i = (set.size()- 1); i>= 0; --i) {
                 if (set.get(i).getIntegerAttribute(SAMConstants.SAM_OPTION_NM)== currStrata[currStrata.length- 1])
                 set.remove(i);
            }
        }
        if (q< currStrataLength- 1)
            System.arraycopy(currStrata, q, currStrata, q+ 1, currStrataLength- (q+ 1));
        currStrata[q]= nm;
        if (currStrataLength== currStrata.length)
            return -1;
        ++currStrataLength;
        return 1;
    }

}
