package barna.astafunk.DP;

import java.util.Comparator;

/**
 * @version 1
 * @autor vitorlc
 * @since 22/02/14
 */

/**
 * HitComparator compares its two hits for order.  Returns a negative integer,
 * zero, or a positive integer as the first argument is less than, equal
 * to, or greater than the second.<p>
 */
public class HitComparator implements Comparator<Hit> {

    /**
     * Compare two hits about bit score.
     * @param o1 Hit object
     * @param o2 Hit object
     * @return a negative integer, zero, or a positive integer as the first argument
     * is less than, equal to, or greater than the second, respectively.
     */
    @Override
    public int compare(Hit o1, Hit o2) {
        return o1.getScore() < o2.getScore()?-1:(o1.getScore() > o2.getScore()? 1 : 0);
    }
}
