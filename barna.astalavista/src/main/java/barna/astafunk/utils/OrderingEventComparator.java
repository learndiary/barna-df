package barna.astafunk.utils;

import barna.model.ASEvent;

import java.util.Comparator;

/**
 * Created by vitor on 6/3/14.
 */

/**
 * This class implements a comparator of AS events regions to sort them.
 */
public class OrderingEventComparator implements Comparator<ASEvent> {

    @Override
    public int compare(ASEvent o1, ASEvent o2) {

        if(o1.getRegionEvent().getStart() < o2.getRegionEvent().getStart()){
            return -1;
        }else if(o1.getRegionEvent().getStart() > o2.getRegionEvent().getStart()){
            return 1;
        }else if(o1.getRegionEvent().getEnd() < o2.getRegionEvent().getEnd()){
            return -1;
        }else if(o1.getRegionEvent().getEnd() > o2.getRegionEvent().getEnd()){
            return 1;
        }
        return 0;
    }
}
