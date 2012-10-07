package barna.model;

import barna.commons.utils.ArrayUtils;

import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Vector;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 10/6/12
 * Time: 12:11 PM
 * To change this template use File | Settings | File Templates.
 */
public class ASEventRegions extends ASEvent implements Comparable {

    DirectedRegion[] reg1= null, reg2= null, reg= null;
    String[] regNredID= null;
    HashMap<String,String> mapRegSymbols= null;
    String regCode = null;
    String regColCode= null;
    int[] reg1Starts, reg1Ends, reg2Starts, reg2Ends;
    static StructureComparator compi= null;

    public static boolean isASVariationWithRegions(ASEvent var) {
        if (var instanceof ASEventRegions)
            return true;
        return false;
    }

    public static boolean isNotASVariationWithRegions(ASEvent var) {
        return !isASVariationWithRegions(var);
    }

    public static StructureComparator getDefaultComparator() {
        if (compi == null) {
            compi = new RegionStructureComparator();
        }

        return compi;
    }

    public int getRegionDegree() {
        return reg1.length+ reg2.length;
    }

    public static String stripOrderSuffix(String domName) {
        if (domName== null|| domName.indexOf('-')< 0)
            return domName;
        int p= domName.lastIndexOf('-');
        try {
            Integer.parseInt(domName.substring(p+1));	// only numbers
        } catch (Exception e) {
            return domName;
        }
        return domName.substring(0, p);
    }
    public ASEventRegions(ASEvent var, DirectedRegion[] regs1, DirectedRegion[] regs2) {

//        ssRegionID3UTR= var.ssRegionID3UTR;
//        ssRegionIDCDS= var.ssRegionIDCDS;
//        ssRegionID5UTR= var.ssRegionID5UTR;
//        trans1= var.trans1;
//        trans2= var.trans2;
        trpts= var.trpts;
        //spliceChain1= var.spliceChain1;	// sorted!
        spliceChains= var.spliceChains;	// sorted!
//        degree= var.degree;
//        asEvents= var.asEvents; // events
        stringRep= var.stringRep;
//        anchors= var.anchors;
//        if (var.attributes!= null)
//            attributes= (HashMap) var.attributes.clone();

        overlap(regs1, regs2);
        //toStringRegionCode();	// why?
    }

    public int compareTo(Object o) {
        return getDefaultComparator().compare(this, o);
    }

    public static class RegionStructureComparator extends ASEvent.StructureComparator {

        @Override
        public int compare(Object arg0, Object arg1) {
            // do that anyway to init ovl chains
            ASEventRegions as1= (ASEventRegions) arg0;
            ASEventRegions as2= (ASEventRegions) arg1;

            int[] domStarts11= as1.getReg1Starts();
            int[] domEnds11= as1.getReg1Ends();
            int[] domStarts12= as1.getReg2Starts();
            int[] domEnds12= as1.getReg2Ends();

            int[] domStarts21= as2.getReg1Starts();
            int[] domEnds21= as2.getReg1Ends();
            int[] domStarts22= as2.getReg2Starts();
            int[] domEnds22= as2.getReg2Ends();

            // schain match between ASs nuked, should be handled in Constructor

            int res= super.compare(arg0, arg1);
            if (res!= 0)
                return res;	// structurally different

            // the ones with many domains first
            if (as1.getRegionDegree()< as2.getRegionDegree())
                return -1;
            if (as1.getRegionDegree()> as2.getRegionDegree())
                return 1;



            // compare
            if (domStarts11.length< domStarts21.length)
                return -1;
            if (domStarts11.length> domStarts21.length)
                return 1;
            if (domStarts12.length< domStarts22.length)
                return -1;
            if (domStarts12.length> domStarts22.length)
                return 1;

            for (int i = 0; i < domStarts11.length; i++) {
                if (domStarts11[i]< domStarts21[i])
                    return -1;
                if (domStarts11[i]> domStarts21[i])
                    return 1;
            }
            for (int i = 0; i < domEnds11.length; i++) {
                if (domEnds11[i]< domEnds21[i])
                    return -1;
                if (domEnds11[i]> domEnds21[i])
                    return 1;
            }
            for (int i = 0; i < domStarts12.length; i++) {
                if (domStarts12[i]< domStarts22[i])
                    return -1;
                if (domStarts12[i]> domStarts22[i])
                    return 1;
            }
            for (int i = 0; i < domEnds12.length; i++) {
                if (domEnds12[i]< domEnds22[i])
                    return -1;
                if (domEnds12[i]> domEnds22[i])
                    return 1;
            }

            // TODO check for domain identity !!!

            return 0; 	// equal
        }
    }

    public static class RegionIdentityComparator extends ASEvent.IdentityComparator {

        @Override
        public int compare(Object arg0, Object arg1) {
            // do that anyway to init ovl chains
            ASEventRegions as1= (ASEventRegions) arg0;
            ASEventRegions as2= (ASEventRegions) arg1;

            int[] domStarts11= as1.getReg1Starts();
            int[] domEnds11= as1.getReg1Ends();
            int[] domStarts12= as1.getReg2Starts();
            int[] domEnds12= as1.getReg2Ends();

            int[] domStarts21= as2.getReg1Starts();
            int[] domEnds21= as2.getReg1Ends();
            int[] domStarts22= as2.getReg2Starts();
            int[] domEnds22= as2.getReg2Ends();

            // schain match between ASs nuked, should be handled in Constructor

            int res= super.compare(arg0, arg1);
            if (res!= 0)
                return res;	// structurally different

            // the ones with many domains first
            if (as1.getRegionDegree()< as2.getRegionDegree())
                return -1;
            if (as1.getRegionDegree()> as2.getRegionDegree())
                return 1;



            // compare
            if (domStarts11.length< domStarts21.length)
                return -1;
            if (domStarts11.length> domStarts21.length)
                return 1;
            if (domStarts12.length< domStarts22.length)
                return -1;
            if (domStarts12.length> domStarts22.length)
                return 1;

            for (int i = 0; i < domStarts11.length; i++) {
                if (domStarts11[i]< domStarts21[i])
                    return -1;
                if (domStarts11[i]> domStarts21[i])
                    return 1;
            }
            for (int i = 0; i < domEnds11.length; i++) {
                if (domEnds11[i]< domEnds21[i])
                    return -1;
                if (domEnds11[i]> domEnds21[i])
                    return 1;
            }
            for (int i = 0; i < domStarts12.length; i++) {
                if (domStarts12[i]< domStarts22[i])
                    return -1;
                if (domStarts12[i]> domStarts22[i])
                    return 1;
            }
            for (int i = 0; i < domEnds12.length; i++) {
                if (domEnds12[i]< domEnds22[i])
                    return -1;
                if (domEnds12[i]> domEnds22[i])
                    return 1;
            }

            // TODO check for domain identity !!!

            return 0; 	// equal
        }
    }

    /**
     * input regions must be sorted
     * @@param regs
     * @@return
     */
    protected void overlap(DirectedRegion[] regs1, DirectedRegion[] regs2) {

        //int[] sp1= SpliceSite.getPositions(spliceChain1);
        //int[] sp2= SpliceSite.getPositions(spliceChain2);
        int[] su= SpliceSite.getPositions(getSpliceUniverse(false));    // TODO true/false=?

        DirectedRegion reg= getRegion(REGION_EVENT)[0];    // TODO which region?
        Vector<Integer> domStarts11= new Vector<Integer>(), domEnds11= new Vector<Integer>(),
                domStarts12= new Vector<Integer>(), domEnds12= new Vector<Integer>();
        Vector<DirectedRegion> dom1V= new Vector<DirectedRegion>();
        Vector<DirectedRegion> dom2V= new Vector<DirectedRegion>();
        for (int i = 0; regs1!= null&& i < regs1.length; i++) {
            if (!regs1[i].overlaps(reg))
                continue;
            domStarts11.add(new Integer(ArrayUtils.convertInsertionPoint(
                    Arrays.binarySearch(su, regs1[i].get5PrimeEdge()))));
            domEnds11.add(new Integer(ArrayUtils.convertInsertionPoint(
                    Arrays.binarySearch(su, regs1[i].get3PrimeEdge()))));
            dom1V.add(regs1[i]);
        }
        for (int i = 0; regs2!= null&& i < regs2.length; i++) {
            if (!regs2[i].overlaps(reg))
                continue;
            domStarts12.add(new Integer(ArrayUtils.convertInsertionPoint(
                    Arrays.binarySearch(su, regs2[i].get5PrimeEdge()))));
            domEnds12.add(new Integer(ArrayUtils.convertInsertionPoint(
                    Arrays.binarySearch(su, regs2[i].get3PrimeEdge()))));
            dom2V.add(regs2[i]);
        }

        if (domStarts11!= null)
            reg1Starts= new int[domStarts11.size()];
        else
            reg1Starts= new int[0];
        for (int i = 0; domStarts11!= null&& i < domStarts11.size(); i++)
            reg1Starts[i]= domStarts11.elementAt(i).intValue();
        if (domEnds11!= null)
            reg1Ends= new int[domStarts11.size()];
        else
            reg1Ends= new int[0];
        for (int i = 0; domEnds11!= null&& i < domEnds11.size(); i++)
            reg1Ends[i]= domEnds11.elementAt(i).intValue();
        if (domStarts12!= null)
            reg2Starts= new int[domStarts12.size()];
        else
            reg2Starts= new int[0];
        for (int i = 0; domStarts12!= null&& i < domStarts12.size(); i++)
            reg2Starts[i]= domStarts12.elementAt(i).intValue();
        if (domEnds12!= null)
            reg2Ends= new int[domEnds12.size()];
        else
            reg2Ends= new int[0];
        for (int i = 0; domEnds12!= null&& i < domEnds12.size(); i++)
            reg2Ends[i]= domEnds12.elementAt(i).intValue();

        reg1= (DirectedRegion[]) ArrayUtils.toField(dom1V);
        if (reg1== null)
            reg1= new DirectedRegion[0];
        reg2= (DirectedRegion[]) ArrayUtils.toField(dom2V);
        if (reg2== null)
            reg2= new DirectedRegion[0];

        Vector v= new Vector();
        v.add(reg1Ends); v.add(reg1);
        ArrayUtils.synchroneousSort(reg1Starts, v);
        orderNeighbors(reg1Starts, reg1Ends, reg1);
        v= new Vector();
        v.add(reg2Ends); v.add(reg2);
        ArrayUtils.synchroneousSort(reg2Starts, v);
        orderNeighbors(reg2Starts, reg2Ends, reg2);
    }

    // pfusch
    private void orderNeighbors(int[] start, int[] end, DirectedRegion[] reg) {
        for (int i = 0; i < start.length- 1; i++) {	// pfusch
            if (start[i]== start[i+1]&&
                    end[i]> end[i+1]) {
                int h= start[i];
                start[i]= start[i+1];
                start[i+1]= h;
                h= end[i];
                end[i]= end[i+1];
                end[i+1]= h;
                DirectedRegion r= reg[i];
                reg[i]= reg[i+1];
                reg[i+1]= r;
            }
        }

    }

    /**
     * domain name with or with-out final numbering
     * @@param domName
     * @@return
     */
    public String getRegionSymbol(String domName) {
        if (mapRegSymbols == null) {	// build up hash
            mapRegSymbols = new HashMap<String,String>();
            char c= 'A';
            for (int i = 0; i < reg1.length; i++) {
                String id= reg1[i].getID();
                String baseID= stripOrderSuffix(id);
                String symbol= mapRegSymbols.remove(baseID);
                if (symbol== null)
                    symbol= new Character(c++).toString();
                mapRegSymbols.put(id, symbol);
                mapRegSymbols.put(baseID, symbol);
            }
            for (int i = 0; i < reg2.length; i++) {
                String id= reg2[i].getID();
                String baseID= stripOrderSuffix(id);
                String symbol= mapRegSymbols.remove(baseID);
                if (symbol== null)
                    symbol= new Character(c++).toString();
                mapRegSymbols.put(id, symbol);
                mapRegSymbols.put(baseID, symbol);
            }

        }

        return mapRegSymbols.get(domName);
    }

    private String toStringRegionCode_addRegion(String baseStr, DirectedRegion[] domains, int[] domStarts, int[] domEnds, int scPos) {

        if (domStarts== null|| domStarts.length== 0)
            return baseStr;

        int i= 0, j= 0;
        while(true) {
            String ins1= "", ins2= "";
            for (; i < domStarts.length; i++) {
                if (domStarts[i]== scPos) {
                    ins1= getRegionSymbol(domains[i++].getID())+"[";
                    break;	// has to close
                }
            }
            for (; j < domEnds.length; j++) {
                if (domEnds[j]== scPos) {
                    ins2= getRegionSymbol(domains[j++].getID())+"]";
                    break;	// has to open another one
                }
            }

            if (ins1.length()> 0&& ins2.length()> 0&& ins1.charAt(0)<= ins2.charAt(0))
                baseStr+= ins1+ ins2;
            else
                baseStr+= ins2+ ins1;



            if (i== domStarts.length&& j== domEnds.length)
                break;
        }
        return baseStr;
    }


    public String toString() {
        return toStringRegionColorCode();
    }

    public String toStringRegionCode() {
        if (regCode== null) {
            String c1= "";
            String c2= "";

            int ltt= 1;
            int p1= 0, p2= 0;
            SpliceSite[] su= getSpliceUniverse(false);  // TODO true/false?
            Comparator compi= new SpliceSite.PositionComparator();
            for (int i = 0; i < su.length; i++) {
                c1= toStringRegionCode_addRegion(c1, reg1, reg1Starts, reg1Ends, i);
                if (getSpliceChain(1)!= null) {
                    int p= Arrays.binarySearch(getSpliceChain(1), su[i], compi);
                    if (p>= 0)
                        c1+= getSpliceChain(1)[p1++].isDonor()?Integer.toString(ltt)+ "^":Integer.toString(ltt)+ "-";
                }
                c2= toStringRegionCode_addRegion(c2, reg2, reg2Starts, reg2Ends, i);
                if (getSpliceChain(2)!= null) {
                    int p= Arrays.binarySearch(getSpliceChain(2), su[i], compi);
                    if (p>= 0)
                        c2+= getSpliceChain(2)[p2++].isDonor()?Integer.toString(ltt)+ "^":Integer.toString(ltt)+ "-";
                }
                ltt++;
            }

            if (getSpliceChain(2)== null|| getSpliceChain(2).length== 0)
                c2+= "0";

            if (reg1.length> 0) {
                if (reg1Starts[reg1Starts.length- 1]== su.length)
                    c1+= getRegionSymbol(reg1[reg1.length- 1].getID())+ "[";
                if (reg1Ends[reg1Ends.length- 1]== su.length)
                    c1+= getRegionSymbol(reg1[reg1.length- 1].getID())+ "]";
            }

            if (reg2.length> 0) {
                if (reg2Starts[reg2Starts.length- 1]== su.length)
                    c2+= getRegionSymbol(reg2[reg2.length- 1].getID())+ "[";
                if (reg2Ends[reg2Ends.length- 1]== su.length)
                    c2+= getRegionSymbol(reg2[reg2.length- 1].getID())+ "]";
            }

            regCode= c1+ " , "+ c2;
        }
        return regCode;

    }


    private String toStringRegionColorCode_addRegion(String baseStr, DirectedRegion[] domains, int[] domStarts, int[] domEnds, int scPos) {

        if (domStarts== null|| domStarts.length== 0)
            return baseStr;

        int i= 0, j= 0;
        while(true) {
            String ins1= "", ins2= "";
            for (; i < domStarts.length; i++) {
                if (domStarts[i]== scPos) {
                    // TODO no gfx here
                    //ins1= Integer.toHexString(SpliceOSigner.getDomainColor(domains[i++].getID()).getRGB())+"[";
                    ins1= domains[i++].getID()+ "[";
                    break;	// has to close
                }
            }
            for (; j < domEnds.length; j++) {
                if (domEnds[j]== scPos) {
                    // TODO no gfx here
                    //ins2= Integer.toHexString(SpliceOSigner.getDomainColor(domains[j++].getID()).getRGB())+"]";
                    ins2= domains[j++].getID()+ "[";
                    break;	// has to open another one
                }
            }

            if (ins1.length()> 0&& ins2.length()> 0&& i== j)
                baseStr+= ins1+ ins2;
            else	// j< oldJ
                baseStr+= ins2+ ins1;



            if (i== domStarts.length&& j== domEnds.length)
                break;
        }
        return baseStr;
    }


    public String toStringRegionColorCode() {
        if (regColCode== null) {
            String c1= "";
            String c2= "";

            int ltt= 1;
            int p1= 0, p2= 0;
            SpliceSite[] su= getSpliceUniverse(false);  // TODO true/false?
            Comparator compi= new SpliceSite.PositionComparator();
            for (int i = 0; i < su.length; i++) {
                c1= toStringRegionColorCode_addRegion(c1, reg1, reg1Starts, reg1Ends, i);
                if (getSpliceChain(1)!= null) {
                    int p= Arrays.binarySearch(getSpliceChain(1), su[i], compi);
                    if (p>= 0)
                        c1+= getSpliceChain(1)[p1++].isDonor()?Integer.toString(ltt)+ "^":Integer.toString(ltt)+ "-";
                }
                c2= toStringRegionColorCode_addRegion(c2, reg2, reg2Starts, reg2Ends, i);
                if (getSpliceChain(2)!= null) {
                    int p= Arrays.binarySearch(getSpliceChain(2), su[i], compi);
                    if (p>= 0)
                        c2+= getSpliceChain(2)[p2++].isDonor()?Integer.toString(ltt)+ "^":Integer.toString(ltt)+ "-";
                }
                ltt++;
            }

            if (getSpliceChain(2)== null|| getSpliceChain(2).length== 0)
                c2+= "0";

            // close unfinished ones
            if (reg1.length> 0) {
                if (reg1Starts[reg1Starts.length- 1]== su.length)
                    // TODO
                    // c1+= Integer.toHexString(SpliceOSigner.getDomainColor(reg1[reg1.length- 1].getID()).getRGB())+ "[";
                    c1+= reg1[reg1.length- 1].getID()+ "[";
                if (reg1Ends[reg1Ends.length- 1]== su.length)
                    // TODO
                    //c1+= Integer.toHexString(SpliceOSigner.getDomainColor(reg1[reg1.length- 1].getID()).getRGB())+ "]";
                    c1+= reg1[reg1.length- 1].getID()+ "]";
            }

            if (reg2.length> 0) {
                if (reg2Starts[reg2Starts.length- 1]== su.length)
                    // TODO
                    //c2+= Integer.toHexString(SpliceOSigner.getDomainColor(reg2[reg2.length- 1].getID()).getRGB())+ "[";
                    c2+= reg2[reg2.length- 1].getID()+ "[";
                if (reg2Ends[reg2Ends.length- 1]== su.length)
                    //c2+= Integer.toHexString(SpliceOSigner.getDomainColor(reg2[reg2.length- 1].getID()).getRGB())+ "]";
                    c2+= reg2[reg2.length- 1].getID()+ "]";
            }

            regColCode= c1+ " , "+ c2;
        }
        return regColCode;

    }

    public DirectedRegion[] getReg1() {
        return reg1;
    }

    public DirectedRegion[] getReg2() {
        return reg2;
    }

    public DirectedRegion[] getRegions() {
        if (reg== null) {
            reg= new DirectedRegion[reg1.length+ reg2.length];
            for (int i = 0; i < reg1.length; i++)
                reg[i]= reg1[i];
            for (int i = 0; i < reg2.length; i++)
                reg[i+reg1.length]= reg2[i];
        }
        return reg;
    }

    public String[] getRegionIDsNonRedundant() {
        if (regNredID == null) {
            DirectedRegion[] r= getRegions();
            HashMap mapReg= new HashMap();
            for (int j = 0; j < r.length; j++) {
                String idStrip= ASEventRegions.stripOrderSuffix(r[j].getID());
                if (mapReg.get(idStrip)!= null)
                    continue;
                mapReg.put(idStrip, idStrip);
            }
            regNredID= (String[]) ArrayUtils.toField(mapReg.keySet());
        }

        return regNredID;

    }

    public String[] getRegionIDs() {
        String[] reg= new String[reg1.length+ reg2.length];
        for (int i = 0; i < reg1.length; i++)
            reg[i]= reg1[i].getID();
        for (int i = 0; i < reg2.length; i++)
            reg[i+reg1.length]= reg2[i].getID();
        return reg;
    }
    public int[] getReg1Ends() {
        return reg1Ends;
    }
    public int[] getReg1Starts() {
        return reg1Starts;
    }
    public int[] getReg2Ends() {
        return reg2Ends;
    }
    public int[] getReg2Starts() {
        return reg2Starts;
    }

}