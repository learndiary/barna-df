package barna.astalavista.gfx;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Random;
import java.util.Vector;

import barna.io.gtf.GTFwrapper;
import barna.model.ASEvent;
import barna.model.AbstractRegion;
import barna.model.DirectedRegion;
import barna.model.SpliceSite;
import barna.model.gff.GFFObject;

import java.util.HashMap;

import javax.swing.JFrame;
import javax.swing.JPanel;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 03/3/06
 * Time: 10:28 PM
 * To change this template use File | Settings | File Templates.
 * @author msammeth
 */
public class SpliceOSigner2 extends JPanel {

    static class IntVector {

        static final int DEFAULT_LOAD_FACTOR= 13;
        int[] vector;
        public int length= -1;
        int incrementSize;

        public IntVector() {
            this(DEFAULT_LOAD_FACTOR);
        }

        public IntVector(int initialCapacity) {
            this(initialCapacity, DEFAULT_LOAD_FACTOR);
        }

        public IntVector(int initialCapacity, int loadFactor) {
            vector= new int[initialCapacity];
            incrementSize= loadFactor;
            length= 0;
        }

        public void add(int x) {
            if (length== vector.length)
                extendVector();
            vector[length++]= x;
        }

        @Override
        protected Object clone() throws CloneNotSupportedException {
            IntVector newV= new IntVector();
            newV.vector= this.vector.clone();
            newV.length= this.length;
            newV.incrementSize= this.incrementSize;
            return newV;
        }

        public IntVector cloneIntVector() {
            try {
                return (IntVector) clone();
            } catch (CloneNotSupportedException e) {
                e.printStackTrace();
            }
            return null;
        }

        void extendVector() {
            int[] newVector= new int[vector.length+ incrementSize];
            for (int i = 0; i < vector.length; i++)
                newVector[i]= vector[i];
            vector= newVector;
        }

        public int size() {
            return length;
        }

        public int remove(int pos) {
            int result= vector[pos];
            --length;
            for (int i = pos; i < length; i++)
                vector[pos]= vector[pos+1];
            return result;
        }

        public int get(int pos) {
            return vector[pos];
        }

        public void set(int pos, int val) {
            if (pos< vector.length)
                vector[pos]= val;
        }
        public int[] toIntArray() {
            int[] result= new int[length];
            for (int i = 0; i < length; i++)
                result[i]= vector[i];
            return result;
        }

        public void insert(int val, int p) {

            if (p< 0)
                p= (p+1)* (-1);

            int[] newA= new int[vector.length+ 1];
            for (int i = 0; i < p; i++)
                newA[i]= vector[i];
            newA[p]= val;
            for (int i = p+1; i < newA.length; i++)
                newA[i]= vector[i-1];

            vector= newA;
            length= newA.length;
        }

        public int getValue(int pos) {
            if (pos< length)
                return vector[pos];
            return 0;
        }

        public void putValue(int pos, int val) {
            if (pos>= length) {
                int rounds= (pos+ 1)- size();
                for (int i = 0; i < rounds; i++)
                    add(0);
            }
            vector[pos]= val;
        }
    }

    static final Dimension DIM_INTRON= new Dimension(30, 1),
            DIM_EXON= new Dimension(20, 10),
            DIM_BORDER= new Dimension(5,2),
            DIM_SPACER= new Dimension(0,5);
    static final Color EXON_COLOR= Color.green.darker().darker(), INTRON_COLOR= Color.black;

    static HashMap<String,String> mapRegSymbols;

    static int getRowDelta() {
        return DIM_EXON.height+ DIM_SPACER.height;
    }

    AbstractRegion[] highlightRegions= null;

    ASEvent event;
    HashMap<Integer, Integer> playground= null;
    HashMap<String, Vector<DirectedRegion>> trptDomMap= null;
    int paintLeftFlank= 0, paintRightFlank= 0;
    HashMap<String, Color> colorTable= null;
    byte firstSStype= SpliceSite.TYPE_NOT_INITED, lastSStype= SpliceSite.TYPE_NOT_INITED;
    String fName= null;
    byte overlapType= -1, overlapMode= -1, overlapInter= -1;

    public static void main(String[] args) {

        // /home/ug/msammeth/workspace/G-Phase/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_splicedESTs_from_UCSC070716_events_k.gtf
        // tst.gtf
        GTFwrapper reader= new GTFwrapper("/home/ug/msammeth/workspace/G-Phase/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_splicedESTs_from_UCSC070716_events_k.gtf");
        reader.setGeneWise(false);
        reader.setReadFeatures(new String[] {ASEvent.GTF_FEATURE_ASEVENT});
        reader.setReadGene(false);
        reader.setReadGTF(true);
        reader.setLimitGTFObs(1);
        reader.setSilent(true);

        while (true) {
            try {
                reader.read();
            } catch (Exception e) {
                e.printStackTrace();
            }
            GFFObject[] gtfs= reader.getGtfObj();
            ASEvent event= new ASEvent(gtfs[0]);
//			if (event.getTranscripts().length< 3)
//				continue;
//			if (event.getSrc().getPos()!= 28368274 || event.getSnk().getPos()!= 2147483647)
//				continue;

            // make some regions
            Random r= new Random();
            HashMap<String, Vector<DirectedRegion>> map= new HashMap<String, Vector<DirectedRegion>>();
            for (int i = 0; i < event.getSpliceChains().length; i++) {
                for (int j = 0; j < event.getSpliceChains()[i].length; j++) {
                    System.out.print(event.getSpliceChains()[i][j]);
                }
                if (event.getSpliceChains()[i].length< 2)
                    continue;
                int lastEnd= event.getSpliceChains()[i][0].getPos();
                int k= 2;
                for (int j = 0; j < k; j++) {
                    int first= lastEnd;
                    int last= event.getSpliceChains()[i][event.getSpliceChains()[i].length-1].getPos();

                    int start= first;
                    int end= last;
                    if (first< last) {
                        start= first+ r.nextInt(last-first)+1;
                        int diff= (last-start)/(k-j);
                        end= start+1;
                        if (diff> 0)
                            end+= r.nextInt(diff);
                        if (end> last)
                            --end;
                    }
//					if (j== 0) {
//						start= first+1;
//						end= first+2;
//					} else {
//						start= last-1;
//						end= last;
//					}
                    lastEnd= end;
                    DirectedRegion reg= new DirectedRegion(start, end, event.getGene().getStrand());
                    reg.setChromosome(event.getGene().getChromosome());

                    String id= event.getTranscripts()[i][0].getTranscriptID();
                    if (event.getTranscripts()[i].length> 1)
                        id= event.getTranscripts()[i][r.nextInt(event.getTranscripts()[i].length-1)].getTranscriptID();
                    Vector<DirectedRegion> v= map.get(id);
                    if (v== null) {
                        v= new Vector<DirectedRegion>();
                        map.put(id, v);
                    }
                    v.add(reg);
                    System.out.print("\t"+reg);
                }
                System.out.println(",");
            }

//			Vector<DirectedRegion> v= new Vector<DirectedRegion>();
//			DirectedRegion reg= new DirectedRegion(14924527,14924547,1);
//			reg.setChromosome("chr10");
//			v.add(reg);
//			reg= new DirectedRegion(14924653,14924669,1);
//			reg.setChromosome("chr10");
//			v.add(reg);
//			reg= new DirectedRegion(14924731,14924747,1);
//			reg.setChromosome("chr10");
//			v.add(reg);
//			map.put(event.getTranscripts()[1][0].getTranscriptID(), v);
//			v= new Vector<DirectedRegion>();
//			reg= new DirectedRegion(14922148,14922162,1);
//			reg.setChromosome("chr10");
//			v.add(reg);
//			reg= new DirectedRegion(14925376,14925377,1);
//			reg.setChromosome("chr10");
//			v.add(reg);
//			map.put(event.getTranscripts()[0][0].getTranscriptID(), v);

            if (map.size()== 0)
                continue;


            System.out.println(event.toString());
            SpliceOSigner2 mySplicOSigner= new SpliceOSigner2(event);
            mySplicOSigner.setTrptDomMap(map);
            JFrame myFrame= new JFrame();
            myFrame.getContentPane().add(mySplicOSigner);
            myFrame.pack();
            myFrame.setVisible(true);

            try {
                System.in.read();
            } catch (IOException e) {
                e.printStackTrace();
            }
            myFrame.dispose();
        }
    }

    public SpliceOSigner2(ASEvent event) {
        this.event= event;
        setOpaque(false);
    }

    HashMap<Integer, Integer> getPlayground() {
        if (playground == null) {

            // modify src/snk if there was no flank
            SpliceSite firstSS= event.getFirstVarSS(), lastSS= event.getLastVarSS();
            if (event.getSrc().getPos()== firstSS.getPos())
                event.getSrc().setPos(firstSS.getPos()-2);
            if (event.getSnk().getPos()== lastSS.getPos())
                event.getSnk().setPos(lastSS.getPos()+2);

            // init mapping
            playground= new HashMap<Integer, Integer>();
            SpliceSite[][] schains= event.getSpliceChains();
            int[] scPos= new int[schains.length];
            for (int i = 0; i < scPos.length; i++)
                scPos[i]= 0;

            int posX= DIM_BORDER.width;
            if (isPaintingLeftFlank()) {
                playground.put(new Integer(Integer.MIN_VALUE), posX);	// event.getSrc().getPos()-1
                if (event.getFirstVarSS().isLeftFlank())
                    posX+= DIM_EXON.width/2;
                else
                    posX+= DIM_INTRON.width/2;
                playground.put(new Integer(event.getSrc().getPos()), posX);
                if (event.getFirstVarSS().isLeftFlank())
                    posX+= DIM_INTRON.width;
                else
                    posX+= DIM_EXON.width;
            }

            while (true) {
                IntVector nextI= new IntVector();
                int nextVal= Integer.MAX_VALUE;
                for (int i = 0; i < schains.length; i++) {
                    if (scPos[i]== schains[i].length)
                        continue;
                    if (schains[i][scPos[i]].getPos()< nextVal) {
                        nextI= new IntVector();
                        nextI.add(i);
                        nextVal= schains[i][scPos[i]].getPos();
                    } else if (schains[i][scPos[i]].getPos()== nextVal)
                        nextI.add(i);
                }

                playground.put(new Integer(schains[nextI.get(0)][scPos[nextI.get(0)]].getPos()), new Integer(posX));

                if (schains[nextI.get(0)][scPos[nextI.get(0)]].isLeftFlank())
                    posX+= DIM_EXON.width;
                else
                    posX+= DIM_INTRON.width;

                for (int i = 0; i < nextI.size(); i++)
                    ++scPos[nextI.get(i)];

                int x= 0;
                for (; x < scPos.length; x++)
                    if (scPos[x]< schains[x].length)
                        break;
                if (x== scPos.length)
                    break;
            }

            if (isPaintingRightFlank()) {
                playground.put(new Integer(event.getSnk().getPos()),  new Integer(posX));
                if (event.getLastVarSS().isLeftFlank())
                    posX+= DIM_INTRON.width/2;	// exon already painted until src
                else
                    posX+= DIM_EXON.width/2;
                playground.put(new Integer(Integer.MAX_VALUE), new Integer(posX));	// event.getSnk().getPos()+1
            }

            // init ovlRegion coordinates
            if (getTrptDomMap()!= null)
                overlap();
        }

        return playground;
    }

    void paintBox(Graphics g, int x, int y, int width, int height, Color c, boolean border) {
        y+= (DIM_EXON.height/2)- (height/2);
        if (border) {
            g.setColor(Color.black);
            g.drawRect(x, y, width, height-1);
        } else {
            g.setColor(c);
            g.fillRect(x, y, width, height);
        }
    }

    private Color getColor(String id) {
        if (getColorTable()== null)
            return Color.red;
        Color c= getColorTable().get(id);
        if (c== null) {
            c= ColorMaster.getRandomColor();
            getColorTable().put(id, c);
        }
        return c;
    }

    protected void paintComponent(Graphics g) {

        HashMap<Integer, Integer> pg= getPlayground();
        int w= 0;
        if (pg.get(new Integer(Integer.MAX_VALUE))!= null)
            w= pg.get(new Integer(Integer.MAX_VALUE))-DIM_BORDER.width-1;
        else
            w= pg.get(new Integer(event.getLastVarSS().getPos())); 	// not - 1, paint it !!

        int start= 0;
        if (!event.getFirstVarSS().isTSS())
            start= DIM_BORDER.width+1;
        //if (event.getLastVarSS().isTES())

        g.setClip(start, 0, w,
                (event.getTranscripts().length* getRowDelta())+ (2* DIM_BORDER.height));

        // ** paint exons **/
        paintExonIntronStructure(g, true, false);

        // paint ovl regions
        int posY= DIM_BORDER.height;
        Vector<DirectedRegion> nonredRegV= new Vector<DirectedRegion>();
        HashMap<DirectedRegion, Integer> rowMap= new HashMap<DirectedRegion, Integer>();
        if (getTrptDomMap()!= null) {
            for (int i = 0; i < event.getTranscripts().length; ++i, posY+= getRowDelta()) {
                for (int j = 0; j < event.getTranscripts()[i].length; j++) {
                    Vector<DirectedRegion> v= getTrptDomMap().get(event.getTranscripts()[i][j].getTranscriptID());
                    for (int k = 0; v!= null&& k < v.size(); k++) {
                        int m= 0;
                        for (; m < nonredRegV.size(); m++) {
                            if (nonredRegV.elementAt(m).getStart()== v.elementAt(k).getStart()
                                    && nonredRegV.elementAt(m).getEnd()== v.elementAt(k).getEnd())
                                break;
                        }
                        //if (m== nonredRegV.size()) {	// redundancy filtering
                        nonredRegV.add(v.elementAt(k));
                        rowMap.put(v.elementAt(k), new Integer(posY));
                        //}
                    }
                }
            }


            DirectedRegion[] nonredRegs= new DirectedRegion[nonredRegV.size()];
            for (int i = 0; i < nonredRegs.length; i++)
                nonredRegs[i]= nonredRegV.elementAt(i);
            Arrays.sort(nonredRegs, DirectedRegion.getDefaultSizeComparator());
            for (int i = nonredRegs.length-1; i >= 0; --i) {
                DirectedRegion reg= nonredRegs[i];
                if (pg.get(reg.get5PrimeEdge())== null|| pg.get(reg.get3PrimeEdge())== null)
                    System.currentTimeMillis();
                int x1= pg.get(reg.get5PrimeEdge()).intValue();
                int x2= pg.get(reg.get3PrimeEdge()).intValue();
                String id= reg.getID();
                if (reg.getAttribute(LaVista.GTF_ATTRIBUTE_GROUP_ID)!= null)
                    id= (String) reg.getAttribute(LaVista.GTF_ATTRIBUTE_GROUP_ID);
                paintBox(g, x1, rowMap.get(reg).intValue(), (x2-x1)+1, DIM_EXON.height-2, getColor(id), false);
            }
        }



        // ** paint introns **/
        paintExonIntronStructure(g, false, false);

        // ** paint exon borders **/
        paintExonIntronStructure(g, false, true);

    }

    private void paintExonIntronStructure(Graphics g, boolean paintExons, boolean paintBorders) {
        SpliceSite[][] schains= event.getSpliceChains();
        int posX= 0, posY= 0, leftX= 0;
        SpliceSite firstSS= event.getFirstVarSS(), lastSS= event.getLastVarSS();
        HashMap<Integer, Integer> pg= getPlayground();

        // paint left flank
        if (isPaintingLeftFlank()) {
            posY= DIM_BORDER.height;
            posX= pg.get(new Integer(Integer.MIN_VALUE)).intValue();	// event.getSrc().getPos()-1
            for (int i = 0; i < schains.length; ++i, posY+= getRowDelta()) {
                if (firstSS.isLeftFlank()) {
                    if (schains[i].length== 0|| !schains[i][0].isTSS())
                        paintBox(g,posX, posY, DIM_EXON.width/2, DIM_EXON.height, EXON_COLOR, (!paintExons));
                } else {
                    if ((!paintExons)&& (schains[i].length== 0|| !schains[i][0].isTSS()))
                        paintBox(g,posX, posY, DIM_INTRON.width/2, DIM_INTRON.height, INTRON_COLOR, false);
                }
            }
            if (firstSS.isLeftFlank())
                leftX= posX+ DIM_EXON.width/2;		// boundary of flanking exon/intron
            else
                leftX= posX+ DIM_INTRON.width/2;
        }

        // paint all transcripts
        posY= DIM_BORDER.height;
        for (int i = 0; i < schains.length; ++i, posY+= getRowDelta()) {
            posX= leftX;
            for (int j = 0; j < schains[i].length; j++) {
                int currX= pg.get(new Integer(schains[i][j].getPos())).intValue();
                if (j== 0) {	// paint flank connection
                    if (isPaintingLeftFlank()) {
                        if (SpliceSite.isLeftFlank(getFirstSStype())) {
                            if ((!paintExons)&& (!schains[i][0].isTSS()))
                                paintBox(g,leftX, posY, (currX- leftX), // +leftX-leftX
                                        DIM_INTRON.height, INTRON_COLOR, false);
                        } else {	// paint downstream
                            if (!schains[i][0].isTSS())
                                paintBox(g,leftX, posY, (currX- leftX),
                                        DIM_EXON.height, EXON_COLOR, (!paintExons));
                        }
                    }
                } else {
                    int lastX= pg.get(new Integer(schains[i][j-1].getPos())).intValue();
                    if (schains[i][j].isLeftFlank()) {
                        if (!paintExons)
                            paintBox(g,lastX, posY,currX- lastX,
                                    DIM_INTRON.height, INTRON_COLOR, false);
                    } else {
                        paintBox(g,lastX, posY, currX- lastX,
                                DIM_EXON.height, EXON_COLOR, (!paintExons));
                    }
                }
            }
        }

        // paint right flank
        if (isPaintingRightFlank()) {
            int rightX= pg.get(new Integer(event.getSnk().getPos()));
            posY= DIM_BORDER.height;
            for (int i = 0; i < schains.length; ++i, posY+= getRowDelta()) {
                int fromPos= leftX;
                if (schains[i].length> 0)
                    fromPos= pg.get(schains[i][schains[i].length-1].getPos()).intValue();
                if (lastSS.isLeftFlank()) {
                    if (schains[i].length== 0|| !schains[i][schains[i].length-1].isTES()) {
                        paintBox(g,fromPos, posY,
                                rightX- fromPos, DIM_EXON.height,
                                EXON_COLOR, (!paintExons));
                        if (!paintExons)
                            paintBox(g,rightX, posY, DIM_INTRON.width/2, DIM_INTRON.height, INTRON_COLOR, false);
                    }
                } else {
                    if (schains[i].length== 0|| !schains[i][schains[i].length-1].isTES()) {
                        if (!paintExons)
                            paintBox(g,fromPos, posY,
                                    rightX- fromPos, DIM_INTRON.height,
                                    INTRON_COLOR, false);
                        paintBox(g,rightX, posY, DIM_EXON.width/2, DIM_EXON.height, EXON_COLOR, (!paintExons));
                    }
                }
            }
        }
    }

    boolean isPaintingRightFlank() {
        if (paintRightFlank == 0) {
            // check last SSs
            SpliceSite[][] schains= event.getSpliceChains();
            boolean hasTES= false, onlyTES= true;
            for (int i = 0; i < schains.length; i++)
                if (schains[i].length> 0) {
                    hasTES|= schains[i][schains[i].length-1].isTES();
                    onlyTES&= schains[i][schains[i].length-1].isTES();
                }
            if ((!hasTES)|| (!onlyTES))
                paintRightFlank= 1;
            else
                paintRightFlank= -1;
        }

        return (paintRightFlank== 1);
    }

    boolean isPaintingLeftFlank() {
        if (paintLeftFlank == 0) {
            // check 1st SSs
            boolean hasTSS= false, onlyTSS= true;
            SpliceSite[][] schains= event.getSpliceChains();
            for (int i = 0; i < schains.length; i++)
                if (schains[i].length> 0) {
                    hasTSS|= schains[i][0].isTSS();
                    onlyTSS&= schains[i][0].isTSS();
                }
            if ((!hasTSS)|| (!onlyTSS))
                paintLeftFlank= 1;
            else
                paintLeftFlank= -1;
        }

        return (paintLeftFlank== 1);
    }

    byte getLastSStype() {
        if (lastSStype == SpliceSite.TYPE_NOT_INITED) {
            SpliceSite[][] schains= event.getSpliceChains();
            for (int i = 0; i < schains.length; i++)
                if (schains[i].length> 0) {
                    lastSStype= schains[i][schains[i].length-1].getType();
                    break;
                }
        }

        return lastSStype;
    }


    byte getFirstSStype() {
        if (firstSStype == SpliceSite.TYPE_NOT_INITED) {
            SpliceSite[][] schains= event.getSpliceChains();
            for (int i = 0; i < schains.length; i++)
                if (schains[i].length> 0) {
                    firstSStype= schains[i][0].getType();
                    break;
                }
        }

        return firstSStype;
    }

    /* (non-Javadoc)
      * @@see javax.swing.JComponent#getPreferredSize()
      */
    public Dimension getPreferredSize2() {

        // x
        Iterator<Integer> iter= getPlayground().values().iterator();
        int max= 0;
        while(iter.hasNext()) {
            int val= iter.next().intValue();
            if (val> max)
                max= val;
        }

        if (isPaintingLeftFlank()) {
            if (SpliceSite.isLeftFlank(getFirstSStype()))
                max+= DIM_INTRON.width+ (DIM_EXON.width/2);
            else
                max+= DIM_EXON.width+ (DIM_INTRON.width/2);
        }

        if (isPaintingRightFlank()) {
            if (SpliceSite.isLeftFlank(getLastSStype()))
                max+= DIM_EXON.width+ (DIM_INTRON.width/2);
            else
                max+= DIM_INTRON.width+ (DIM_EXON.width/2);
        }

        max+= 2*DIM_BORDER.width;

        // y
        int may= (2*DIM_BORDER.height)
                + (event.getSpliceChains().length* getRowDelta())
                - DIM_SPACER.height;	// for before/after first/last no spacer

        return new Dimension(max, may);
    }

    /* (non-Javadoc)
      * @@see javax.swing.JComponent#getPreferredSize()
      */
    public Dimension getPreferredSize() {

        int may= (2* DIM_BORDER.height)+ (getRowDelta()* event.getTranscripts().length);
        int max= getPlayground().get(new Integer(event.getLastVarSS().getPos()));
        if (getPlayground().get(new Integer(Integer.MAX_VALUE))!= null)
            max= getPlayground().get(new Integer(Integer.MAX_VALUE));
        max+= (2* DIM_BORDER.width);

//		max= 600;
        return new Dimension(max, may);
    }

    public HashMap<String, Vector<DirectedRegion>> getTrptDomMap() {
        return trptDomMap;
    }

    public void setTrptDomMap(HashMap<String, Vector<DirectedRegion>> trptDomMap)  {

        // compose super-regions
        HashMap<String, Vector<DirectedRegion>> remap= null;
        if (overlapMode== LaVista.OVL_MODUS_ALL) {
            remap= new HashMap<String, Vector<DirectedRegion>>();
            Iterator<String> iter= trptDomMap.keySet().iterator();
            while (iter.hasNext()) {
                String tid= iter.next();
                Vector<DirectedRegion> v= trptDomMap.get(tid);
                HashMap<String, Vector<DirectedRegion>> tmpMap= new HashMap<String, Vector<DirectedRegion>>();
                for (int i = 0; i < v.size(); i++) {
                    Vector<DirectedRegion> vv= tmpMap.get(v.elementAt(i).getAttribute(LaVista.GTF_ATTRIBUTE_GROUP_ID));
                    if (vv== null) {
                        vv= new Vector<DirectedRegion>();
                        tmpMap.put((String) v.elementAt(i).getAttribute(LaVista.GTF_ATTRIBUTE_GROUP_ID), vv);
                    }
                    vv.add(v.remove(i--));	// empty original vector
                }

                Iterator<String> iter2= tmpMap.keySet().iterator();
                while (iter2.hasNext()) {
                    Vector<DirectedRegion> vv= tmpMap.get(iter2.next());
                    if (vv.size()== 1) {	// leave original
                        //vv.elementAt(0).setID((String) vv.elementAt(0).getAttribute(LaVista.GTF_ATTRIBUTE_GROUP_ID));
                        v.add(vv.elementAt(0));
                    } else {	// >1
                        int min= Integer.MAX_VALUE, max= Integer.MIN_VALUE;
                        for (int i = 0; i < vv.size(); i++) {
                            if (vv.elementAt(i).get5PrimeEdge()< min)
                                min= vv.elementAt(i).get5PrimeEdge();
                            if (vv.elementAt(i).get3PrimeEdge()> max)
                                max= vv.elementAt(i).get3PrimeEdge();
                        }
                        DirectedRegion reg= new DirectedRegion(min, max, vv.elementAt(0).getStrand());
                        reg.setChromosome(vv.elementAt(0).getChromosome());
                        reg.setID(vv.elementAt(0).getID());	// ID for remapping
                        reg.addAttribute(LaVista.GTF_ATTRIBUTE_GROUP_ID, vv.elementAt(0).getAttribute(LaVista.GTF_ATTRIBUTE_GROUP_ID));
                        v.add(reg);
                        remap.put(tid+"_"+reg.getAttribute(LaVista.GTF_ATTRIBUTE_GROUP_ID), vv);
                    }
                }
            }
        }



        // get event overlapping domains
        DirectedRegion[] varRegs= new DirectedRegion[] {event.getRegionEvent()};
        if (overlapType== LaVista.OVL_OPTION_EXONIC)
            varRegs= event.getExonicRegions();
        else if (overlapType== LaVista.OVL_OPTION_VARIABLE)
            varRegs= event.getVariableRegions();
        else if (overlapType== LaVista.OVL_OPTION_VAR_EXONIC)
            varRegs= event.getVariableExonicRegions();
        HashMap<String, Vector<DirectedRegion>> map= new HashMap<String, Vector<DirectedRegion>>();
        for (int i = 0; i < event.getTranscripts().length; i++) {
            for (int j = 0; j < event.getTranscripts()[i].length; j++) {
                String id= event.getTranscripts()[i][j].getTranscriptID();
                Vector<DirectedRegion> v= trptDomMap.get(id);
                Vector<DirectedRegion> vv= map.get(id);
                if (vv== null)
                    vv= new Vector<DirectedRegion>();

                HashMap<String, Vector<DirectedRegion>> tmpMap= new HashMap<String, Vector<DirectedRegion>>();
                for (int k = 0; v!= null&& k < v.size(); k++) {
                    int x= 0;
                    for (; x < varRegs.length; x++)
                        if ((getOverlapInter()== LaVista.OVL_INTERSECT_OVERLAP  && v.elementAt(k).overlaps(varRegs[x]))
                                || (getOverlapInter()== LaVista.OVL_INTERSECT_INCLUDE  && varRegs[x].contains(v.elementAt(k))))
                            break;

                    if (x< varRegs.length) {
                        if (remap== null|| remap.get(id+"_"+v.elementAt(k).getAttribute(LaVista.GTF_ATTRIBUTE_GROUP_ID))== null)
                            vv.add((DirectedRegion) v.elementAt(k).clone());	// clone for later modification
                        else {
                            Vector<DirectedRegion> vvv= remap.get(id+"_"+v.elementAt(k).getAttribute(LaVista.GTF_ATTRIBUTE_GROUP_ID));	// remap to origninal regions
                            for (int m = 0; m < vvv.size(); m++)
                                vv.add((DirectedRegion) vvv.elementAt(m).clone());
                        }
                    } else {
                        Vector<DirectedRegion> vx= tmpMap.get(v.elementAt(k).getAttribute(LaVista.GTF_ATTRIBUTE_GROUP_ID));	// collect non-included regions
                        if (vx== null) {
                            vx= new Vector<DirectedRegion>();
                            tmpMap.put((String) v.elementAt(k).getAttribute(LaVista.GTF_ATTRIBUTE_GROUP_ID), vx);
                        }
                        vx.add(v.elementAt(k));
                    }
                }

                // now only if mandatory, throw out all regions that are not completely contained
                if (overlapMode== LaVista.OVL_MODUS_MANDATORY) {	// not single
                    for (int k = 0; k < vv.size(); k++) {
                        Vector<DirectedRegion> vx= tmpMap.remove(vv.elementAt(k).getAttribute(LaVista.GTF_ATTRIBUTE_GROUP_ID));
                        if (vx!= null)
                            vv.remove(i--);
                    }
                }

                // only add non-empty
                if (vv.size()> 0)
                    map.put(id, vv);
            }
        }


        this.trptDomMap = map;
    }

    public String getOvlStructure() {

        return "";
    }

    // cannot overlap regions of length 1 (start== end)
    void overlap() {

        SpliceSite[] su= event.getSpliceUniverseWithFlanks();

        // get var sites (realPos)
        Iterator<Integer> it1= getPlayground().keySet().iterator();
        int[] vPos= new int[getPlayground().keySet().size()];
        int c= 0;
        while (it1.hasNext())
            vPos[c++]= it1.next().intValue();
        Arrays.sort(vPos);

        // get starts/ ends of regions
        Iterator <Vector<DirectedRegion>> iter= getTrptDomMap().values().iterator();
        HashMap<Integer,Integer> startEndMap= new HashMap<Integer, Integer>();
        while (iter.hasNext()) {
            Vector<DirectedRegion> v= iter.next();
            for (int i = 0; i < v.size(); i++) {
                Integer start= new Integer(v.elementAt(i).getStart());
                Integer end= new Integer(v.elementAt(i).getEnd());
                if (start.equals(end)) {	// regions need an extension for hashing different start and end
                    SpliceSite dummy= new SpliceSite(start, SpliceSite.TYPE_NOT_INITED, event.getGene());
                    int p= Arrays.binarySearch(su, dummy, SpliceSite.getDefaultPositionComparator());
                    if (p>= 0) {
                        if (su[p].isRightFlank()) {
                            --start;
                            v.elementAt(i).setStart(start);
                        } else {
                            ++end;
                            v.elementAt(i).setStart(end);
                        }
                    } else {	// arbitrarily extend start or end
                        --start;
                        v.elementAt(i).setStart(start);
                    }
                }
                startEndMap.put(start,start);
                startEndMap.put(end, end);
            }
        }
        int[] startEnds= new int[startEndMap.size()];
        c= 0;
        it1= startEndMap.keySet().iterator();
        while (it1.hasNext())
            startEnds[c++]= it1.next().intValue();
        Arrays.sort(startEnds);

        // overlap
        HashMap<Integer, Integer> pg= getPlayground();
        for (int i = 0; i < vPos.length-1; i++) {
            // segment
            int st= Arrays.binarySearch(startEnds, vPos[i]);
            int nd= Arrays.binarySearch(startEnds, vPos[i+1]);
            int segs= 0;
            if (st< 0)
                st= -(st+1);
            if (nd< 0)
                nd= -(nd+1);
            if (st> nd|| st>= startEnds.length)
                continue;

            segs+= (nd-st)+ 1;

            // map
            int mpStart= pg.get(new Integer(vPos[i])).intValue();
            int mpEnd= pg.get(new Integer(vPos[i+1])).intValue();
            int len= mpEnd- mpStart;	//-1= +1 -2 (border) for exons, +1 for introns, so nothing
            int segLen= len/ segs;
            int rest= len- (segLen* segs);
            int restDiv= 0;
            if (rest> 0) {
                restDiv= (nd- st)/ rest;
                if ((nd- st)% rest!= 0)
                    ++restDiv;
            }
            int offset= 0;
            if (startEnds[st]!= vPos[i])
                offset+= segLen;
            for (int j = st; j < nd; ++j,offset+=segLen) {
                if (j> st&& restDiv> 0&& (((j-st)% restDiv)== 0))
                    ++offset;	// divide rest
                if (pg.get(new Integer(startEnds[j]))== null) { // mod playground only if a mapping is not already contained in there
                    if (i== 0&& j== st&& pg.get(Integer.MIN_VALUE)!= null)
                        pg.put(new Integer(startEnds[j]), pg.get(Integer.MIN_VALUE));	// stretch before src, after snk
                    else if (i== (vPos.length-2)&& j== (nd-1)&& pg.get(Integer.MAX_VALUE)!= null)
                        pg.put(new Integer(startEnds[j]), pg.get(Integer.MAX_VALUE));	// extend to infinity, if crossing sink and if there is infinity coord
                    else
                        pg.put(new Integer(startEnds[j]), new Integer(mpStart+offset));
                }
            }
        }

    }

    public boolean hasOverlappingRegions() {
        return (trptDomMap!= null&& trptDomMap.size()!= 0);
    }

    public String getRegionSymbol(String domName) {
        if (mapRegSymbols == null&& trptDomMap!= null) {	// build up hash
            mapRegSymbols = new HashMap<String,String>();
            char c= 'A';
            Iterator<Vector<DirectedRegion>> iter= trptDomMap.values().iterator();
            while (iter.hasNext()) {
                Vector<DirectedRegion> v= iter.next();
                for (int i = 0; i < v.size(); i++) {
                    String id= v.elementAt(i).getID();
                    String baseID= stripOrderSuffix(id);
                    String symbol= mapRegSymbols.remove(baseID);
                    if (symbol== null)
                        symbol= new Character(c++).toString();
                    if (c>'Z')	// overflow
                        c= 'A';
                    mapRegSymbols.put(id, symbol);
                    mapRegSymbols.put(baseID, symbol);
                }
            }
        }

        return mapRegSymbols.get(domName);
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

    public HashMap<String, Color> getColorTable() {
        return colorTable;
    }

    public void setColorTable(HashMap<String, Color> colorTable) {
        this.colorTable = colorTable;
    }

    public String getFName() {
        return fName;
    }

    public void setFName(String name) {
        fName = name;
    }

    public byte getOverlapType() {
        return overlapType;
    }

    public void setOverlapType(byte overlapType) {
        this.overlapType = overlapType;
    }

    public byte getOverlapMode() {
        return overlapMode;
    }

    public void setOverlapMode(byte overlapMode) {
        this.overlapMode = overlapMode;
    }

    public byte getOverlapInter() {
        return overlapInter;
    }

    public void setOverlapInter(byte overlapInter) {
        this.overlapInter = overlapInter;
    }


}
