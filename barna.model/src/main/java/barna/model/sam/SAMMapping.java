package barna.model.sam;

import barna.commons.ByteArrayCharSequence;
import barna.model.Mapping;
//import com.sun.xml.internal.fastinfoset.algorithm.BooleanEncodingAlgorithm;
import net.sf.samtools.*;

import java.util.*;

/**
 * @author Emilio Palumbo (emiliopalumbo@gmail.com)
 */
public class SAMMapping implements Mapping{

    private String readName;
    private String referenceName;

    ArrayList<Integer[]> blocks;
    int currentBlock = 0;
    private int alignmentStart;
    private int alignmentEnd;
    private int readLength;
    private int mappingQuality;
    private byte strandFlag;
    private String cigarString;
    private HashMap<String,String> alternates;

    public SAMMapping() {
    }

    public SAMMapping(SAMRecord r) {

        readName = r.getReadName();
        referenceName = r.getHeader().getSequence(r.getReferenceIndex()).getSequenceName();
        alignmentStart = r.getAlignmentStart()-1;
        alignmentEnd = r.getAlignmentEnd();
        readLength = r.getReadLength();
        mappingQuality = r.getMappingQuality();
        strandFlag = r.getReadNegativeStrandFlag()?(byte)-1:(byte)1;
        cigarString = r.getCigarString();

        initBlocks(r.getCigar());

        initAlternates(r.getAttributes());
    }

    public SAMMapping(SAMRecord r, String suffix) {

        this(r);
        this.readName+=suffix;
    }

    @Override
    public String getName() {
        return readName;
    }

    @Override
    public String getChromosome() {
        return referenceName;
    }

    @Override
    public int getStart() {
        return alignmentStart;
    }

    @Override
    public int getEnd() {
        return alignmentEnd;
    }

    @Override
    public int getLength() {
        return readLength;
    }

    @Override
    public int getScore() {
        return mappingQuality;
    }

    @Override
    public byte getStrand() {
        return strandFlag;
    }

    @Override
    public int getBlockCount() {
        if (blocks != null)
            return blocks.size();
        else return -1;
    }

    private void initBlocks(Cigar c) {
        int bStart =0, bLength = 0;
        blocks = new ArrayList<Integer[]>();
        if (blocks.size() == 0) {
//            for (CigarElement e : c.getCigarElements()) {
//                if (!e.getOperator().equals(CigarOperator.N)) {
//                    bLength+=e.getLength();
//                } else {
//                    blocks.add(new Integer[]{bStart,bLength});
//                    bStart+=bLength;//+e.getLength();
//                    bLength = 0;
//                }
//            }
//            blocks.add(new Integer[]{bStart,bLength});
//            bStart+=bLength;

            for (CigarElement e : c.getCigarElements()) {
                if (!e.getOperator().equals(CigarOperator.N)) {
                    if (e.getOperator().equals(CigarOperator.M)) {
                        bLength+=e.getLength();
                    }
                } else {
                    blocks.add(new Integer[]{bStart, bLength});
                    bStart+=bLength+e.getLength();
                    bLength = 0;
                }
            }
            blocks.add(new Integer[]{bStart,bLength});
            bStart+=bLength;

//            int start = blocksList.get(0).getReferenceStart();
//
//            for(AlignmentBlock block : blocksList) {
//                blocks.add(new Integer[]{block.getReferenceStart()-start,block.getLength()});
//            }
        }
    }

    private void initAlternates(List<SAMRecord.SAMTagAndValue> attributes) {
        for (SAMRecord.SAMTagAndValue attr : attributes) {
            if (attr.tag.equals("XA")) {
                if (alternates==null)
                    alternates=new HashMap<String,String>();
                String[] alts = attr.value.toString().split(";");
                for (String s : alts) {
                    String[] fields=s.split(",");
                    String key = fields[0]+","+fields[1]+","+fields[2];
                    if (this.getString().equals(key))
                        continue;
                    if (!alternates.containsKey(key)) {
                        alternates.put(key,fields[3]);
                    }
                }
            }
        }
    }

    @Override
    public int getNextBlockStart() {
        try {
            return blocks.get(currentBlock)[0];
        }   catch (Exception e){
            return -1;
        }
    }

    @Override
    public int getNextBlockSize() {
        try {
            return blocks.get(currentBlock++)[1];
        }   catch (Exception e){
            return -1;
        }
    }

    public List<SAMMapping> getAlternates() {
        if (alternates==null)
            return null;

        ArrayList<SAMMapping> list = new ArrayList<SAMMapping>();

        for (String alt : alternates.keySet()) {
            String[] s = alt.split(",");
            SAMMapping mapping = new SAMMapping();

            Cigar c = TextCigarCodec.getSingleton().decode(s[2]);

            mapping.readName = readName;
            mapping.referenceName = s[0];
            mapping.alignmentStart = Integer.parseInt(s[1].substring(1))-1;
            mapping.alignmentEnd = mapping.alignmentStart+c.getReferenceLength();
            mapping.readLength = c.getReferenceLength();
            mapping.mappingQuality = Integer.parseInt(alternates.get(alt));//Integer.parseInt(s[3]);
            mapping.strandFlag = s[1].substring(0,1).equals("+")?(byte)1:-1;

            mapping.initBlocks(c);

            list.add(mapping);
        }

        return list;
    }

    public boolean hasAlternates() {
        return (alternates!=null);
    }

    public String getString() {
        return this.getChromosome()+","+(this.getStrand()>0?"+":"-")+(this.getStart()+1)+","+this.cigarString;
    }

    public static class SAMIdComparator implements Comparator<SAMMapping> {
        public int compare(SAMMapping o1, SAMMapping o2) {
            return o1.getName().compareTo(o2.getName());
        }
    }

    @Override
    public boolean equals(Mapping otherMapping) {
        SAMMapping sam = null;
        try {
            sam = (SAMMapping)otherMapping;
            if (!this.referenceName.equals(sam.getChromosome()))
                return false;
            if (!this.readName.equals(sam.getName()))
                return false;
            if (this.alignmentStart!=sam.getStart())
                return false;
            if (this.alignmentEnd!=sam.getEnd())
                return false;
            if (this.hasAlternates()) {
                for (String s : alternates.keySet()) {
                    if (!sam.alternates.containsKey(s))
                        return false;
                }
            }
            return true;
        } catch (Exception e) {
            return false;
        }
    }
}
