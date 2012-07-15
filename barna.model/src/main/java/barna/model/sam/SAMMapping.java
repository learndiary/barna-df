package barna.model.sam;

import barna.commons.ByteArrayCharSequence;
import barna.model.Mapping;
//import com.sun.xml.internal.fastinfoset.algorithm.BooleanEncodingAlgorithm;
import net.sf.samtools.*;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

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
    private List<String> alternates;

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

        initBlocks(r.getCigar());

        initAlternates(r.getAttributes());
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
                    alternates=new ArrayList<String>();
                String[] alts = attr.value.toString().split(";");
                for (String s : alts) {
                    alternates.add(s);
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

        for (String alt : alternates) {
            String[] s = alt.split(",");
            SAMMapping mapping = new SAMMapping();

            Cigar c = TextCigarCodec.getSingleton().decode(s[2]);

            mapping.readName = readName;
            mapping.referenceName = s[0];
            mapping.alignmentStart = Integer.parseInt(s[1].substring(1))-1;
            mapping.alignmentEnd = mapping.alignmentStart+c.getReferenceLength();
            mapping.readLength = c.getReferenceLength();
            mapping.mappingQuality = Integer.parseInt(s[3]);
            mapping.strandFlag = s[1].substring(0,1).equals("+")?(byte)1:-1;

            mapping.initBlocks(c);

            list.add(mapping);
        }

        return list;
    }

    public Boolean hasAlternates() {
        return (alternates!=null);
    }

    public static class SAMIdComparator implements Comparator<SAMMapping> {
        public int compare(SAMMapping o1, SAMMapping o2) {
            return o1.getName().compareTo(o2.getName());
        }
    }
}
