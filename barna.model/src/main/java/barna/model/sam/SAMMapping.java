package barna.model.sam;

import barna.model.Mapping;
import net.sf.samtools.*;

import java.util.ArrayList;
import java.util.Comparator;

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
    private int length;
    private int mappingQuality;
    private byte strandFlag;
    private String cigarString;

    public SAMMapping(SAMRecord r) {

        readName = r.getReadName();
        referenceName = r.getHeader().getSequence(r.getReferenceIndex()).getSequenceName();
        alignmentStart = r.getAlignmentStart()-1;
        alignmentEnd = r.getAlignmentEnd();
        length = 0;
        mappingQuality = r.getMappingQuality();
        strandFlag = r.getReadNegativeStrandFlag()?(byte)-1:(byte)1;
        cigarString = r.getCigarString();

        initBlocks(r.getCigar());
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
        if (length == 0) {
            Cigar c = TextCigarCodec.getSingleton().decode(cigarString);
            for (CigarElement e : c.getCigarElements()) {
                if (e.getOperator().equals(CigarOperator.M)||e.getOperator().equals(CigarOperator.D))
                    length+=e.getLength();
            }

        }
        return length;
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

    public String getString() {
        return this.getChromosome()+","+(this.getStrand()>0?"+":"-")+(this.getStart()+1)+","+this.cigarString;
    }

    public static class SAMIdComparator implements Comparator<SAMMapping> {
        public int compare(SAMMapping o1, SAMMapping o2) {
            return o1.getName().compareTo(o2.getName());
        }
    }
}