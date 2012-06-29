package barna.model.sam;

import barna.commons.ByteArrayCharSequence;
import barna.model.Mapping;
import net.sf.samtools.*;

import java.util.ArrayList;

/**
 * @author Emilio Palumbo (emiliopalumbo@gmail.com)
 */
public class SAMMapping implements Mapping{

    private CharSequence readName;
    private CharSequence referenceName;

    ArrayList<Integer[]> blocks;
    int currentBlock = 0;
    private int alignmentStart;
    private int alignmentEnd;
    private int readLength;
    private int mappingQuality;
    private byte strandFlag;

    public SAMMapping(SAMRecord r) {

        readName = r.getReadName();
        referenceName = r.getReferenceName();
        alignmentStart = r.getAlignmentStart();
        alignmentEnd = r.getAlignmentEnd();
        readLength = r.getReadLength();
        mappingQuality = r.getMappingQuality();
        strandFlag = r.getReadNegativeStrandFlag()?(byte)-1:(byte)1;


        initBlocks(r.getCigar());
    }

    @Override
    public CharSequence getName() {
        return readName;
    }

    @Override
    public CharSequence getChromosome() {
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
        int bStart =getStart(), bLength = 0;
        blocks = new ArrayList<Integer[]>();
        if (blocks.size() == 0) {
            for (CigarElement e : c.getCigarElements()) {
                if (!e.getOperator().equals(CigarOperator.N)) {
                    bLength+=e.getLength();
                } else {
                    blocks.add(new Integer[]{bStart,bLength});
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
            return blocks.get(++currentBlock)[0];
        }   catch (Exception e){
            return -1;
        }
    }

    @Override
    public int getNextBlockSize() {
        try {
            return blocks.get(currentBlock)[1];
        }   catch (Exception e){
            return -1;
        }
    }
}
