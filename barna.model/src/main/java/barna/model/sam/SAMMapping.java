package barna.model.sam;

import barna.commons.ByteArrayCharSequence;
import barna.model.Mapping;
import net.sf.samtools.*;

import java.util.ArrayList;

/**
 * @author Emilio Palumbo (emiliopalumbo@gmail.com)
 */
public class SAMMapping implements Mapping{

    SAMRecord mapping = null;
    ArrayList<Integer[]> blocks;
    int currentBlock = 0;

    public SAMMapping(SAMRecord r) {
        mapping = r;
        initBlocks();
    }

    @Override
    public CharSequence getName() {
        return mapping.getReadName();
    }

    @Override
    public CharSequence getChromosome() {
        return mapping.getReferenceName();
    }

    @Override
    public int getStart() {
        return mapping.getAlignmentStart();
    }

    @Override
    public int getEnd() {
        return mapping.getAlignmentEnd();
    }

    @Override
    public int getLength() {
        return mapping.getReadLength();
    }

    @Override
    public int getScore() {
        return mapping.getMappingQuality();
    }

    @Override
    public byte getStrand() {
        if (mapping.getReadNegativeStrandFlag())
            return -1;
        return 1;
    }

    @Override
    public int getBlockCount() {
        if (blocks != null)
            return blocks.size();
        else return -1;
    }

    private void initBlocks() {
        int bStart =getStart(), bLength = 0;
        blocks = new ArrayList<Integer[]>();
        if (blocks.size() == 0) {
            for (CigarElement e : mapping.getCigar().getCigarElements()) {
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
