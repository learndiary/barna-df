package barna.model.sam;

import barna.model.Mapping;
import net.sf.samtools.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;

/**
 * @author Emilio Palumbo (emiliopalumbo@gmail.com)
 */
public class SAMMapping implements Mapping{

    public static String SAM_OPTION_NH= "NH";
    public static String SAM_OPTION_XT= "XT";

    private String readName;
    private String referenceName;

    ArrayList<Integer[]> blocks;
    int currentBlock = 0;
    private int alignmentStart;
    private int alignmentEnd;
    private int mateAlignmentStart;
    private String mateReferenceName;
    private int length;
    private int mappingQuality;
    private byte strandFlag;
    private byte[] sequence;
    private Cigar cigar;

    public int getHits() {
        return hits;
    }

    private int hits;
    private boolean primary;
    private int insertSize;
    private boolean paired;
    private boolean properlyPaired;
    private Character xt;

    public SAMMapping(SAMRecord r) {

        readName = r.getReadName();
        referenceName = r.getHeader().getSequence(r.getReferenceIndex()).getSequenceName();
        alignmentStart = r.getAlignmentStart()-1;
        alignmentEnd = r.getAlignmentEnd();
        mateAlignmentStart = r.getMateAlignmentStart()-1;
        mateReferenceName = r.getMateReferenceName();
        insertSize = r.getInferredInsertSize();
        length = 0;
        mappingQuality = r.getMappingQuality();
        strandFlag = r.getReadNegativeStrandFlag()?(byte)-1:(byte)1;
        cigar = TextCigarCodec.getSingleton().decode(r.getCigarString());
        sequence = r.getReadBases();
        hits = r.getIntegerAttribute(SAM_OPTION_NH)!=null ? r.getIntegerAttribute(SAM_OPTION_NH) : -1;
        xt = r.getCharacterAttribute(SAM_OPTION_XT);
        primary = !r.getNotPrimaryAlignmentFlag();
        paired = r.getReadPairedFlag();
        properlyPaired = paired ? r.getProperPairFlag() : false;
        initBlocks();
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

    public int getMateStart() {
        return mateAlignmentStart;
    }

    public String getMateChromosome() {
        if (mateReferenceName.equals("="))
                return referenceName;
        return mateReferenceName;
    }

    public boolean isPrimary() {
        return primary;
    }

    public int getInsertSize() {
        return insertSize;
    }

    public boolean isProperlyPaired() {
        return properlyPaired;
    }

    public boolean isPaired() {
        return paired;
    }

    public boolean isUnique() {
        if (xt != null)
            return xt == 'U';
        return hits == 1;
    }

    @Override
    public int getLength() {
        if (length == 0) {
            for (CigarElement e : cigar.getCigarElements()) {
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

    private void initBlocks() {
        int bStart =0, bLength = 0;
        blocks = new ArrayList<Integer[]>();
        if (blocks.size() == 0) {
            for (CigarElement e : cigar.getCigarElements()) {
                if (!e.getOperator().equals(CigarOperator.N)) {
                    if (e.getOperator().equals(CigarOperator.M)||e.getOperator().equals(CigarOperator.D)) {
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
            if (currentBlock == blocks.size())
                currentBlock = 0;
            return blocks.get(currentBlock)[0];
        }   catch (Exception e){
            return -1;
        }
    }

    @Override
    public int getNextBlockSize() {
        try {
            if (currentBlock == blocks.size())
                currentBlock = 0;
            return blocks.get(currentBlock++)[1];
        }   catch (Exception e){
            return -1;
        }
    }

    @Override
    public CharSequence getSequence() {
        return new String(sequence);
    }

    @Override
    public CharSequence getCigar() {
        return cigar.toString();
    }

    @Override
    public double getCount(boolean weighted) {
        return (weighted && this.hits > 0 ? 1.0/(double)this.hits : 1.0);
    }

    public String getString() {
        return this.getChromosome()+","+(this.getStrand()>0?"+":"-")+(this.getStart()+1)+","+this.cigar.toString();
    }

    public boolean isMateOf(SAMMapping mapping) {
        if (!this.getChromosome().equals(mapping.getMateChromosome()))
            return false;
        if (!this.getMateChromosome().equals(mapping.getChromosome()))
            return false;
        if (this.getStart() != mapping.getMateStart())
            return false;
        if (this.getMateStart() != mapping.getStart())
            return false;
        if (this.getInsertSize() != -mapping.getInsertSize())
            return false;

        return true;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        SAMMapping that = (SAMMapping) o;

        if (alignmentEnd != that.alignmentEnd) return false;
        if (alignmentStart != that.alignmentStart) return false;
        if (hits != that.hits) return false;
        if (length != that.length) return false;
        if (mappingQuality != that.mappingQuality) return false;
        if (mateAlignmentStart != that.mateAlignmentStart) return false;
        if (strandFlag != that.strandFlag) return false;
        if (!cigar.equals(that.cigar)) return false;
        if (!mateReferenceName.equals(that.mateReferenceName)) return false;
        if (!readName.equals(that.readName)) return false;
        if (!referenceName.equals(that.referenceName)) return false;
        if (!Arrays.equals(sequence, that.sequence)) return false;

        return true;
    }

    public static class SAMIdComparator implements Comparator<SAMMapping> {
        public int compare(SAMMapping o1, SAMMapping o2) {
            return o1.getName().compareTo(o2.getName());
        }
    }
}
