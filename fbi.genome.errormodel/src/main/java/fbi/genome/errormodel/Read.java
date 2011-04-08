package fbi.genome.errormodel;

/**
 * Read information. Mutable class so we do not have to create a new object for each read.
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class Read {
   private CharSequence name;
   private CharSequence sequence;
   private int length;
   private int[] qualities;


    public CharSequence getName() {
        return name;
    }

    public void setName(CharSequence name) {
        this.name = name;
    }

    public CharSequence getSequence() {
        return sequence;
    }

    public void setSequence(CharSequence sequence) {
        this.sequence = sequence;
        setLength(sequence.length());

        // init the quality array
        if(qualities == null || qualities.length < getLength()) qualities = new int[getLength()];
    }

    public int getLength() {
        return length;
    }

    public void setLength(int length) {
        this.length = length;
    }

    public int[] getQualities() {
        return qualities;
    }

    public void setQualities(int[] qualities) {
        this.qualities = qualities;
    }
}
