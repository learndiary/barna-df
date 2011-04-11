package fbi.genome.errormodel;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
abstract class Distribution {
    /**
     * Number of quality values
     */
    protected int size = 0;
    /**
     * The values
     */
    protected int[] values;

    /**
     * Number of reads
     */
    protected int reads;

    public Distribution(int size) {
        this.size = size;
        this.values = new int[size];
    }

    public abstract void addRead(Read read);

    public double getValue(int position){
        return values[position]/(double)reads;
    }

    public String toString(){
        StringBuffer b = new StringBuffer();
        for (int i = 0; i < size; i++) {
            b.append(i).append("\t").append(getValue(i)).append("\n");
        }
        return b.toString();
    }


}
