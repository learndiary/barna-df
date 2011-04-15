package fbi.genome.errormodel;

/**
 * Distribution of a Character to its quality values
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class CharacterQualityDistribution extends Distribution{
    /**
     * If true, we trim the end of reads to cu away all qualities scores == 2
     */
    private static final boolean TRIM_ILLUMINA_BAD_ENDS = false;
    private char character;

    public CharacterQualityDistribution(char character, int size) {
        super(size);
        this.character = Character.toUpperCase(character);
    }

    public void addRead(Read read){
        // trim ends

        int finalPositon =read.getLength()-1;
        int[] q = read.getQualities();
        if(TRIM_ILLUMINA_BAD_ENDS){
            for(;finalPositon > 0;finalPositon--){
                if(q[finalPositon] != 2) break;
            }
        }


        //for (int i = 0; i < q.length; i++) {
        for (int i = 0; i <= finalPositon; i++) {
            char readChar = Character.toUpperCase(read.getSequence().charAt(i));
            if(readChar == character){
                reads++;
                values[q[i]]++;
            }
        }
    }
}
