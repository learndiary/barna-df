package fbi.genome.errormodel;

import java.util.regex.Pattern;

/**
 * Read information. Mutable class so we do not have to create a new object for each read.
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class Read {

   private static Pattern MISMATCH_PATTERN= Pattern.compile("([A,C,G,T,N,a,c,g,t,n]\\d{1,3}+)");
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
/*
    public void parseMissmatches(ByteArrayCharSequence line){
            Matcher m= MISMATCH_PATTERN.matcher(line);
            if (!m.find())
                return;

            mm[0]= base.subSequence(pos[3]+ 1+ m.start(1)+1, pos[3]+ 1+ m.end(1)).parseInt();
            //char c1= base.charAt(pos[0]+ mm[0]);
            mm[2]= base.charAt(pos[0]+ mm[0]);	// subst!
            //char c2= base.charAt(pos[3]+ 1+ m.start(1));
            mm[1]= base.charAt(pos[3]+ 1+ m.start(1));

            // DEBUG
            String s= m.group(1);
            char c1= base.charAt(pos[0]+ mm[0]);
            char c2= base.charAt(pos[0]+ mm[0]);

            return m;
        }

    }
*/
}
