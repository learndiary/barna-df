package fbi.genome.errormodel;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

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
   private List<Mapping> mappings;



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

    public Mapping addMapping(){
        if (mappings == null) mappings = new ArrayList<Mapping>();
        Mapping mapping = new Mapping();
        mappings.add(mapping);
        return mapping;
    }

    public List<Mapping> getMappings(){
        return mappings;
    }

    public void reset(){
        if(mappings != null) mappings.clear();
        length = 0;
        name = null;
        sequence = null;
    }


    class Mapping{
        private List<Missmatch> missmatches;
        private int quality;
        private String name;

        public void addMissmatch(int position, char genomicCharacter){
            if(missmatches == null) missmatches = new ArrayList<Missmatch>();
            missmatches.add(new Missmatch(position, genomicCharacter, getSequence().charAt(position-1)));
        }

        public List<Missmatch> getMissmatches() {
            return missmatches;
        }

        public List<Character> getMissmatches(int position) {
            if(missmatches == null || missmatches.size() == 0) return null;
            ArrayList<Character> cc = new ArrayList<Character>();
            for (Missmatch missmatch : missmatches) {
                if(missmatch.getPosition() == position) cc.add(missmatch.getGenomicCharacter());
            }
            return cc.size() == 0 ? null : cc;
        }

        public void setQuality(int quality) {
            this.quality = quality;
        }

        public int getQuality() {
            return quality;
        }

        public void setName(String name) {
            this.name = name;
        }

        public String getName() {
            return name;
        }
    }

    class Missmatch{
        private int position;
        private char genomicCharacter;
        private char readCharacter;

        Missmatch(int position, char genomicCharacter, char readCharacter) {
            this.position = position;
            this.genomicCharacter = genomicCharacter;
            this.readCharacter = readCharacter;
        }

        public int getPosition() {
            return position;
        }

        public void setPosition(int position) {
            this.position = position;
        }

        public char getGenomicCharacter() {
            return genomicCharacter;
        }

        public void setGenomicCharacter(char genomicCharacter) {
            this.genomicCharacter = genomicCharacter;
        }

        public char getReadCharacter() {
            return readCharacter;
        }

        public void setReadCharacter(char readCharacter) {
            this.readCharacter = readCharacter;
        }
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
