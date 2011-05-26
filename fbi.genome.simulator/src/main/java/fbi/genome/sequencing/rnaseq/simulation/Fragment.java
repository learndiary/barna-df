package fbi.genome.sequencing.rnaseq.simulation;

/**
 * Model fragments based on the global id ({@code <chromosome>:<position>@<transcriptID>})
 * start and end coordinates.
 */
public class Fragment  {
    /**
     * The global id
     */
	private CharSequence id;
    /**
     * The start
     */
    private int start;

    /**
     * The end
     */
    private int end;

    /**
     * Create a new Fragment
     *
     * @param id the id
     * @param start the start
     * @param end the end
     */
    public Fragment(final CharSequence id, final int start, final int end) {
        if(id == null) throw new NullPointerException("Fragment with NULL id not permitted");
        if(id.length() == 0) throw new IllegalArgumentException("Fragment with empty ID not permitted");
        if(start > end) throw new IllegalArgumentException("Fragment start > end! " + start + "->" + end);
        this.id = id;
        this.start = start;
        this.end = end;
    }

    /**
     * Get the ID
     *
     * @return id the global fragment id
     */
    public CharSequence getId() {
        return id;
    }

    /**
     * Get the start
     * @return start the start
     */
    public int getStart() {
        return start;
    }

    /**
     * Get the end
     * @return end the end
     */
    public int getEnd() {
        return end;
    }

    /**
     * The length
     *
     * @return length the fragment length
     */
    public int length(){
        return end - start + 1;
    }

    @Override
    public String toString() {
        return start+"\t"+end+"\t"+id;
    }
}
