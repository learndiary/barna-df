/*
 * Copyright (c) 2010, Micha Sammeth
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * The names of its contributors may be not used to endorse or promote
 *       products derived from this software without specific prior written
 *       permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL MICHA SAMMETH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package barna.flux.simulator.fragmentation;

/**
 * Model fragments based on the global id ({@code <chromosome>:<position>@<transcriptID>})
 * start and end coordinates.
 */
public class Fragment {
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
     * Number of duplications created during amplification
     */
    private int duplicates;

    /**
     * Create a new Fragment
     *
     * @param id    the id
     * @param start the start
     * @param end   the end
     */
    public Fragment(final CharSequence id, final int start, final int end) {
        if (id == null) {
            throw new NullPointerException("Fragment with NULL id not permitted");
        }
        if (id.length() == 0) {
            throw new IllegalArgumentException("Fragment with empty ID not permitted");
        }
        if (start > end) {
            throw new IllegalArgumentException("Fragment start > end! " + start + "->" + end);
        }
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
     *
     * @return start the start
     */
    public int getStart() {
        return start;
    }

    /**
     * Get the end
     *
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
    public int length() {
        return end - start + 1;
    }

    /**
     * Get the number of duplicates created during amplification
     *
     * @return duplicates the number of duplicates
     */
    public int getDuplicates() {
        return duplicates;
    }

    /**
     * Set the number of duplicates
     *
     * @param duplicates number of duplicates
     */
    public void setDuplicates(final int duplicates) {
        if(duplicates < 0) throw new IllegalArgumentException("Number of Duplicates must be >= 0");
        this.duplicates = duplicates;
    }

    @Override
    public String toString() {
        return start + "\t" + end + "\t" + id + "\t" + duplicates;
    }
}
