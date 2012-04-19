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

package barna.flux.simulator;

import barna.commons.ByteArrayCharSequence;
import barna.commons.log.Log;
import barna.io.FileHelper;

import java.io.*;
import java.util.Map;

/**
 * Manage the profiler file
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class ProfilerFile {
    /**
     * Coverage attributes column
     */
    public static final int PRO_COL_NR_COV = 10;
    /**
     * Sequence column
     */
    public static final int PRO_COL_NR_SEQ = 8;
    /**
     * Molecules column
     */
    public static final int PRO_COL_NR_MOL = 4;
    /**
     * Nr column
     */
    public static final int PRO_COL_NR_FRG = 6;

    /**
     * CDS identifier in the profile file
     */
    public static final ByteArrayCharSequence PRO_FILE_CDS = new ByteArrayCharSequence("CDS");
    /**
     * NC identifier in th profile file
     */
    public static final ByteArrayCharSequence PRO_FILE_NC = new ByteArrayCharSequence("NC");
    /**
     * Profile file entry separator
     */
    public static final String PRO_FILE_SEP = "\t";
    /**
     * Profile file newline character
     */
    public static final String PRO_FILE_CR = "\n";

    /**
     * This appends counts, relative and absolute, at the given column. If the column already exists, the columns and
     * everything after the column will be deleted.
     * <p>
     * The map that contains the absolute counts must contain all counts, we use the sum of all counts
     * to compute the relative value. The keys are created from the global ids, that is: {@code <chromosome>:<position>@<transcriptid>}.
     * For example: {@code chr1:9788211-9892649W@NM_177547}
     * </p>
     *
     * @param proFile the profile
     * @param colNr   the column to write to (must be {@literal > 3})
     * @param map  hash mapping profile entry id to a <code>Number</code>
     * @param writeFrac additionally write fraction of sum-of-values for each line
     * @return success true if successfully updated
     */
    public static boolean appendProfile(File proFile, int colNr, Map<CharSequence, Number> map, boolean writeFrac) {
        if (colNr <= 3) {
            throw new IllegalArgumentException("You can not append to column <= 3!");
        }
        // allow empty string for no-pcr
//      if (map.size()== 0)
//          throw new IllegalArgumentException("Empty map");

        BufferedReader buffy = null;
        BufferedWriter wright = null;
        
        //DecimalFormat df = null, df2= null;
        
        // debug
/*        Iterator<Number> iter= map.values().iterator();
        Number n= iter.next();
        // get 3 significant digits for floats
        if (n instanceof Double|| n instanceof Float) {
        	double min= n.doubleValue()== 0? Double.MAX_VALUE: n.doubleValue();
        	while (iter.hasNext()) {
        		double dn= iter.next().doubleValue();
        		if (dn> 0&& dn< min)
        			min= dn;
        	}
        	StringBuilder sb= new StringBuilder("####0.");
        	if (min!= Double.MAX_VALUE) {
	        	double base= 1/ 10;
	        	for (int i= 0; i< 3;base/= 10) {
	        		sb.append("0");
	        		if (base< min)
	        			++i;
	        	}
        	}
        	
        	df= new DecimalFormat(sb.toString());
        }
        // 3 significant digits for fraction
        if (writeFrac) {
        	iter= map.values().iterator();
        	double sum= 0d, min= Double.MAX_VALUE;
        	while (iter.hasNext()) {
        		double dn= iter.next().doubleValue();
        		if (dn< min&& dn> 0)
        			min= dn;
        		sum+= dn;
        	}
        	
        	StringBuilder sb= new StringBuilder("####0.");
        	if (min!= Double.MAX_VALUE) {
	        	double base= 1/ 10;
	        	for (int i= 0; i< 3;base/= 10) {
	        		sb.append("0");
	        		if (base< min)
	        			++i;
	        	}
        	}
        	
        	df2= new DecimalFormat(sb.toString());
        }
*/        
        try {
            Log.progressStart("Updating .pro file ");

            /**
             * Count fragments
             */
            double total = 0;
            if (writeFrac) {
	            for (final Object aLong : map.values()) {
	                total += ((Number) aLong).doubleValue();
	            }
            }

            /*
            Create temp file and init reader/writer
             */
            File tmpF = FileHelper.createTempFile("updatePro", ".tmp");
            buffy = new BufferedReader(new FileReader(proFile));
            wright = new BufferedWriter(new FileWriter(tmpF));

            // line tokens
            String[] token;
            long bytesRead = 0;
            long bytesTotal = proFile.length();
            // limit floats to 3 significant decimals
            String nullStr = Double.toString(0d) + PRO_FILE_SEP + Long.toString(0);

            ByteArrayCharSequence key = new ByteArrayCharSequence(100);
            for (String s = null; (s = buffy.readLine()) != null; ) {
                // progress report
                bytesRead += s.length() + PRO_FILE_CR.length();
                Log.progress(bytesRead, bytesTotal);

                // split
                token = s.split(PRO_FILE_SEP);

                // if we can simply append, just write the full line
                // otherwise write the tokens up to the column we want to append
                if (token.length == colNr) {
                    wright.write(s);
                    wright.write(PRO_FILE_SEP);
                } else {
                    for (int i = 0; i < colNr; i++) {
                        wright.write(token[i]);
                        wright.write(PRO_FILE_SEP);
                    }
                }


                // check if the entry is in the map and we have counts otherwise append
                // the null string
                //String id = token[0] + "@" + token[1];
                key.clear();
                key.append(token[0] + "@" + token[1]);
                if (map.containsKey(key)) {
                    Number absCnt = (Number) map.get(key);
                    if (writeFrac&& total> 0) {
	                    double relFreq = absCnt.doubleValue() / (double) total;
	                    wright.write(Double.toString(relFreq));	// df2.format(relFreq)
	                    wright.write(PRO_FILE_SEP);
                    }
                    if (absCnt instanceof Double|| absCnt instanceof Float)
                    	wright.write(absCnt.toString());	// df.format(absCnt.doubleValue())
                    else
                    	wright.write(absCnt.toString());
                } else {
                    if (writeFrac&& total> 0) {
                        wright.write(nullStr);
	                    wright.write(PRO_FILE_SEP);
                    }
                    wright.write(nullStr);
                }
                // write final newline
                wright.write(PRO_FILE_CR);
            }


            /*
            Check this for issue #58 and make sure we close all streams before moving the files
             */
            buffy.close();
            wright.close();


            if (!proFile.delete()) {
                Log.warn("Unable to remove original .pro file, I try the move anyways");
            }
            if (!FileHelper.move(tmpF, proFile, null)) {
                throw new RuntimeException("Unable to move new .pro file to its location!");
            }
            Log.progressFinish("OK", false);
            return true;

        } catch (Exception e) {
            Log.progressFailed("Error");
            Log.error("Error while updating .pro file: " + e.getMessage(), e);
            return false;
        } finally {
            if (buffy != null) {
                try {
                    buffy.close();
                } catch (IOException ignore) {
                }
            }
            if (wright != null) {
                try {
                    wright.close();
                } catch (IOException ignore) {
                }
            }

        }
    }


    /**
     * Writes initial profile consisting of
     * <p/>
     * <pre>
     *     1. chromosome and global position
     *     2. transcript ID
     *     3. type (CDS or NC)
     *     4. length
     * </pre>
     * <p/>
     * Expression values are not written
     *
     * @param profile the profile
     * @param target  the target file
     * @throws IOException in case of an error
     */
    public static void writeProfile(Profiler profile, File target) throws IOException {
        BufferedWriter writer = null;
        try {
            writer = new BufferedWriter(new FileWriter(target));
            for (int i = 0; i < profile.size(); i++) {
                writer.write(
                        profile.getLociId(i) + ProfilerFile.PRO_FILE_SEP +
                                profile.getId(i) + ProfilerFile.PRO_FILE_SEP +
                                (profile.isCds(i) ? ProfilerFile.PRO_FILE_CDS : ProfilerFile.PRO_FILE_NC) + ProfilerFile.PRO_FILE_SEP +
                                Integer.toString(profile.getLength(i)) +
                                "\n");
            }
            writer.flush();
        } finally {
            if (writer != null) {
                try {
                    writer.close();
                } catch (IOException ignore) {
                }
            }
        }
    }

}
