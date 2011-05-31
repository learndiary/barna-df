package fbi.genome.sequencing.rnaseq.simulation;

import fbi.commons.ByteArrayCharSequence;
import fbi.commons.Log;
import fbi.commons.file.FileHelper;

import java.io.*;
import java.util.Map;

/**
 * Manage the profiler file
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class ProfilerFile {
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
     * @param counts  map from profile entry id to the absolute count
     * @return success true if successfully updated
     */
    public static boolean appendProfile(File proFile, int colNr, Map<ByteArrayCharSequence, Long> counts) {
        if (colNr <= 3) {
            throw new IllegalArgumentException("You can not append to column <= 3!");
        }
        BufferedReader buffy = null;
        BufferedWriter wright = null;

        try {
            Log.progressStart("Updating .pro file ");

            /**
             * Count fragments
             */
            long total = 0;
            for (final Long aLong : counts.values()) {
                total += aLong;
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
                key.reset();
                key.append(token[0] + "@" + token[1]);
                if (total > 0 && counts.containsKey(key)) {
                    long absCnt = counts.get(key);
                    double relFreq = absCnt / (double) total;
                    wright.write(Double.toString(relFreq));
                    wright.write(PRO_FILE_SEP);
                    wright.write(Long.toString(absCnt));
                } else {
                    wright.write(nullStr);
                }
                // write final newline
                wright.write(PRO_FILE_CR);
            }

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
