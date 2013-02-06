package barna.flux.capacitor.profile;

import barna.commons.system.OSChecker;
import barna.commons.utils.StringUtils;
import barna.model.commons.Coverage;

import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Emilio Palumbo
 */
public class CoverageStats {

    /**
     * The coverage profile of systematic biases.
     */
    private Coverage coverage = null;

    public CoverageStats() {
        this(0);
    }

    public CoverageStats(int exonicLength) {
         this.coverage = new Coverage(exonicLength);
    }

    public Coverage getCoverage() {
        return coverage;
    }

    public void setCoverage(Coverage coverage) {
        this.coverage = coverage;
    }

    /**
     * Writes coverage statistics of a transcript to disk.
     *
     * @param geneID       locus identifier
     * @param transcriptID transcript identifier
     * @param cds          flag to indicate whether transcript has an annotated ORF
     * @param length       (processed) length of the transcript
     * @param nrReads      number of Mappings
     */
    public boolean writeCoverageStats(BufferedWriter coverageWriter,String geneID, String transcriptID,
                                   boolean cds, int length, long nrReads
    ) {

        try {
            List<String> line = new ArrayList<String>();
            line.add(geneID);
            line.add(transcriptID);
            line.add(cds ? "CDS" : "NC");
            line.add(Integer.toString(length));
            line.add(Long.toString(nrReads));
            line.add(Float.toString(coverage.getFractionCovered()));
            line.add(Long.toString(coverage.getChiSquare(true, false)));
            line.add(Float.toString((float) coverage.getCV(true, true))+OSChecker.NEW_LINE);

            coverageWriter.write(
                    StringUtils.join("\t", line)
            );
            return true;

        } catch (Exception e) {
            return false;
        }
    }
}
