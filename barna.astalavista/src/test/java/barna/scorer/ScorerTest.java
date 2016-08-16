package barna.scorer;

import barna.commons.Execute;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;

/**
 * Created by micha on 15/8/16.
 */
public class ScorerTest {

    @BeforeClass
    public static void initExecuter() {
        Execute.initialize(2);
    }

    @AfterClass
    public static void shutdownExecuter() {
        Execute.shutdown();
    }

    @Test
    public void testChr19() throws Exception {

        File gtf= new File(getClass().getResource("/gencode.v19.chr19.annotation.gtf.gz").getFile());
        File chrDir= new File(getClass().getResource("/./").getFile());
        File vcf= new File(getClass().getResource("/chr19_markers_CEU-only_0.9.vcf").getFile());
        File temp= File.createTempFile(getClass().getName(), "sites");

        ScorerSettings settings= new ScorerSettings();
        settings.set(ScorerSettings.IN_FILE, gtf);
        settings.set(ScorerSettings.VARIANT_FILE, vcf);
        settings.set(ScorerSettings.CHR_SEQ, chrDir);
        settings.set(ScorerSettings.SITES_FILE, temp);

        Scorer.runIt(settings, null);

    }


}
