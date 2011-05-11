package fbi.genome.sequencing.rnaseq.simulation;

import fbi.commons.ByteArrayCharSequence;
import fbi.commons.Log;
import fbi.commons.file.FileHelper;
import fbi.genome.model.constants.Constants;

import java.io.*;
import java.util.Hashtable;
import java.util.Iterator;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class ProfilerFile {
    public static final int PRO_COL_NR_SEQ= 8;
    public static final int PRO_COL_NR_MOL= 4;
    public static final int PRO_COL_NR_FRG= 6;
    public static final String PRO_FILE_SEP= "\t";
    public static final ByteArrayCharSequence PRO_FILE_CDS= new ByteArrayCharSequence("CDS");
    public static final ByteArrayCharSequence PRO_FILE_NC= new ByteArrayCharSequence("NC");
    public static String PRO_FILE_CR= "\n";
    public static String PRO_FILE_TAB= "\t";

    public static boolean appendProfile(FluxSimulatorSettings settings, int colNr, Hashtable<CharSequence,Long> mapFrags) {
        try {
            Log.progressStart("Updating .pro file ");

            long total= 0;
            Iterator<Long> iter= mapFrags.values().iterator();
            while(iter.hasNext())
                total+= iter.next();

            File proFile = settings.get(FluxSimulatorSettings.PRO_FILE);
            BufferedReader buffy= new BufferedReader(new FileReader(proFile));
            File tmpF= new File(System.getProperty(Constants.PROPERTY_TMPDIR)+ File.separator+ proFile.getName()+ Fragmenter.TMP_SFX);
            BufferedWriter wright= new BufferedWriter(new FileWriter(tmpF));
            String[] token = new String[0];
            long bytesRead= 0, bytesTotal= proFile.length();
            int perc= 0, lineCtr= 0;
            String nullStr= Double.toString(0d)+ PRO_FILE_TAB+Long.toString(0);
            for (String s= null; (s= buffy.readLine())!= null;++lineCtr) {

                bytesRead+= s.length()+ PRO_FILE_CR.length();
                Log.progress(bytesRead, bytesTotal);
                token= s.split(PRO_FILE_SEP);
                if (token.length== colNr) {
                    wright.write(s);
                    wright.write(PRO_FILE_SEP);
                } else
                    for (int i = 0; i < colNr; i++) {
                        wright.write(token[i]);
                        wright.write(PRO_FILE_SEP);
                    }
                String id= token[0]+ "@"+ token[1];
                if (mapFrags.containsKey(id)) {
                    long absCnt= mapFrags.get(id);
                    double relFreq= absCnt/ (double) total;
                    if (Double.isNaN(relFreq))
                        System.currentTimeMillis();
                    wright.write(Double.toString(relFreq));
                    wright.write(PRO_FILE_SEP);
                    wright.write(Long.toString(absCnt));
                } else
                    wright.write(nullStr);

                wright.write(PRO_FILE_CR);
                if (lineCtr%1000== 0)
                    wright.flush();
            }
            buffy.close();
            wright.flush();
            wright.close();

            proFile.delete();
            FileHelper.move(tmpF, proFile, null);

            Log.progressFinish("OK", false);
            return true;

        } catch (Exception e) {
            e.printStackTrace();
            return false;
        }
    }
}
