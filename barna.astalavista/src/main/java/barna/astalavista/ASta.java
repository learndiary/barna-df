package barna.astalavista;

import barna.commons.Execute;
import barna.commons.cli.jsap.JSAPParameters;
import barna.commons.log.Log;
import barna.model.Gene;
import barna.model.Transcript;
import barna.model.commons.MyFile;
import barna.model.splicegraph.SplicingGraph;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;

import java.text.DecimalFormat;
import java.util.Date;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 2/5/13
 * Time: 3:01 PM
 */
public class ASta extends AStalavista {

    /**
     * Counter for introns with invalid lengths/sequence attributes.
     */
    static long invalidIntrons= 0;

    /**
     * Counter for number of introns analyzed.
     */
    static long totalIntrons= 0;

    private int evBefore;

    @Override
    public String getName() {
        return "asta";
    }

    @Override
    public String getDescription() {
        return "AStalavista event retriever";
    }

    @Override
    public String getLongDescription() {
        return "The AStalavista event retriever decomposes an input annotation systematically into AS events.";
    }

    protected AStalavistaSettings getSettings() {
        if (settings== null)
            return new AStaSettings();
        return settings;    // else
    }

    public static void main(String[] args) {

        Execute.initialize(2);
        ASta myAsta= new ASta();

        // construct to register parameters in JSAP
        List<Parameter> parameter = JSAPParameters.getJSAPParameter(new AStaSettings());
        JSAP jsap = JSAPParameters.registerParameters(parameter);

        // parse
        try{
            JSAPResult toolParameter = jsap.parse(args);
            if (!myAsta.validateParameter(toolParameter)){
                System.exit(-1);
            }
        } catch (Exception e) {
            Log.error("Parameter error : " + e.getMessage(), e);
            e.printStackTrace();
            System.exit(-1);
        }

        Future<Void> captain= Execute.getExecutor().submit(myAsta);
        try {
            captain.get();
        } catch (InterruptedException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (ExecutionException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        Execute.shutdown();
    }




    /**
     * Summarizes the settings before the run
     * @deprecated uses deprecated member variables,
     */
    @Override
    protected void printSettings(AStalavistaSettings settings) {
        super.printSettings(settings);
        if (!settings.get(AStaSettings.EVENTS).isEmpty()) {
            Log.message("# EVENTS #");
            Log.message("# output\t"+ settings.get(AStaSettings.EVENTS_FILE));
            Log.message("# dimension\t" + settings.get(AStaSettings.EVENTS_DIMENSION));
            Log.message("# as_events\t"+ (SplicingGraph.retrieveASEvents? "true": "false"));
            Log.message("# ds_events\t"+ (SplicingGraph.retrieveDSEvents? "true": "false"));
            Log.message("# vs_events\t"+ (SplicingGraph.retrieveVSEvents? "true": "false"));
            Log.message("# internalOnly\t" + SplicingGraph.onlyInternal);
            if (!SplicingGraph.onlyInternal)
                Log.message("# edgeConfidenceLevel " + Transcript.getEdgeConfidenceLevel());
            Log.message("# canonicalSS\t"+ SplicingGraph.canonicalSS);
            Log.message("# acceptableIntrons\t"+ settings.get(AStaSettings.EVENTS_ATR).contains(AStaSettings.EventOptions.IOK));
            if (SplicingGraph.acceptableIntrons)
                Log.message("# intronConfidenceLevel " + SplicingGraph.intronConfidenceLevel);
        }

    }

    private void getEventExctractor() {
        // init event writer thread
        // TODO make real instance, not static
        EventExtractor.writerThread= new WriterThread();
        if (settings.get(AStaSettings.EVENTS_FILE)!= null)
            EventExtractor.writerThread.outputFname= settings.get(AStaSettings.EVENTS_FILE).getAbsolutePath();
        else
            EventExtractor.writerThread.outputFname= settings.get(AStalavistaSettings.IN_FILE).getAbsolutePath()+"_astalavista.gtf.gz";

        if (EventExtractor.writerThread.outputFname!= null&& new MyFile(EventExtractor.writerThread.outputFname).exists()) {
            // Confirm o..+"\n by typing \'yes\':"
            Log.warn("Overwriting output file " + EventExtractor.writerThread.outputFname + ".");
        }
        EventExtractor.writerThread.start();
    }

    @Override
    protected void callBegin() throws Exception {
        super.callBegin();
        getEventExctractor();
    }

    @Override
    protected void callLoop(Gene g) throws Exception {
        if (g.getTranscriptCount() == 1) {
            return;
        }

        // sets types of events to be extracted, k, etc..
        EventExtractor extractor= new EventExtractor(g, (AStaSettings) settings);
        Thread extractorThread= new Thread(extractor);
        extractorThread.start();
        extractorThread.join();
    }

    @Override
    protected void callFinish() throws Exception {
        try {
            EventExtractor.writerThread.setKill(true);
            EventExtractor.writerThread.interrupt();
            EventExtractor.writerThread.join();
        } catch (InterruptedException e) {
            ; // :)
        }

        Log.info("found "+ EventExtractor.counter+" events.");

        if (settings.get(AStaSettings.EVENTS_ATR).contains(AStaSettings.EventOptions.IOK)) {
            DecimalFormat df = new DecimalFormat("#.##");
            System.err.println("discarded " + invalidIntrons + " introns, " +
                    "found " + (totalIntrons - invalidIntrons) + " valid ones when checking splice sites: " +
                    "ratio (invalid/total) = " + df.format(((double) invalidIntrons) / totalIntrons));
        }
        Log.message("AStalavista.");

    }

    @Override
    protected void callBatchBegin(Gene g) {
        evBefore = EventExtractor.counter;
    }

    @Override
    protected void callBatchFinish(Gene g) {
        // prepare for next batch
        int div = (int) (EventExtractor.cumulGC +
                EventExtractor.cumulGF +
                EventExtractor.cumulEV) / 1000;
        int frac = 0;
        if (div > 0)
            frac = (EventExtractor.counter - evBefore) / div;
        else
            frac = (EventExtractor.counter - evBefore);

        Log.debug("[" + new Date(System.currentTimeMillis()) + "] " + g.getChromosome() +
                " graph construct " + (EventExtractor.cumulGC / 1000) +
                //" sec, fuzzy flanks "+(cumulGF/1000) +
                " sec, contraction " + (EventExtractor.cumulGT / 1000) +
                " sec, extract events " + (EventExtractor.cumulEV / 1000) +
                " sec, found " + (EventExtractor.counter - evBefore) +
                " events, " + frac + " ev/sec.");
        System.gc();
        EventExtractor.writerThread.interrupt();

    }

    @Override
    public boolean validateParameter(JSAPResult args) {
        if (!validateParameter(new AStaSettings(), args))
            return false;
        return true;
    }

}
