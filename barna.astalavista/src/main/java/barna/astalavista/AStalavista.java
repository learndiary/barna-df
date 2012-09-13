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

package barna.astalavista;

import barna.commons.Execute;
import barna.commons.cli.jsap.JSAPParameters;
import barna.commons.launcher.FluxTool;
import barna.commons.log.Log;
import barna.io.GeneAheadReaderThread;
import barna.io.gtf.GTFwrapper;
import barna.model.ASEvent;
import barna.model.Graph;
import barna.model.Species;
import barna.model.Transcript;
import barna.model.commons.MyFile;
import barna.model.constants.Constants;
import barna.model.splicegraph.SplicingGraph;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;

import java.io.File;
import java.io.OutputStreamWriter;
import java.lang.reflect.Method;
import java.rmi.ServerError;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Executor;
import java.util.concurrent.Future;


public class AStalavista implements FluxTool<Void>{

    public static void main(String[] args) {

        Execute.initialize(2);
        JSAP jsap = new JSAP();
        AStalavista myAsta= new AStalavista();
        List<Parameter> parameter = myAsta.getParameter();
        if(parameter != null){
            for (Parameter p : parameter) {
                try {
                    jsap.registerParameter(p);
                } catch (JSAPException e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
            }
        }
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
     * Input with the transcriptome annotation in GTF.
     */
    private static File inputFile;
    public static int counter= 0;
    private static int readAheadLimit= -1;
    static boolean onlyInternal= true;
    public static boolean DEBUG= false;
    static long invalidIntrons= 0, totalIntrons= 0;
    static boolean acceptableIntrons= false; // schmu-buh, immer true sobald intronConfidence gesetzt

    @Override
    public Void call() throws Exception {

        // TODO parse arguments
        //inputFile= new MyFile(parseArguments(SplicingGraph.writerThread, args).getAbsolutePath());

        //
        // /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716.gtf
        // /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716.gtf
        // /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_splicedESTs_from_UCSC070716.gtf
        // /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_splicedESTs_from_UCSC070716_chr11.gtf
        // /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_splicedESTs_from_UCSC070716_chr6.gtf
        //
        // /home/ug/msammeth/annotations/mm8_0602_RefSeq_fromUCSC_070807.gtf
        // /home/ug/msammeth/annotations/mm8_0602_RefSeq_fromUCSC_070807_mRNAs_fromUCSC070919.gtf
        // /home/ug/msammeth/annotations/mm8_0602_RefSeq_fromUCSC_070807_mRNAs_fromUCSC070919_splicedESTs_fromUCSC070919.gtf
        boolean output= false, output2= true;
//		if (rusc)
//			outputFname= "delme.asta";

        SplicingGraph.writerThread.start();

        // init and start threads
        long t0= System.currentTimeMillis();
        if (output2) {
            // writerThread
            // outputStats(SplicingGraph.writerThread, new OutputStreamWriter(System.err));
            // TODO write settings
            //Date ti= new Date(t0);
            //System.out.println("["+ti+"]  started, k= "+EventExtractorThread.n+" species "+EventExtractorThread.species+", input file "+inputFile.getAbsolutePath()+", output file= "+outputFname);
        }
        //GTFChrReader reader= new GTFChrReader(file.getAbsolutePath());
        //ChromosomeReaderThread readerThread= new ChromosomeReaderThread(reader);
        GTFwrapper reader= new GTFwrapper(inputFile.getAbsolutePath());
        if (readAheadLimit> 0)
            reader.setReadAheadLimit(readAheadLimit);
        reader.setNoIDs(null);
        //reader.sweepToChromosome("chr17");
        GeneAheadReaderThread readerThread= new GeneAheadReaderThread(reader);
        readerThread.setOutput(output);
        readerThread.setOutput2(output2);
        readerThread.start();
        try {
            readerThread.join();
            readerThread.getDownstreamThread().join();
        } catch (InterruptedException e1) {
            ;	// :)
        }

        System.err.println("took "+((System.currentTimeMillis()- t0)/1000)+" sec.");
        try {
            SplicingGraph.writerThread.setKill(true);
            SplicingGraph.writerThread.interrupt();
            SplicingGraph.writerThread.join();
        } catch (InterruptedException e) {
            // TODO Auto-generated catch block
        }
        System.err.println("found "+counter+" events.");
        if (acceptableIntrons) {
            DecimalFormat df = new DecimalFormat("#.##");
            System.err.println("discarded " + invalidIntrons + " introns, " +
                    "found " + (totalIntrons - invalidIntrons) + " valid ones when checking splice sites: " +
                    "ratio (invalid/total) = " + df.format(((double) invalidIntrons) / totalIntrons));
        }


        return null;
    }

    @Override
    public String getName() {
        return "astalavista";
    }

    @Override
    public String getDescription() {
        return "The AStalavista event retriever";
    }

    @Override
    public List<Parameter> getParameter() {

        // converts parameter file parameters to CLI parameters

        ArrayList<Parameter> parameters = new ArrayList<Parameter>();
        AStalavistaSettings settings= new AStalavistaSettings();
        Collection<barna.commons.parameters.Parameter> pars=
                settings.getParameters().values();
        for (barna.commons.parameters.Parameter parameter : pars) {

            Class c= parameter.getType();
            Parameter p= null;
            if (c.equals(Boolean.class)) {
                p= JSAPParameters.switchParameter(
                        parameter.getLongOption(),
                        parameter.getShortOption())
                        .defaultValue(parameter.getDefault().toString())
                        .type(c)
                        .help(parameter.getDescription())
                        .get();
            } else {
                p= JSAPParameters.flaggedParameter(
                    parameter.getLongOption(),
                    parameter.getShortOption())
                    .type(c)
                    .help(parameter.getDescription())
                    .valueName(parameter.getName())
                    .get();
            }
            // TODO required() not implemented
            if (parameter.getLongOption()!= null|| parameter.getShortOption()!= 0)
                parameters.add(p);
        }

       return parameters;
    }

    @Override
    public boolean validateParameter(JSAPResult args) {

        // output help
        if (args.userSpecified(AStalavistaSettings.PRINT_PARAMETERS.getName())) {
            AStalavistaSettings settings= new AStalavistaSettings();
            settings.write(System.out);
            return false;
        }

        if(args.getFile(AStalavistaSettings.GEN_DIR.getName())!= null)
            acceptableIntrons= true;

        if (args.userSpecified(AStalavistaSettings.EXT_EVENTS.getLongOption()))
            acceptableIntrons= true;
        else
            acceptableIntrons= false;

        // check for necessary components
        /* if (args.userSpecified(AStalavistaSettings.REF_FILE.getShortOption()))
            inputFile= args.getFile(AStalavistaSettings.REF_FILE.getShortOption());
        else */
        if (args.userSpecified(AStalavistaSettings.REF_FILE.getLongOption()))
            inputFile= args.getFile(AStalavistaSettings.REF_FILE.getLongOption());
        if (inputFile== null|| !inputFile.exists()) {
            Log.error("Hey, you forgot to specify a valid input file! ");
            if (inputFile!= null)
                Log.error("Cannot find: "+inputFile.getAbsolutePath());
            return false;
        }

        if ((SplicingGraph.canonicalSS|| acceptableIntrons)&& Graph.overrideSequenceDirPath== null) {
            Log.error("You want me to check introns for valid/canonical splice sites, but you did not provide a valid sequence directory");
            return false;
        }

        SplicingGraph.writerThread= new SplicingGraph.WriterThread();
        if (SplicingGraph.writerThread.outputFname== null&& (!SplicingGraph.writeStdOut)) {
            SplicingGraph.writerThread.outputFname= inputFile.getAbsolutePath()+"_astalavista.gtf.gz";
        }
        if (SplicingGraph.writerThread.outputFname!= null&& new MyFile(SplicingGraph.writerThread.outputFname).exists()) {
            // Confirm o..+"\n by typing \'yes\':"
            System.err.println("Overwriting output file "+ SplicingGraph.writerThread.outputFname+".");
        }

        //setFile(args.getFile(""))

        return true;
    }
}
