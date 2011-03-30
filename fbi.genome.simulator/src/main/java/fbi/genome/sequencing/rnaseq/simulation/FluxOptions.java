package fbi.genome.sequencing.rnaseq.simulation;

import fbi.commons.Log;
import fbi.commons.file.FileHelper;
import fbi.commons.options.Options;
import fbi.genome.model.constants.Constants;

import java.io.*;

/**
 * Wrap around command line options
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class FluxOptions {
    /**
     * Show help
     */
    private boolean help;
    /**
     * Expression mode
     */
    private boolean expression;
    /**
     * Library mode
     */
    private boolean library;
    /**
     * Sequence mode
     */
    private boolean sequence;
    /**
     * extract splice junctions
     */
    private boolean extractSpliceJunctions;
    /**
     * Install demo data
     */
    private boolean install;
    /**
     * PAR file
     */
    private File file;
    /**
     * IModel file
     */
    private File fileIM;
    /**
     * Genome directory
     */
    private File fileGenome;
    /**
     * The EFlanks
     */
    private int[] eflanks;
    /**
     * Store the options
     */
    private Options options;

    public FluxOptions() {
        // create options
        options = new Options(this);
        options.addOption("setHelp", "request command line options", null, "help");
        options.addParameter("setFile", "specify parameter file (PAR file)", 'p', "par", "parameter");
        options.addOption("setDoExpr", "simulate expression", 'x', "expr", "express");
        options.addOption("setDoLib", "simulate library construction", 'l', "lib", "library");
        options.addOption("setDoSeq", "simulate sequencing", 's', "seq", "sequence");
        options.addParameter("setDoSJ", "extract splice junctions (GTF file)", 'j', "sj", "junctions");
        options.addParameter("set5flank", "exonic flank 5' of intron (-sj)", null, "5flank");
        options.addParameter("set3flank", "exonic flank 3' of intron (-sj)", null, "3flank");
        options.addParameter("setGenomeDir", "set the path to the directory with genomic sequences (-sj)", 'g', "genome");
        options.addParameter("setImodel", "specify the intron model (-sj)", null, "imodel");
        options.addOption("setDoInstall", "install the demonstration (.par) files", null, "install");
    }

    /**
     * Parse the command line parameters and return true if everything is fine.
     * Return false otherwise
     *
     * @param args the command line args
     * @return valid returns true if all command line arguments are valid
     */
    public boolean parseParameters(String[] args) throws Exception {
        return options.parse(args);
    }

    /**
     * Check if everything is valid and apply the options to the given simulator.
     *
     * Note that this method DOES call System.exit eventually!
     *
     * @param sim the simulator
     */
    public void apply(FluxSimulator sim){
        if(help){
            if (options != null){
                options.printUsage(System.err);
                System.exit(0);
            }
        }

        // do install
        if(install){
            install();
            System.exit(0);
        }


        // check consistency of parameters
        if (extractSpliceJunctions && (expression|| library|| sequence)) {
            Log.error("Cannot mix splice junction extraction with an RNAseq experiment " +
                    "(expression, library construction and sequencing).");
            System.exit(-1);
        }

        // apply the parameters
        sim.eflanks = getEFlanks();
        sim.file = file;
        sim.fileGenome = fileGenome;
        sim.fileIM = fileIM;
        sim.doExpr = expression;
        sim.doLib = library;
        sim.doSeq = sequence;
        sim.doSJ = extractSpliceJunctions;
    }

    private static void install() {
        String wrapper= System.getProperty(Constants.PROPERTY_KEY_WRAPPER_BASE);
        try {
            Log.info("installing..");
            File demoDir= new File(wrapper+ File.separator+ ".."+ File.separator+ "demo");
            if ((!demoDir.exists())|| (!demoDir.isDirectory())) {
                    Log.error("didn't find \'demo\' folder.");
            }
            String[] files= demoDir.list();
            String demoPath= demoDir.getAbsolutePath();
            try {
                demoPath= demoDir.getCanonicalPath();
            } catch (Exception exx) {
                ; // :)
            }
            for (int i = 0; i < files.length; i++) {
                if (!files[i].endsWith(FluxSimulatorSettings.DEF_SFX_PAR))
                    continue;
                File tmpFile= File.createTempFile("simulator", "install");
                File parFile= new File(demoPath+ File.separator+ files[i]);
                if (!FileHelper.copy(parFile, tmpFile)) {
                    Log.error("couldn't copy to "+ System.getProperty("java.io.tmpdir"));
                }
                BufferedReader buffy= new BufferedReader(new FileReader(tmpFile));
                BufferedWriter writer= new BufferedWriter(new FileWriter(parFile));
                int ctr= 0;
                for (String s = null; (s= buffy.readLine())!= null;++ctr) {
                    if (!(s.startsWith(FluxSimulatorSettings.PAR_FRG_FNAME)||
                            s.startsWith(FluxSimulatorSettings.PAR_PRO_FNAME)||
                            s.startsWith(FluxSimulatorSettings.PAR_SEQ_FNAME)||
                            s.startsWith(FluxSimulatorSettings.PAR_REF_FNAME)||
                            s.startsWith(FluxSimulatorSettings.PAR_ERR_FNAME))) {
                        writer.write(s);
                        writer.write(System.getProperty("line.separator"));
                        continue;
                    }
                    String[] ss= s.split("\\s");
                    writer.write(ss[0]);
                    writer.write("\t");
                    int p1= ss[1].lastIndexOf("\\"), p2= ss[1].lastIndexOf('/');
                    int p= Math.max(p1, p2);
                    // IzPack variable substitution, eg ${INSTALL_PATH}${FILE_SEPARATOR}testRun.pro
                    int p3= ss[1].lastIndexOf("}");
                    p= Math.max(p, p3);
                    writer.write(demoPath+ File.separator+ ss[1].substring(p+1));
                    writer.write(System.getProperty("line.separator"));
                }
                buffy.close();
                writer.flush();
                writer.close();
                tmpFile.delete();
                Log.info("\twrote "+ctr+" lines to "+parFile.getName());
            }
            Log.info("[OK]");
        } catch (Exception e) {
            Log.error(e.getMessage());
        }
    }




    /**
     * Set the file
     *
     * @param fName file
     */
    public void setFile(String fName) {
        file= new File(fName);
    }

    /**
     * Set the IModel
     * @param fName model file
     */
    public void setImodel(String fName) {
        fileIM= new File(fName);
    }

    /**
     * Extract junktions
     * @param fName
     */
    public void setDoSJ(String fName) {
        extractSpliceJunctions= true;
        file= new File(fName);
    }

    /**
     * Install
     */
    public void setDoInstall(){
        install= true;
    }

    /**
     * Activate expression mode
     */
    public void setDoExpr(){
        expression= true;
    }

    /**
     * Activate library mode
     */
    public void setDoLib(){
        library = true;
    }

    /**
     * Activate sequence mode
     */
    public void setDoSeq(){
        sequence = true;
    }

    /**
     * Activate help
     */
    public void setHelp(){
        help= true;
    }

    /**
     * Set the genome directory
     *
     * @param path path to the genome directory
     */
    public void setGenomeDir(String path) {
        fbi.genome.model.Graph.overrideSequenceDirPath= path;
        this.fileGenome= new File(path);
    }

    public int[] getEFlanks() {
        if (eflanks == null) {
            eflanks = new int[2];
            eflanks[0]= -1;
            eflanks[1]= -1;
        }
        return eflanks;
    }
    public void set5flank(String length) {
        int x= -1;
        try {
            x= Integer.parseInt(length);
            if (x<= 0)
                throw new RuntimeException("Illegal Value");
        } catch (Exception e) {
            Log.error("Not a valid length for 5' exon flank: " + length);
        }
        getEFlanks()[0]= x;
    }
    public void set3flank(String length) {
        int x= -1;
        try {
            x= Integer.parseInt(length);
            if (x<= 0)
                throw new RuntimeException("Illegal Value");
        } catch (Exception e) {
            Log.error("Not a valid length for 3' exon flank: " + length);
        }
        getEFlanks()[1]= x;
    }

}
