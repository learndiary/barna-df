package fbi.genome.errormodel;

import fbi.commons.Log;
import fbi.commons.StringConstants;
import fbi.commons.io.Serializer;
import fbi.commons.options.HelpPrinter;
import fbi.commons.tools.Qualities;
import fbi.genome.sequencing.rnaseq.simulation.FluxTool;
import org.cyclopsgroup.jcli.ArgumentProcessor;
import org.cyclopsgroup.jcli.annotation.Cli;
import org.cyclopsgroup.jcli.annotation.Option;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
@Cli(name = "model", description = "create markov error model")
public class MarkovErrorModel implements FluxTool{

    private File file;
    private int readLength;
    private File output;
    private File input;
    private int limit = -1;
    private boolean printQualityDistribution;
    private boolean printReadQualityDistribution;

    public File getFile() {
        return file;
    }
    @Option(name = "f", longName = "file",description = ".map input file")
    public void setFile(File file) {
        this.file = file;
    }

    public int getReadLength() {
        return readLength;
    }
    @Option(name = "l", longName="length", description = "read length")
    public void setReadLength(int readLength) {
        this.readLength = readLength;
    }

    public File getOutput() {
        return output;
    }
    @Option(name = "o", longName = "output", description = "output file name")
    public void setOutput(File output) {
        this.output = output;
    }

    public File getInput() {
        return input;
    }
    @Option(name = "i", longName = "input", description = "input file name")
    public void setInput(File input) {
        this.input = input;
    }

    public int getLimit() {
        return limit;
    }
    @Option(name = "s", longName = "limit", description = "read limit number of sequences to create the model")
    public void setLimit(int limit) {
        this.limit = limit;
    }

    public boolean isPrintQualityDistribution() {
        return printQualityDistribution;
    }

    @Option(name = "q", longName = "qualityDistribution", description = "print the quality distribution")
    public void setPrintQualityDistribution(boolean printQualityDistribution) {
        this.printQualityDistribution = printQualityDistribution;
    }

    public boolean isPrintReadQualityDistribution() {
        return printReadQualityDistribution;
    }

    @Option(name = "r", longName = "readDistribution", description = "print the read quality distribution")
    public void setPrintReadQualityDistribution(boolean printReadQualityDistribution) {
        this.printReadQualityDistribution = printReadQualityDistribution;
    }

    public boolean validateParameters(HelpPrinter printer, ArgumentProcessor toolArguments) {
        if(getInput() != null){
            if(!getInput().exists()){
                printer.out.print("Model input file "+ getInput() + " not found!");
                printer.print(toolArguments);
                return false;
            }
            return true;
        }

        if(getFile() == null){
            printer.out.println("No input file specified!\n");
            printer.print(toolArguments);
            return false;
        }else if(!getFile().exists()){
            printer.out.println(getFile().getAbsolutePath() + " does not exist!\n");
            printer.print(toolArguments);
            return false;
        }
        if(readLength == 0){
            printer.out.println("Please specify a read length!\n");
            printer.print(toolArguments);
            return false;
        }
        if(getOutput() == null){
            printer.out.println("Please specify an output file!\n");
            printer.print(toolArguments);
            return false;
        }
        return true;
    }

    public Object call() throws Exception {

        if(getInput() != null){
            Log.message("Reading model from " + getInput().getAbsolutePath());
            FileInputStream in = new FileInputStream(getInput());
            QualityTransitions trans = (QualityTransitions) Serializer.load(in);
            in.close();
            return trans;
        }

        int numStates = Qualities.PHRED_RANGE[1];
        Log.progressStart("Creating Markov Model");
        MapFileReader reader = new MapFileReader(getFile(), Qualities.Technology.Phred);
        QualityTransitions trans = new QualityTransitions(numStates, readLength);

        ReadQualityDistribution readQuals = new ReadQualityDistribution(numStates);
        QualityDistribution qualityDistribution = new QualityDistribution(numStates);

        int c = 0;
        Read read = null;
        while((read = reader.parseNext(false)) != null && (limit < 0 || c < limit )){
            trans.addRead(read);
            readQuals.addRead(read);
            qualityDistribution.addRead(read);
            Log.progress(c++,limit);
        }
        Log.progressFinish(StringConstants.OK, true);

        Log.message("Writing model to " + getOutput().getAbsolutePath());
        FileOutputStream out = new FileOutputStream(getOutput());
        Serializer.save(trans, out);
        out.close();

        if (isPrintQualityDistribution()){
            System.out.println("QUALITY DISTRIBUTION");
            System.out.println(qualityDistribution.toString());
        }
        if(isPrintReadQualityDistribution()){
            System.out.println("READ LENGTH");
            System.out.println(readQuals.toString());
        }

        return null;
    }
}