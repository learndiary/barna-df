package fbi.genome.errormodel;

import fbi.commons.Log;
import fbi.commons.StringConstants;
import fbi.commons.options.HelpPrinter;
import fbi.commons.tools.Qualities;
import fbi.genome.sequencing.rnaseq.simulation.FluxTool;
import org.cyclopsgroup.jcli.ArgumentProcessor;
import org.cyclopsgroup.jcli.annotation.Cli;
import org.cyclopsgroup.jcli.annotation.Option;

import java.io.File;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
@Cli(name = "model", description = "create markov error model")
public class MarkovErrorModel implements FluxTool{

    private File file;
    private int readLength;

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

    public boolean validateParameters(HelpPrinter printer, ArgumentProcessor toolArguments) {
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
        return true;
    }

    public Object call() throws Exception {
        Log.progressStart("Creating Markov Model");
        MapFileReader reader = new MapFileReader(getFile(), Qualities.Technology.Phred);
        int numStates = Qualities.PHRED_RANGE[1];
        QualityTransitions trans = new QualityTransitions(numStates, readLength);
        int c = 0;
        Read read = null;
        int limit = 100000;
        while((read = reader.parseNext()) != null && c < limit ){
            trans.addRead(read);
            Log.progress(c++,limit);
        }
        Log.progressFinish(StringConstants.OK, true);
        return null;
    }
}
