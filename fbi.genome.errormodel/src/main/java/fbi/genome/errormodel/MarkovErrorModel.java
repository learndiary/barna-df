package fbi.genome.errormodel;

import fbi.commons.Log;
import fbi.commons.options.HelpPrinter;
import fbi.commons.tools.Qualities;
import fbi.genome.sequencing.rnaseq.simulation.FluxTool;
import org.cyclopsgroup.jcli.ArgumentProcessor;
import org.cyclopsgroup.jcli.annotation.Cli;
import org.cyclopsgroup.jcli.annotation.Option;

import java.io.File;
import java.util.Random;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
@Cli(name = "model", description = "create markov error model")
public class MarkovErrorModel implements FluxTool{

    private File file;

    public File getFile() {
        return file;
    }
    @Option(name = "f", longName = "file",description = ".map input file")
    public void setFile(File file) {
        this.file = file;
    }

    public boolean validateParameters(HelpPrinter printer, ArgumentProcessor toolArguments) {
        if(getFile() == null){
            printer.out.println("No input file specified!");
            printer.print(toolArguments);
            return false;
        }else if(!getFile().exists()){
            printer.out.println(getFile().getAbsolutePath() + " does not exist!");
            printer.print(toolArguments);
            return false;
        }
        return true;
    }

    public Object call() throws Exception {
        Log.info("Create Markov Model");
        MapFileReader reader = new MapFileReader(getFile(), Qualities.Technology.Phred);
        int numStates = Qualities.PHRED_RANGE[1];
        QualityTransitions trans = new QualityTransitions(numStates);
        int c = 0;
        Read read = null;
        Random random = new Random();

        while((read = reader.parseNext()) != null && c < 100000 ){

            trans.addRead(read);

        }
        return null;
    }
}
