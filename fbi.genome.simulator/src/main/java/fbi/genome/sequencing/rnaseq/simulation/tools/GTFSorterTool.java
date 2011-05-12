package fbi.genome.sequencing.rnaseq.simulation.tools;

import fbi.commons.options.HelpPrinter;
import fbi.genome.io.gff.GFFSorter;
import fbi.genome.sequencing.rnaseq.simulation.FluxTool;
import org.cyclopsgroup.jcli.ArgumentProcessor;
import org.cyclopsgroup.jcli.annotation.Cli;
import org.cyclopsgroup.jcli.annotation.Option;

import java.io.File;

/**
 * Sort GTF Files
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
@Cli(name = "sortGtf", description = "Sort a GTF file")
public class GTFSorterTool implements FluxTool {

    private File intFile;
    private File outFile;

    public File getIntFile() {
        return intFile;
    }
    @Option(name = "i", longName = "input", description = "GTF input file")
    public void setIntFile(final File intFile) {
        this.intFile = intFile;
    }

    public File getOutFile() {
        return outFile;
    }
    @Option(name = "o", longName = "output", description = "GTF output file")
    public void setOutFile(final File outFile) {
        this.outFile = outFile;
    }

    @Override
    public boolean validateParameters(final HelpPrinter printer, final ArgumentProcessor toolArguments) {
        if(getIntFile() == null ){
            printer.out.println("Please specify an input file");
            return false;
        }

        if(getOutFile() == null){
            printer.out.println("Please specify an output file");
            return false;
        }
        return true;
    }

    @Override
    public Object call() throws Exception {
        GFFSorter.sort(intFile, outFile);
        return null;
    }
}
