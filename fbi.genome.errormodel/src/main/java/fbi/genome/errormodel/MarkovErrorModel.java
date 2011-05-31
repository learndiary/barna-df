package fbi.genome.errormodel;

import fbi.commons.Log;
import fbi.commons.StringUtils;
import fbi.commons.flux.FluxTool;
import fbi.commons.io.Serializer;
import fbi.commons.options.HelpPrinter;
import fbi.commons.tools.Qualities;
import fbi.commons.tools.XYLinePlotter;
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
public class MarkovErrorModel implements FluxTool {

    private File file;
    private int readLength;
    private File output;
    private File input;
    private int limit = -1;
    private boolean printQualityDistribution;
    private boolean printReadQualityDistribution;
    private String plotDir;
    private Qualities.Technology technology;

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
    @Option(name = "l", longName="length", description = "read length", required = false)
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

    public String getPlotDir() {
        return plotDir;
    }
    @Option(name = "p", longName = "plots", description = "Directory to plot results")
    public void setPlotDir(String plotDir) {
        this.plotDir = plotDir;
    }
    @Option(name = "", longName = "tech", description = "Technology [phred|solexa|illumina13|illumina18]", required = false)
    public void setTechnology(String technology){
        if(technology != null && technology.length()>0){
            technology = technology.trim();
            if (technology.equalsIgnoreCase("phred")){
                this.technology = Qualities.Technology.Phred;
            }else if(technology.equalsIgnoreCase("solexa")){
                this.technology = Qualities.Technology.Solexa;
            }else if(technology.equalsIgnoreCase("illumina13")){
                this.technology = Qualities.Technology.Illumina13;
            }else if(technology.equalsIgnoreCase("illumina18")){
                this.technology = Qualities.Technology.Illumina18;
            }
        }
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
        if(getOutput() == null){
            printer.out.println("Please specify an output file!\n");
            printer.print(toolArguments);
            return false;
        }

        if(technology == null){
            printer.out.println("No technology specified, please specify the technology used to create the quality scores!\n");
            printer.print(toolArguments);

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



        if(getReadLength() == 0){
            Log.message("No read length specified (-l). Trying to figure the read length");
            MapFileReader reader = new MapFileReader(getFile(), technology);
            int reads = 0;
            // check the first 10 reads
            while(reads++ < 10){
                Read read = reader.parseNext(false);
                if(read == null) break;
                if(readLength == 0){
                    readLength = read.getLength();
                }else{
                    if(readLength != read.getLength()){
                        Log.error("Looks like your reads are of different length ... sorry ... I can not handle that!");
                        return null;
                    }
                }
            }
            reader.close();

            if(readLength == 0){
                Log.error("Unable to figure out the read length");
                return null;
            }
            Log.message("Read length " + readLength);

        }

        int numStates = Qualities.PHRED_RANGE[1];
        switch (technology){
            case Phred:numStates = Qualities.PHRED_RANGE[1];break;
            case Solexa:numStates = Qualities.SOLEXA_RANGE[1];break;
            case Illumina13:numStates = Qualities.ILLUMINA_13_RANGE[1];break;
            case Illumina18:numStates = Qualities.ILLUMINA_18_RANGE[1];break;
        }




        Log.progressStart("Creating Markov Model");
        MapFileReader reader = new MapFileReader(getFile(), technology);
        QualityTransitions trans = new QualityTransitions(numStates, readLength);

        ReadQualityDistribution readQuals = new ReadQualityDistribution(numStates);
        QualityDistribution qualityDistribution = new QualityDistribution(numStates);
        ReadLengthToQualityDistribution lengthDist = new ReadLengthToQualityDistribution(readLength);
        CrossTalkModel crossTalkQuality = new CrossTalkModel(numStates, true);
        CrossTalkModel crossTalkPosition = new CrossTalkModel(readLength, false);

        // charachter quality distributions
        CharacterQualityDistribution dA = new CharacterQualityDistribution('A',numStates);
        CharacterQualityDistribution dC = new CharacterQualityDistribution('C',numStates);
        CharacterQualityDistribution dG = new CharacterQualityDistribution('G',numStates);
        CharacterQualityDistribution dT = new CharacterQualityDistribution('T',numStates);
        CharacterQualityDistribution dN = new CharacterQualityDistribution('N',numStates);

        int c = 0;
        Read read = null;

        while((read = reader.parseNext(false)) != null && (limit < 0 || c < limit )){
            // TEST only seq that mapped to something NOT a chromosome
            if(read.getMappings() == null || read.getMappings().size() == 0){
                continue;
            }

            //for (Read.Mapping mapping : read.getMappings()) {
            //    if(!mapping.getName().startsWith("chr")) continue;
            //}

            trans.addRead(read);
            readQuals.addRead(read);
            qualityDistribution.addRead(read);
            lengthDist.addRead(read);
            crossTalkQuality.addRead(read);
            crossTalkPosition.addRead(read);

            dA.addRead(read);
            dC.addRead(read);
            dG.addRead(read);
            dT.addRead(read);
            dN.addRead(read);


            Log.progress(c++,limit);
        }
        Log.progressFinish(StringUtils.OK, true);

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

        System.out.println("QUALITY PER POSITION");
        System.out.println(lengthDist.toString());



        System.out.println(crossTalkQuality.toString());
        System.out.println(crossTalkQuality.printStateTable());


        if(plotDir != null){
            File dir = new File(plotDir);
            dir.mkdirs();

            /// plot crosstalk quality
            new XYLinePlotter("CrossTalk A", "Quality", "Rate")
            .dataset("A->C", crossTalkQuality.getDistribution('A', 'C'))
            .dataset("A->G", crossTalkQuality.getDistribution('A', 'G'))
            .dataset("A->T", crossTalkQuality.getDistribution('A', 'T'))
            .dataset("A->N", crossTalkQuality.getDistribution('A', 'N'))
            .pdf(new File(dir, "crosstalk-A.pdf"), 500, 500);

            new XYLinePlotter("CrossTalk C", "Quality", "Rate")
            .dataset("C->A", crossTalkQuality.getDistribution('C', 'A'))
            .dataset("C->G", crossTalkQuality.getDistribution('C', 'G'))
            .dataset("C->T", crossTalkQuality.getDistribution('C', 'T'))
            .dataset("C->N", crossTalkQuality.getDistribution('C', 'N'))
            .pdf(new File(dir, "crosstalk-C.pdf"), 500, 500);

            new XYLinePlotter("CrossTalk G", "Quality", "Rate")
            .dataset("G->A", crossTalkQuality.getDistribution('G', 'A'))
            .dataset("G->C", crossTalkQuality.getDistribution('G', 'C'))
            .dataset("G->T", crossTalkQuality.getDistribution('G', 'T'))
            .dataset("G->N", crossTalkQuality.getDistribution('G', 'N'))
            .pdf(new File(dir, "crosstalk-G.pdf"), 500, 500);

            new XYLinePlotter("CrossTalk T", "Quality", "Rate")
            .dataset("T->A", crossTalkQuality.getDistribution('T', 'A'))
            .dataset("T->C", crossTalkQuality.getDistribution('T', 'C'))
            .dataset("T->G", crossTalkQuality.getDistribution('T', 'G'))
            .dataset("T->N", crossTalkQuality.getDistribution('T', 'N'))
            .pdf(new File(dir, "crosstalk-T.pdf"), 500, 500);



            /// plot crosstalk position
            new XYLinePlotter("CrossTalk A", "Position", "Rate")
            .dataset("A->C", crossTalkPosition.getDistribution('A', 'C'))
            .dataset("A->G", crossTalkPosition.getDistribution('A', 'G'))
            .dataset("A->T", crossTalkPosition.getDistribution('A', 'T'))
            .dataset("A->N", crossTalkPosition.getDistribution('A', 'N'))
            .pdf(new File(dir, "crosstalkPosition-A.pdf"), 500, 500);

            new XYLinePlotter("CrossTalk C", "Position", "Rate")
            .dataset("C->A", crossTalkPosition.getDistribution('C', 'A'))
            .dataset("C->G", crossTalkPosition.getDistribution('C', 'G'))
            .dataset("C->T", crossTalkPosition.getDistribution('C', 'T'))
            .dataset("C->N", crossTalkPosition.getDistribution('C', 'N'))
            .pdf(new File(dir, "crosstalkPosition-C.pdf"), 500, 500);

            new XYLinePlotter("CrossTalk G", "Position", "Rate")
            .dataset("G->A", crossTalkPosition.getDistribution('G', 'A'))
            .dataset("G->C", crossTalkPosition.getDistribution('G', 'C'))
            .dataset("G->T", crossTalkPosition.getDistribution('G', 'T'))
            .dataset("G->N", crossTalkPosition.getDistribution('G', 'N'))
            .pdf(new File(dir, "crosstalkPosition-G.pdf"), 500, 500);

            new XYLinePlotter("CrossTalk T", "Position", "Rate")
            .dataset("T->A", crossTalkPosition.getDistribution('T', 'A'))
            .dataset("T->C", crossTalkPosition.getDistribution('T', 'C'))
            .dataset("T->G", crossTalkPosition.getDistribution('T', 'G'))
            .dataset("T->N", crossTalkPosition.getDistribution('T', 'N'))
            .pdf(new File(dir, "crosstalkPosition-T.pdf"), 500, 500);



            /// plot character distributions
            File characters2QualPDF = new File(dir, "characters2Qualities.pdf");
            XYLinePlotter c2qPlotter = new XYLinePlotter("Characters 2 Qualities", "Quality", "% Characters");

            c2qPlotter.dataset("A", dA.getDistribution());
            c2qPlotter.dataset("C", dC.getDistribution());
            c2qPlotter.dataset("G", dG.getDistribution());
            c2qPlotter.dataset("T", dT.getDistribution());
            c2qPlotter.dataset("N", dN.getDistribution());

            c2qPlotter.yBounds(0, 1.0);

            c2qPlotter.pdf(characters2QualPDF, 500, 500);



            // print read length 2 quality
            new XYLinePlotter("Read Length 2 Quality", "Position", "Quality")
                    .dataset("", lengthDist.getDistribution())
                    .pdf(new File(dir, "length2Quality.pdf"),500,500);

            new XYLinePlotter("Quality Distribution", "Rate", "Quality")
                    .dataset("", readQuals.getDistribution())
                    .xBounds(0, 35)
                    .pdf(new File(dir, "averageQuality.pdf"), 500, 500);


        }

        return null;
    }
}