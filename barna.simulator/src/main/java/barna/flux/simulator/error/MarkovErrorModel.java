/*
 * This file is part of the Flux Library.
 *
 * The code of the Flux Library may be freely distributed and modified under the terms of the
 * European Union Public Licence (EUPL) published on the web site <http://www.osor.eu/eupl/european-union-public-licence-eupl-v.1.1>.
 * Copyright for the code is held jointly by the individual authors, who should be listed
 * in @author doc comments. According to Article 5 and Article 11 of the EUPL, publications that
 * include results produced by the Flux Library are liable to reference the Work,
 * see the Flux Library homepage <http://flux.sammeth.net> for more information.
 */

package barna.flux.simulator.error;

import barna.commons.cli.jsap.JSAPParameters;
import barna.commons.launcher.FluxTool;
import barna.commons.log.Log;
import barna.commons.utils.StringUtils;
import barna.io.FileHelper;
import barna.model.Qualities;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;
import com.thoughtworks.xstream.XStream;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * Error model driver class to create error models from GEM mapping files.
 * <p>
 * For easy access to a stored model, use the {@link #loadErrorModel(java.io.File)}
 * method.
 * </p>
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */

public class MarkovErrorModel implements FluxTool {
    /**
     * The mapping file
     */
    private File file;
    /**
     * The read length
     */
    private int readLength;
    /**
     * The output file
     */
    private File output;
    /**
     * Limit number of reads that will be taken into account
     */
    private int limit = -1;
    /**
     * The technology used
     */
    private Qualities.Technology technology;

    /**
     * Get the mapping file
     *
     * @return file the mapping file
     */
    public File getFile() {
        return file;
    }

    /**
     * Set the mapping file
     *
     * @param file the mapping file
     */
    public void setFile(File file) {
        this.file = file;
    }

    /**
     * Get the readlength. Returns -1 if not set explicitly. In this case it is
     * guessed.
     *
     * @return readLength the read length
     */
    public int getReadLength() {
        return readLength;
    }

    /**
     * Set the read length
     *
     * @param readLength the read length
     */
    public void setReadLength(int readLength) {
        if (readLength <= 0) {
            throw new IllegalArgumentException("Read length <= 0 not permitted!");
        }
        this.readLength = readLength;
    }

    /**
     * Get the output file name
     *
     * @return output the output filename
     */
    public File getOutput() {
        return output;
    }

    /**
     * Set the output filename
     *
     * @param output the output filename
     */
    public void setOutput(File output) {
        this.output = output;
    }

    /**
     * Returns the current sequence limit or -1
     *
     * @return limit the limit or -1
     */
    public int getLimit() {
        return limit;
    }

    /**
     * Set the current sequence limit
     *
     * @param limit the limit
     */
    public void setLimit(int limit) {
        if (limit <= 0) {
            throw new IllegalArgumentException("Limit <= 0 not permitted!");
        }
        this.limit = limit;
    }

    /**
     * Set the technology. This is used to properly translate quality characters
     *
     * @param technology the technology
     */
    public void setTechnology(String technology) {
        if (technology != null && technology.length() > 0) {
            technology = technology.trim();
            if (technology.equalsIgnoreCase("phred")) {
                this.technology = Qualities.Technology.Phred;
            } else if (technology.equalsIgnoreCase("solexa")) {
                this.technology = Qualities.Technology.Solexa;
            } else if (technology.equalsIgnoreCase("illumina13")) {
                this.technology = Qualities.Technology.Illumina13;
            } else if (technology.equalsIgnoreCase("illumina18")) {
                this.technology = Qualities.Technology.Illumina18;
            }
        }
    }

    @Override
    public String getName() {
        return "errormodel";
    }

    @Override
    public String getDescription() {
        return "Create an error model";
    }

    @Override
    public List<Parameter> getParameter() {
        ArrayList<Parameter> parameters = new ArrayList<Parameter>();
        parameters.add(JSAPParameters.flaggedParameter("file", 'f').type(File.class).help("gem mapping input file").valueName("map").required().get());
        parameters.add(JSAPParameters.flaggedParameter("output", 'o').type(File.class).help("output file name").valueName("output").required().get());
        parameters.add(JSAPParameters.flaggedParameter("tech", 't').help("Technology [phred|solexa|illumina13|illumina18]").valueName("tech").required().get());
        parameters.add(JSAPParameters.flaggedParameter("limit", 's').help("read limit number of sequences to create the model").get());
        parameters.add(JSAPParameters.flaggedParameter("length", 'l').help("read length").get());
        return parameters;

    }

    @Override
    public boolean validateParameter(JSAPResult args) {
        setFile(args.getFile("file"));
        setOutput(args.getFile("output"));
        setTechnology(args.getString("tech"));
        if(args.userSpecified("limit"))setLimit(args.getInt("limit"));
        if(args.userSpecified("length"))setReadLength(args.getInt("length"));

        if (getFile() == null) {
            Log.error("No input file specified!\n");
            return false;
        } else if (!getFile().exists()) {
            Log.error(getFile().getAbsolutePath() + " does not exist!\n");
            return false;
        }
        if (getOutput() == null) {
            Log.error("Please specify an output file!\n");
            return false;
        }

        if (technology == null) {
            Log.error("No technology specified, please specify the technology used to create the quality scores!\n");
            return false;
        }

        return true;
    }

    public Object call() throws Exception {
        if (getReadLength() == 0) {
            Log.message("No read length specified (-l). Trying to figure the read length");
            MapFileReader reader = new MapFileReader(getFile(), technology);
            int reads = 0;
            // check the first 10 reads
            while (reads++ < 10) {
                Read read = reader.parseNext(false);
                if (read == null) break;
                if (readLength == 0) {
                    readLength = read.getLength();
                } else {
                    if (readLength != read.getLength()) {
                        Log.error("Looks like your reads are of different length ... sorry ... I can not handle that!");
                        return null;
                    }
                }
            }
            reader.close();

            if (readLength == 0) {
                Log.error("Unable to figure out the read length");
                return null;
            }
            Log.message("Read length " + readLength);

        }

        int numStates = Qualities.PHRED_RANGE[1];
        switch (technology) {
            case Phred:
                numStates = Qualities.PHRED_RANGE[1];
                break;
            case Solexa:
                numStates = Qualities.SOLEXA_RANGE[1];
                break;
            case Illumina13:
                numStates = Qualities.ILLUMINA_13_RANGE[1];
                break;
            case Illumina18:
                numStates = Qualities.ILLUMINA_18_RANGE[1];
                break;
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
        CharacterQualityDistribution dA = new CharacterQualityDistribution('A', numStates);
        CharacterQualityDistribution dC = new CharacterQualityDistribution('C', numStates);
        CharacterQualityDistribution dG = new CharacterQualityDistribution('G', numStates);
        CharacterQualityDistribution dT = new CharacterQualityDistribution('T', numStates);
        CharacterQualityDistribution dN = new CharacterQualityDistribution('N', numStates);

        int c = 0;
        Read read = null;
        int sum = limit;
        if (sum <= 0) {
            long ll = FileHelper.countLines(getFile());
            if(ll > Integer.MAX_VALUE) throw new RuntimeException("Unable to cache " + ll + " elements... value > Integer.MAX_VALUE");
            sum = (int) ll;
        }
        while ((read = reader.parseNext(false)) != null && (limit < 0 || c < limit)) {
            if (read.getMappings() == null || read.getMappings().size() == 0) {
                continue;
            }

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
            Log.progress(c++, sum);
        }
        Log.progressFinish(StringUtils.OK, true);

        Log.info("Writing model to " + getOutput().getAbsolutePath());
        OutputStream out = new GZIPOutputStream(new FileOutputStream(getOutput()));
        // prepare model
        QualityErrorModel qualityErrorModel = new QualityErrorModel(technology, readLength, trans, crossTalkQuality);
        XStream ss = createXStream();
        ss.toXML(qualityErrorModel, out);
        out.close();


        // try to read just to make sure that
        // the model is valid
        try {
            Log.info("Validating model");
            QualityErrorModel restoredModel = loadErrorModel(getOutput());
            if (restoredModel == null) {
                throw new RuntimeException("Unable to load created model. Something went wrong while saving !");
            }
        } catch (Exception e) {
            throw new RuntimeException("Unable to load created model. Something went wrong while saving !", e);
        }

        Log.info("\tError model stats");
        Log.info("\t\tAverage mutations " + crossTalkQuality.getAverageMutations());
        Log.info("\t\tAverage quality " + qualityDistribution.getAverageQuality());

        return null;
    }

    /**
     * Load an error model from file
     *
     * @return model the model
     * @throws IOException in case of any errors
     */
    public static QualityErrorModel loadErrorModel(File file) throws IOException {
        return loadErrorModel(file.getAbsolutePath(), new FileInputStream(file));
    }

    /**
     * Load an error model from file
     *
     * @return input the input stream
     * @throws IOException in case of any errors
     */
    public static QualityErrorModel loadErrorModel(String name, InputStream inputStream) throws IOException {
        Log.info("Reading error model " + (name == null ? "" : name));
        GZIPInputStream gz = new GZIPInputStream(inputStream);
        XStream xx = createXStream();
        QualityErrorModel trans = (QualityErrorModel) xx.fromXML(gz);
        gz.close();
        return trans;
    }

    static XStream createXStream() {
        XStream stream = new XStream();
        stream.alias("fbi.genome.sequencing.rnaseq.simulation.error.QualityErrorModel", QualityErrorModel.class);
        stream.alias("QualityErrorModel", QualityErrorModel.class);
        return stream;
    }
}