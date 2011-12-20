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

import barna.commons.ByteArrayCharSequence;
import barna.commons.io.IOHandler;
import barna.commons.io.IOHandlerFactory;
import barna.model.Qualities;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Read .map files
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class MapFileReader {

    /**
     * Pattern no match missmatches in the .map file
     */
    private static final Pattern MISSMATCH_PATTERN = Pattern.compile("([A,C,G,T,N,a,c,g,t,n])(\\d{1,3}+)");
    /**
     * Match the quality string
     */
    private static final Pattern QUALITY_PATTERN = Pattern.compile("\\@(\\d+)/(\\d+)");
    /**
     * The input file
     */
    private File file;
    /**
     * The technology
     */
    private Qualities.Technology qualityTechnology;
    /**
     * Current read
     */
    private Read read;
    /**
     * The input stream
     */
    private BufferedInputStream stream;
    /**
     * The io handler
     */
    private IOHandler io;

    public MapFileReader(File file, Qualities.Technology qualityTechnology) {
        this.file = file;
        this.qualityTechnology = qualityTechnology;
        this.read = new Read();

    }

    Read parseNext(boolean skip) throws IOException {
        if (stream == null) {
            stream = new BufferedInputStream(new FileInputStream(file));
            io = IOHandlerFactory.createDefaultHandler();
            io.addStream(stream);
        }

        ByteArrayCharSequence line = io.readLine(stream);
        if (line == null) return null;
        if (skip) return read;


        read.reset();

        // parse the line and
        String[] split = line.toString().split("\t");
        read.setName(split[0]);
        read.setSequence(split[1]);
        // todo check what happens with files without any qualities
        String quals = split[2];
        String mappings = split[3];
        if (!mappings.equals("0:0:0")) {

            // extract all best mappings
            String[] splitedMappings = mappings.split(":");

            for (int i = 0; i < splitedMappings.length; i++) {
                String splitedMapping = splitedMappings[i];
                int n = Integer.parseInt(splitedMapping);
                if (n > 0) {
                    // this is the best mapping with index n mappings and i missmatches

                    String missmatches = split[4];
                    String[] splittedMatches = missmatches.split(",");
                    for (int j = 0; j < n && n <= splittedMatches.length; j++) {
                        Read.Mapping mapping = read.addMapping();
                        String splittedMatch = splittedMatches[j];
                        Matcher matcher = MISSMATCH_PATTERN.matcher(splittedMatch);
                        int missmatchStart = -1;
                        while (matcher.find()) {
                            if (missmatchStart < 0) missmatchStart = matcher.start();
                            mapping.addMissmatch(Integer.parseInt(matcher.group(2)), matcher.group(1).charAt(0));
                        }

                        Matcher qualityMatcher = QUALITY_PATTERN.matcher(splittedMatch);
                        if (qualityMatcher.find()) {
                            if (missmatchStart < 0) missmatchStart = qualityMatcher.start();
                            int avgQuality = Integer.parseInt(qualityMatcher.group(1));
                            int num = Integer.parseInt(qualityMatcher.group(2));
                            mapping.setQuality(avgQuality);
                        }

                        if (missmatchStart >= 0) {
                            // get the name
                            mapping.setName(splittedMatch.substring(0, missmatchStart));
                        }

                    }

                    break; // only best
                }

            }
        }


        for (int i = 0; i < read.getLength(); i++) {
            char q = quals.charAt(i);
            read.getQualities()[i] = Qualities.quality(qualityTechnology, q);
        }

        return read;
    }

    public void close() {
        if (io != null) {
            io.close();
        }
    }


}
