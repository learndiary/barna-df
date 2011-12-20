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

package barna.commons.launcher;

import barna.commons.log.Log;

/**
 * Static helper methods for command line tools. This class is bound to the {@link Log} implementation
 * to check for quiet state and log levels.
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class CommandLine {
    /**
     * Yes options for user confirmation
     */
    public static final String[] user_yes = new String[]{"yes", "y", "si", "yo", "ja", "ok"};
    /**
     * No options for user confirmation
     */
    public static final String[] user_no = new String[]{"no", "n", "nein", "nope", "nix", "noe"};


    /**
     * Ask the user to confirm the given question. If logging is disabled (shutup mode), this will always return true
     *
     * @param message the question (null permitted)
     * @return confirmed true if the user confirmed with yes
     */
    public static boolean confirm(String message) {
        if (Log.getLogLevel() == Log.Level.NONE || !Log.isInteractive()) {
            return true;
        }

        while (true) {
            StringBuilder sb = new StringBuilder();
            int in;

            if (message != null) {
                System.err.print(message + " ");
                System.err.flush();
            }

            try {
                while ((in = System.in.read()) != '\n') {
                    sb.append((char) in);
                }
            } catch (Exception ignore) {
                // ignore exceptions
            }

            String s = sb.toString().toLowerCase().trim();
            for (String user_ye : user_yes) {
                if (s.equals(user_ye)) {
                    return true;
                }
            }
            for (String anUser_no : user_no) {
                if (s.equals(anUser_no)) {
                    return false;
                }
            }
        }
    }

}
